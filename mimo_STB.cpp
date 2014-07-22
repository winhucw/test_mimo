#include <itpp/itcomm.h>
using namespace itpp;
using namespace std;

int main(int argc,char **argv) {
    // -- modulation and channel parameters (taken from command line input) --
    int nC;                    // type of constellation  (1=QPSK, 2=16-QAM, 3=64-QAM)
    int nRx;                   // number of receive antennas
    int nTx;                   // number of transmit antennas
    
    if (argc != 4) {
        cout << "Usage: ./mimo nTx nRx nC Tc" << endl << "Example: ./mimo 2 2 1" << endl;
        exit(1);
    }
    else {
        sscanf(argv[1], "%i", &nTx);
        sscanf(argv[2], "%i", &nRx);
        sscanf(argv[3], "%i", &nC);
    }
    
    // -- initialize MIMO channel with uniform QAM per complex dimension and Gray coding --
    ND_UQAM chan;
    chan.set_M(1, 1 << (2*nC));
    
    AWGN_Channel awgn;

    // -- STC --
    STC stc("Alamouti_2xN", 1 << (2*nC));
    nTx = stc.get_nb_emission_antenna();
    const int channel_uses = stc.get_channel_uses();
    const int symb_block = stc.get_nb_symbols_per_block();
    
    SISO siso;
    cout << "Check : "chan.bits_per_symbol()(0) << "   " << chan.get_symbols()(0) << "   " << chan.get_bits2symbols()(0) << endl;
    siso.set_constellation(chan.bits_per_symbol()(0), chan.get_symbols()(0), chan.get_bits2symbols()(0));
    siso.set_demapper_method("Alamouti_maxlogMAP");
    siso.set_st_block_code(stc.get_nb_symbols_per_block(), stc.get_1st_gen_matrix(), stc.get_2nd_gen_matrix(), nRx);
    
    cout << nTx << " TX antennas, " << nRx << " RX antennas, " << (1 << (2*nC)) << "- QAM " << endl;
    // -- simulation control parameters --
    const vec EbN0db = "5:5:20";           // SNR range
    const int Nbits = 1000;           // number of bits to ever simulate per SNR point
    
    // -- Channel code parameters --
    Convolutional_Code code;
    ivec generator(3);
    generator(0) = 0133;                            // use rate 1/3 code
    generator(1) = 0165;
    generator(2) = 0171;
    double rate = 1.0 / 3.0;
    code.set_generator_polynomials(generator, 7);   // The encoder has constraint length 7.
    bvec dummy;
    code.encode_tail(randb(Nbits), dummy);
    const int Nc = length(dummy);      // find out how long the coded blocks are
    const int Nctx = (int)(2 * nC * ceil(double(Nc) / double(2 * nC)));
    
    // -- ofdm modulation --
    const int Nfft = 1024; // 1024
    const int Ncp = 128;  // 128, 1/8 of the FFT size
    OFDM ofdm;
    ofdm.set_parameters(Nfft,Ncp);

    RNG_randomize();
    
    cout.setf(ios::fixed, ios::floatfield);
    cout.setf(ios::showpoint);
    cout.precision(5);
    
    // ================== Run simulation =======================
    for (int nsnr = 0; nsnr < length(EbN0db); nsnr++) {
        cout << " Run " << nsnr << " time " << endl;
        const double Eb = 1.0;      // transmitted energy per information bit
        const double N0 = inv_dB(-EbN0db(nsnr));
        const double sigma2 = N0;   // Variance of each scalar complex noise sample
        const double Es = rate * 2 * nC * Eb; // Energy per complex scalar symbol
        BERC berc;     // counter for coded BER
        BERC bercu;
        
        awgn.set_noise(N0);
        siso.set_noise(N0 / 2);
            
        // -- generate and encode random data --
        bvec inputbits = randb(Nbits);
        bvec txbits = code.encode_tail(inputbits);
        cout << " Length of bit stream after convolution = " << Nc << endl;
        
        // -- mapping --
        int Nmap = Nc / (2 * nC);  // length of complex symbols after mapping
        if (Nc != Nctx) {
            Nmap++;
            txbits = concat(txbits, randb(Nctx - Nc));
        }
        cout << " Length of complex symbols after mapping = " << Nmap << endl;
        cvec inputsymbol(Nmap);
        for (int k = 0; k < Nmap; k++) {
            bvec bitstmp = txbits(k * (2 * nC), (k + 1) * (2 * nC) - 1);
            inputsymbol.set_subvector(k, chan.modulate_bits(bitstmp));
        }
//        cout << " inputsymbol \n" << inputsymbol(0, 5) << '\n' << inputsymbol(Nmap - 5, Nmap - 1) << endl;
        
        // -- STC --
        int Ns = Nmap / symb_block * channel_uses; // length of complex symbols after stc for each antenna
        cmat symbol_stc = stc.encode(inputsymbol);
        cout << " Check Ns = " << Ns << " , symbol_stc.rows() = " << symbol_stc.rows() << " , symbol_stc.cols() = " << symbol_stc.cols() << endl;
        
        // -- IFFT --
        int add;
        if (Ns % Nfft != 0) {
            add = (int)(Ns / Nfft + 1) * Nfft - Ns;
            cvec addzeros = zeros_c(symb_block);
            for (int i = 0; i < add; i++)
                symbol_stc.append_row(addzeros);
        }
        int Nf = symbol_stc.rows(); // length of complex symbols, multiple of fft
//        cout << " symbol_stc lowest rows : \n" << symbol_stc.get_rows(2045,2047) << endl;
//        cout << " symbol_stc rightest columns \n " << symbol_stc.get_col(1) << endl;
        cmat txsymbol(Nf / Nfft * (Nfft + Ncp), symb_block); // complex transmitted data after ifft and cp; will pass through channel
        for (int k = 0; k < nTx; k++) {
            txsymbol.set_col(k, ofdm.modulate(symbol_stc.get_col(k)));
        }
//        cout << " txsymbol #rows = " << txsymbol.rows() << " , txsymbol # cols = " << txsymbol.cols() << endl;
//        cout << " txsymbol first col : \n" << txsymbol.get_col(0) << endl;
//        cout << " txsymbol first row : \n" << txsymbol.get_row(0) << endl;
  
        // -- generate channel and data ----
        const int tx_duration = txsymbol.rows();  // multiple of (fft + cp)
        cmat H(nTx * nRx, tx_duration) ;   // channel matrix
        H.ones();
        H *= sqrt(Es);
        
        siso.set_impulse_response(H);
        cmat Y(tx_duration, nRx);            // received data
        
        for (int i = 0; i < tx_duration / channel_uses; i++ ) {
            Y.set_submatrix(i * channel_uses, 0, txsymbol(i * channel_uses, (i+1) * channel_uses - 1, 0, nTx - 1) * reshape(H.get_col(i), nTx, nRx));
        }
        Y = awgn(Y);
        
        // -- FFT and demodulate --
        cmat rxsymbol(Nf, symb_block);  // complex received data after fft
        for (int k = 0; k < nTx; k++) {
            rxsymbol.set_col(k, ofdm.demodulate(Y.get_col(k)));
        }
        rxsymbol.del_rows(Nf - add, Nf - 1); // extra adds before to meet multiple of fft
//        cout << " rxsymbol #row = " << rxsymbol.rows() << "rxsymbol #col = " << rxsymbol.cols() << endl;
//        cout << " rxsymbol first col : \n" << rxsymbol.get_col(0) << endl;
//        cout << " rxsymbol first row : \n" << rxsymbol.get_row(0) << endl;
        
        vec dem_extridata(Nctx);
        bvec dem_data(Nc);
        vec dem_apriodata(Nctx);
        siso.demapper(dem_extridata, rxsymbol, dem_apriodata);
        for (int i = 0; i < Nctx; i++) {
            dem_data(i) = (dem_extridata(i) > 0);
        }
//        cout << " # dem_extridata = " << dem_extridata.size() << endl;
        
//        cout << " decode data : \n " << dec_data << " \n # decode data = " << dec_data.size() << endl;
        cout << "txbits \n" << txbits(0,9) << endl;
        cout << "dem data \n" << dem_data(0,9) << endl;
        cout << "Uncoded : # txbits = " << txbits.size() << " , # dem_data = " << dem_data.size() << endl;
        bercu.count(txbits, dem_data);
        cout << " --> BER for uncoded = " << bercu.get_errorrate() << endl;
        
        bvec dec_data;
        dem_extridata.del(Nctx - (Nctx - Nc), Nctx - 1);
        code.decode_tail(dem_extridata, dec_data);
        
        cout << "input bits \n" << inputbits(0,9) << endl;
        cout << "dec bits \n" << dec_data(0,9) << endl;
        berc.count(inputbits, dec_data);
        cout << "Coded : # inputbits = " << inputbits.size() << " , # dec_data = " << dec_data.size() << endl;
        cout << " --> BER for coded = " << berc.get_errorrate() << endl;
        cout.flush();
    }
    return 0;
}
