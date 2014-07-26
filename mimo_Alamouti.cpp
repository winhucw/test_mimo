// MIMO + OFDM
// C.W
//
#include <itpp/itcomm.h>
using namespace itpp;
using namespace std;

int main(int argc,char **argv) {
    // -- modulation and channel parameters (taken from command line input) --
    int nC;                    // type of constellation  (1=QPSK, 2=16-QAM, 3=64-QAM)
    int nRx;                   // number of receive antennas
    int nTx;                   // number of transmit antennas
//    int Tc = 512;              // symbol durations, multiple of T, T<=coherence_time<=tx_duration, T is the ST code duration
    int nb_iter = 5;
    string map_metric = "maxlogMAP";
    string code_name = "Alamouti_2xN";
    string demapper_method = "Alamouti_maxlogMAP";
    ivec gen = "0133 0165 0171";
    int constraint_len = 7;
    double rate = 1.0 / 3.0;
    const vec EbN0db = "-5:1:5";           // SNR range
    const int Nbits = int(1e6);           // number of bits to ever simulate per SNR point
    
    if (argc != 4) {
        cout << "Usage: ./mimo nTx nRx nC" << endl << "Example: ./mimo 2 2 1" << endl;
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
    STC stc(code_name, 1 << (2*nC));
    nTx = stc.get_nb_emission_antenna();
    const int channel_uses = stc.get_channel_uses();
    const int symb_block = stc.get_nb_symbols_per_block();
    
    // -- Convolutional code parameters --
    Convolutional_Code code;
    code.set_generator_polynomials(gen, constraint_len);
    bvec dummy;
    code.encode(randb(Nbits), dummy);
    const int Nc = length(dummy);
    const int Nctx = (int)(2 * nC * symb_block * ceil(double(Nc) / double(2 * nC * symb_block)));
    
    // -- SISO demapper --
    SISO siso;
    siso.set_map_metric(map_metric);
    siso.set_generators(gen, constraint_len);
    siso.set_constellation(chan.bits_per_symbol()(0), chan.get_symbols()(0), chan.get_bits2symbols()(0));
    siso.set_demapper_method(demapper_method);
    siso.set_st_block_code(stc.get_nb_symbols_per_block(), stc.get_1st_gen_matrix(), stc.get_2nd_gen_matrix(), nRx);
    
    // -- SISO decoder --
    vec dem_extridata(Nctx);
    vec dem_apriodata(Nctx);
    dem_apriodata.zeros();
    
    vec nsc_intrinsic_coded(Nctx);
	vec nsc_apriori_data(Nbits);
	nsc_apriori_data.zeros();
	vec nsc_extrinsic_coded(Nctx);
	vec nsc_extrinsic_data(Nbits);
    
    
    // -- ofdm modulation --
    const int Nfft = 1024; // 1024
    const int Ncp = 128;  // 128, 1/8 of the FFT size
    OFDM ofdm;
    ofdm.set_parameters(Nfft,Ncp);

    RNG_randomize();
    
    BPSK bpsk;
    BERC berc;
    vec ber(EbN0db.length());
    bvec rec_bits(Nbits);
    
    cout.setf(ios::fixed, ios::floatfield);
    cout.setf(ios::showpoint);
    cout.precision(5);
    
    // ================== Run simulation =======================
    for (int nsnr = 0; nsnr < length(EbN0db); nsnr++) {
        cout << " --------- Run " << nsnr << " time -------- " << endl;
        const double Eb = 1.0;      // transmitted energy per information bit
        const double N0 = inv_dB(-EbN0db(nsnr));
        const double sigma2 = N0;   // Variance of each scalar complex noise sample
        const double Es = rate * 2 * nC * Eb; // Energy per complex scalar symbol
        
        
        awgn.set_noise(N0);
        siso.set_noise(N0 / 2);
            
        // -- generate and encode random data --
        bvec inputbits = randb(Nbits);
        bvec txbits = code.encode_tail(inputbits);

        // -- mapping --
        if (Nc != Nctx) {
            txbits = concat(txbits, randb(Nctx - Nc));
        }
        const int Nmap = Nctx / (2 * nC);  // length of complex symbols after mapping
        cvec inputsymbol(Nmap);
        for (int k = 0; k < Nmap; k++) {
            bvec bitstmp = txbits(k * (2 * nC), (k + 1) * (2 * nC) - 1);
            inputsymbol.set_subvector(k, chan.modulate_bits(bitstmp));
        }
        
        // -- STC --
        const int Ns = Nmap / symb_block * channel_uses; // length of complex symbols after stc for each antenna
        cmat symbol_stc = stc.encode(inputsymbol);

        
        // -- IFFT --
        const int Nf = (Ns / Nfft + 1) * Nfft; // length of complex symbols, multiple of fft
        int add;
        if (Ns % Nfft != 0) {
            add = Nf - Ns;
            cvec addzeros = zeros_c(symb_block);
            for (int i = 0; i < add; i++)
                symbol_stc.append_row(addzeros);
        }
        cmat txsymbol(Nf / Nfft * (Nfft + Ncp), symb_block); // complex transmitted data after ifft and cp; will pass through channel
        for (int k = 0; k < nTx; k++) {
            txsymbol.set_col(k, ofdm.modulate(symbol_stc.get_col(k)));
        }
  
        // -- generate channel and data ----
        const int tx_duration = txsymbol.rows();
        cmat H(nTx * nRx, tx_duration) ;
        H.ones();
        H *= sqrt(Es);
        
        siso.set_impulse_response(H);
        cmat Y(tx_duration, nRx);
        
        for (int i = 0; i < tx_duration / channel_uses; i++ ) {  // Y = Hx + n
            Y.set_submatrix(i * channel_uses, 0, txsymbol(i * channel_uses, (i+1) * channel_uses - 1, 0, nTx - 1) * reshape(H.get_col(i), nTx, nRx));
        }
        Y = awgn(Y);
        
        // -- FFT and demodulate --
        cmat rxsymbol(Nf, symb_block);  // complex received data after fft
        for (int k = 0; k < nTx; k++) {
            rxsymbol.set_col(k, ofdm.demodulate(Y.get_col(k)));
        }
        rxsymbol.del_rows(Nf - add, Nf - 1); // extra adds before to meet multiple of fft

        siso.demapper(dem_extridata, rxsymbol, dem_apriodata);
        nsc_intrinsic_coded = dem_extridata;
            
        siso.nsc(nsc_extrinsic_coded, nsc_extrinsic_data, nsc_intrinsic_coded, nsc_apriori_data);
        rec_bits = bpsk.demodulate_bits(-nsc_extrinsic_data);
            
        berc.clear();
        berc.count(inputbits, rec_bits);
        ber(nsnr) = berc.get_errorrate();
        cout << ber(nsnr) << endl;
    }
    cout.flush();
    return 0;
}
