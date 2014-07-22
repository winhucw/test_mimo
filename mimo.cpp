#include <itpp/itcomm.h>
using namespace itpp;
using namespace std;

int main(int argc,char **argv) {
    // -- modulation and channel parameters (taken from command line input) --
    int nC;                    // type of constellation  (1=QPSK, 2=16-QAM, 3=64-QAM)
    int nRx;                   // number of receive antennas
    int nTx;                   // number of transmit antennas
    int Tc;                    // coherence time
    
    if (argc != 5) {
        cout << "Usage: ./mimo nTx nRx nC Tc" << endl << "Example: ./mimo 2 2 1 18 " << endl;
        exit(1);
    }
    else {
        sscanf(argv[1], "%i", &nTx);
        sscanf(argv[2], "%i", &nRx);
        sscanf(argv[3], "%i", &nC);
        sscanf(argv[4], "%i", &Tc);
    }
    
    // -- STC --
    string code_name = "Alamouti_2xN";
    STC st_block_code(code_name, 1 << (2*nC));
    nTx = st_block_code.get_nb_emission_antenna();
    int channel_uses = st_block_code.get_channel_uses();
    const int out_symb_block = st_block_code.get_nb_emission_antenna();
    
    cout << nTx << " TX antennas, " << nRx << " RX antennas, " << (1 << nC) << "-PAM per dimension, coherence time " << Tc << endl;
    // -- simulation control parameters --
    const vec EbN0db = "5:20:15";           // SNR range
    const int Nmethods = 2;                 // number of demodulators(ML,ZF) to try
    const int Nbits = 10;           // number of bits to ever simulate per SNR point
    
    // -- Channel code parameters --
    Convolutional_Code code;
    ivec generator(3);
    generator(0) = 0133;                            // use rate 1/3 code
    generator(1) = 0165;
    generator(2) = 0171;
    double rate = 1.0 / 3.0;
    code.set_generator_polynomials(generator, 7);   // The encoder has constraint length 7.
    
    // -- initialize MIMO channel with uniform QAM per complex dimension and Gray coding --
    ND_UQAM chan;
    chan.set_M(nTx, 1 << (2*nC));
    //cout << chan << endl;
    
    // -- ofdm modulation --
    int Nfft = 1024;
    int Ncp = 128;  // 1/8 of the FFT size
    OFDM ofdm;
    ofdm.set_parameters(Nfft,Ncp);

    RNG_randomize();
    
    ivec Contflag = ones_i(Nmethods);   // flag to determine whether to run a given demodulator
    if (pow(2.0, nC*2.0*nTx) > 256) {   // ML decoder is too complex
        Contflag(1) = 0;
    }
    if (nTx > nRx) {
        Contflag(0) = 0;                  // ZF not for underdetermined systems
    }
    cout << "Running methods: " << Contflag << endl;
    
    cout.setf(ios::fixed, ios::floatfield);
    cout.setf(ios::showpoint);
    cout.precision(5);
    
    // ================== Run simulation =======================
    for (int nsnr = 0; nsnr < length(EbN0db); nsnr++) {
        const double Eb = 1.0;      // transmitted energy per information bit
        const double N0 = inv_dB(-EbN0db(nsnr));
        const double sigma2 = N0;   // Variance of each scalar complex noise sample
        const double Es = rate * 2 * nC * Eb; // Energy per complex scalar symbol
        const int Nbitspvec = 2 * nC * nTx;        // Number of bits per channel vector
        Array<BERC> berc(Nmethods);     // counter for coded BER
            
        // -- generate and encode random data --
        bvec inputbits = randb(Nbits);
        bvec txbits;
        code.encode_tail(inputbits, txbits);
        const int Nc = length(txbits);
            
        // -- mapping and IFFT --
        const int Nmap = Nc / (2 * nC);  // original length of complex symbols
        cvec inputfft(Nmap);
        cvec txfft;
        cvec outputfft;
        for (int k = 0; k < Nmap / nTx ; k++) {
            bvec bitstmp = txbits(k * Nbitspvec, (k + 1) * Nbitspvec - 1);
            cvec symbols = chan.modulate_bits(bitstmp);
            inputfft.set_subvector(k * nTx, symbols);
        }
        int nf = Nmap;
        if (Nmap % Nfft != 0) {
            nf = (Nmap / Nfft + 1) * Nfft;  // multiple of Nfft
            int addzeros = nf - Nmap;
            inputfft = concat(inputfft, zeros_c(addzeros));
        }
        cout << "----- The " << nsnr << " time, Nf (should be a multiple of 1024) " << nf << endl;
        ofdm.modulate(inputfft, txfft);
            
        // -- generate channel and data ----
        const int nfc = txfft.size(); // multiple of (Nfft + Ncp)
        cout << "----- txfft.size (should be a multiple of 1152) " << nfc << endl;
        const int Nvec = nfc / nTx;
        cvec rxfft(nfc);
        Array<cvec> Y(Nvec);            // received data
        Array<cmat> H(Nvec / Tc + 1);   // channel matrix (new matrix for each coherence interval)
            
        for (int k = 0; k < Nvec; k++) {
            if (k % Tc == 0) {
                H(k / Tc) = sqrt(Es) * randn_c(nRx, nTx);
            }
            cvec x = txfft(k * nTx, (k + 1) * nTx - 1);
            cvec e = sqrt(sigma2) * randn_c(nRx);
            Y(k) = H(k / Tc) * x + e;
            rxfft.set_subvector(k * nTx, Y(k)) ;
        }
            
        // -- FFT and demodulate --
        ofdm.demodulate(rxfft, outputfft);
        cout << "----- outputfft.size (should be a multiple of 1024) " << outputfft.size() << endl;
            
        Array<QLLRvec> LLRin(Nmethods);
        for (int i = 0; i < Nmethods; i++) {
            LLRin(i) = zeros_i(nf);
        }
            
        QLLRvec llr_apr = zeros_i(Nbitspvec);
        QLLRvec llr_apost = zeros_i(Nbitspvec);
        for (int k = 0; k < Nmap / nTx; k++) {
            Y(k) = outputfft(k * nTx, (k+1) * nTx -1);
            // zero forcing demodulation
            if (Contflag(0)) {
                chan.demodulate_soft_bits(Y(k), H(k / Tc), sigma2, llr_apr, llr_apost, ND_UQAM::ZF_LOGMAP);
                LLRin(0).set_subvector(k * Nbitspvec, llr_apost);
            }
            
            // ML demodulation
            if (Contflag(1)) {
                chan.demodulate_soft_bits(Y(k), H(k / Tc), sigma2, llr_apr, llr_apost);
                LLRin(1).set_subvector(k * Nbitspvec, llr_apost);
            }
        }

            // -- decode and count errors --
            for (int i = 0; i < Nmethods; i++) {
                bvec decoded_bits;
                if (Contflag(i)) {
                    // QLLR values must be converted to real numbers (for convolutional decoder)
                    vec llr = chan.get_llrcalc().to_double(LLRin(i).left(nf));
                    code.decode_tail(llr, decoded_bits);
                    berc(i).count(inputbits(0, Nbits - 1), decoded_bits(0, Nbits - 1));  // coded BER
                }
            }
        
        cout << "-----------------------------------------------------" << endl;
        cout << "Eb/N0: " << EbN0db(nsnr) << " dB. Simulated " << Nbits << " bits." << endl;
        cout << "Coded BER:   " << berc(0).get_errorrate()  << " (ZF);     " << berc(1).get_errorrate()  << " (ML)" << endl;
        cout.flush();
    }
    return 0;
}
