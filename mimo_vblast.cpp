// 0728~
// MIMO + vblast
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
    int channel_uses = 1;
    string map_metric = "maxlogMAP";
    string code_name = "V-BLAST_MxN";
    string demapper_method = "zfPIC";  // zfPIC, Hassibi_maxlogMAP
    ivec gen = "0133 0165 0171";
    int constraint_len = 7;
    double rate = 1.0 / 3.0;
    const vec EbN0db = "-5:2:15";           // SNR range
    const int Nbits = int(1e3);           // number of bits to ever simulate per SNR point
    
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
    STC stc(code_name, 1 << (2*nC), nTx, channel_uses);
    nTx = stc.get_nb_emission_antenna();
    channel_uses = stc.get_channel_uses();
    const int symb_block = stc.get_nb_symbols_per_block();
    
    // -- Convolutional code parameters --
    Convolutional_Code code;
    code.set_generator_polynomials(gen, constraint_len);
    bvec dummy, dummy2;
   // code.encode(randb(Nbits), dummy);
   // const int Nc = length(dummy);
   // const int Nctx = (int)(2 * nC * symb_block * ceil(double(Nc) / double(2 * nC * symb_block)));
    code.encode(randb(Nbits/nTx), dummy2);
    const int Nc = length(dummy2);
    const int Nctx = (int)(2 * nC * symb_block * ceil(double(Nc) / double(2 * nC * symb_block)));
    
    // -- SISO demapper --
    SISO siso;
    siso.set_map_metric(map_metric);
    siso.set_generators(gen, constraint_len);
    siso.set_constellation(chan.bits_per_symbol()(0), chan.get_symbols()(0), chan.get_bits2symbols()(0));
    siso.set_demapper_method(demapper_method);
    siso.set_st_block_code(stc.get_nb_symbols_per_block(), stc.get_1st_gen_matrix(), stc.get_2nd_gen_matrix(), nRx);
    
    RNG_randomize();
    
    BPSK bpsk;
    BERC berc;
    vec ber(EbN0db.length());
    
    cout.setf(ios::fixed, ios::floatfield);
    cout.setf(ios::showpoint);
    cout.precision(5);
    
    // ================== Run simulation =======================
    for (int nsnr = 0; nsnr < length(EbN0db); nsnr++) {
 //       cout << " --------- Run " << nsnr << " time -------- " << endl;
        const double Eb = 1.0;      // transmitted energy per information bit
        const double N0 = inv_dB(-EbN0db(nsnr));
        const double sigma2 = N0;   // Variance of each scalar complex noise sample
        const double Es = rate * 2 * nC * Eb; // Energy per complex scalar symbol
        
        
        awgn.set_noise(N0);
        siso.set_noise(N0 / 2);
        
        // -- generate and encode random data --
        bvec inputbits = randb(Nbits);
        bmat bits_tx = reshape(inputbits, nTx, Nbits / nTx).transpose();
        bmat txbits(Nctx, nTx);
        for (int i = 0; i < nTx; i++) {
           txbits.set_col(i,code.encode(bits_tx.get_col(i)));
        }
        
        // -- mapping --
        if (Nc != Nctx) {
            int add = Nctx - Nc;
            for (int i = 0; i < add; i++)
                txbits.append_row(randb(nTx));
        }
        const int Nmap = Nctx / (2 * nC) * nTx;  // length of complex symbols after mapping
        const int Nm = Nctx / (2 * nC); // length of complex symbols per antenna after mapping
        cvec inputsymbol(Nmap);
        for (int k = 0; k < Nm; k++) {
            for (int i = 0; i < nTx; i++) {
                bvec bitstmp = txbits.get_col(i)(k*2*nC,(k+1)*2*nC-1);
                inputsymbol.set_subvector(k*nTx+i, chan.modulate_bits(bitstmp));
            }
        }
        // -- STC --
        const int Ns = Nmap / symb_block * channel_uses; // length of complex symbols after stc for each antenna
        cmat symbol_stc = stc.encode(inputsymbol);
        
        
        // -- generate channel and data ----
        const int tx_duration = symbol_stc.rows();
        cmat H(nTx * nRx, tx_duration);
    //    H.ones();
    //    H *= sqrt(Es);
        H = randn_c(nTx * nRx, tx_duration);
        H *= sqrt(Es);
        cmat Y(tx_duration, nRx); 
        
        for (int i = 0; i < tx_duration / channel_uses; i++ ) {  // Y = Hx + n  
          Y.set_submatrix(i * channel_uses, 0, symbol_stc(i * channel_uses, (i+1) * channel_uses - 1, 0, nTx - 1) * reshape(H.get_col(i), nTx, nRx));
        }
        Y = awgn(Y);
        
        siso.set_impulse_response(H);


        // -- SISO decoder --
        vec dem_extridata(Nctx * nTx);
        vec dem_apriodata(Nctx * nTx);
        dem_apriodata.zeros();
        
        //vec nsc_intrinsic_coded(Nctx * nTx);
        mat nsc_intrinsic_coded(Nctx, nTx);
        vec nsc_apriori_data(Nbits / nTx);
        nsc_apriori_data.zeros();
        vec nsc_extrinsic_coded(Nctx);
        vec nsc_extrinsic_data(Nbits / nTx);
        mat data(Nbits / nTx, nTx);
        vec rec_data(Nbits);
        bvec rec_bits(Nbits);

        siso.demapper(dem_extridata, Y, dem_apriodata);
        
        for (int i = 0; i < Nm; i++) {
            nsc_intrinsic_coded.set_submatrix(i*2*nC, 0, reshape(dem_extridata(i*2*nC*nTx, (i+1)*2*nC*nTx -1), 2*nC, nTx));
        }
        
        for (int k = 0; k < nTx; k++) {
            siso.nsc(nsc_extrinsic_coded, nsc_extrinsic_data, nsc_intrinsic_coded.get_col(k), nsc_apriori_data);
            data.set_col(k, nsc_extrinsic_data);
        }
        rec_data = rvectorize(data);
        rec_bits = bpsk.demodulate_bits(-rec_data);
        
        berc.clear();
        berc.count(inputbits, rec_bits);
        ber(nsnr) = berc.get_errorrate();
        cout << ber(nsnr) << endl;
    }
    cout.flush();
    return 0;
}
