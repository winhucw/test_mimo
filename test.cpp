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
    STC stc(code_name, 1 << (2*nC));
    nTx = stc.get_nb_emission_antenna();
    const int channel_uses = stc.get_channel_uses();
    const int symb_block = stc.get_nb_symbols_per_block();
    
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
    chan.set_M(1, 1 << (2*nC));
    //cout << chan << endl;
    
    ND_UQAM chan2;
    chan2.set_M(nTx, 1<< (2*nC));

    const int Nbits = 8;
    bvec inputbits = "1,0,1,0,1,1,1,1";
    bvec txbits;
    code.encode_tail(inputbits,txbits);
    int Nc = txbits.size();
    cout << " input bit stream: " << inputbits << " , length = " << Nbits << endl;
    cout << " encoded bit stream: " << txbits << " , length = " << Nc << endl;
    
    cvec inputsymbol(Nc/(2*nC));
    cout << " modulated symbol for one Tx: " << endl;
    for (int i = 0; i < Nc / (2*nC); i++) {
        bvec bitstmp = txbits(i*2*nC, (i+1)*2*nC-1);
        cvec symbol = chan.modulate_bits(bitstmp);
        cout << " the " << i << " th symbol " << symbol << endl;
        inputsymbol.set_subvector(i,symbol);
    }
    cout << " length = " << inputsymbol.size() << endl;
    
    cmat symbol_stc = stc.encode(inputsymbol);
    int Ns = symbol_stc.size() / nTx / symb_block * channel_uses;
    cout << " symbol_stc \n" << symbol_stc << endl;
    cout << " estimate = symbol_stc size = " << Ns << " , row = " << symbol_stc.get_row(1).size() << ", column = " << symbol_stc.get_col(1).size() << endl;
    
    if (Ns % 30 != 0) {
        cvec z = zeros_c(nTx);
        for (int i = 0; i< 4; i++) {
            symbol_stc.append_row(z);
        }

        cout << " new symbol_stc \n" << symbol_stc << endl;
        cout << " new symbol_stc row = " << symbol_stc.get_row(1).size() << ", column = " << symbol_stc.get_col(1).size() << endl;
    }
    return 0;
}
