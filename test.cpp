#include <itpp/itcomm.h>
using namespace itpp;
using namespace std;

int main(int argc,char **argv) {
    TDL_Channel chan;
    double norm = 0.1;
    chan.set_norm_doppler(norm);
    int sample = 10;
    cmat ch_coeff;
    chan.generate(sample, ch_coeff);
    cout << ch_coeff << endl;
    return 0;
}
