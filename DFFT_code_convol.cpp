#include <iostream>     
#include <complex>      
#include <cmath>
#include <iomanip>
#include <vector>
 
using namespace std;
 
double PI = acos(0) * 2;
typedef complex<double> xd;
typedef vector<double> dvec;
typedef vector<xd> xvec;
const xd J(0, 1); // sqrt(-1)
 
xvec dft(const dvec &input)
{
    double N = input.size();
    xvec X(N);
 
    for (double k = 0; k < N; ++k) {
        for (double n = 0; n < N; ++n) {
            X[k] += (double)input[n] * exp(-2. * J * PI * n * k / N);
        }
    }
 
    return X;
}
 
dvec idft(const xvec &input)
{
    double N = input.size();
    xvec x(N);
    dvec out(N);
 
    for (double k = 0; k < N; ++k) {
        for (double n = 0; n < N; ++n) {
            x[k] += input[n] * exp(2. * J * PI * n * k / N);
        }
        out[k] = x[k].real() / N;
    }
 
    return out;
}

dvec convolve(const dvec &a, const dvec &b) 
{
    // calculate degree of resulting polynomial
    size_t N = 2 * a.size() - 1;
     
    // extend size and pad with 0
    dvec acof(N, 0), bcof(N, 0);
    copy(a.begin(), a.end(), acof.begin());
    copy(b.begin(), b.end(), bcof.begin());
     
    xvec apv, bpv, cpv(N);
     
    // evaluation
    apv = dft(acof);
    bpv = dft(bcof);
     
    // point-wise multiplcation
    for (size_t i = 0; i < N; ++i) {
        cpv[i] = apv[i] * bpv[i];
    }
     
    // interpolation
    return idft(cpv);
}
int main()
{
    cout << fixed << setprecision(2);
 
    dvec input = { 1, 2, 3, 4 };
     
    // convert from time to frequency domain
    xvec freqdom = dft(input);
 
    for (const auto &f : freqdom) {
        cout << f << endl;
    }
    cout << endl;
     
    // convert from frequency to time domain
    dvec timedom = idft(freqdom);
     
    for (const auto &t : timedom) {
        cout << t << ' ';
    }
    cout << endl;
}
