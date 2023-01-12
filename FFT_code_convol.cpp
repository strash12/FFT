#include <iostream>     
#include <complex>      
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>
 
using namespace std;
 
double PI = acos(0) * 2;
typedef complex<double> xd;
typedef vector<double> dvec;
typedef vector<xd> xvec;
const xd J(0, 1); // sqrt(-1)
 
inline xd omega(const double &p, const double &q)
{
    return exp((2. * PI * J * q) / p);
}
 
xvec _fft(xvec &f)
{
    double N = f.size();
     
    if (N == 1) return f;
     
    xvec fe, fo;
    fe.reserve(N / 2);
    fo.reserve(N / 2);
     
    for (int i = 0; i < N; i += 2) {
        fe.push_back(f[i]);     // even
        fo.push_back(f[i + 1]); // odd
    }
     
    fe = _fft(fe);
    fo = _fft(fo);
     
    for (int m = 0; m < N / 2; ++m) {
        xd omfo = omega(N, -m) * fo[m];
        f[m]         = fe[m] + omfo;
        f[m + N / 2] = fe[m] - omfo;
    }
 
    return f;
}
 
xvec fft(const dvec &x)
{
    xvec f(x.size());
     
    for (size_t i = 0; i < x.size(); ++i) { 
        f[i] = xd(x[i], 0);
    }
     
    return _fft(f);
}
 
xvec _ifft(xvec &x)
{
    double N = x.size();
     
    if (N == 1) return x;
     
    xvec xe, xo;
    xe.reserve(N / 2);
    xo.reserve(N / 2);
     
    for (int i = 0; i < N; i += 2) {
        xe.push_back(x[i]);     // even
        xo.push_back(x[i + 1]); // odd
    }
     
    xe = _ifft(xe);
    xo = _ifft(xo);
     
    for (int m = 0; m < N / 2; ++m) {
        xd iomxo = omega(N, m) * xo[m];
        x[m]         = xe[m] + iomxo;
        x[m + N / 2] = xe[m] - iomxo;
    }
     
    return x;
}
 
dvec ifft(xvec f)
{
    double N = f.size();
     
    xvec xcomplex = _ifft(f);
    dvec x(N);
     
    for (int i = 0; i < N; ++i) {
        x[i] = xcomplex[i].real() / N;
    }
     
    return x;
}
 
// vector convolution
dvec convolve(const dvec &a, const dvec &b) 
{
    // calculate degree of resulting polynomial
    size_t N = 2 * a.size() - 1;
     
    // extend size to match result
    dvec acof(N), bcof(N);
    copy(a.begin(), a.end(), acof.begin());
    copy(b.begin(), b.end(), bcof.begin());
     
    xvec apv, bpv, cpv(N);
     
    // evaluation
    apv = fft(acof);
    bpv = fft(bcof);
     
    // point-wise multiplcation
    for (size_t i = 0; i < N; ++i) {
        cpv[i] = apv[i] * bpv[i];
    }
     
    for (const auto &t : cpv)  cout << t << ' ';
    cout << endl;
     
    // interpolation
    return ifft(cpv);
}
int main()
{
    cout << fixed << setprecision(2);
 
    dvec input = { 1,6,3,8,9,5,4,2 };
     
    // convert from time to frequency domain
    xvec freqdom = fft(input);
 
    for (const auto &f : freqdom) {
        cout << f << endl;
    }
    cout << endl;
     
    // convert from frequency to time domain
    auto timedom = ifft(freqdom);
     
    for (const auto &t : timedom) {
        cout << t << ' ';
    }
    cout << endl;
}
