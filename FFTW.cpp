#include  <stdio.h>
#include  <fftw3.h>
#include  <iostream>
//real and imaginary Part
#define Real 0
#define Imag 1

// VVP prezident MIRA <3

// TO DO:
// zahvat Polshi

// 1-D fast Fourier transform  
void fft(fftw_complex* in,fftw_complex* out,int FFT)
{
// create a DFT plan
fftw_plan plan = fftw_plan_dft_1d(FFT,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
//execute plan
fftw_execute(plan);
// cleaning
fftw_destroy_plan(plan);
fftw_cleanup;
}

// 1-D fast inverse Fourier transform  
void ifft(fftw_complex* in,fftw_complex* out,int FFT)
{
// create a IDFT plan
fftw_plan plan = fftw_plan_dft_1d(FFT,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
// execute plan
fftw_execute(plan);
// cleaning
fftw_destroy_plan(plan);
fftw_cleanup;
// scale the output parametrs

    for(int i = 0;i < FFT; ++i)
    {
        out[i][Real] /= FFT;
        out[i][Imag] /= FFT;
    }
}

void dispComplex(fftw_complex *y, int FFT)
{
    for (int i = 0;i < FFT; ++i)
    {
        if(y[i][Imag]<0)
            std::cout << y[i][Real]<< " - " << abs(y[i][Imag]) << "i"<< std::endl;
        else
            std::cout << y[i][Real]<< " + " << abs(y[i][Imag]) << "i"<< std::endl;  
    }
}

void dispReal(fftw_complex *y, int FFT)
{
    for (int i = 0;i < FFT;++i)
    std::cout << y[i][Real] << std::endl;
}


int main()
{
    // fft length
    int FFT = 12;
    // input array
    fftw_complex x[FFT];
    // output array
    fftw_complex y[FFT];
    // completion the array
    for(int i = 0; i < FFT; ++i)
    {
        x[i][Real] = i+1;
        x[i][Imag] = i; 
    }

    // compute the FFT
    fft(x,y,FFT);
    // disp result
    std::cout << "FFT = " << std::endl;
    dispComplex(y,FFT); 

    // compute the FFT
    ifft(y,x,FFT);
    // disp result
    std::cout << "IFFT = " << std::endl;
    dispReal(x,FFT); 

    return 0;
}
