/**********************************
This is the FFTW solver. FFTW is a package 
developed by MIT to solver problems by FFT.

Author : Chao Yan
Date   : May 21st, 2016
**********************************/

#include <fftw3.h>
#include <iostream>

#define numOfXGrid 4
#define numOfYGrid 2

#define REAL 0
#define IMAG 1

using namespace std;

class FFTW_solver {
public:
    FFTW_solver();
    ~FFTW_solver();
    void solve();

private:

    //int numOfXGrid;
    //int numOfYGrid;

    // These are the arrays in real space;
    double **v;
    // temp is used in executing the plans;
    double *temp_Velocity;
    
    // These are the arrays in fourier space
    fftw_complex **V;
    // temp is used in executing the plans;
    fftw_complex *temp_U;
    
    // FFT transform 
    fftw_plan plan_r2c;
    fftw_plan plan_c2r;

};
