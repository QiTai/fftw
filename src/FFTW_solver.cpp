/*********************************
This is the CPP file of FFTW solver

Author : Chao Yan
Date   : May 21st, 2016
**********************************/

#include "FFTW_solver.h"

FFTW_solver::FFTW_solver() {

    /*=====================================
    Initializing arrays in real space
    =====================================*/
    v = new double*[numOfXGrid];

    temp_Velocity = new double[numOfXGrid * numOfYGrid];

    for (int i = 0; i < numOfXGrid; i++) {
        v[i] = new double[numOfYGrid]; 
    }

//    v[0][0] = 1;
//    v[0][1] = 2;
//    v[1][0] = 3;
//    v[1][1] = 4;
//    v[2][0] = 5;
//    v[2][1] = 6;
//    v[3][0] = 7;
//    v[3][1] = 8;
//
//    for (int i = 0; i < numOfXGrid * numOfYGrid; i++) {
//        temp_Velocity[i] = 0;
//    }
//    
//    /*=======================================
//        initializing multiple threads
//      =====================================*/
//    //if (fftw_init_threads()) {
//    //    fftw_plan_with_nthreads(THREADS);
//    //    cout << "Using " << THREADS << "threads" << endl;
//    //} else {
//    //    cout << "Using multiple threads failed" << endl;
//    //    exit(0);
//    //}
//
//    /*=====================================
//      initializing arrays in fourier space
//      ===================================*/
//
//    V = (fftw_complex**)fftw_malloc(sizeof(fftw_complex*) * numOfXGrid);
    temp_U = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * numOfXGrid * (numOfYGrid / 2 + 1));

    for (int i = 0; i < numOfXGrid; i++) {
        V[i] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (numOfYGrid / 2 + 1));
        for (int j = 0; j < numOfYGrid; j++) {
            V[i][j][REAL] = 0;
            V[i][j][IMAG] = 0;
        }
    }

    for (int i = 0; i < numOfXGrid * (numOfYGrid / 2 + 1); i++) {
        temp_U[i][REAL] = 0;
        temp_U[i][IMAG] = 0;    
    }

    /*=========================================
      Initializing plans
      =======================================*/
    plan_r2c = fftw_plan_dft_r2c_2d(numOfXGrid, numOfYGrid, temp_Velocity, temp_U, FFTW_ESTIMATE);    
    plan_c2r = fftw_plan_dft_c2r_2d(numOfXGrid, numOfYGrid, temp_U, temp_Velocity, FFTW_ESTIMATE);

    // test for complex multiplication
   // fftw_complex *temp_R = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * numOfXGrid * (numOfYGrid / 2 + 1));

   // for (int i = 0; i < numOfXGrid * (numOfYGrid / 2 + 1); i++) {
   //     temp_R[i][REAL] = 0;
   //     temp_R[i][IMAG] = 0;
   // }

   // temp_U[0][REAL] = 1;
   // temp_U[0][IMAG] = 2;
   // temp_U[1][REAL] = 3;
   // temp_U[1][IMAG] = 4;
   // temp_U[2][REAL] = 5;
   // temp_U[2][IMAG] = 6;
   // temp_U[3][REAL] = 7;
   // temp_U[3][IMAG] = 8;
   // temp_U[4][REAL] = 9;
   // temp_U[4][IMAG] = 10;
   // temp_U[5][REAL] = 11;
   // temp_U[5][IMAG] = 12;
   // temp_U[6][REAL] = 13;
   // temp_U[6][IMAG] = 14;
   // temp_U[7][REAL] = 15;
   // temp_U[7][IMAG] = 16;

   // temp_R = temp_U * temp_U;

   // for (int i = 0; i < numOfXGrid * (numOfYGrid / 2 + 1); i++) {
   //     cout << "temp_R[" << i << "] : " << temp_R[i][REAL] << " + " << temp_R[i][IMAG] << "i" << endl;
   // }

}

FFTW_solver::~FFTW_solver() {
}

void FFTW_solver::solve() {
    

    /*==================================================
      set the double* in fftw_plan_dft_r2c_2d and execute
      ================================================*/
    for (int i = 0; i < numOfXGrid; i++) {
        for (int j = 0; j < numOfYGrid; j++) {
            temp_Velocity[i * numOfYGrid + j] = v[i][j];
        }
    }
    
    fftw_execute(plan_r2c);

    /*===================================================
      set the complex* in fftw_plan_dft_c2r_2d and execute
      =================================================*/
    for (int i = 0; i < numOfXGrid; i++) {
        for (int j = 0; j < numOfYGrid; j++) {
            V[i][j][REAL] = temp_U[i * (numOfYGrid / 2 + 1) + j][REAL];
            V[i][j][IMAG] = temp_U[i * (numOfYGrid / 2 + 1) + j][IMAG];
            cout << "V[" << i << "][" << j << "] : " << V[i][j][REAL] << " + " << V[i][j][IMAG] << "i" << endl;
        }
    }

    fftw_execute(plan_c2r);


    for (int i = 0; i < numOfXGrid; i++) {
        for (int j = 0; j < numOfYGrid; j++) {
            v[i][j] = temp_Velocity[i * numOfYGrid + j] / (numOfXGrid * numOfYGrid);
            cout << "v[" << i << "][" << j << "] : " << v[i][j] << endl;
        }    
    }

    /*======================================
      destory the fftw paln
      ====================================*/
    fftw_destroy_plan(plan_r2c);
    fftw_destroy_plan(plan_c2r);

    /*=====================================
        optional : deallocate all of other persistent data, and reset FFTW to the pristine state
      ===================================*/
    fftw_cleanup();

}
