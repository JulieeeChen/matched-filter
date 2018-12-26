//============================================================================
// Name        : MF6.cpp
// Author      : Sohrab
// Version     : 1
// Copyright   : Hi!
// Description : Matched Filter in C++, Ansi-style
//============================================================================

#include <iostream>
#include <string>
#include <cmath>
#include <math.h>
#include <ctime>
#include <complex>
#include <vector>
#include <string>
#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include <thrust/complex.h>
// For the CUDA runtime routines (prefixed with "cuda_")
#include <cuda_runtime.h>

// references
#define date_ref "1015"
#define obj_ref "1"
#define run_ref 2

// 0: down sampling,  1: averaging, 2: nothing
#define average 2

//internal distance
#define int_dst 2.615

// TX relative position to TX's starting point
#define Tx_pos_x 0.41
#define Tx_pos_y -0.028
#define Tx_pos_z -0.012

//starting point of using samples
#define N_spl_fr 1000
#define N_lfreq_spl  0
#define N_hfreq_spl  0
#define N_mfreq_spl (N_spl_fr/2)-N_lfreq_spl-N_hfreq_spl
#define N_mfreq_spl_slow 2*((N_spl_fr/2)-N_lfreq_spl-N_hfreq_spl)

// Number of frames for each axis
#define N_x_stg 20 //1667
#define N_z_stg 20

//constants
#define Ts 6e-3
#define Rs 5e5
#define lambda 4.983873e-3

//step size between each two frames considered
#define dlx 0.005 * 0.006
#define dlz lambda/2
#define linxmax -dlx/2-(N_x_stg-1)*dlx
#define linzmax -dlz/2-(N_z_stg-1)*dlz

// environment dimensions
#define xmin -5
#define xmax 5


#define ymin 0
#define ymax 10


#define zmin -1.5
#define zmax 1.5


//resolution
#define res 0.07

//scientific values for some constants
#define sci_fac 1e8

#define c_sci 2.9979
#define fc_sci 609
#define As1_sci 1.5001e4
#define As2_sci 7.5005e3


#define file_size 100020000
#define size1 3000

#define BLOCK_WIDTH 8

#define Beat_I(uu, vv, nn) Beat_I[uu*N_z_stg*size1 + vv*size1 + nn]
#define Beat_R(uu, vv, nn) Beat_R[uu*N_z_stg*size1 + vv*size1 + nn]
#define cell_MF(xx, yy, zz) cell_MF[xx*Ny*Nz + yy*Nz + zz]
#define deviceCellMF(xx, yy, zz) deviceCellMF[xx*Ny*Nz + yy*Nz + zz]
using namespace std;

/****************** FUNCTIONS ******************/

struct indices {
	int kx;
	int ky;
	int kz;
};

indices idxfinder(int n1, int n2, int n3, int k) {

	k = k % (n1 * n2 * n3);

	indices I;
	I.kx = k % n1;
	I.ky = ((int) floor(k / n1)) % n2;
	I.kz = ((int) floor(k / (n1 * n2))) % n3;

	return I;
}

/****************************************************/

/************* KERNEL CALL *************************/

__global__ void matchedFilterKernel(float* Beat_R, float* Beat_I, thrust::complex<float>* cell_MF, int Nx, int Ny, int Nz) {

    #define MF_x_axis(xx) (xx*res + xmin)
    #define MF_y_axis(yy) (yy*res + ymin)
    #define MF_z_axis(zz) (zz*res + zmin)
    #define u_axis(uu) (-dlx/2 - uu*dlx)
    #define v_axis(vv) (-dlz/2 - vv*dlz)

    const float pi = acosf(-1);
	const thrust::complex<double> i(0, 1);
	const thrust::complex<float> i_float(0, 1);

    int xx, yy, zz;
    xx = blockIdx.x * blockDim.x + threadIdx.x;
    yy = blockIdx.y * blockDim.y + threadIdx.y;
    zz = blockIdx.z * blockDim.z + threadIdx.z;

    if(xx < Nx && yy < Ny && zz < Nz) {
        float cell_z = MF_z_axis(zz);
        float cell_y = MF_y_axis(yy);
        float cell_x = MF_x_axis(xx);

        thrust::complex<float> cell_sum = 0;

        // for(int nn = 0; nn < size1; nn++)  // 3000
        //             Beat[nn] = Beat_R(uu, vv, nn) + i_float * Beat_I(uu, vv, nn);
        // __shared__ complex<float> Beat[size1]


        float cell_dist_t = sqrtf(
                (cell_x - Tx_pos_x) * (cell_x - Tx_pos_x)
                        + (cell_y - Tx_pos_y) * (cell_y - Tx_pos_y)
                        + (cell_z - Tx_pos_z) * (cell_z - Tx_pos_z));


        for (int uu = 0; uu < N_x_stg; uu++) { // N_x_stg  

            float x_diff = (cell_x - u_axis(uu)) * (cell_x - u_axis(uu));

            for (int vv = 0; vv < N_z_stg; vv++) { // 2d receiver 1667*20

                float temp_tau = (cell_dist_t + int_dst * 2 + sqrtf( x_diff +
                        (cell_z - v_axis(vv)) * (cell_z - v_axis(vv)) + cell_y * cell_y) ) / c_sci;

                thrust::complex<float> temp_sig = exp(-i_float * (float) fmod((float)2.0 * pi * fc_sci * temp_tau, 2*pi) );
                thrust::complex<float> Beat[size1];


                thrust::complex<float> cell_sig_fst_temp[N_mfreq_spl];
                thrust::complex<float> cell_sig_slow_temp[N_mfreq_spl_slow];

                for(int nn = 0; nn < size1; nn++)  // 3000
                    Beat[nn] = Beat_R(uu, vv, nn) + i_float * Beat_I(uu, vv, nn);



                for (int nn = 0; nn < N_mfreq_spl; nn++) { // for each fixed receiver and object location, 3000 samples
                    cell_sig_fst_temp[nn] = temp_sig * exp(-i_float * (float) fmod((float)(2.0 * pi *
                            As1_sci * (N_lfreq_spl / Rs + nn / Rs) * temp_tau), 2*pi));
                    cell_sum += cell_sig_fst_temp[nn] * (Beat_R(uu, vv, nn) + i_float * Beat_I(uu, vv, nn)); //Beat[nn];
                }


                for (int nn = N_mfreq_spl; nn < 2*N_mfreq_spl; nn++) {
                    cell_sum += cell_sig_fst_temp[2*N_mfreq_spl-1-nn] * Beat[nn];
                }



                for (int nn = 0; nn < N_mfreq_spl_slow; nn++) {
                    cell_sig_slow_temp[nn] = temp_sig * exp(-i_float * (float) fmod((float)(2.0 * pi
                            * As2_sci * (N_lfreq_spl * 2 / Rs + nn / Rs) * temp_tau), 2*pi) );
                    cell_sum += cell_sig_slow_temp[nn] * Beat[nn+2*N_mfreq_spl];
                }

                for (int nn = N_mfreq_spl_slow; nn < 2* N_mfreq_spl_slow; nn++) {
                    cell_sum += cell_sig_slow_temp[2*N_mfreq_spl_slow-1-nn] * Beat[nn+2*N_mfreq_spl];
                }


            }



        }

        cell_MF(xx, yy, zz) = cell_sum;
    }


    #undef MF_x_axis
    #undef MF_y_axis
    #undef MF_z_axis
    #undef u_axis
    #undef v_axis

}


/**************************************************/




int
main(void)
{
    cudaError_t err = cudaSuccess;
    /************* LARGE ARRAY DECLRATATIONS AND NX, NY, NZ************/
    int Nx = 143; // (int) floor((xmax-xmin)/res)+1; //143
    int Ny = 10; //(int) floor((ymax-ymin)/res)+1; //143
    int Nz = 10; //(int) floor((zmax-zmin)/res)+1; //43

    // complex<float> cell_sig_fst[N_x_stg][N_z_stg][N_mfreq_spl];
    // complex<float> cell_sig_slow[N_x_stg][N_z_stg][N_mfreq_spl_slow];

    // Allocate host memory

    float* Beat_R = (float *)malloc(N_x_stg * N_z_stg * size1 * sizeof(float)); //[N_x_stg][N_z_stg][size1] = {};
    float* Beat_I = (float *)malloc(N_x_stg * N_z_stg * size1 * sizeof(float)); //[N_x_stg][N_z_stg][size1] = {};
    thrust::complex<float>* cell_MF = (thrust::complex<float>*)malloc(Nx * Ny * Nz * sizeof(thrust::complex<float>)); //[Nx][Ny][Nz] 143 * 143 *43
    
    // Verify that allocations succeeded
    if (Beat_R == NULL || Beat_I == NULL || cell_MF == NULL )
    {
        fprintf(stderr, "Failed to allocate host vectors!\n");
        exit(EXIT_FAILURE);
    }
    /**************************************************/

	clock_t begin = clock();
    clock_t end;
	// srand (time(NULL));
	// for (int ii = 0; ii < N_z_stg; ii++)
	// 	for (int jj = 0; jj < N_x_stg; jj++)
	// 		for (int kk = 0; kk < size1; kk++) {
	// 			Beat_R(jj, ii, kk) = (rand()%10 + 1)/10;
	// 			Beat_I(jj, ii, kk) = (rand()%10 + 1)/10;
    // 		}
    

	/*********** READ THE .BIN FILES ************/
	FILE *fp = fopen("/home/synrg-gpu1/Desktop/MF6/testReal.bin","rb");


	for (int ii = 0; ii < N_z_stg; ii++){
		for (int jj = 0; jj < N_x_stg; jj++) {

			float b[size1];
			fseek(fp, (ii*N_x_stg + jj)*size1*4, SEEK_SET);
			fread(b, sizeof *b, size1, fp);
			for(int kk = 0; kk < size1; kk++) {
				Beat_R(jj, ii, kk) = 1;
				//if (ii == 0 && jj == 1 && kk < 500) cout << b[kk] << endl;
			}
		}
	}

	fclose(fp);

	cout << "Successfully read the file in " << (double) (clock() - begin) / CLOCKS_PER_SEC << " seconds!" << endl;

	FILE *fp2 = fopen("/home/synrg-gpu1/Desktop/MF6/testImag.bin","rb");


	for (int ii = 0; ii < N_z_stg; ii++){
		for (int jj = 0; jj < N_x_stg; jj++) {

			float b[size1];
			fseek(fp2, (ii*N_x_stg + jj)*size1*4, SEEK_SET);
			fread(b, sizeof *b, size1, fp2);

			for(int kk = 0; kk < size1; kk++) {
				Beat_I(jj, ii, kk) = 0;

			}
		}
	}

	fclose(fp2);
	cout << "Successfully read the files in " << (double) (clock() - begin) / CLOCKS_PER_SEC << " seconds!" << endl;
	// cout << Beat_I(149, 14, 149)<< endl << endl;

	/******************** END OF READ FILE *********************/

	//some constants
	const float pi = acos(-1);
	const thrust::complex<double> i(0, 1);
    const thrust::complex<float> i_float(0, 1);


    for (int i = 0; i < 1; i++){
        for (int j = 0; j < 1; j++){
            for (int k = 0; k < 10; k++) {
                cout << cell_MF(k, j, i) << " ";
            }
            std::endl( std::cout );
        }
       std::endl( std::cout );
  }

    float* deviceBeatI;
    float* deviceBeatR;
    thrust::complex<float>* deviceCellMF;

    clock_t begin_mem = clock();
    
    // Allocate GPU memory
    err = cudaMalloc((void **) &deviceBeatR , N_z_stg * N_x_stg * size1 * sizeof(float));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate deviceBeatR (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **) &deviceBeatI , N_z_stg * N_x_stg * size1 * sizeof(float));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate deviceBeatI (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    err = cudaMalloc((void **) &deviceCellMF , Nx * Ny * Nz * sizeof(thrust::complex<float>));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate deviceCellMF (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }


    printf("Copy input data from the host memory to the CUDA device\n");

    err = cudaMemcpy(deviceBeatR, Beat_R, N_z_stg * N_x_stg * size1 * sizeof(float), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy deviceBeatR from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(deviceBeatI, Beat_I, N_z_stg * N_x_stg * size1 * sizeof(float), cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy deviceBeatI from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    cout << "hi!" << endl;

    
    dim3 DimGrid(ceil(Nx * 1.0 / BLOCK_WIDTH), ceil(Ny * 1.0 / BLOCK_WIDTH), ceil(Nz * 1.0 /BLOCK_WIDTH));
    dim3 DimBlock(BLOCK_WIDTH, BLOCK_WIDTH, BLOCK_WIDTH);

    cout << "Allocating & copying memory DONE! Time taken:" << (double) (clock() - begin_mem) / CLOCKS_PER_SEC;
    
    matchedFilterKernel<<<DimGrid, DimBlock>>>(deviceBeatR, deviceBeatI, deviceCellMF, Nx, Ny, Nz);
    cudaDeviceSynchronize();
    err = cudaGetLastError();

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to launch matchedFilterKernel  (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    begin_mem = clock();

    printf("Copy output data from the CUDA device to the host memory\n");

    err = cudaMemcpy(cell_MF, deviceCellMF, Nx * Ny * Nz * sizeof(thrust::complex<float>), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy deviceCellMF from device to host (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }


    cout << "Copying memory back DONE! Time taken:" << (double) (clock() - begin_mem) / CLOCKS_PER_SEC;


	cout << "Hi! \n";
	cout << Nx << endl;
	cout << Ny << endl;
	cout << Nz << endl;
	cout << N_x_stg << endl;
	cout << N_z_stg << endl;
	cout << N_mfreq_spl_slow << endl;

	end = clock();
    cout << "DONE! Time taken:" << (double) (end - begin) / CLOCKS_PER_SEC;

    for (int i = 0; i < 1; i++){
        for (int j = 0; j < 5; j++){
            for (int k = 0; k < 5; k++) {
                cout << cell_MF(k, j, i) << " ";
            }
            std::endl( std::cout );
        }
        std::endl( std::cout );
    }

    cudaFree(deviceBeatR);
    cudaFree(deviceBeatI);
    cudaFree(deviceCellMF);

    free(Beat_R);
    free(Beat_I);
    free(cell_MF);

	return 0;
}
