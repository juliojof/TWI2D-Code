#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "PropagationFunctions_GPU_Wrap.h"
#include "PropagationFunctions_GPU.h"
#include <cuda_runtime.h>

dim3 GPU_nBlocks;
dim3 GPU_nThreads;

float *GPU_wav1;
float *GPU_wav2;
float *GPU_lap;
float *GPU_sigma;
float *GPU_sigmaInv;
float *GPU_ws;
float *GPU_source;
float *GPU_seismogram;

int GPU_nx;
int GPU_ny;
int GPU_nz;
int GPU_nt;

int GPU_dx;
int GPU_dy;
int GPU_dz;
int GPU_dt;

int GPU_nrec;

void GPU_initEnvironment(int nx, int nz, int nt, int nrec) {
    GPU_nThreads.x = 32;
    GPU_nThreads.y = 1;
    GPU_nThreads.z = 16;
    GPU_nBlocks.x  = nx/GPU_nThreads.x;
    GPU_nBlocks.y  = 1;
    GPU_nBlocks.z  = nz/GPU_nThreads.z;

    GPU_nx = nx;
    GPU_nz = nz;
    GPU_nt = nt;
    GPU_nrec = nrec;
}

void GPU_alocArrays_acoustic() {
    size_t sz = GPU_nx * GPU_nz * sizeof(float);
    // GPU_alocMem((void**) &GPU_wav1, sz);
    GPU_alocMem((void**) &GPU_wav2, sz);
    GPU_alocMem((void**) &GPU_ws, sz);
    GPU_alocMem((void**) &GPU_lap, sz);
    GPU_alocMem((void**) &GPU_sigma, sz);
    GPU_alocMem((void**) &GPU_sigmaInv, sz);
    GPU_alocMem((void**) &GPU_source, GPU_nt*sizeof(float));
    GPU_alocMem((void**) &GPU_seismogram, GPU_nrec*GPU_nt*sizeof(float));
}

void GPU_freeArrays_acoustic() {
    cudaFree(GPU_wav1);
    cudaFree(GPU_wav2);
    cudaFree(GPU_ws);
    cudaFree(GPU_lap);
    cudaFree(GPU_sigma);
    cudaFree(GPU_sigmaInv);
    cudaFree(GPU_source);
    cudaFree(GPU_seismogram);
}

void GPU_getSourceWavelet_CPU2GPU(float *source) {
    cudaMemcpy(GPU_source, source, GPU_nt*sizeof(float), cudaMemcpyHostToDevice);
}
void GPU_getRecordedData_GPU2CPU(float *seismogram) {
    cudaMemcpy(seismogram, GPU_seismogram, GPU_nrec*GPU_nt*sizeof(float), cudaMemcpyDeviceToHost);
}

void GPU_copyPresentWavefield_CPU2GPU(float *wavefield) {
    cudaMemcpy(GPU_wav2, wavefield, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyHostToDevice);
}
void GPU_copyLaplacian_GPU2CPU(float *wavefield) {
    cudaMemcpy(GPU_wav2, wavefield, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyHostToDevice);
}
void GPU_laplacian_10thOrder_2D_Acoustic_Iso() {
    laplacian_10thOrder_2D_Acoustic_Iso_GPU <<< GPU_nBlocks, GPU_nThreads >>> (GPU_lap, GPU_wav2, GPU_nx, GPU_nz);
}