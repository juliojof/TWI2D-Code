#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>


#define IC_REGULAR_IMAGE        0
#define IC_REGULAR_ILUMIN       1
#define IC_BORNSCT_IMAGE        2
#define IC_BORNSCT_IMAGE2       3
#define IC_BORNSCT_ILUMIN       4

#define IC_SOU_WAVE             100
#define IC_REC_WAVE             101


// 2nd derivative not staggered
#define C5_10    (+0.0003174603f)
#define C4_10    (-0.004960318f)
#define C3_10    (+0.03968254f)
#define C2_10    (-0.2380952f)
#define C1_10    (+1.666667f)
#define C0_10    (+1.463611f)

// 1st derivative staggered - 10th order
#define B5_10    (+0.000118679f)
#define B4_10    (-0.00176566f)
#define B3_10    (+0.0138428f)
#define B2_10    (-0.0897217f)
#define B1_10    (+1.21124f)

// 1st derivative staggered - 20th order
#define B10_20    (-0.000000037237585f)
#define B09_20    (+0.000000883780612f)
#define B08_20    (-0.0000102165f)
#define B07_20    (+0.0000770772f)
#define B06_20    (-0.000430613f)
#define B05_20    (+0.00192978f)
#define B04_20    (-0.00744345f)
#define B03_20    (+0.0270942f)
#define B02_20    (-0.112892f)
#define B01_20    (+1.24181f)

// 1st derivative staggered - 40th order
#define B05_40    (+0.005662415831386f)
#define B04_40    (-0.014040751058073f)
#define B03_40    (+0.037232819311337f)
#define B02_40    (-0.126407610083589f)
#define B01_40    (+1.257424773937471f)

// 1st derivative staggered - 80th order
#define B05_80    (+0.0569362f)
#define B04_80    (-0.438909f)
#define B03_80    (+1.8755f)
#define B02_80    (-3.45354f)
#define B01_80    (+1.26531f)

// 1st derivative not staggered
#define A5_10_f  (+1/630)
#define A4_10_f  (-5/252)
#define A3_10_f  (+5/42)
#define A2_10_f  (-10/21)
#define A1_10_f  (+5/3)
#define A5_10    (+0.00158730)
#define A4_10    (-0.01984127)
#define A3_10    (+0.11904762)
#define A2_10    (-0.47619048)
#define A1_10    (+1.66666667)

#define OL 5
#define id2(ix,iz)                         ( (ix) + (iz)*nx )
#define id2t(ix,iz,it)                     ( (ix) + (iz)*nx + (it)*nx*nz )
#define id3_eic_lx(lx,ix,iz)               ( (lx)+lx0 + (ix)*nlx + (iz)*nx*nlx)
#define id3_eic_lx_lz(lx,lz,ix,iz)         ( (lx)+lx0 + ((lz)+lz0)*nlx + (ix)*nlx*nlz + (iz)*nx*nlx*nlz)
#define id3_eic_lx_transp(lx,iz,ix)        ( (lx)+lx0 + (iz)*nlx + (ix)*nz*nlx)
#define id3_eic_th(ix,iz,itheta)           ( (ix) + (iz)*nx + ((itheta)+ltheta0)*nx*nz )
#define id3_eic_th_transp(itheta,iz,ix)    ( (itheta)+ltheta0 + (iz)*ntheta + (ix)*nz*ntheta )

#define jd2(jx,jz)                         ( (jx) + (jz)*(bdx+2*OL) )

#define tx threadIdx.x
#define ty threadIdx.y
#define tz threadIdx.z
#define bx blockIdx.x
#define by blockIdx.y
#define bz blockIdx.z

#define bdx blockDim.x
#define bdy blockDim.y
#define bdz blockDim.z
#define gdx gridDim.x
#define gdy gridDim.y
#define gdz gridDim.z

// Imaging condition
__global__ void GPU_DEV_imagingCondition_transmission(float *vel, float *image, float *ws, float *wr, int nx, int nz, int nt, \
                                                      float dx, float dz, float dt, int lx0, int lz0, int nxb, int nzb, int it);
__global__ void GPU_DEV_imagingCondition(float *image, const float *u1, const float *u2, int nx, int nz);
__global__ void GPU_DEV_eimagingCondition_lx(float *eimage, float *ws, float *wr, int nx, int nz, int lx0, int nxb, int nzb);
__global__ void GPU_DEV_eimagingCondition_lx_lz(float *eimage, float *ws, float *wr, int nx, int nz, int lx0, int lz0, int nxb, int nzb);

// Data recording
__global__ void GPU_DEV_recordData(int mode, float *seismogramMultiplexed, int *recIndexes, float *wav, int nx, int nz, int nrec);

// Source injection
__global__ void GPU_DEV_injectSource_Acoustic(float *u, float *source, int it, int ix0, int iz0, int nx, int nz, float dx, float dz);

// Compute explicit time-step
void GPU_LNCH_timeStep_Acoustic_abs2(float *GPU_wav2, float *GPU_wav1, float *GPU_lap);
__global__ void GPU_DEV_timeStep_Acoustic_abs2(float *u1, const float *u2, const float *lap, const float *vp2dt2dx2, \
                                               const float *sigma, const float *sigmaInv, int nx, int nz);

// Compute the 2D laplacian to 10th order
void GPU_LNCH_laplacian_10thOrder_2D_Acoustic_Iso(float *GPU_wav2, float *GPU_lap);
__global__ void GPU_DEV_laplacian_10thOrder_2D_Acoustic_Iso(float *lap, const float *u2, int nx, int nz);
__global__ void GPU_DEV_laplacian_10thOrder_2D_Acoustic_Iso_opt0(float *lap, const float *u2, int nx, int nz);

// Born Scattering Kernel Functions
__global__ void GPU_OPER_wavefieldTimeSecondDerivative_part1(float *secDer, float *u2, float *u1, float invDt2, int nx);
__global__ void GPU_OPER_wavefieldTimeSecondDerivative_part2(float *secDer, float *u2, float invDt2, int nx);
__global__ void GPU_DEV_scatter_Acoustic_lx_opt(float *eref, float *ws, float *sct, int nx, int nz, int lx0, int nxb, int nzb);

// Demoveout operator
__global__ void GPU_DEV_Contraction_lx(int add, float *eimage, float *eimage_contracted,  \
                                  int nx, float dx, int nz, float dz, int lx0, \
                                       int ixMin, int ixMax, int izMin, int izMax, float contraction_factor);
__global__ void GPU_DEV_Demoveout(int add, int adj, float *eimage, float *eimage_deMO,  \
                                  int nx, float dx, int nz, float dz, int lx0, \
                                  int ixMin, int ixMax, int izMin, int izMax, \
                                  float contraction_factor, float lx_null);
__global__ void GPU_DEV_Demoveout_Hicks(int add, int adj, float *eimage, float *eimage_deMO,  \
                                  int nx, float dx, int nz, float dz, int lx0, \
                                  int ixMin, int ixMax, int izMin, int izMax, \
                                  float contraction_factor, float lx_null);
__global__ void GPU_DEV_Demoveout_Hicks_old(int add, float *eimage, float *eimage_deMO,  \
                                  int nx, float dx, int nz, float dz, int lx0, \
                                  int ixMin, int ixMax, int izMin, int izMax, \
                                  float contraction_factor, float lx_null);

// 
__global__ void GPU_DEV_lx2theta(int adj, int add, float *eimage, float *eimage_theta, \
                                 int nx, float dx, int nz, float dz, int nxb, int nzb, \
                                 int lx0, int ntheta, float dtheta, float otheta);

__device__ float fsignf_GPU(float x);


// Hicks interp
__device__ void  DEV_interSinc1D_8P(float* coef, float alpha);
__device__ void  DEV_interSinc2D_8P(float *coef, float alpha_x, float alpha_y);
__device__ float DEV_sinc(float x);
__device__ float DEV_kaiser(float x, float r, float b);
__device__ float DEV_mBFZ(float x);


dim3 GPU_nBlocks;
dim3 GPU_nThreads;

float *GPU_vp;
float *GPU_wav1;
float *GPU_wav2;
float *GPU_secDer;
float *GPU_secDer_sct;
float *GPU_sct1;
float *GPU_sct2;
float *GPU_ws;
float *GPU_wc;
float *GPU_lap;
float *GPU_sigma;
float *GPU_sigmaInv;

float *GPU_sou_inc;
float *GPU_sou_sct;
float *GPU_rec_inc;
float *GPU_rec_sct;

float *GPU_eref;
float *GPU_eimage;
float *GPU_image;
float *GPU_ilumin;
float *GPU_imageSouSide;
float *GPU_imageRecSide;
float *GPU_iluminSouSide;
float *GPU_iluminRecSide;


float *GPU_source;
float *GPU_seismogram;
float *GPU_seismogramMultiplexed;
int   *GPU_recIndexes;

int GPU_nx;
int GPU_ny;
int GPU_nz;
int GPU_nt;

int GPU_dx;
int GPU_dy;
int GPU_dz;
int GPU_dt;

int    GPU_EIC_jt;
float  GPU_EIC_dt;
int    GPU_EIC_nt;

float GPU_invDt2;

int GPU_nrec;

int GPU_lx0;
int GPU_lz0;
int GPU_nlx;
int GPU_nlz;

int GPU_nxb;
int GPU_nzb;

int GPU_izs, GPU_ixs, GPU_it;

int GPU_idev;


extern "C" void GPU_OPER_initEnvironment(int idev, int nx, int nz, int nt, float dx, float dz, float dt, int nrec, int lx0, int lz0, int nxb, int nzb) {

    GPU_idev = idev;
    cudaSetDevice(GPU_idev);
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, GPU_idev);
    int nDev;
    cudaGetDeviceCount(&nDev);
    // printf("\n\n Cuda device properties: nDev=%d    sharedMemPerBlock=%zd    totalGlobalMem=%zd\n\n", nDev, \
    prop.sharedMemPerBlock, prop.totalGlobalMem);

    GPU_nThreads.x = 32;
    GPU_nThreads.y = 1;
    GPU_nThreads.z = 8;
    GPU_nBlocks.x  = nx/GPU_nThreads.x;
    GPU_nBlocks.y  = 1;
    GPU_nBlocks.z  = nz/GPU_nThreads.z;

    GPU_dx = dx;
    GPU_dz = dz;
    GPU_dt = dt;
    GPU_nx = nx;
    GPU_nz = nz;
    GPU_nt = nt;
    GPU_nrec = nrec;
    GPU_nxb = nxb;
    GPU_nzb = nzb;
    GPU_lx0 = lx0;
    GPU_lz0 = lz0;
    GPU_nlx = 2*lx0+1;
    GPU_nlz = 2*lz0+1;

    printf("\n\n  Sizes for InitEnv:  GPU_dt=%f  GPU_nt=%d  \n\n", GPU_dt, GPU_nt);

    GPU_invDt2 = 1.0f/(dt*dt);
}

extern "C" void GPU_OPER_alocMem(void** dev_a, size_t size) {
    cudaMalloc(dev_a, size);
    cudaMemset(*dev_a, 0, size);
}

extern "C" void GPU_OPER_alocArrays_acoustic() {
    size_t sz = GPU_nx * GPU_nz * sizeof(float);
    cudaSetDevice(GPU_idev);
    GPU_OPER_alocMem((void**) &GPU_vp                     , sz);
    GPU_OPER_alocMem((void**) &GPU_wav1                   , sz);
    GPU_OPER_alocMem((void**) &GPU_wav2                   , sz);
    GPU_OPER_alocMem((void**) &GPU_secDer                 , sz);
    GPU_OPER_alocMem((void**) &GPU_ws                     , sz);
    GPU_OPER_alocMem((void**) &GPU_lap                    , sz);
    GPU_OPER_alocMem((void**) &GPU_sigma                  , sz);
    GPU_OPER_alocMem((void**) &GPU_sigmaInv               , sz);
    GPU_OPER_alocMem((void**) &GPU_image                  , sz);
    GPU_OPER_alocMem((void**) &GPU_ilumin                 , sz);
    GPU_OPER_alocMem((void**) &GPU_eimage                 , sz*GPU_nlx*GPU_nlz);
    GPU_OPER_alocMem((void**) &GPU_source                 , GPU_nt         *sizeof(float));
    GPU_OPER_alocMem((void**) &GPU_seismogram             , GPU_nrec*GPU_nt*sizeof(float));
    GPU_OPER_alocMem((void**) &GPU_seismogramMultiplexed  , GPU_nrec*GPU_nt*sizeof(float));
    GPU_OPER_alocMem((void**) &GPU_recIndexes             , GPU_nrec*2     *sizeof(int  ));
}
extern "C" void GPU_OPER_freeArrays_acoustic() {
    cudaSetDevice(GPU_idev);
    cudaFree(GPU_vp);
    cudaFree(GPU_wav1);
    cudaFree(GPU_wav2);
    cudaFree(GPU_secDer);
    cudaFree(GPU_ws);
    cudaFree(GPU_lap);
    cudaFree(GPU_sigma);
    cudaFree(GPU_sigmaInv);
    cudaFree(GPU_image);
    cudaFree(GPU_ilumin);
    cudaFree(GPU_eimage);
    cudaFree(GPU_source);
    cudaFree(GPU_seismogram);
    cudaFree(GPU_seismogramMultiplexed);
    cudaFree(GPU_recIndexes);
}

extern "C" void GPU_OPER_alocArrays_acoustic_ExtBornSct() {
    size_t sz = GPU_nx * GPU_nz * sizeof(float);
    cudaSetDevice(GPU_idev);
    GPU_OPER_alocMem((void**) &GPU_secDer_sct             , sz);
    GPU_OPER_alocMem((void**) &GPU_sct1                   , sz);
    GPU_OPER_alocMem((void**) &GPU_sct2                   , sz);
    GPU_OPER_alocMem((void**) &GPU_wc                     , sz);
    GPU_OPER_alocMem((void**) &GPU_eref                   , sz*GPU_nlx*GPU_nlz);

    GPU_OPER_alocMem((void**) &GPU_imageSouSide           , sz);
    GPU_OPER_alocMem((void**) &GPU_imageRecSide           , sz);
    GPU_OPER_alocMem((void**) &GPU_iluminSouSide          , sz);
    GPU_OPER_alocMem((void**) &GPU_iluminRecSide          , sz);
}
extern "C" void GPU_OPER_freeArrays_acoustic_ExtBornSct() {
    cudaSetDevice(GPU_idev);
    cudaFree(GPU_secDer_sct);
    cudaFree(GPU_sct1);
    cudaFree(GPU_sct2);
    cudaFree(GPU_wc);
    cudaFree(GPU_eref);

    cudaFree(GPU_imageSouSide);
    cudaFree(GPU_imageRecSide);
    cudaFree(GPU_iluminSouSide);
    cudaFree(GPU_iluminRecSide);
}


extern "C" void GPU_OPER_alocArrays_transmissionEIC(int jt_EIC) {
    GPU_EIC_jt = jt_EIC;
    if(jt_EIC<=0)    return;
    GPU_EIC_dt =  GPU_dt * GPU_EIC_jt;
    GPU_EIC_nt = ((GPU_dt*(GPU_nt-1))/GPU_EIC_dt) + 1;
    size_t sz = GPU_nx * GPU_nz * GPU_EIC_nt * sizeof(float);
    printf("\n\n  Sizes for IC:  GPU_EIC_jt=%d  GPU_EIC_dt=%f  GPU_EIC_nt=%d   GPU_dt=%f  GPU_nt=%d  \n\n", GPU_EIC_jt, GPU_EIC_dt, GPU_EIC_nt, GPU_dt, GPU_nt);
    cudaSetDevice(GPU_idev);
    GPU_OPER_alocMem((void**) &GPU_sou_inc   , sz);
    GPU_OPER_alocMem((void**) &GPU_sou_sct   , sz);
    GPU_OPER_alocMem((void**) &GPU_rec_inc   , sz);
    GPU_OPER_alocMem((void**) &GPU_rec_sct   , sz);
}
extern "C" void GPU_OPER_freeArrays_transmissionEIC() {
    if(GPU_EIC_jt<=0)    return;
    cudaSetDevice(GPU_idev);
    cudaFree(GPU_sou_inc);
    cudaFree(GPU_sou_sct);
    cudaFree(GPU_rec_inc);
    cudaFree(GPU_rec_sct);
}

extern "C" void GPU_OPER_getVelP_CPU2GPU(float *vp) {
    int ix, iz;
    float *transp = (float *) calloc(GPU_nx*GPU_nz, sizeof(float));
    for(ix=0; ix<GPU_nx; ix++)
        for(iz=0; iz<GPU_nz; iz++)    transp[ix+iz*GPU_nx] = vp[iz+ix*GPU_nz];

    cudaSetDevice(GPU_idev);
    cudaMemcpy(GPU_vp, transp, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyHostToDevice);

    free(transp);
}
extern "C" void GPU_OPER_getSigmas_CPU2GPU(float *sigma, float *sigmaInv) {
    int ix, iz;
    float *transp = (float *) calloc(GPU_nx*GPU_nz, sizeof(float));
    
    for(ix=0; ix<GPU_nx; ix++)
        for(iz=0; iz<GPU_nz; iz++)    transp[ix+iz*GPU_nx] = sigma[iz+ix*GPU_nz];
    cudaSetDevice(GPU_idev);
    cudaMemcpy(GPU_sigma   , transp, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyHostToDevice);

    for(ix=0; ix<GPU_nx; ix++)
        for(iz=0; iz<GPU_nz; iz++)    transp[ix+iz*GPU_nx] = sigmaInv[iz+ix*GPU_nz];
    cudaSetDevice(GPU_idev);
    cudaMemcpy(GPU_sigmaInv, transp, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyHostToDevice);

    free(transp);
}
extern "C" void GPU_OPER_getSourceWavelet_CPU2GPU(float *source) {
    cudaSetDevice(GPU_idev);
    cudaMemcpy(GPU_source, source, GPU_nt*sizeof(float), cudaMemcpyHostToDevice);
}
extern "C" void GPU_OPER_getRecordedData_GPU2CPU(float *seismogram) {
    float *seisTemp = (float *) calloc(GPU_nt*GPU_nrec, sizeof(float));

    cudaSetDevice(GPU_idev);
    cudaMemcpy(seisTemp, GPU_seismogramMultiplexed, GPU_nrec*GPU_nt*sizeof(float), cudaMemcpyDeviceToHost);
    
    int ir, it;
    for(ir=0; ir<GPU_nrec; ir++)
        for(it=0; it<GPU_nt; it++)    seismogram[it+ir*GPU_nt] = seisTemp[ir+it*GPU_nrec];

    free(seisTemp);
}
extern "C" void GPU_OPER_getRecordedData_CPU2GPU(float *seismogram) {
    float *seisTemp = (float *) calloc(GPU_nt*GPU_nrec, sizeof(float));
    int ir, it;
    for(ir=0; ir<GPU_nrec; ir++)
        for(it=0; it<GPU_nt; it++)    seisTemp[ir+it*GPU_nrec] = seismogram[it+ir*GPU_nt];

    cudaSetDevice(GPU_idev);
    cudaMemcpy(GPU_seismogramMultiplexed, seisTemp, GPU_nrec*GPU_nt*sizeof(float), cudaMemcpyHostToDevice);
    free(seisTemp);
}

extern "C" void GPU_OPER_getSourceReceiverIndexes_CPU2GPU(int ixs, int izs, int *ixr, int *izr, int nrec) {
 
    GPU_ixs = ixs;
    GPU_izs = izs;

    int ir;
    int *recIndexes = (int *) calloc(2*nrec, sizeof(int));
    for(ir=0; ir<nrec; ir++) {
       recIndexes[2*ir + 0] = izr[ir];
       recIndexes[2*ir + 1] = ixr[ir];
    }
    cudaSetDevice(GPU_idev);
    cudaMemcpy(GPU_recIndexes, recIndexes, GPU_nrec*2*sizeof(int), cudaMemcpyHostToDevice);
    free(recIndexes);
}
extern "C" void GPU_OPER_zeroWavefieldArrays() {
    cudaSetDevice(GPU_idev);
    cudaMemset(GPU_wav2, 0, GPU_nx*GPU_nz*sizeof(float));
    cudaMemset(GPU_wav1, 0, GPU_nx*GPU_nz*sizeof(float));
    cudaMemset(GPU_lap,  0, GPU_nx*GPU_nz*sizeof(float));
}
extern "C" void GPU_OPER_zeroScatteredWavefieldArrays() {
    cudaSetDevice(GPU_idev);
    cudaMemset(GPU_sct2, 0, GPU_nx*GPU_nz*sizeof(float));
    cudaMemset(GPU_sct1, 0, GPU_nx*GPU_nz*sizeof(float));
}
extern "C" void GPU_OPER_copyEImage_GPU2CPU(float *eimage) {
    int lx, lz, ix, iz;
    float *transp = (float *) calloc(GPU_nlx*GPU_nlz*GPU_nlz*GPU_nx*GPU_nz, sizeof(float));
    cudaSetDevice(GPU_idev);
    cudaMemcpy(transp, GPU_eimage, GPU_nlx*GPU_nlz*GPU_nx*GPU_nz*sizeof(float), cudaMemcpyDeviceToHost);
    for(ix=0; ix<GPU_nx; ix++)
        for(iz=0; iz<GPU_nz; iz++)
            for(lx=0; lx<GPU_nlx; lx++)    
                for(lz=0; lz<GPU_nlz; lz++) { 
                    int i  = lx + lz*GPU_nlx + iz*GPU_nlx*GPU_nlz + ix*GPU_nlx*GPU_nlz*GPU_nz;
                    int i_ = lx + lz*GPU_nlx + ix*GPU_nlx*GPU_nlz + iz*GPU_nlx*GPU_nlz*GPU_nx;
                    eimage[i] = transp[i_];
                }
    free(transp);
}
extern "C" void GPU_OPER_copyEImage_CPU2GPU(float *eimage) {
    int lx, lz, ix, iz;
    float *transp = (float *) calloc(GPU_nlx*GPU_nlz*GPU_nx*GPU_nz, sizeof(float));
    for(ix=0; ix<GPU_nx; ix++)
        for(iz=0; iz<GPU_nz; iz++)
            for(lx=0; lx<GPU_nlx; lx++)    
                for(lz=0; lz<GPU_nlz; lz++) {
                    int i  = lx + lz*GPU_nlx + iz*GPU_nlx*GPU_nlz + ix*GPU_nlx*GPU_nlz*GPU_nz;
                    int i_ = lx + lz*GPU_nlx + ix*GPU_nlx*GPU_nlz + iz*GPU_nlx*GPU_nlz*GPU_nx;
                    transp[i_] = eimage[i];
                }
    cudaSetDevice(GPU_idev);
    cudaMemcpy(GPU_eimage, transp, GPU_nlx*GPU_nlz*GPU_nx*GPU_nz*sizeof(float), cudaMemcpyHostToDevice);
    free(transp);
}
extern "C" void GPU_OPER_copyImage_GPU2CPU(float *image) {
    int ix, iz;
    float *transp = (float *) calloc(GPU_nx*GPU_nz, sizeof(float));
    cudaSetDevice(GPU_idev);
    cudaMemcpy(transp, GPU_image, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyHostToDevice);
    for(ix=0; ix<GPU_nx; ix++)
        for(iz=0; iz<GPU_nz; iz++)
            image[iz+ix*GPU_nz] = transp[ix+iz*GPU_nx];
    free(transp);
}
extern "C" void GPU_OPER_copyImage_CPU2GPU(float *image) {
    int ix, iz;
    float *transp = (float *) calloc(GPU_nx*GPU_nz, sizeof(float));
    for(ix=0; ix<GPU_nx; ix++)
        for(iz=0; iz<GPU_nz; iz++)
            transp[ix+iz*GPU_nx] = image[iz+ix*GPU_nz];
    cudaSetDevice(GPU_idev);
    cudaMemcpy(GPU_image, transp, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyHostToDevice);
    free(transp);
}
extern "C" void GPU_OPER_copyIlumin_GPU2CPU(float *ilumin) {
    int ix, iz;
    float *transp = (float *) calloc(GPU_nx*GPU_nz, sizeof(float));
    cudaSetDevice(GPU_idev);
    cudaMemcpy(transp, GPU_ilumin, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyDeviceToHost);
    for(ix=0; ix<GPU_nx; ix++)
        for(iz=0; iz<GPU_nz; iz++)
            ilumin[iz+ix*GPU_nz] = transp[ix+iz*GPU_nx];
    free(transp);
}
extern "C" void GPU_OPER_copyImage_BRNSCT_GPU2CPU(float *imageSouSide, float *imageRecSide, int add) {
    int ix, iz;
    float *transpSouSide = (float *) calloc(GPU_nx*GPU_nz, sizeof(float));
    float *transpRecSide = (float *) calloc(GPU_nx*GPU_nz, sizeof(float));
    cudaSetDevice(GPU_idev);
    cudaMemcpy(transpSouSide, GPU_imageSouSide, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(transpRecSide, GPU_imageRecSide, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyDeviceToHost);
    for(ix=0; ix<GPU_nx; ix++)
        for(iz=0; iz<GPU_nz; iz++) {
            if(add==0) {
                imageSouSide[iz+ix*GPU_nz] = transpSouSide[ix+iz*GPU_nx];
                imageRecSide[iz+ix*GPU_nz] = transpRecSide[ix+iz*GPU_nx];
            }
            else {
                imageSouSide[iz+ix*GPU_nz] += transpSouSide[ix+iz*GPU_nx];
                imageRecSide[iz+ix*GPU_nz] += transpRecSide[ix+iz*GPU_nx];
            }
        }
    free(transpSouSide);
    free(transpRecSide);
}
extern "C" void GPU_OPER_copyIlumin_BRNSCT_GPU2CPU(float *iluminSouSide, float *iluminRecSide, int add) {
    int ix, iz;
    float *transpSouSide = (float *) calloc(GPU_nx*GPU_nz, sizeof(float));
    float *transpRecSide = (float *) calloc(GPU_nx*GPU_nz, sizeof(float));
    cudaSetDevice(GPU_idev);
    cudaMemcpy(transpSouSide, GPU_iluminSouSide, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(transpRecSide, GPU_iluminRecSide, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyDeviceToHost);
    for(ix=0; ix<GPU_nx; ix++)
        for(iz=0; iz<GPU_nz; iz++) {
            if(add==0) {
                iluminSouSide[iz+ix*GPU_nz] = transpSouSide[ix+iz*GPU_nx];
                iluminRecSide[iz+ix*GPU_nz] = transpRecSide[ix+iz*GPU_nx];
            }
            else {
                iluminSouSide[iz+ix*GPU_nz] += transpSouSide[ix+iz*GPU_nx];
                iluminRecSide[iz+ix*GPU_nz] += transpRecSide[ix+iz*GPU_nx];
            }
        }
    free(transpSouSide);
    free(transpRecSide);
}

extern "C" void GPU_OPER_copyPresentWavefieldWS_CPU2GPU(float *wavefield) {
    cudaSetDevice(GPU_idev);
    cudaMemcpy(GPU_ws, wavefield, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyHostToDevice);
}
extern "C" void GPU_OPER_copyPresentWavefieldWC_CPU2GPU(float *wavefield) {
    cudaSetDevice(GPU_idev);
    cudaMemcpy(GPU_wc, wavefield, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyHostToDevice);
}
extern "C" void GPU_OPER_copyPresentWavefieldWS_GPU2CPU(float *wavefield) {
    cudaSetDevice(GPU_idev);
    cudaMemcpy(wavefield, GPU_ws, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyDeviceToHost);
}
extern "C" void GPU_OPER_copyPresentWavefieldWC_GPU2CPU(float *wavefield) {
    cudaSetDevice(GPU_idev);
    cudaMemcpy(wavefield, GPU_wc, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyDeviceToHost);
}
extern "C" void GPU_OPER_copyPresentWavefield_CPU2GPU(float *wavefield) {
    cudaSetDevice(GPU_idev);
    cudaMemcpy(GPU_wav2, wavefield, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyHostToDevice);
}
extern "C" void GPU_OPER_copyPresentWavefield_GPU2CPU(float *wavefield)  {
    cudaSetDevice(GPU_idev);
    cudaMemcpy(wavefield, GPU_wav2, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyDeviceToHost);
}
extern "C" void GPU_OPER_copyPresentScatteredWavefield_GPU2CPU(float *wavefield)  {
    cudaSetDevice(GPU_idev);
    cudaMemcpy(wavefield, GPU_sct2, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyDeviceToHost);
}
extern "C" void GPU_OPER_copyPresentWavefield_secDer_GPU2CPU(float *wavefield)  {
    cudaSetDevice(GPU_idev);
    cudaMemcpy(wavefield, GPU_secDer    , GPU_nx*GPU_nz*sizeof(float), cudaMemcpyDeviceToHost);
}
extern "C" void GPU_OPER_copyPresentScatteredWavefield_secDer_GPU2CPU(float *wavefield)  {
    cudaSetDevice(GPU_idev);
    cudaMemcpy(wavefield, GPU_secDer_sct, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyDeviceToHost);
}
extern "C" void GPU_OPER_copyLaplacian_GPU2CPU(float *lap) {
    int ix, iz;
    float *transp = (float *) calloc(GPU_nx*GPU_nz, sizeof(float));
    cudaSetDevice(GPU_idev);
    cudaMemcpy(transp, GPU_lap, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyDeviceToHost);
    for(ix=0; ix<GPU_nx; ix++)
        for(iz=0; iz<GPU_nz; iz++)    lap[iz+ix*GPU_nz] = transp[ix+iz*GPU_nx];
    free(transp);
}
extern "C" void GPU_OPER_setTimeStep(int it) {
    GPU_it = it;
}

extern "C" void GPU_OPER_zeroSecDerIncArray() {
    cudaMemset(GPU_secDer    , 0, GPU_nx*GPU_nz);
}
extern "C" void GPU_OPER_zeroSecDerSctArray() {
    cudaMemset(GPU_secDer_sct, 0, GPU_nx*GPU_nz);
}

extern "C" void GPU_OPER_wavefieldTimeSecondDerivative_part1(int mode){
    float *secDer, *u2, *u1;
    cudaSetDevice(GPU_idev);
    if(mode==0) {   
        secDer = GPU_secDer;
        u2 = GPU_wav2;
        u1 = GPU_wav1;
    }
    if(mode==1) {
        secDer = GPU_secDer_sct;
        u2 = GPU_sct2;
        u1 = GPU_sct1;
    }
    GPU_OPER_wavefieldTimeSecondDerivative_part1 <<< GPU_nBlocks, GPU_nThreads >>> (secDer, u2, u1, GPU_invDt2, GPU_nx);
    cudaDeviceSynchronize();
}
__global__ void GPU_OPER_wavefieldTimeSecondDerivative_part1(float *secDer, float *u2, float *u1, float invDt2, int nx) {
    int ix, iz, i;
    ix = tx + bx*bdx;
    iz = tz + bz*bdz;
    i = id2(ix,iz);
    secDer[i] = (-2.0f*u2[i] + u1[i]) * invDt2;
    __syncthreads();
}

extern "C" void GPU_OPER_wavefieldTimeSecondDerivative_part2(int mode){
    float *secDer, *u2;
    cudaSetDevice(GPU_idev);
    if(mode==0) {   
        secDer = GPU_secDer;
        u2 = GPU_wav2;
    }
    if(mode==1) {   
        secDer = GPU_secDer_sct;
        u2 = GPU_sct2;
    }
    GPU_OPER_wavefieldTimeSecondDerivative_part2 <<< GPU_nBlocks, GPU_nThreads >>> (secDer, u2, GPU_invDt2, GPU_nx);
    cudaDeviceSynchronize();
}
__global__ void GPU_OPER_wavefieldTimeSecondDerivative_part2(float *secDer, float *u2, float invDt2, int nx) {
    int ix, iz, i;
    ix = tx + bx*bdx;
    iz = tz + bz*bdz;
    i = id2(ix,iz);
    secDer[i] += u2[i] * invDt2;
    __syncthreads();
}

extern "C" void GPU_WRAP_scatter_Acoustic_lx_opt() {
    dim3 blocks, threads;
    cudaSetDevice(GPU_idev);
    threads.x = 32;
    threads.z = 8;
    blocks.x = GPU_nx-2*GPU_nxb;
    blocks.z = (GPU_nz-2*GPU_nzb)/8;
    GPU_DEV_scatter_Acoustic_lx_opt <<< blocks, threads >>> (GPU_eimage, GPU_secDer, GPU_sct1, GPU_nx, GPU_nz, GPU_lx0, GPU_nxb, GPU_nzb);
    cudaDeviceSynchronize();
return;
}
__global__ void GPU_DEV_scatter_Acoustic_lx_opt(float *eref, float *ws, float *sct, int nx, int nz, int lx0, int nxb, int nzb) {
    int ix, iz, lx, nlx=2*lx0+1;
    ix =      bx     + nxb;
    iz = tz + bz*bdz + nzb;
    int lxIni = -lx0  + tx;
    int lxFin = +lx0;
    for(lx=lxIni; lx<=lxFin; lx+=bdx)
        if( (ix+lx)<(nx-2*nxb)  &&  (ix-lx)>nxb ) {
            // sct[id2(ix+lx,iz)] += ws[id2(ix-lx,iz)] * eref[id3_eic_lx(lx,ix,iz)];
            float val = ws[id2(ix-lx,iz)] * eref[id3_eic_lx(lx,ix,iz)];
            atomicAdd(&(sct[id2(ix+lx,iz)]), val);
        }
    
    __syncthreads();
}

// >>>>>>>>>>>>>>>>>>>>>>>  <<<<<<<<<<<<<<<<<<<<<<<<
// >>>>>>>>> Wavefield imaging functions <<<<<<<<<<<
// >>>>>>>>>>>>>>>>>>>>>>>  <<<<<<<<<<<<<<<<<<<<<<<<
// Compute the 2D explicit time-step
extern "C" void GPU_WRAP_imagingCondition(int mode) {
    cudaSetDevice(GPU_idev);
    if(mode==IC_REGULAR_IMAGE) {
        GPU_DEV_imagingCondition <<< GPU_nBlocks, GPU_nThreads >>> (GPU_image, GPU_ws, GPU_wav2, GPU_nx, GPU_nz);
    }
    if(mode==IC_REGULAR_ILUMIN) {
        GPU_DEV_imagingCondition <<< GPU_nBlocks, GPU_nThreads >>> (GPU_ilumin, GPU_wav2, GPU_wav2, GPU_nx, GPU_nz);
    }
    if(mode==IC_BORNSCT_IMAGE) {
        GPU_DEV_imagingCondition <<< GPU_nBlocks, GPU_nThreads >>> (GPU_imageSouSide, GPU_sct2, GPU_ws, GPU_nx, GPU_nz);
        GPU_DEV_imagingCondition <<< GPU_nBlocks, GPU_nThreads >>> (GPU_imageRecSide, GPU_wav2, GPU_wc, GPU_nx, GPU_nz);
    }
    if(mode==IC_BORNSCT_IMAGE2) {
        GPU_DEV_imagingCondition <<< GPU_nBlocks, GPU_nThreads >>> (GPU_imageSouSide, GPU_sct2  , GPU_ws, GPU_nx, GPU_nz);
        GPU_DEV_imagingCondition <<< GPU_nBlocks, GPU_nThreads >>> (GPU_imageRecSide, GPU_secDer, GPU_wc, GPU_nx, GPU_nz);
    }
    if(mode==IC_BORNSCT_ILUMIN) {
        GPU_DEV_imagingCondition <<< GPU_nBlocks, GPU_nThreads >>> (GPU_iluminSouSide, GPU_wav2, GPU_wav2, GPU_nx, GPU_nz);
        GPU_DEV_imagingCondition <<< GPU_nBlocks, GPU_nThreads >>> (GPU_iluminRecSide, GPU_sct2, GPU_sct2, GPU_nx, GPU_nz);
    }
    cudaDeviceSynchronize();
}
__global__ void GPU_DEV_imagingCondition(float *image, const float *u1, const float *u2, int nx, int nz) {
    int ix, iz, i;
    ix = tx + bx*bdx;
    iz = tz + bz*bdz;
    i = id2(ix,iz);
    image[i] += u2[i]*u1[i];
    __syncthreads();
}



extern "C" void GPU_OPER_EIC_storeWavefields_GPU2GPU(int it, int side, int mode)  {
    if( GPU_EIC_jt<=0  ||  (it%GPU_EIC_jt)!=0 )    return;

    cudaSetDevice(GPU_idev);
    size_t shift = (it/GPU_EIC_jt) * GPU_nx*GPU_nz;

    if(mode==IC_BORNSCT_IMAGE)
    {
        if(side==IC_SOU_WAVE) {
            cudaMemcpy(GPU_sou_inc+shift, GPU_wav2, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyDeviceToDevice);
            cudaMemcpy(GPU_sou_sct+shift, GPU_sct2, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyDeviceToDevice);
        }
        if(side==IC_REC_WAVE) {
            cudaMemcpy(GPU_rec_inc+shift, GPU_wav2, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyDeviceToDevice);
            cudaMemcpy(GPU_rec_sct+shift, GPU_sct2, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyDeviceToDevice);
        }
    }
    if(mode==IC_BORNSCT_IMAGE2)
    {
        if(side==IC_SOU_WAVE) {
            cudaMemcpy(GPU_sou_inc+shift, GPU_secDer    , GPU_nx*GPU_nz*sizeof(float), cudaMemcpyDeviceToDevice);
            cudaMemcpy(GPU_sou_sct+shift, GPU_secDer_sct, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyDeviceToDevice);
        }
        if(side==IC_REC_WAVE) {
            cudaMemcpy(GPU_rec_inc+shift, GPU_wav2, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyDeviceToDevice);
            cudaMemcpy(GPU_rec_sct+shift, GPU_sct2, GPU_nx*GPU_nz*sizeof(float), cudaMemcpyDeviceToDevice);
        }
    }

}
extern "C" void GPU_WRAP_imagingCondition_transmission() {
    
    int it;
    dim3 blocks, threads;
    cudaSetDevice(GPU_idev);    

    threads.x = 32;
    threads.z = 8;
    blocks.x = GPU_nx-2*GPU_nxb;
    blocks.z = (GPU_nz-2*GPU_nzb)/8;
    
    for(it=0; it<GPU_EIC_nt; it++)
    {
        GPU_DEV_imagingCondition_transmission <<< blocks, threads >>> (GPU_vp, GPU_imageSouSide, GPU_sou_inc, GPU_rec_sct, \
                                                                       GPU_nx, GPU_nz, GPU_EIC_nt, GPU_dx, GPU_dz, GPU_EIC_dt, \
                                                                       GPU_lx0, GPU_lz0, GPU_nxb, GPU_nzb, it);
        GPU_DEV_imagingCondition_transmission <<< blocks, threads >>> (GPU_vp, GPU_imageRecSide, GPU_sou_sct, GPU_rec_inc, \
                                                                       GPU_nx, GPU_nz, GPU_EIC_nt, GPU_dx, GPU_dz, GPU_EIC_dt, \
                                                                       GPU_lx0, GPU_lz0, GPU_nxb, GPU_nzb, it);  
    }
    cudaDeviceSynchronize();
return;
}

__global__ void GPU_DEV_imagingCondition_transmission(float *vel, float *image, float *ws, float *wr, int nx, int nz, int nt, \
                                                       float dx, float dz, float dt, int lx0, int lz0, int nxb, int nzb, int it) {
    
    int ix, iz, lx, lz, its0, its1, itr0, itr1, i;
    float  t, ts, tr, tau, lambda, lambda_x, lambda_z, dt_;

    ix =      bx     + nxb;
    iz = tz + bz*bdz + nzb;

    int lxIni = -lx0  + tx;
    int lxFin = +lx0;
    int lzIni = -lz0;
    int lzFin = +lz0;

    t = it*dt;
    dt_ = 1.0f/dt;
    
    for(lx=lxIni; lx<=lxFin; lx+=bdx)
        for(lz=lzIni; lz<=lzFin; lz++)
            if( nxb<=(ix+lx) && (ix+lx)<(nx-nxb)  &&  nxb<=(ix-lx) && (ix-lx)<(nx-nxb) )
            {
                i = id2(ix,iz);
                lambda_x = lx*dx;
                lambda_z = lz*dz;
                lambda = sqrtf(lambda_x*lambda_x + lambda_z*lambda_z);
                tau = lambda/vel[i];
                ts = t+tau;
                tr = t-tau;
                its0 = (int)(ts*dt_);
                its1 = its0 + 1;
                itr0 = (int)(tr*dt_);
                itr1 = itr0 + 1;
                if(its0>=0 && its1<nt  &&  itr0>=0 && itr1<nt)
                {
                    float ps0 = its1 - ts*dt_;
                    float ps1 = ts*dt_ - its0;
                    float pr0 = itr1 - tr*dt_;
                    float pr1 = tr*dt_ - itr0;
                    float pstot = ps0+ps1;
                    float prtot = pr0+pr1;
                    ps0 /= pstot;
                    ps1 /= pstot;
                    pr0 /= prtot;
                    pr1 /= prtot;
                    image[i] += (ps0*ws[id2t(ix-lx,iz-lz,its0)]+ps1*ws[id2t(ix-lx,iz-lz,its1)]) * (pr0*wr[id2t(ix+lx,iz+lz,itr0)]+pr1*wr[id2t(ix+lx,iz+lz,itr1)])
                              + (pr0*ws[id2t(ix-lx,iz-lz,itr0)]+pr1*ws[id2t(ix-lx,iz-lz,itr1)]) * (ps0*wr[id2t(ix+lx,iz+lz,its0)]+ps1*wr[id2t(ix+lx,iz+lz,its1)]);
                }
            }
        __syncthreads();
}


extern "C" void GPU_WRAP_eimagingCondition_lx() {
    dim3 blocks, threads;
    cudaSetDevice(GPU_idev);
    // cudaDeviceProp prop;
    threads.x = 32;
    threads.z = 8;
    blocks.x = GPU_nx-2*GPU_nxb;
    blocks.z = (GPU_nz-2*GPU_nzb)/8;
    size_t szShared = (threads.x + 2*GPU_lx0) * threads.z;
    GPU_DEV_eimagingCondition_lx <<< blocks, threads, szShared >>> (GPU_eimage, GPU_ws, GPU_wav2, GPU_nx, GPU_nz, GPU_lx0, GPU_nxb, GPU_nzb);
    cudaDeviceSynchronize();
return;
}
__global__ void GPU_DEV_eimagingCondition_lx(float *eimage, float *ws, float *wr, int nx, int nz, int lx0, int nxb, int nzb) {
    int ix, iz, lx, nlx=2*lx0+1;

    ix =      bx     + nxb;
    iz = tz + bz*bdz + nzb;

    int lxIni = -lx0  + tx;
    int lxFin = lx0;
    /*
    for(lx=lxIni; lx<=lxFin; lx+=bdx) {
        slice[lx+lx0 + tz*nlx          ] = ws[id2(ix-lx,iz)];
        slice[lx+lx0 + tz*nlx + nlx*bdz] = wr[id2(ix+lx,iz)];
    }
    */
    for(lx=lxIni; lx<=lxFin; lx+=bdx)
        // if( (ix+lx)<(nx-2*nxb)  &&  (ix-lx)>nxb ) {
        if( nxb<=(ix+lx) && (ix+lx)<(nx-nxb)  &&  nxb<=(ix-lx) && (ix-lx)<(nx-nxb) ) {
            atomicAdd(&(eimage[id3_eic_lx(lx,ix,iz)]), -(ws[id2(ix-lx,iz)]*wr[id2(ix+lx,iz)]) );
            // eimage[id3_eic_lx(lx,ix,iz)] -= ws[id2(ix-lx,iz)] * wr[id2(ix+lx,iz)];
        }
    
    __syncthreads();

}

//*
extern "C" void GPU_WRAP_eimagingCondition_lx_lz() {
    dim3 blocks, threads;
    cudaSetDevice(GPU_idev);
    threads.x = 32;
    threads.z = 8;
    blocks.x = GPU_nx-2*GPU_nxb;
    blocks.z = (GPU_nz-2*GPU_nzb)/8;
    GPU_DEV_eimagingCondition_lx_lz <<< blocks, threads >>> (GPU_eimage, GPU_ws, GPU_wav2, GPU_nx, GPU_nz, GPU_lx0, GPU_lz0, GPU_nxb, GPU_nzb);
    cudaDeviceSynchronize();
return;
}
__global__ void GPU_DEV_eimagingCondition_lx_lz(float *eimage, float *ws, float *wr, int nx, int nz, \
                                                    int lx0, int lz0, int nxb, int nzb) {
    
    int ix, iz, lx, lz;
    int nlx =2*lx0+1;
    int nlz =2*lz0+1;

    ix =      bx     + nxb;
    iz = tz + bz*bdz + nzb;

    int lxIni = -lx0  + tx;
    int lxFin = +lx0;
    int lzIni = -lz0;
    int lzFin = +lz0;
    
    for(lx=lxIni; lx<=lxFin; lx+=bdx)
        for(lz=lzIni; lz<=lzFin; lz++)
            // if( (ix+lx)<(nx-2*nxb)  &&  (ix-lx)>nxb  &&  (iz+lz)<(nz-2*nzb)  &&  (iz-lz)>nzb)
            if( nxb<=(ix+lx) && (ix+lx)<(nx-nxb)  &&  nxb<=(ix-lx) && (ix-lx)<(nx-nxb)  && \
                nzb<=(iz+lz) && (iz+lz)<(nz-nzb)  &&  nzb<=(iz-lz) && (iz-lz)<(nz-nzb) ) 
            {    
                eimage[id3_eic_lx_lz(lx,lz,ix,iz)] -= ws[id2(ix-lx,iz-lz)] * wr[id2(ix+lx,iz+lz)];
            }
    
    __syncthreads();

}
//*/

/*
void GPU_WRAP_eimagingCondition_lx_lz_lt() {
    dim3 blocks, threads;
    cudaSetDevice(GPU_idev);
    threads.x = 32;
    threads.z = 8;
    blocks.x = GPU_nx-2*GPU_nxb;
    blocks.z = (GPU_nz-2*GPU_nzb)/8;
    GPU_DEV_eimagingCondition_lx_lz_tau <<< blocks, threads >>> (GPU_eimage, GPU_ws, GPU_wav2, GPU_nx, GPU_nz, GPU_lx0, GPU_lz0, GPU_lt0, GPU_nxb, GPU_nzb);
    cudaDeviceSynchronize();
return;
}
__global__ void GPU_DEV_eimagingCondition_lx_lz_tau(float *eimage, float *ws, float *wr, int nx, int nz, \
                                                    int lx0, int lz0, int lt0, int nxb, int nzb) {
    
    int ix, iz, lx, lz, tau;
    int nlx =2*lx0+1;
    int nlz =2*lz0+1;
    int nlt=2*lt0+1;

    ix =      bx     + nxb;
    iz = tz + bz*bdz + nzb;

    int lxIni = -lx0  + tx;
    int lxFin = +lx0;
    int lzIni = -lz0;
    int lzFin = +lz0;
    int ltIni = -lt0;
    int ltFin = +lt0;
    
    for(lx=lxIni; lx<=lxFin; lx+=bdx)
        for(lz=lzIni; lz<=lzFin; lz++)
            for(lt=ltIni; lt<=ltFin; lt++)
                if( (ix+lx)<(nx-2*nxb)  &&  (ix-lx)>nxb )
                    eimage[id3_eic_lx_lz_lt(lx,lz,tau,ix,iz)] -= ws[id2t(ix-lx,iz-lz,tau)] * wr[id2t(ix+lx,iz+lz,-tau)];
    
    __syncthreads();

}
//*/


// >>>>>>>>>>>>>>>>>>>>>>>  <<<<<<<<<<<<<<<<<<<<<<<<
// >>>>>>>> Wavefield propagation functions <<<<<<<<
// >>>>>>>>>>>>>>>>>>>>>>>  <<<<<<<<<<<<<<<<<<<<<<<<

extern "C" void GPU_WRAP_recordInjectData(int mode) {
    dim3 threads, blocks;
    threads.x = 32;
    threads.y = 1;
    threads.z = 1;
    blocks.x  = 1 + GPU_nrec/threads.x;
    blocks.y  = 1;
    blocks.z  = 1;
    cudaSetDevice(GPU_idev);
    float *wavefieldArray;
    if(mode==0)    wavefieldArray = GPU_wav2;
    if(mode==1)    wavefieldArray = GPU_wav1;
    if(mode==2)    wavefieldArray = GPU_sct2;
    GPU_DEV_recordData <<< blocks, threads >>> (mode, GPU_seismogramMultiplexed+(GPU_it*GPU_nrec), \
                                                GPU_recIndexes, wavefieldArray, GPU_nx, GPU_nz, GPU_nrec);
    cudaDeviceSynchronize();
}
__global__ void GPU_DEV_recordData(int mode, float *seismogramMultiplexed, int *recIndexes, float *wav, int nx, int nz, int nrec) {
    
    int ir = tx + bx*bdx;
    if(ir<nrec)
    {
        int izr = recIndexes[2*ir + 0];
        int ixr = recIndexes[2*ir + 1];
        int i = id2(ixr,izr);
        if(mode==0 || mode==2)    seismogramMultiplexed[ir]  = wav[i];
        if(mode==1)               wav[i] -= seismogramMultiplexed[ir];
    }
    __syncthreads();

return;
}


extern "C" void GPU_WRAP_injectSource_Acoustic()  {
    cudaSetDevice(GPU_idev);
    GPU_DEV_injectSource_Acoustic <<< 1,1 >>>(GPU_wav1, GPU_source, GPU_it, GPU_ixs, GPU_izs, GPU_nx, GPU_nz, GPU_dx, GPU_dz);
    cudaDeviceSynchronize();
return;
}
__global__ void GPU_DEV_injectSource_Acoustic(float *u, float *source, int it, int ix0, int iz0, int nx, int nz, float dx, float dz)
{
    if(bx==0 && by==0 && bz==0 && tx==0 && ty==0 && tz==0)
    {
        int i = id2(ix0,iz0);
        u[i] -= source[it]/(dx*dz);
    }
    __syncthreads();
return;
}

extern "C" void GPU_WRAP_propagate_Acoustic_abs2() {
    GPU_LNCH_laplacian_10thOrder_2D_Acoustic_Iso(GPU_wav2, GPU_lap);
    GPU_LNCH_timeStep_Acoustic_abs2(GPU_wav2, GPU_wav1, GPU_lap);
    float *tmp = GPU_wav2;
    GPU_wav2   = GPU_wav1;
    GPU_wav1   = tmp;

return;
}
extern "C" void GPU_WRAP_propagate_Acoustic_Scattered_abs2() {
    GPU_LNCH_laplacian_10thOrder_2D_Acoustic_Iso(GPU_sct2, GPU_lap);
    GPU_LNCH_timeStep_Acoustic_abs2(GPU_sct2, GPU_sct1, GPU_lap);
    float *tmp = GPU_sct2;
    GPU_sct2   = GPU_sct1;
    GPU_sct1   = tmp;
return;
}
// Compute the 2D explicit time-step
void GPU_LNCH_timeStep_Acoustic_abs2(float *u2, float *u1, float *lap) {
    cudaSetDevice(GPU_idev);
    GPU_DEV_timeStep_Acoustic_abs2 <<< GPU_nBlocks, GPU_nThreads >>> (u1, u2, lap, GPU_vp, GPU_sigma, GPU_sigmaInv, GPU_nx, GPU_nz);
    cudaDeviceSynchronize();
}
__global__ void GPU_DEV_timeStep_Acoustic_abs2(float *u1, const float *u2, const float *lap, const float *vp2dt2dx2, \
                                               const float *sigma, const float *sigmaInv, int nx, int nz) {
    int ix, iz, i;
    ix = tx + bx*bdx;
    iz = tz + bz*bdz;
    i = id2(ix,iz);
    u1[i] = (2.0f*u2[i] - (1.0f-sigma[i])*u1[i] + vp2dt2dx2[i]*lap[i]) * sigmaInv[i];
    __syncthreads();
}

// Compute the 2D laplacian to 10th order
// Division by dx^2 is embedded in the coefficients C0 through C5.
void GPU_LNCH_laplacian_10thOrder_2D_Acoustic_Iso(float *wavefield, float *lap) {
    cudaSetDevice(GPU_idev);
    size_t sizeForShared = (GPU_nThreads.x+2*OL) * (GPU_nThreads.z+2*OL) * sizeof(float);
    GPU_DEV_laplacian_10thOrder_2D_Acoustic_Iso_opt0 <<< GPU_nBlocks, GPU_nThreads, sizeForShared >>> (lap, wavefield, GPU_nx, GPU_nz);
    cudaDeviceSynchronize();
}
__global__ void GPU_DEV_laplacian_10thOrder_2D_Acoustic_Iso(float *lap, const float *u2, int nx, int nz) {
    int ix, iz, i;
    ix = tx + bx*bdx;
    iz = tz + bz*bdz;
    i = id2(ix,iz);
    lap[i] =      C5_10 * (u2[id2(ix  ,iz+5)]+u2[id2(ix  ,iz-5)]
                          +u2[id2(ix+5,iz  )]+u2[id2(ix-5,iz  )])
           +      C4_10 * (u2[id2(ix  ,iz+4)]+u2[id2(ix  ,iz-4)]
                          +u2[id2(ix+4,iz  )]+u2[id2(ix-4,iz  )])
           +      C3_10 * (u2[id2(ix  ,iz+3)]+u2[id2(ix  ,iz-3)]
                          +u2[id2(ix+3,iz  )]+u2[id2(ix-3,iz  )])
           +      C2_10 * (u2[id2(ix  ,iz+2)]+u2[id2(ix  ,iz-2)]
                          +u2[id2(ix+2,iz  )]+u2[id2(ix-2,iz  )])
           +      C1_10 * (u2[id2(ix  ,iz+1)]+u2[id2(ix  ,iz-1)]
                          +u2[id2(ix+1,iz  )]+u2[id2(ix-1,iz  )])
           - 4.0f*C0_10 *  u2[id2(ix  ,iz  )];
    
    __syncthreads();
}
__global__ void GPU_DEV_laplacian_10thOrder_2D_Acoustic_Iso_opt0(float *lap, const float *u2, int nx, int nz)
{
    int ix, iz, i, jx, jz;
    ix = tx + bx*bdx;
    iz = tz + bz*bdz;
    i = id2(ix,iz);

    jx = tx + OL;
    jz = tz + OL;

    extern __shared__ float ushared[];

    __syncthreads();

    
    if(tx<OL) {
        if(ix>=OL  )    ushared[jd2(jx-OL, jz   )] = u2[id2(ix-OL,iz   )];
        else            ushared[jd2(jx-OL, jz   )] = 0.0f;
    }
    if(tx>=bdx-OL) {
        if(ix<nx-OL)    ushared[jd2(jx+OL, jz   )] = u2[id2(ix+OL,iz   )];
        else            ushared[jd2(jx+OL, jz   )] = 0.0f;
    }       
    if(tz<OL) {
        if(iz>=OL  )    ushared[jd2(jx   , jz-OL)] = u2[id2(ix   ,iz-OL)];
        else            ushared[jd2(jx   , jz-OL)] = 0.0f;
    }
    if(tz>=bdz-OL) {
        if(iz<nz-OL)    ushared[jd2(jx   , jz+OL)] = u2[id2(ix   ,iz+OL)];
        else            ushared[jd2(jx   , jz+OL)] = 0.0f;
    }
    ushared[jd2(jx,jz)] = u2[id2(ix,iz)];

    __syncthreads();

    
    lap[i] =        C5_10 * ( ushared[jd2(jx  ,jz+5)] + ushared[jd2(jx  ,jz-5)]
                            + ushared[jd2(jx+5,jz  )] + ushared[jd2(jx-5,jz  )])
           +        C4_10 * ( ushared[jd2(jx  ,jz+4)] + ushared[jd2(jx  ,jz-4)]
                            + ushared[jd2(jx+4,jz  )] + ushared[jd2(jx-4,jz  )])
           +        C3_10 * ( ushared[jd2(jx  ,jz+3)] + ushared[jd2(jx  ,jz-3)]
                            + ushared[jd2(jx+3,jz  )] + ushared[jd2(jx-3,jz  )])
           +        C2_10 * ( ushared[jd2(jx  ,jz+2)] + ushared[jd2(jx  ,jz-2)]
                            + ushared[jd2(jx+2,jz  )] + ushared[jd2(jx-2,jz  )])
           +        C1_10 * ( ushared[jd2(jx  ,jz+1)] + ushared[jd2(jx  ,jz-1)]
                            + ushared[jd2(jx+1,jz  )] + ushared[jd2(jx-1,jz  )])
           - 4.0f * C0_10 *   ushared[jd2(jx  ,jz  )];

    __syncthreads();
}



// ******************************************************************
// ******************************************************************
// ********************* Math Functions in GPU **********************
// ******************************************************************
// ******************************************************************
//*
extern "C" void GPU_WRAP_lx2theta(int adj, int gpuDev, int add, float *eimage, \
                       float *eimage_theta, int nx, float dx, int nz, \
                       float dz, int nxb, int nzb, int lx0, int ntheta, \
                       float dtheta, float otheta) {

    cudaSetDevice(gpuDev);
    
    float *eimage_gpu, *eimage_theta_gpu;

    int nlx = 2*lx0 + 1;
    GPU_OPER_alocMem((void**) &eimage_gpu      , nx*nz*nlx    * sizeof(float));
    GPU_OPER_alocMem((void**) &eimage_theta_gpu, nx*nz*ntheta * sizeof(float));

    if(adj==0)
        cudaMemcpy(eimage_gpu      , eimage      , nx*nz*nlx   *sizeof(float), cudaMemcpyHostToDevice);
    else
        cudaMemcpy(eimage_theta_gpu, eimage_theta, nx*nz*ntheta*sizeof(float), cudaMemcpyHostToDevice);

    dim3 blocks, threads;
    threads.x = 32;
    threads.y = 1;
    threads.z = 8;
    blocks.x  = nx/threads.x;
    blocks.y  = 1;
    blocks.z  = nz/threads.z;

    GPU_DEV_lx2theta <<< blocks, threads >>> (adj, 0, eimage_gpu, eimage_theta_gpu, \
                                              nx, dx, nz, dz, nxb, nzb, lx0, ntheta, \
                                              dtheta, otheta);

    if(adj==0)
       cudaMemcpy(eimage_theta, eimage_theta_gpu, nx*nz*ntheta*sizeof(float), cudaMemcpyDeviceToHost);
    else
       cudaMemcpy(eimage      , eimage_gpu      , nx*nz*nlx   *sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(eimage_gpu);
    cudaFree(eimage_theta_gpu);

    cudaDeviceSynchronize();
}
__global__ void GPU_DEV_lx2theta(int adj, int add, float *eimage, float *eimage_theta, \
                                 int nx, float dx, int nz, float dz, int nxb, int nzb, \
                                 int lx0, int ntheta, float dtheta, float otheta) {

    int ix, iz, lx, itheta, ltheta0=ntheta/2;
    int nlx=2*lx0+1;

    ix = tx + bx*bdx;
    iz = tz + bz*bdz;

    if(add == 0)
    {
        if(adj==0) {
            for(itheta=-ltheta0; itheta<=ltheta0; itheta++) {
                int   i_th  = id3_eic_th_transp(itheta,iz,ix);
                eimage_theta[i_th] = 0.0f;
            }
        }
        else {
            for(lx=-lx0; lx<=lx0; lx++) {
                int i = id3_eic_lx_transp(lx,iz,ix);
                eimage[i] = 0.0;
            }
        }
    }

    __syncthreads();

    // Slant Stack to get Angle Gather
    int   iz_0, iz_1, i_0, i_1;
    float len, deltaZ, z, z_0, z_1;
    float pz_0, pz_1, ptot;
    for(itheta=-ltheta0; itheta<=ltheta0; itheta++)
    {
        int   i_th  = id3_eic_th_transp(itheta,iz,ix);
        float theta = (itheta*dtheta) * M_PI/180.0f;
        for(lx=-lx0; lx<=lx0; lx++)
        {
            len    = lx*dx;
            deltaZ = len * tanf(theta);
            z      = dz*iz + deltaZ;
            iz_0   = floorf(z/dz);
            if(iz_0<=nzb)         continue;
            if(iz_0>=nz-nzb-1)    continue;
            iz_1   = iz_0 + 1;
            z_0    = iz_0*dz;
            z_1    = iz_1*dz;
            pz_0   = z_1 - z;
            pz_1   = z - z_0;
            ptot   = pz_0 + pz_1;
            i_0    = id3_eic_lx_transp(lx,iz_0,ix);
            i_1    = id3_eic_lx_transp(lx,iz_1,ix);
            if(adj==0) {
                atomicAdd(&(eimage_theta[i_th]), (pz_0*eimage[i_0] + pz_1*eimage[i_1])/ptot);
                // eimage_theta[i_th] += (pz_0*eimage[i_0] + pz_1*eimage[i_1])/ptot;
            }
            else {
                atomicAdd(&(eimage[i_0]), pz_0*eimage_theta[i_th]/(ptot*ntheta));
                atomicAdd(&(eimage[i_1]), pz_1*eimage_theta[i_th]/(ptot*ntheta));
                // eimage[i_0] += pz_0*eimage_theta[i_th]/(ptot*ntheta);
                // eimage[i_1] += pz_1*eimage_theta[i_th]/(ptot*ntheta);
            }
        }
    }
    __syncthreads();
return;
}
//*
extern "C" void GPU_WRAP_Contraction_lx(int gpuDev, int add, float *eimageInp, float *eimageOut, \
                                        int nx, float dx, int nz, float dz, int lx0, int nxb, int nzb, \
                                        float xMin, float xMax, float zMin, float zMax, float contraction_factor) {

    cudaSetDevice(gpuDev);

    int ixMin = nxb + xMin/dx;
    int ixMax = nxb + xMax/dx;
    int izMin = nzb + zMin/dz;
    int izMax = nzb + zMax/dz;
    
    float *eimageInp_gpu, *eimageOut_gpu;

    int nlx = 2*lx0 + 1;
    GPU_OPER_alocMem((void**) &eimageInp_gpu, nx*nz*nlx * sizeof(float));
    GPU_OPER_alocMem((void**) &eimageOut_gpu, nx*nz*nlx * sizeof(float));
    cudaMemcpy(eimageInp_gpu, eimageInp, nx*nz*nlx*sizeof(float), cudaMemcpyHostToDevice);

    dim3 blocks, threads;
    threads.x = 32;
    threads.y = 1;
    threads.z = 8;
    blocks.x  = nx/threads.x;
    blocks.y  = 1;
    blocks.z  = nz/threads.z;

    GPU_DEV_Contraction_lx <<< blocks, threads >>> (0, eimageInp_gpu, eimageOut_gpu, nx, dx, nz, dz, lx0, \
                                                    ixMin, ixMax, izMin, izMax, contraction_factor);

    cudaMemcpy(eimageOut, eimageOut_gpu, nx*nz*nlx*sizeof(float), cudaMemcpyDeviceToHost);
    cudaFree(eimageInp_gpu);
    cudaFree(eimageOut_gpu);

    cudaDeviceSynchronize();
}
__global__ void GPU_DEV_Contraction_lx(int add, float *eimage, float *eimage_contracted,  \
                                       int nx, float dx, int nz, float dz, int lx0, \
                                       int ixMin, int ixMax, int izMin, int izMax, float contraction_factor) {

    int ix, iz, lx;
    int nlx=2*lx0+1;

    ix = tx + bx*bdx;
    iz = tz + bz*bdz;

    if(add == 0)
    {
        for(lx=-lx0; lx<=lx0; lx++) {
            int   i_th  = id3_eic_lx_transp(lx,iz,ix);
            eimage_contracted[i_th] = 0.0f;
        }
    }

    __syncthreads();

    // ODCIG contraction
    float lambda, lambda_;
    int   lx_0;
    int   i;
    if(ixMin<=ix && ix<=ixMax  &&  izMin<=iz && iz<=izMax)
    {
        for(lx=-lx0; lx<=lx0; lx++)
        {
            // Contraction
            lambda  = lx * dx;
            lambda_ = lambda * contraction_factor;
            lx_0    = floorf(lambda_/dx);

            float alpha = (lambda_ - dx*lx_0)/dx;
            float coef[8];
            DEV_interSinc1D_8P(coef, alpha);

            i = id3_eic_lx_transp(lx, iz, ix);
            int i_;
            for(i_=0; i_<8; i_++)
            {
                int ll = lx_0-3+i_;
                int ii =  id3_eic_lx_transp(ll, iz, ix);
                if(-lx0<=ll  &&  ll<=lx0)
                    eimage_contracted[ii] += contraction_factor * coef[i_] * eimage[i];
            }
        } // loop over lx
    }
    __syncthreads();
return;
}
//*/
extern "C" void GPU_WRAP_Demoveout(int gpuDev, int add, int adj, float *eimageInp, float *eimageOut, \
                        int nx, float dx, int nz, float dz, int lx0, \
                        int nxb, int nzb, \
                        float xMin, float xMax, float zMin, float zMax, \
                        float contraction_factor, float lx_null) {

    cudaSetDevice(gpuDev);
    
    float *eimageInp_gpu, *eimageOut_gpu;

    int nlx = 2*lx0 + 1;
    GPU_OPER_alocMem((void**) &eimageInp_gpu, nx*nz*nlx * sizeof(float));
    GPU_OPER_alocMem((void**) &eimageOut_gpu, nx*nz*nlx * sizeof(float));
    cudaMemcpy(eimageInp_gpu, eimageInp, nx*nz*nlx*sizeof(float), cudaMemcpyHostToDevice);

    int ixMin = nxb + xMin/dx;
    int ixMax = nxb + xMax/dx;
    int izMin = nzb + zMin/dz;
    int izMax = nzb + zMax/dz;
    // int ixMin = 0;
    // int ixMax = nx;
    // int izMin = 0;
    // int izMax = nz;

    dim3 blocks, threads;
    threads.x = 32;
    threads.y = 1;
    threads.z = 8;
    blocks.x  = nx/threads.x;
    blocks.y  = 1;
    blocks.z  = nz/threads.z;

    GPU_DEV_Demoveout <<< blocks, threads >>> (0, adj, eimageInp_gpu, eimageOut_gpu, nx, dx, nz, dz, lx0, ixMin, \
                                               ixMax, izMin, izMax, contraction_factor, lx_null);

    cudaMemcpy(eimageOut, eimageOut_gpu, nx*nz*nlx*sizeof(float), cudaMemcpyDeviceToHost);
    cudaFree(eimageInp_gpu);
    cudaFree(eimageOut_gpu);

    cudaDeviceSynchronize();
}

__global__ void GPU_DEV_Demoveout(int add, int adj, float *eimage, float *eimage_deMO,  \
                                  int nx, float dx, int nz, float dz, int lx0, \
                                  int ixMin, int ixMax, int izMin, int izMax, \
                                  float contraction_factor, float lx_null) {

    int ix, iz, lx, iz_, lx_, il;
    int nlx=2*lx0+1;

    ix = tx + bx*bdx;
    iz = tz + bz*bdz;

    float thetaMin = -65.0 * M_PI/180.0;
    float thetaMax = +65.0 * M_PI/180.0;
    float dtheta = 2.0 * M_PI/180.0;
    float theta;
    int   nl = 11;
    float L  = (nl-1)*dx;
    int   nl0 = (nl-1)/2;
    float dl = L/(nl-1);
    float l, Delta_lx, Delta_z;
    float stack_value, plx_0, plx_1, pz_0, pz_1;
    float lx_pos, z_pos, lx_move, z_move, move;

    if(add == 0)
    {
        for(lx=-lx0; lx<=lx0; lx++) {
            int   i_th  = id3_eic_lx_transp(lx,iz,ix);
            eimage_deMO[i_th] = 0.0f;
        }
    }

    __syncthreads();

    // ODCIG contraction
    int   i_z0_lx0, i_z0_lx1, i_z1_lx0, i_z1_lx1;
    int  cond;
    if(ixMin<=ix  &&  ix<=ixMax  &&  izMin<=iz  &&  iz<=izMax)
    {
    for(lx=-lx0+(nl/2); lx<=lx0-(nl/2); lx++)
    {
        // DeMoveout
        lx_pos = lx*dx; 
        z_pos  = iz*dz;
        for(theta=thetaMin; theta<=thetaMax; theta+=dtheta)
        {
            // Slant-Stack
            stack_value = 0.0f;
            int nn = 0;
            for(il=0; il<nl; il++)
            {
                l = (il-nl0)*dl;
                Delta_lx = + l * cosf(theta);
                Delta_z  = - l * sinf(theta);

                lx_   = floorf((lx_pos+Delta_lx)/dx);
                plx_1 = ((lx_pos+Delta_lx) - lx_*dx)/dx;
                plx_0 = ((lx_+1)*dx - (lx_pos+Delta_lx))/dx;

                iz_   = floorf((z_pos+Delta_z)/dz);
                pz_1  = ((z_pos+Delta_z) - iz_*dz)/dz;
                pz_0  = ((iz_+1)*dz-(z_pos+Delta_z))/dz;

                i_z0_lx0 = id3_eic_lx_transp(lx_  , iz_  , ix);
                i_z1_lx0 = id3_eic_lx_transp(lx_  , iz_+1, ix);
                i_z0_lx1 = id3_eic_lx_transp(lx_+1, iz_  , ix);
                i_z1_lx1 = id3_eic_lx_transp(lx_+1, iz_+1, ix);

                // cond = (<=lx_ && lx_<lx0-1) && (0<=iz_ && iz_<nz-1);
                cond = ( (-lx0)<=lx_ && lx_<=(+lx0) )  &&  ( 0<=iz_  &&  iz_<(nz-1) );
                if( !cond )    continue;
                
                nn++;

                float ptot = plx_0*pz_0 + plx_0*pz_1 + plx_1*pz_0 + plx_1*pz_1;

                stack_value += (plx_0*pz_0 * eimage[i_z0_lx0]
                             +  plx_0*pz_1 * eimage[i_z1_lx0]
                             +  plx_1*pz_0 * eimage[i_z0_lx1]
                             +  plx_1*pz_1 * eimage[i_z1_lx1])/ptot;
            }
            stack_value /= nn;

            // Slant-Spread
            if ( (abs(lx)*dx)>=lx_null )
                move = (1.0-contraction_factor) * lx_null;
            else
                move = (1.0-contraction_factor) * abs(lx)*dx;

            if(adj)    move *= -1.0;
            // else if(fabsf(theta)>0.0)
            //     move = (1.0-contraction_factor) * abs(lx)*dx / fabsf(sinf(theta));
        
            if(fsignf_GPU(theta)==fsignf_GPU(lx)) move *= -1.0;
            if(fsignf_GPU(theta)!=fsignf_GPU(lx)) move *= +1.0;
            lx_move = lx_pos + move * sinf(theta);
            z_move =   z_pos + move * cosf(theta);
            //  z_move =  z_pos + move;

            for(il=0; il<nl; il++)
            {
                l = (il-nl0)*dl;
                Delta_lx = + l * cosf(theta);
                Delta_z  = - l * sinf(theta);
                lx_   = floorf((lx_move+Delta_lx)/dx);
                plx_1 = ((lx_move+Delta_lx) - lx_*dx)/dx;
                plx_0 = 1.0 - plx_1;

                iz_   = floorf((z_move+Delta_z)/dz);
                pz_1  = ((z_move+Delta_z) - iz_*dz)/dz;
                pz_0  = 1.0 - pz_1;

                i_z0_lx0 = id3_eic_lx_transp(lx_  , iz_  , ix);
                i_z1_lx0 = id3_eic_lx_transp(lx_  , iz_+1, ix);
                i_z0_lx1 = id3_eic_lx_transp(lx_+1, iz_  , ix);
                i_z1_lx1 = id3_eic_lx_transp(lx_+1, iz_+1, ix);
                cond = ( (-lx0)<=lx_ && lx_<=(lx0-1) )  &&  ( 0<=iz_  &&  iz_<nz-1 );
                if( !cond )    continue;

                float ptot = plx_0*pz_0 + plx_0*pz_1 + plx_1*pz_0 + plx_1*pz_1;

                atomicAdd(&(eimage_deMO[i_z0_lx0]), plx_0*pz_0 * stack_value/ptot);
                atomicAdd(&(eimage_deMO[i_z1_lx0]), plx_0*pz_1 * stack_value/ptot);
                atomicAdd(&(eimage_deMO[i_z0_lx1]), plx_1*pz_0 * stack_value/ptot);
                atomicAdd(&(eimage_deMO[i_z1_lx1]), plx_1*pz_1 * stack_value/ptot);
            }
        } // loop over theta
    } // loop over lx
    }
    __syncthreads();
return;
}

__global__ void GPU_DEV_Demoveout_Hicks(int add, int adj, float *eimage, float *eimage_deMO,  \
                                  int nx, float dx, int nz, float dz, int lx0, \
                                  int ixMin, int ixMax, int izMin, int izMax, \
                                  float contraction_factor, float lx_null) {

    int ix, iz, lx, iz_, lx_, il;
    int nlx=2*lx0+1;

    ix = tx + bx*bdx;
    iz = tz + bz*bdz;

    float thetaMin = -55.0 * M_PI/180.0;
    float thetaMax = +55.0 * M_PI/180.0;
    float dtheta   =   1.0 * M_PI/180.0;
    float theta;
    int   nl = 11;
    float L  = (nl-1)*dx;
    int   nl0 = (nl-1)/2;
    float dl = L/(nl-1);
    float l, Delta_lx, Delta_z;
    float stack_value, alpha_lx, alpha_z;
    float lx_pos, z_pos, lx_move, z_move, move;

    int diam = 8;
    int r = 4;
    int r_ = r-1;
    float coef[8*8];

    if(add == 0)
    {
        for(lx=-lx0; lx<=lx0; lx++) {
            int   i_th  = id3_eic_lx_transp(lx,iz,ix);
            eimage_deMO[i_th] = 0.0f;
        }
    }

    __syncthreads();

    // ODCIG contraction
    if(ixMin<=ix  &&  ix<=ixMax  &&  izMin<=iz  &&  iz<=izMax)
    {
    for(lx=-lx0; lx<=lx0; lx++)
    {
        // DeMoveout
        lx_pos = lx*dx; 
        z_pos  = iz*dz;
        for(theta=thetaMin; theta<=thetaMax; theta+=dtheta)
        {
            // Slant-Stack
            stack_value = 0.0f;
            for(il=0; il<nl; il++)
            {
                l = (il-nl0)*dl;
                Delta_lx = + l * cosf(theta);
                Delta_z  = - l * sinf(theta);

                lx_      = floorf((lx_pos+Delta_lx)/dx);
                alpha_lx = ((lx_pos+Delta_lx) - lx_*dx)/dx;
                iz_      = floorf((z_pos+Delta_z)/dz);
                alpha_z  = ((z_pos+Delta_z) - iz_*dz)/dz;

                DEV_interSinc2D_8P(coef, alpha_lx, alpha_z);

                int ih_lx, ih_z;
                for(ih_z=-r_; ih_z<=r; ih_z++) {
                    for(ih_lx=-r_; ih_lx<=r; ih_lx++) {
                        int icoef = (ih_lx+r_) + (ih_z+r_)*diam;
                        int jlx = (lx_+ih_lx);
                        int jz  = (iz_+ih_z );
                        if(lx0>=jlx && jlx<=lx0 && 0<=jz && jz <nz)
                        {
                            int j_z_lx = id3_eic_lx_transp(jlx, jz, ix);
                            stack_value += coef[icoef] * eimage[j_z_lx];                            
                        }
                    }
                }
            }

            // Slant-Spread
            if ( (abs(lx)*dx)>=lx_null )
                move = (1.0-contraction_factor) * lx_null;
            else
                move = (1.0-contraction_factor) * abs(lx)*dx;

            if(adj)    move *= -1.0;
            // else if(fabsf(theta)>0.0)
            //     move = (1.0-contraction_factor) * abs(lx)*dx / fabsf(sinf(theta));
        
            if(fsignf_GPU(theta)==fsignf_GPU(lx)) move *= -1.0;
            if(fsignf_GPU(theta)!=fsignf_GPU(lx)) move *= +1.0;
            lx_move = lx_pos + move * sinf(theta);
            z_move =   z_pos + move * cosf(theta);

            for(il=0; il<nl; il++)
            {
                l = (il-nl0)*dl;
                Delta_lx = + l * cosf(theta);
                Delta_z  = - l * sinf(theta);
                lx_      = floorf((lx_move+Delta_lx)/dx);
                alpha_lx = ((lx_move+Delta_lx) - lx_*dx)/dx;

                iz_      = floorf((z_move+Delta_z)/dz);
                alpha_z  = ((z_move+Delta_z) - iz_*dz)/dz;

                DEV_interSinc2D_8P(coef, alpha_lx, alpha_z);

                int ih_lx, ih_z;
                for(ih_z=-r_; ih_z<=r; ih_z++) {
                    for(ih_lx=-r_; ih_lx<=r; ih_lx++) {
                        int icoef = (ih_lx+r_) + (ih_z+r_)*diam;
                        int jlx = (lx_+ih_lx);
                        int jz  = (iz_+ih_z );
                        if(lx0>=jlx && jlx<=lx0 && 0<=jz && jz <nz)
                        {
                            int j_z_lx = id3_eic_lx_transp(jlx, jz, ix);
                            atomicAdd(&(eimage_deMO[j_z_lx]),stack_value * coef[icoef]);
                        }
                    }
                }

            }
        } // loop over theta
    } // loop over lx
    }
    __syncthreads();
return;
}

__global__ void GPU_DEV_Demoveout_Hicks_old(int add, float *eimage, float *eimage_deMO,  \
                                        int nx, float dx, int nz, float dz, int lx0, \
                                        int ixMin, int ixMax, int izMin, int izMax, \
                                        float contraction_factor, float lx_null) {

    int ix, iz, lx, iz_, lx_, il;
    int nlx=2*lx0+1;

    ix = tx + bx*bdx;
    iz = tz + bz*bdz;

    float  thetaMin = -55.0 * M_PI/180.0;
    float  thetaMax = +55.0 * M_PI/180.0;
    float  dtheta = 1.0 * M_PI/180.0;
    float  theta;
    int    nl = 11;
    float  L  = (nl-1)*dx;
    int    nl0 = (nl-1)/2;
    float  dl = L/(nl-1);
    float  l, Delta_lx, Delta_z;
    float  stack_value;
    float  lx_pos, z_pos, lx_move, z_move, move;
    int    diam = 8, r=diam/2, r_=r-1;
    float  coef[8*8];

    if(add == 0)
    {
        for(lx=-lx0; lx<=lx0; lx++) {
            int   i_th  = id3_eic_lx_transp(lx,iz,ix);
            eimage_deMO[i_th] = 0.0f;
        }
    }

    __syncthreads();

    // ODCIG contraction
    int  cond;
    if(ixMin<=ix  &&  ix<=ixMax  &&  izMin<=iz  &&  iz<=izMax)
    {
    for(lx=-lx0; lx<=lx0; lx++)
    {
        // DeMoveout
        lx_pos = lx*dx; 
        z_pos  = iz*dz;
        for(theta=thetaMin; theta<=thetaMax; theta+=dtheta)
        {
            // Slant-Stack
            stack_value = 0.0f;
            for(il=0; il<nl; il++)
            {
                l = (il-nl0)*dl;
                Delta_lx = + l * cosf(theta);
                Delta_z  = - l * sinf(theta);

                lx_   = floorf((lx_pos+Delta_lx)/dx);
                float alpha_x = ((lx_pos+Delta_lx) - lx_*dx)/dx;

                iz_   = floorf((z_pos+Delta_z)/dz);
                float alpha_z  = ((z_pos+Delta_z) - iz_*dz)/dz;
                
                DEV_interSinc2D_8P(coef, alpha_x, alpha_z);
                int iix, iiz;
                for(iix=-r_; iix<=r; iix++) {
                    for(iiz=-r_; iiz<=r; iiz++) {
                        int ii_ = (iix+r_) + (iiz+r_)*diam;
                        int i_z_lx = id3_eic_lx_transp(lx_+iix, iz_+iiz, ix);
                        cond = ((-lx0)<=lx_+iix && lx_+iix<=(lx0))  &&  (0<=iz_+iiz  &&  iz_+iiz<=nz-1);
                        if( !cond )    continue;
                        stack_value += coef[ii_] * eimage[i_z_lx];
                    }
                }
            }

            // Slant-Spread
            if ( (abs(lx)*dx)>=lx_null )
                move = (1.0-contraction_factor) * lx_null;
            else
                move = (1.0-contraction_factor) * abs(lx)*dx;
            // else if(fabsf(theta)>0.0)
            //     move = (1.0-contraction_factor) * abs(lx)*dx / fabsf(sinf(theta));
        
            if(fsignf_GPU(theta)==fsignf_GPU(lx)) move *= -1.0;
            if(fsignf_GPU(theta)!=fsignf_GPU(lx)) move *= +1.0;
            lx_move = lx_pos + move * sinf(theta);
            z_move  =  z_pos + move * cosf(theta);
            //  z_move =  z_pos + move;

            for(il=0; il<nl; il++)
            {
                l = (il-nl0)*dl;
                Delta_lx = + l * cosf(theta);
                Delta_z  = - l * sinf(theta);

                lx_      = floorf((lx_move+Delta_lx)/dx);
                float alpha_x  = ((lx_move+Delta_lx) - lx_*dx)/dx;

                iz_      = floorf((z_move+Delta_z)/dz);
                float alpha_z  = ((z_move+Delta_z) - iz_*dz)/dz;

                DEV_interSinc2D_8P(coef, alpha_x, alpha_z);
                int iix, iiz;
                for(iix=-r_; iix<=r; iix++) {
                    for(iiz=-r_; iiz<=r; iiz++) {
                        int ii_ = (iix+r_) + (iiz+r_)*diam;
                        int i_z_lx = id3_eic_lx_transp(lx_+iix, iz_+iiz, ix);
                        cond = ((-lx0)<=lx_+iix && lx_+iix<=(lx0))  &&  (0<=iz_+iiz  &&  iz_+iiz<=nz-1);
                        if( !cond )    continue;
                        atomicAdd(&(eimage_deMO[i_z_lx]), coef[ii_] * stack_value);
                    }
                }
            }

        } // loop over theta
    } // loop over lx
    }
    __syncthreads();
return;
}

__device__ float fsignf_GPU(float x) {
    return (x/fabsf(x));
}


// Hicks interpolation
__device__ void DEV_interSinc1D_8P(float* coef, float alpha) {
    int r = 4;
    int diam = 2*r;
    float b = 6.31f;
    float x, sum = 0.0f;
    int i;
    int r_ = r-1;
    
    for(i=0; i<diam; i++)    coef[i] = 0.0;

    for(i=0; i<diam; i++)
    {
        x = 1.0*(i-r_) - alpha;
        coef[i] = DEV_sinc(x) * DEV_kaiser(x,r,b);
        sum += coef[i];
    }
    for(i=0; i<diam; i++) {
        coef[i] /= sum;
    }
return;
}
__device__ void DEV_interSinc2D_8P(float *coef, float alpha_x, float alpha_y) {
    int r = 4;
    int diam = 2*r;
    float b = 6.31f;
    float x, y, sum = 0.0f;
    float coef1[8];
    float coef2[8];
    int i1, i2, i;
    int r_ = r-1;

    for(i=0; i<diam*diam; i++)    coef[i] = 0.0;

    for(i=0; i<diam; i++) {
        x = 1.0f*(i-r_) - alpha_x;
        coef1[i] = DEV_sinc(x) * DEV_kaiser(x,r,b);
    }
    for(i=0; i<diam; i++) {
        y = 1.0f*(i-r_) - alpha_y;
        coef2[i] = DEV_sinc(y) * DEV_kaiser(y,r,b);
    }
    for(i2=0; i2<diam; i2++) {
        for(i1=0; i1<diam; i1++) {
            i = i1 + i2*diam;
            coef[i] = coef1[i1]*coef2[i2];
            sum += coef[i];
        }
    }
    for(i=0; i<diam*diam; i++)    coef[i] /= sum;

return;
}
// Sinc function
__device__ float DEV_sinc(float x) {
    if(x==0.0f)
        return 1.0f;
    else
        return (sinf(M_PI*x)/(M_PI*x));
}
// Kaiser window
__device__ float DEV_kaiser(float x, float r, float b) {
    float window;
    if(fabsf(x)<=r)
        window = DEV_mBFZ(b*sqrtf(1.0 - (x/r)*(x/r))) / DEV_mBFZ(b);
    else 
        window = 0.0;
return window;
}
// Modified Bessel First Kind Zero Order
/*
__device__ float DEV_mBFZ(float x) {
    double add;
    double out  = 0.0;
    double samp = 1.0;
    double arg = 0.25f * x*x;
    long long int kfat  = 1;
    long long int kfat1 = 1;
    long long int k;
    k = 0;
    do {
        // out += ((samp/(kfat))/kfat1);
        add = samp/(kfat*kfat1);
        out += add;
        samp *= arg;
        k++;
        kfat  *= k;
        kfat1 *= (k+1);
    } while(add>0.000001*arg);

return ((float) out);
}
//*/
__device__ float DEV_mBFZ(float x) {
    float add;
    float out  = 0.0;
    float samp = 1.0;
    float arg = 0.25f * x*x;
    long long int kfat  = 1;
    long long int kfat1 = 1;
    long long int k;
    k = 0;
    do {
        // out += ((samp/(kfat))/kfat1);
        add = samp/(kfat*kfat1);
        out += add;
        samp *= arg;
        k++;
        kfat  *= k;
        kfat1 *= (k+1);
    } while(add>0.00001*arg);

return ((float) out);
}

//*/
    