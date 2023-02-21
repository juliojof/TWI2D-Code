#include "cufft.h"

#ifndef PROP_GPU_H
#define PROP_GPU_H

#define IC_REGULAR_IMAGE          0
#define IC_REGULAR_IMAGE2         1
#define IC_REGULAR_IMAGE3         2
#define IC_REGULAR_ILUMIN         3
#define IC_BORNSCT_IMAGE          4
#define IC_BORNSCT_IMAGE2         5
#define IC_BORNSCT_IMAGE3         6
#define IC_BORNSCT_ILUMIN         7
#define IC_BORNSCT_ILUMIN_INC     8
#define IC_BORNSCT_ILUMIN_SCT     9
#define IC_REGULAR_EXPREF         10
#define IC_REGULAR_EXPREF_SECDER  11

#define IC_SOU_WAVE             100
#define IC_REC_WAVE             101

#define MODE_CIGS_FULL          1
#define MODE_CIGS_ORTH          2
#define MODE_CIGS_TIME          3

#define min(a,b) ((a)<=(b)) ? (a) : (b)
#define max(a,b) ((a)>=(b)) ? (a) : (b)

void GPU_OPER_initAll(int idev, int nx, int nz, int nt, float dx, float dz, float dt, int nrec);

void GPU_OPER_initEnvironment(int idev, int nx, int nz, int nt, float dx, float dz, float dt, \
                              int nsou, int nrec, int lx0, int lz0, int lt0, int nxb, int nzb, int MODE_CIGS);

void GPU_OPER_initTAU(int tau0, float dtau);

void GPU_OPER_alocArrays_acoustic();
void GPU_OPER_freeArrays_acoustic();
void GPU_OPER_alocArrays_ExtExpRef(int jt);
void GPU_OPER_freeArrays_ExtExpRef();
void GPU_OPER_init_ExtExpRef_Docig(int jt, int lambda0, int ndip, float dipMin, float ddip);
void GPU_OPER_freeArrays_ExtExpRef_Docig();
void GPU_OPER_zeroArrays_acoustic();
void GPU_OPER_alocArrays_acoustic_BornSct();
void GPU_OPER_freeArrays_acoustic_BornSct();
void GPU_OPER_zeroArrays_BornSct();
void GPU_OPER_alocArrays_acoustic_ExtBornSct();
void GPU_OPER_freeArrays_acoustic_ExtBornSct();
void GPU_OPER_zeroArrays_ExtBornSct();
void GPU_OPER_setStore_inc_sct(int store_inc, int store_sct);
void GPU_OPER_alocArrays_EIC(int jt_EIC, int store_inc, int store_sct);
void GPU_OPER_freeArrays_EIC();
void GPU_OPER_alocArrays_EIC_recSecDer(int jt_EIC, int store_inc, int store_sct);
void GPU_OPER_freeArrays_EIC_recSecDer();

void GPU_OPER_getVelP_CPU2GPU(float *vp);
void GPU_OPER_getSigmas_CPU2GPU(float *sigma, float *sigmaInv);

void GPU_OPER_getSourceWavelet_CPU2GPU(float *source);
void GPU_OPER_getSourceWavelet_CPU2GPU_old(float *source);
void GPU_OPER_getRecordedData_GPU2CPU(float *seismogram);
void GPU_OPER_getRecordedData_CPU2GPU(float *seismogram);
void GPU_OPER_getSource_CPU2GPU(float *source);
void GPU_OPER_getSourceReceiverIndexes_CPU2GPU(int *ixs, int *izs, int *ixr, int *izr);

void GPU_OPER_zeroIluminArray();
void GPU_OPER_zeroWavefieldArrays();
void GPU_OPER_zeroScatteredWavefieldArrays();

float* GPU_WRAP_convolution(int idev, float *a, float *b, int na, int nb);
float* GPU_WRAP_SeismogramCoding(int idev, float *a, float *b, int ns, int na, int nb);
void   GPU_WRAP_SeismogramResamp(int idev, float *traces, float *traces_resamp, int ns, \
                                 int nt, float dt, int nt_resamp, float dt_resamp);

void GPU_OPER_copyTLCIG_GPU2CPU(float *TLcig);
void GPU_OPER_copyTLCIG_CPU2GPU(float *TLcig);
void GPU_OPER_copyTCIG_GPU2CPU(float *Tcig);
void GPU_OPER_copyTCIG_CPU2GPU(float *Tcig);
void GPU_OPER_copyOCIGS_GPU2CPU(float *Hocig, float *Vocig);
void GPU_OPER_copyOCIGS_CPU2GPU(float *Hocig, float *Vocig);
void GPU_OPER_copyResidualOCIGS_CPU2GPU(float *Hocig, float *Vocig);
void GPU_OPER_copyOCIGS_phaseShift_CPU2GPU(float *Hocig, float *Vocig);
void GPU_OPER_copyResidualOCIGS_phaseShift_CPU2GPU(float *Hocig, float *Vocig);

void GPU_OPER_copyDOCIG_ExtExpRef_CPU2GPU(float *Docig_bck_org, float *Docig_res_org, float *Docig_bck_psh, float *Docig_res_psh);

void GPU_OPER_copyEImage_GPU2CPU(float *eimage);
void GPU_OPER_copyEImage_CPU2GPU(float *eimage);
void GPU_OPER_copyImage_GPU2CPU(float *image);
void GPU_OPER_copyImage_CPU2GPU(float *image);
void GPU_OPER_copyRef_GPU2CPU(float *ref);
void GPU_OPER_copyRef_CPU2GPU(float *ref);
void GPU_OPER_copyIlumin_GPU2CPU(float *ilumin);
void GPU_OPER_copyImage_BRNSCT_GPU2CPU(float *imageSouSide, float *imageRecSide, int add);
void GPU_OPER_copyIlumin_BRNSCT_GPU2CPU(float *iluminSouSide, float *iluminRecSide, int add);

void GPU_OPER_copyPresentWavefieldWS_CPU2GPU(float *wavefield);
void GPU_OPER_copyPresentWavefieldWC_CPU2GPU(float *wavefield);
void GPU_OPER_copyPresentWavefieldWS_GPU2CPU(float *wavefield);
void GPU_OPER_copyPresentWavefieldWC_GPU2CPU(float *wavefield);

void GPU_OPER_copyPresentWavefield_CPU2GPU(float *wavefield);
void GPU_OPER_copyPresentWavefield_GPU2CPU(float *wavefield);
void GPU_OPER_copyPresentScatteredWavefield_GPU2CPU(float *wavefield);
void GPU_OPER_copyPresentWavefield_secDer_GPU2CPU(float *wavefield);
void GPU_OPER_copyPresentScatteredWavefield_secDer_GPU2CPU(float *wavefield);

void GPU_WRAP_imagingCondition_transmission();
void GPU_OPER_EIC_storeWavefields_GPU2GPU(int it, int side, int mode);

void GPU_OPER_copyLaplacian_GPU2CPU(float *lap);
void GPU_OPER_copyPresentWavefieldWS_CPU2GPU_opt(float *wavefield);
void GPU_OPER_copyPresentWavefield_GPU2CPU_opt(float *wavefield);

void GPU_OPER_alocMem(void** dev_a, size_t size);

void GPU_OPER_setTimeStep(int it);

// Born scattering functions
void GPU_OPER_zeroSecDerIncArray();
void GPU_OPER_zeroSecDerSctArray();
void GPU_OPER_wavefieldTimeSecondDerivative_part1(int mode);
void GPU_OPER_wavefieldTimeSecondDerivative_part2(int mode);
void GPU_WRAP_scatter_Acoustic_OCIGS_opt2(int useVocig);
void GPU_WRAP_scatter_Acoustic_OCIGS_opt(int useVocig);
void GPU_WRAP_scatter_Acoustic();
void GPU_WRAP_scatter_Acoustic_lx_opt();
void GPU_WRAP_scatter_Acoustic_lz_opt();
void GPU_WRAP_scatter_Acoustic_lx_lz_opt();

void GPU_WRAP_Contraction_theta(int gpuDev, int add, float *adcigInp, float *adcigOut, \
                                int nx, float dx, int nz, float dz, int ltheta0, float dtheta, int nxb, int nzb, \
                                float xMin, float xMax, float zMin, float zMax, float stretching_factor, float thetaMax);
void GPU_WRAP_Contraction_lx(int gpuDev, int add, float *eimageInp, float *eimageOut, \
                             int nx, float dx, int nz, float dz, int lx0, int nxb, int nzb, \
                             float xMin, float xMax, float zMin, float zMax, float contraction_factor);
void GPU_WRAP_Contraction_lx_lz(int gpuDev, int add, float *eimageInp, float *eimageOut, \
                                int nx, float dx, int nz, float dz, int lx0, int lz0, \
                                float dipMin, float dipMax, float ddip, float lambdaMax, float dlambda, \
                                int nxb, int nzb, float xMin, float xMax, float zMin, float zMax, float contraction_factor);
void GPU_WRAP_adcig_theta2tanTheta(int adj, int gpuDev, int add, float *adcigTheta, float *adcigTanth, \
                                   int nx, float dx, int nz, float dz, \
                                   int ltheta0, float dtheta, int ltanth0, float dtanth, \
                                   int nxb, int nzb, float xMin, float xMax, float zMin, float zMax);
void GPU_WRAP_Demoveout_AD(int gpuDev, int add, int adj, float *adcigInp, float *adcigOut, \
                           int nx, float dx, int nz, float dz, int lx0, int ltheta0, float dtheta, \
                           int nxb, int nzb, float xMin, float xMax, float zMin, float zMax, \
                           float contraction_factor, float lx_null);
void GPU_WRAP_Demoveout(int gpuDev, int add, int adj, float *eimageInp, float *eimageOut, \
                        int nx, float dx, int nz, float dz, int lx0, \
                        int nxb, int nzb, \
                        float xMin, float xMax, float zMin, float zMax, \
                        float contraction_factor, float lx_null);
void GPU_WRAP_Demoveout_Tcig(int gpuDev, int add, int adj, float *eimageInp, float *eimageOut, \
                             int nx, float dx, int nz, float dz, int tau0, \
                             int nxb, int nzb, float xMin, float xMax, float zMin, float zMax, \
                             float contraction_factor, float tau_null);
void GPU_WRAP_Demoveout_Docig(int gpuDev, int add, int adj, float *eimageInp, float *eimageOut, \
                              long long int nx, float dx, long long int nz, float dz, \
                              long long int nxb, long long int nzb, \
                              long long int ld0, float dip, \
                              float xMin, float xMax, float zMin, float zMax, \
                              float contraction_factor, float lx_null);
void GPU_WRAP_Demoveout_Gocig(int gpuDev, int add, int adj, float *eimageInp, float *eimageOut, \
                              long long int nx, float dx, long long int nz, float dz, long long int nxb, long long int nzb, \
                              long long int lx0, long long int lz0, float lambdaMax, float dipMin, float dipMax, float ddip, \
                              float xMin, float xMax, float zMin, float zMax, \
                              float contraction_factor, float lx_null);
void GPU_WRAP_Demoveout_lx_lz(int gpuDev, int add, int adj, float *eimageInp, float *eimageOut, \
                              int nx, float dx, int nz, float dz, int lx0, int lz0, \
                              float dipMin, float dipMax, float ddip, float lambdaMax, float dlambda, \
                              int nxb, int nzb, float xMin, float xMax, float zMin, float zMax, \
                              float contraction_factor, float lx_null);

void GPU_WRAP_lxlzlt2LambdaLtDip(int adj, int gpuDev, int add, float *eimage_lxlzlt, \
                                 float *eimage_LambdaLtDip, int nx, float dx, int nz, \
                                 float dz, int nxb, int nzb, int lx0, int lz0, int lt0, \
                                 float lambdaMax, float dlambda, int nlambda, \
                                 float dipMin, float dipMax, float ddip, int ndip);

void GPU_WRAP_lxlz2LambdaDip_old(int adj, int gpuDev, int add, float *eimage_lxlz, \
                             float *eimage_LambdaDip, int nx, float dx, int nz, \
                             float dz, int nxb, int nzb, int nx_, int nz_, int lx0, int lz0, \
                             float lambdaMax, float dlambda, int nlambda, \
                             float dipMin, float dipMax, float ddip, int ndip);

void GPU_WRAP_lxlz2LambdaDip(int adj, int gpuDev, int add, float *eimage_lxlz, \
                             float *eimage_LambdaDip, long long int nx, float dx, long long int nz, \
                             float dz, long long int nxb, long long int nzb, \
                             long long int nx_, float ox_, long long int nz_, float oz_, long long int lx0, long long int lz0, \
                             float lambdaMax, float dlambda, long long int nlambda, \
                             float dipMin, float dipMax, float ddip, long long int ndip);

void GPU_WRAP_lx2theta(int adj, int gpuDev, int add, float *eimage, \
                       float *eimage_theta, int nx, float dx, int nz, \
                       float dz, int nxb, int nzb, int lx0, int ntheta, \
                       float dtheta, float otheta);

void GPU_WRAP_lz2theta(int adj, int gpuDev, int add, float *eimage, \
                       float *eimage_theta, int nx, float dx, int nz, \
                       float dz, int nxb, int nzb, int lz0, int ntheta, \
                       float dtheta, float otheta);

//ODCIG Phase-Shift
void GPU_WRAP_HocigPhaseShift(int idev, float *Hocig_out, float *Hocig_inp, float *vel, \
                              long long int nx, long long int nz, long long int nt, \
                              long long int nxb, long long int nzb, \
                              float dx, float dz, float dt, int lx0, long long int ntstep);

void GPU_WRAP_DocigPhaseShift(int idev, float *Docig_out, float *Docig_inp, float *vel, float dip, \
                              long long int nx, long long int nz, long long int nt, \
                              long long int nxb, long long int nzb, \
                              float dx, float dz, float dt, int lx0, long long int ntstep);

void GPU_WRAP_HocigPhaseShift_kernel(int idev, cuComplex *Hocig_inp, cuComplex *Hocig_psh, float *vel, \
                                    long long int lx0, long long int nt, long long int nz, long long int nx, \
                                    long long int nzb, float dx, float dt, float dz, long long int ix);

void GPU_WRAP_DocigPhaseShift_kernel(int idev, cuComplex *Hocig_inp, cuComplex *Hocig_psh, float *vel, float dip, \
                                    long long int lx0, long long int nt, long long int nz, long long int nx, \
                                    long long int nzb, long long int nxb, float dx, float dt, float dz, \
                                    long long int ix, long long int ntstep);


void GPU_WRAP_Hocig_t_lx_2_w_klx(int idev, cuComplex *Hocig_w_klx_z, float *Hocig, long long int nlx, \
                                 long long int nt, long long int nz, int direction);

void GPU_WRAP_fixOriginOfHocigFFT(int idev, cufftComplex *Hocig, long long int lx0, \
                                  long int nt, long long int nz, \
                                  long long int nzb, float dx, float dt);

// Transformations VHocig-Docig
void GPU_WRAP_VHocig2Docig(int idev, float *Docig_out, float *Hocig_inp, float *Vocig_inp, \
                           long long int nx, long long int nz, \
                           float dx, float dz, long long int ndip, float ddip, float dipMin, \
                           long long int lx0, long long int lv0, long long int ld0, int iproc, int useVocig);
void GPU_WRAP_Docig2VHocig(int idev, float *Docig_inp, float *Hocig_out, float *Vocig_out, \
                           long long int nx, long long int nz, \
                           float dx, float dz, long long int ndip, float ddip, float dipMin, \
                           long long int lx0, long long int lv0, long long int ld0, int iproc, int useVocig);
void GPU_WRAP_VHocig2Docig_kernel(int idev, int adj, cuComplex *Docig_fft, cuComplex *Hocig_fft, \
                                  long long int lx0, long long int ld0, \
                                  long long int nx, long long int nz, float dx, float dz, \
                                  long long int ndip, float ddip, float dipMin, int mode);
void GPU_WRAP_VHocig_xz2KxKz(int idev, cuComplex *a, float *b, long long int nx, long long int nz, long long int lp0, int direction);
void GPU_WRAP_Docig_xz2KxKz(int idev, cuComplex *complex, float *real, long long int nx, long long int nz, long long int lp0, long long int ndip, int direction);

// Misc
// void GPU_WRAP_transp_dim3(int idev, float *mou, float *min, long long int n1, long long int n2, long long int n3, int plane);
// void GPU_WRAP_applycuFFT1D_inverseC2R(int idev, float *out, cufftComplex *in, long long int n1, long long int batch);
// void GPU_WRAP_applycuFFT1D_forwardR2C(int idev, cufftComplex *ou, float *inp, long long int n1, long long int batch);
// void GPU_WRAP_applycuFFT2D_forwardR2C(int idev, cufftComplex *ou, float *inp, long long int n1, long long int n2, long long int batch);
// void GPU_WRAP_applycuFFT2D_inverseC2R(int idev, float *out, cufftComplex *in, long long int n1, long long int n2, long long int batch);
// void GPU_WRAP_applyFFTW2D_inverseC2R(int idev, float *out, cufftComplex *in, long long int n1, long long int n2);
// void GPU_WRAP_applyFFTW2D_forwardR2C(int idev, cufftComplex *ou, float *inp, long long int n1, long long int n2);
// void GPU_WRAP_complex2real(float *ou, cufftComplex *in, int n, float coef);

// Imaging conditions
void GPU_WRAP_imagingCondition(int mode);
void GPU_WRAP_imagingCondition_smart(int mode, int it);
void GPU_WRAP_imagingCondition_ExtExpRef(int mode);
void GPU_WRAP_eimagingCondition_OCIGS_smart(int useVocig, int it);
void GPU_WRAP_eimagingCondition_OCIGS(int useVocig);
void GPU_WRAP_eimagingCondition_lx();
void GPU_WRAP_eimagingCondition_lt();
void GPU_WRAP_scatteringCondition_TAU(int side);
void GPU_WRAP_eimagingCondition_TAU();
void GPU_WRAP_eimagingCondition_lx_lz();
void GPU_WRAP_eimagingCondition_lx_lz_lt();

// Demoveout Operator

// >>>>>>>>>>>>>>>>>>>>>>>  <<<<<<<<<<<<<<<<<<<<<<<<
// >>>>>>>> Wavefield propagation functions <<<<<<<<
// >>>>>>>>>>>>>>>>>>>>>>>  <<<<<<<<<<<<<<<<<<<<<<<<
void GPU_WRAP_injectBackgroundExtRef_TimeFunction_Acoustic(int useVocig);
void GPU_WRAP_injectResidualExtRef_TimeFunction_Acoustic(int useVocig);

int  GPU_WRAP_injectBackgroundExtRef_Acoustic(int useVocig);
void GPU_WRAP_injectResidualExtRef_Acoustic(int useVocig);
int  GPU_WRAP_injectExtRef_Acoustic(int useVocig);
int  GPU_WRAP_injectExtRef_Acoustic_Docig();

void GPU_WRAP_recordInjectData(int mode);

void GPU_WRAP_injectSource_Acoustic();
void GPU_WRAP_injectSource_Acoustic_old();

void GPU_WRAP_propagate_Acoustic_abs2();
void GPU_WRAP_propagate_Acoustic_Scattered_abs2();


// Envelope, Hilbert
void GPU_WRAP_envelope(int idev, float *out, float *traces, long long int nt, float dt, long long int ntraces);
void GPU_WRAP_hilbert(int idev, float *traces_out, float *traces, int ntraces, int nt, float dt);

// Smoothing operator
void GPU_WRAP_smooth2D(int idev, float *out, float *inp, int n1, int n2, int r1, int r2);




// 2D wavenumber bandpass
void GPU_WRAP_wavefield_bandpass_kx_kz_iso_GPU2GPU(int idev, long long int nx, long long int nz, float dx, float dz, \
                                                   float lambda1, float lambda2, float lambda3, float lambda4);

void GPU_WRAP_wavefield_bandpass_kx_kz_GPU2GPU(int idev, long long int nx, long long int nz, float dx, float dz, \
                                               float lambda1, float lambda2, float lambda3, float lambda4);

void GPU_WRAP_bandpass_kx_kz_GPU2GPU(int idev, float *out, float *inp, long long int nx, long long int nz, float dx, float dz, \
                                     float kx1, float kx2, float kx3, float kx4, float kz1, float kz2, float kz3, float kz4);

void GPU_WRAP_2D_fourier(int idev, cuComplex *array_kx_kz, float *array_x_z, long long int nx, long long int nz, long long int batch, int direction);

void GPU_WRAP_bandpass_kx_kz_GPU2GPU_kernel(int idev, cuComplex *array_kx_kz, long long int nx, long long int nz, float dx, float dz, \
                                            float kx1, float kx2, float kx3, float kx4, float kz1, float kz2, float kz3, float kz4);

void GPU_WRAP_bandpass_kx_kz_iso_GPU2GPU_kernel(int idev, cuComplex *array_kx_kz, long long int nx, long long int nz, float dx, float dz, \
                                                float k1, float k2, float k3, float k4);


#endif