#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>
#include "omp.h"
#include <time.h>
#include <mpi.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "PropagationFunctions_GPU.h"
#include "cufft.h"
#include "shot_3D.h"

int EIMSH_CIGS_ORTH;
int EIMSH_CIGS_VCIG;

int    nxb, nzb;
float  dxPSF;
float  dzPSF;
int    jxPSF;
int    jzPSF;

int gpuDev = 1;
int ngpu = 4;
int procGPU;
int firstGPU = 4;

size_t AllocatedMem_CPU;

#include "DataTypes.h"

#include "Play.h"
#include "Solver_TWI2D.h"
#include "ObjFuncExplorer.h"
#include "utl.h"

/////////////// Main Function //////////////
int main(int argc, char **argv) {


// Initialize the MPI environment
MPI_Init(NULL, NULL);

// Get the number of processes
int nproc;
MPI_Comm_size(MPI_COMM_WORLD, &nproc);

// Get the rank of the process
int iproc;
MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

procGPU = iproc%ngpu + firstGPU;


MPI_Barrier(MPI_COMM_WORLD);
AllocatedMem_CPU = 0;

int ret;

// Reading inputFileName from command line
memset(inputFileName,  0, sizeof inputFileName);
ret = readCommandLineArg( (void *) inputFileName, argc, argv, "parmsFile", "mandatory", TYPE_STRING);
if(ret==1) { printf("\n\n  ERROR: Execution parameters [ parmsFile= ] was not informed upon program call. Program will stop.  \n\n"); exit(-1); }
// ret = readCommandLineArg( (void *) &ngpu, argc, argv, "ngpu", "mandatory", TYPE_INT);
// if(ret==1) { printf("\n\n  ERROR: Execution parameters [ ngpu= ] was not informed upon program call. Program will stop.  \n\n"); exit(-1); }
// ret = readCommandLineArg( (void *) &firstGPU, argc, argv, "firstGPU", "mandatory", TYPE_INT);
// if(ret==1) { printf("\n\n  ERROR: Execution parameters [ firstGPU= ] was not informed upon program call. Program will stop.  \n\n"); exit(-1); }

// ********* Marmousi successfull experiment **************
Parms parms;
ret = readParms(&parms);
if(ret != 0) {
    printf("\n\n  **********  ATTENTION:  ********** \n\n      Error while reading input parameters! Program will abort.    \n\n\n");
    exit(-1);
}

if(iproc==0)    printf("\n\n  inputFileName=%s  \n\n", inputFileName);

// Output directories
initExecPathName(parms.executionDirectory);


nxb = 50;
nzb = 50;
int nx, nz, nx_, nz_, nt, jt, jt_EIC;
float  LX=parms.LX, LZ=parms.LZ, LT=parms.LT;
float  dx=parms.dx, dz=parms.dz, dt=parms.dt;
float  ox=parms.ox, oz=parms.oz, ot;


nx_= floor(0.5 + LX/dx) + 2*nxb;
nz_= floor(0.5 + LZ/dz) + 2*nzb;
nx = floor(0.5 + ((1.0*nx_)/32.0)) * 32 - 2*nxb;
nz = floor(0.5 + ((1.0*nz_)/ 8.0)) *  8 - 2*nzb;
nt = floor(0.5 + LT/dt) + 1;
jt     = parms.jt;
jt_EIC = parms.jt_EIC;

LX = (nx-1)*dx;
LZ = (nz-1)*dz;

time_t time_now, time_start, time_begin, time_elapsed, time_begin_iter;
int    ix, iz, it;

float perc = 0.0;
float offMax = parms.offMax;
int   acquisitionGeom = 0;
float dxShot = parms.dxShot;

nx += 2*nxb;
nz += 2*nzb;
ox -= nxb*dx;
oz -= nzb*dz;
ot  = 0.0;

MSH msh;
init_msh(&msh, nx, nz, nt, dx, dz, dt, ox, oz, ot, jt, jt_EIC, nxb, nzb);

Eimage_MSH eimsh;
init_eimsh(&eimsh, &parms, &msh);

EIProc eiproc;
init_eiproc(&eiproc, &parms);

InvParms invparms;
init_InvParms(&invparms, &parms);

GRDProc grdProc;
init_grdProc(&grdProc, &parms);

float agc_time_window = parms.agc_time_window;

// Arrays
float  *vp;

// Focusing parameter
int outputs;
if(iproc==0)    outputs=1;
else            outputs=0;
float muFocus = parms.muFocus, epsFocus = parms.epsFocus, lambdaNull = parms.lambdaNull;

time_begin = time(NULL);
if(iproc==0) {
    printf("\n Total elapsed time: %3ld \n", time(NULL)-time_begin);
    printf("\n\n  parms.filePathVel = %s "     , parms.filePathVel);
}


// Reading velocity model
vp = readInputVel(parms.filePathVel, "vel_init", parms.nxInput_vel, parms.nzInput_vel, \
                  parms.deltaXInput_vel, parms.deltaZInput_vel, parms.OrigXInput_vel, parms.OrigZInput_vel, \
                  nx, nz, dx, dz, dt, ox, oz, nxb, nzb, parms.velScalar);


if(iproc==0) { printf("\n\n Done all readings \n\n"); fflush(stdout); }

// Initializing dataset
MPI_Barrier(MPI_COMM_WORLD);    
if(iproc==0)  { printf("\n\n Initializing datasets \n\n"); fflush(stdout); }
int    nshot=floor(0.5+LX/dxShot)+1, ishot;
float  dx_rec = dx;
Dataset dataset;
initDataset(&dataset, nshot, dx_rec, nt, nx, nz, dt, dx, dz, offMax, parms.f1, parms.f2, parms.f3, parms.f4);
for(ishot=0; ishot<nshot; ishot++) {
    readShot(parms.filePathInputData, parms.filePathInputSource, &((dataset.shot)[ishot]));
}
MPI_Barrier(MPI_COMM_WORLD);

// Saving data to file
/*
FILE *fSeis_Hdr, *fSeis_Bin;
if(iproc==0)
{
    initFilesSmart3d("seismogram_input", &fSeis_Hdr, &fSeis_Bin, msh.nt, msh.dt, msh.ot, (dataset.shot)[0].nrec, 1, 0, dataset.nshot/1, 1, 0);
    for(ishot=0; ishot<nshot; ishot+=1) {
        writeSamples((dataset.shot)[ishot].seismogram, fSeis_Bin, msh.nt*(dataset.shot)[0].nrec);
    }
    fclose(fSeis_Bin);
    fclose(fSeis_Hdr);
}
//*/
    MPI_Barrier(MPI_COMM_WORLD);

// Save source to file
if(iproc==0) outputSmart1d("source", (dataset.shot)[0].source, nt, dt, ot);

MPI_Barrier(MPI_COMM_WORLD);

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>> Solver TWI - Begin <<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
int RUN_MODE = parms.RUN_MODE;
if(RUN_MODE == RUN_TWI_REG) 
{
    Dataset dataset_residuals;
    initDataset(&dataset_residuals, nshot, dx_rec, nt, nx, nz, dt, dx, dz, offMax, parms.f1, parms.f2, parms.f3, parms.f4);
    if(iproc==0) { printf("\n\n Starting Solver_TWI2D. \n\n"); fflush(stdout); }
    Solver_TWI2D(&dataset, &dataset_residuals, vp, &invparms, &msh, &eimsh, &eiproc, &grdProc, \
                 acquisitionGeom, offMax, perc, parms.agc_time_window, iproc, nproc);
}
else if(RUN_MODE == RUN_TWI_OBJFUNCEXP) 
{
    int lref, lmis;
    float dref, dmis;
    dref = parms.OF_vel_ref_step;
    dmis = parms.OF_vel_mis_step;
    lref = parms.OF_vel_ref_mag/dref;
    lmis = parms.OF_vel_mis_mag/dmis;
    

    if(iproc==0) { printf("\n\n Starting ObjFuncExplorer_TWI2D. \n\n"); fflush(stdout); }
    Dataset dataset_ref, dataset_mis, dataset_dif;
    initDataset(&dataset_ref, nshot, dx_rec, nt, nx, nz, dt, dx, dz, offMax, parms.f1, parms.f2, parms.f3, parms.f4);
    initDataset(&dataset_mis, nshot, dx_rec, nt, nx, nz, dt, dx, dz, offMax, parms.f1, parms.f2, parms.f3, parms.f4);
    initDataset(&dataset_dif, nshot, dx_rec, nt, nx, nz, dt, dx, dz, offMax, parms.f1, parms.f2, parms.f3, parms.f4);
    ObjFuncExplorer_TWI2D(&dataset, &dataset_ref, &dataset_mis, &dataset_dif, vp, \
                          &msh, &eimsh, &eiproc, lref, lmis, dref, dmis, \
                          acquisitionGeom, offMax, perc, parms.agc_time_window, iproc, nproc);
}

MPI_Barrier(MPI_COMM_WORLD);
if(iproc==0) {
    printf("\n Finished Solver_TWI2D: Total elapsed time: %3ld \n", time(NULL)-time_begin);
}
MPI_Barrier(MPI_COMM_WORLD);
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>> Solver TWI - End <<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

free(vp);  

// Finalizing datasets
finalizeDataset(&dataset);

if(iproc==0) {
    printf("\n\n The program is being concluded normally \n\n");
    printf("\n Total elapsed time: %3ld \n", time(NULL)-time_begin);
    fflush(stdout);
}

// Finalize the MPI environment.
MPI_Finalize();

return 0;    
}

#include "shot_3D.c"
#include "Solver_TWI2D.c"
#include "ObjFuncExplorer.c"
#include "utl.c"
#include "ModMig.c"
#include "FocusingOperator.c"
#include "shot.c"
#include "ForwardAdjoint.c"
#include "PropagationFunctions.c"
#include "ForwardAdjoint_Elastic.c"
#include "PropagationFunctions_Elastic.c"
#include "Model.c"
#include "IO.c"
#include "ConjugateGradient_PSF.c"
#include "Filter.c"
#include "Math.c"
#include "CG_theta2lx.c"
#include "CG_theta2lx_GPU.c"
#include "CG_Docig2VHocig.c"
#include "CG_tanTheta2theta.c"
#include "ParmRead.c"
#include "CIGProcessing.c"
// #include "Residuals.c"


float *pullOffsetPlane(float *eimage, int nx, int nz, int lx0, int lx) {
    
    int ix, iz;
    int nlx = 2*lx0+1;

    float *oplane = CPU_zaloc1F(nx*nz);

    for(ix=0; ix<nx; ix++)
        for(iz=0; iz<nz; iz++)
        {
            int i_oplane = id2(iz,ix);
            int i_eimage = id3_eic_lx(lx,iz,ix);
            oplane[i_oplane] = eimage[i_eimage];
        }

return oplane;
}
void pushOffsetPlane(float *eimage, float *plane, int nx, int nz, int lx0, int lx) {
    
    int ix, iz;
    int nlx = 2*lx0+1;

    for(ix=0; ix<nx; ix++)
        for(iz=0; iz<nz; iz++)
        {
            int i_oplane = id2(iz,ix);
            int i_eimage = id3_eic_lx(lx,iz,ix);
            eimage[i_eimage] = plane[i_oplane];
        }

return;
}
float *pullThetaPlane(float *eimage_theta, int nx, int nz, int ltheta0, int itheta) {
    
    int ix, iz;
    int ntheta = 2*ltheta0+1;

    float *thetaPlane = CPU_zaloc1F(nx*nz);

    for(ix=0; ix<nx; ix++)
        for(iz=0; iz<nz; iz++)
        {
            int i_plane        = id2(iz,ix);
            int i_eimage_theta = id3_eic_th(itheta,iz,ix);
            thetaPlane[i_plane] = eimage_theta[i_eimage_theta];
        }

return thetaPlane;
}
void pushThetaPlane(float *eimage_theta, float *thetaPlane, int nx, int nz, int ltheta0, int itheta) {
    
    int ix, iz;
    int ntheta = 2*ltheta0+1;

    for(ix=0; ix<nx; ix++)
        for(iz=0; iz<nz; iz++)
        {
            int i_plane        = id2(iz,ix);
            int i_eimage_theta = id3_eic_th(itheta,iz,ix);
            eimage_theta[i_eimage_theta] = thetaPlane[i_plane];
        }

return;
}
void addWhiteNoise2EImagePSF(float *eimage_psf, int nx, int nz, int lx0, float epsilon) {
    
    int   ix, iz;
    int   lx;
    int   nlx = 2*lx0+1;
    float *rms = CPU_zaloc1F(nlx);
    float avgRMS_pos = 0.0f;
    float avgRMS_neg = 0.0f;

    for(lx=-lx0; lx<=lx0; lx++)
    {
        for(ix=0; ix<nx; ix++)
            for(iz=0; iz<nz; iz++)            
            {
                int i_eimage = id3_eic_lx(lx,iz,ix);
                rms[lx+lx0] += (eimage_psf[i_eimage]*eimage_psf[i_eimage])/(nx*nz);
            }
        rms[lx+lx0] = sqrtf(rms[lx+lx0]);
    }
    for(lx=-lx0; lx<=lx0; lx++) {
        if(lx<=-lx0/2)    avgRMS_neg += rms[lx+lx0];
        if(lx>=+lx0/2)    avgRMS_pos += rms[lx+lx0];
    }
    for(lx=-lx0; lx<=lx0; lx++)
        for(ix=0; ix<nx; ix++)
            for(iz=0; iz<nz; iz++)            
            {
                int i_eimage = id3_eic_lx(lx,iz,ix);
                if(lx<=-lx0/2)    eimage_psf[i_eimage] += epsilon*avgRMS_neg;
                if(lx>=+lx0/2)    eimage_psf[i_eimage] += epsilon*avgRMS_pos;
            }
return;
}

void accumulate(float *out, const float *inp, size_t n) {
    size_t i;
    for(i=0; i<n; i++)    
        out[i] += inp[i];
return;
}






void filterECIG(float *eimage, EIProc *eiproc, Eimage_MSH *emsh, MSH *msh)
{ 
    if(eiproc->applyBandass2EImage)
        lxApplyBandpass(eimage, msh->nx, msh->nz, msh->dz, emsh->lx0, eiproc->lambda1, eiproc->lambda2, eiproc->lambda3, eiproc->lambda4);
    MPI_Barrier(MPI_COMM_WORLD);

    if(eiproc->applyXTaper)
        lxTapering_x(eimage, eimage, msh->nx, msh->nz, msh->dx, emsh->lx0, eiproc->x1, eiproc->x2, eiproc->x3, eiproc->x4);
    MPI_Barrier(MPI_COMM_WORLD);

    if(eiproc->applyZTaper) {
        lxTapering(  eimage, eimage, msh->nx, msh->nz, msh->dz, emsh->lx0, emsh->lz0, eiproc->z0, eiproc->z1);
        lxTapering(  eimage, eimage, msh->nx, msh->nz, msh->dz, emsh->lx0, emsh->lz0, eiproc->z3, eiproc->z2);
    }
    // outputSmart3d("Hocig_rightAfter_ZTaper", eimage, 2*emsh->lx0+1, msh->dx, -emsh->lx0*msh->dx,  msh->nz, msh->dz, msh->oz,  msh->nx, msh->dx, msh->ox);
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(eiproc->applyLXTaper)
        lxTapering_lx_varDepth(eimage, eimage, msh->nx, msh->nz, msh->dx, emsh->lx0, \
                               eiproc->lx1_i, eiproc->lx2_i, eiproc->lx3_i, eiproc->lx4_i, \
                               eiproc->lx1_f, eiproc->lx2_f, eiproc->lx3_f, eiproc->lx4_f, \
                               eiproc->z_i, eiproc->z_f);
    MPI_Barrier(MPI_COMM_WORLD);

return;
}

void normalizeHocigByEnvelope(float *Hocig, float *Hocig_normalizer, long long int lx0, long long int nx, long long int nz, float dz)
{
    // Compute to envelope
    long long int nlx = 2*lx0+1;
    float *Hocig_transp            = transp_dim3(Hocig           , nlx, nz, nx, 12);
    float *Hocig_normalizer_transp = transp_dim3(Hocig_normalizer, nlx, nz, nx, 12);
    // float *Hocig_envelope = transform2Envelope(Hocig_transp, nz, dz, nx*nlx);
    float *Hocig_envelope = CPU_zaloc1F(nx*nz*nlx);
    // GPU_WRAP_envelope(procGPU, Hocig_envelope, Hocig_normalizer_transp, nz, dz, nx*nlx);
    takeAbsluteValueArrays(Hocig_envelope, Hocig_normalizer_transp, nlx*nx*nz);
    
    // Smooth and normalize envelope
    int ix;
    for(ix=0; ix<nx; ix++) {
        long long int shiftMem = ix*nlx*nz;
        int rz  = 250.0/dz;
        int rlx = lx0;
        // modelSmooth(Hocig_envelope+shiftMem, Hocig_envelope+shiftMem, nlx, nz, rlx, rz);
        GPU_WRAP_smooth2D(procGPU, Hocig_envelope+shiftMem, Hocig_envelope+shiftMem, nz, nlx, rz, rlx);
    }
    normalizeArray(Hocig_envelope, nx*nz*nlx);
    outputSmart3d("Hocig_envelope", Hocig_envelope, nz, dz, 0.0, nlx, dz, -lx0*dz, nx, dz, 0.0);

    // Divide input Hocig by envelope
    float noiseLevel = 0.1;
    divideArraysWhiteNoise(Hocig_transp, Hocig_transp, Hocig_envelope, nx*nz*nlx, noiseLevel);
    float *Hocig_normalized = transp_dim3(Hocig_transp, nz, nlx, nx, 12);

    outputSmart3d("Hocig_normalized", Hocig_normalized, nlx, dz, -lx0*dz, nz, dz, 0.0, nx, dz, 0.0);
    outputSmart3d("Hocig_normalizer", Hocig_normalizer, nlx, dz, -lx0*dz, nz, dz, 0.0, nx, dz, 0.0);

    // Free memory
    free(Hocig_transp);
    free(Hocig_normalizer_transp);
    free(Hocig_envelope);

    memcpy(Hocig, Hocig_normalized, nx*nz*nlx*sizeof(float));

    free(Hocig_normalized);

return;
}
void normalizeVocigByEnvelope(float *Vocig, float *Vocig_normalizer, long long int lv0, long long int nx, long long int nz, float dx)
{
    // Compute envelope
    long long int nlv = 2*lv0+1;
    float *Vocig_transp   = transp_dim3(Vocig_normalizer, nlv, nz, nx, 12);
    float *Vocig_transp2  = transp_dim3(Vocig_transp    , nz, nlv, nx, 13);
    // float *Vocig_envelope = transform2Envelope(Vocig_transp2, nx, dx, nz*nlv);
    float *Vocig_envelope = CPU_zaloc1F(nx*nz*nlv);
    // GPU_WRAP_envelope(procGPU, Vocig_envelope, Vocig_transp2, nx, dx, nz*nlv);
    takeAbsluteValueArrays(Vocig_envelope, Vocig_transp2, nlv*nx*nz);

    // Smooth and normalize envelope
    int iz;
    for(iz=0; iz<nz; iz++) {
        long long int shiftMem = iz*nlv*nx;
        int rx  = 250.0/dx;
        int rlv = lv0;
        // modelSmooth(Vocig_envelope+shiftMem, Vocig_envelope+shiftMem, nlv, nx, rlv, rx);
        GPU_WRAP_smooth2D(procGPU, Vocig_envelope+shiftMem, Vocig_envelope+shiftMem, nx, nlv, rx, rlv);
    }
    normalizeArray(Vocig_envelope, nx*nz*nlv);

    // Divide input Vocig by envelope
    float noiseLevel = 0.1;
    float *Vocig_envelope_untransp  = transp_dim3(Vocig_envelope         , nx, nlv, nz, 13);
    float *Vocig_envelope_untransp2 = transp_dim3(Vocig_envelope_untransp, nz, nlv, nx, 12);
    divideArraysWhiteNoise(Vocig, Vocig, Vocig_envelope_untransp2, nx*nz*nlv, noiseLevel);
    

    // float *Vocig_normalized_transp2 = CPU_zaloc1F(nx*nz*nlv);
    // divideArraysWhiteNoise(Vocig_normalized_transp2, Vocig_transp2, Vocig_envelope, nx*nz*nlv, noiseLevel);
    // float *Vocig_normalized_transp  = transp_dim3(Vocig_normalized_transp2, nx, nlv, nz, 13);
    // float *Vocig_normalized = transp_dim3(Vocig_normalized_transp, nz, nlv, nx, 12);

    // Free memory
    free(Vocig_transp);
    free(Vocig_transp2);
    free(Vocig_envelope);
    free(Vocig_envelope_untransp);
    free(Vocig_envelope_untransp2);
    // free(Vocig_normalized_transp);
    // free(Vocig_normalized_transp2);
return;
}

float* filterECIG_Angle(float *eimage, EIProc *eiproc, int nx, int nz, int nxb, int nzb, \
                        float dx, float dz, float ox, float oz, int ltheta0, float dtheta, int lx0, \
                        int iproc, int procGPU)
{
    if(!eiproc->applyTaperAngle)    return NULL;

    int      nlx = 2*lx0 + 1;
    int   ntheta = 2*ltheta0 + 1;
    float otheta = -ltheta0*dtheta;

    float *eimage_theta        = CPU_zaloc1F(nx*nz*ntheta);

    GPU_WRAP_lx2theta(0, procGPU, 0, eimage, eimage_theta, nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);
        
    // Applying cosine to ADCIG (Jacobian of the Offset to Angle transformation)
    // applyCosADCIG(eimage_theta, eimage_theta, nx, nz, ntheta, dtheta);

    // Angle tapering
    // printf("\n  eiproc->theta1=%f  eiproc->theta2=%f  eiproc->theta3=%f  eiproc->theta4=%f  \n", \
            eiproc->theta1, eiproc->theta2, eiproc->theta3, eiproc->theta4);
    taperADCIG(eimage_theta, eimage_theta, nx, nz, ntheta, dtheta, eiproc->theta2,  eiproc->theta1, 1.0);
    taperADCIG(eimage_theta, eimage_theta, nx, nz, ntheta, dtheta, eiproc->theta3,  eiproc->theta4, 1.0);
    taperADCIG(eimage_theta, eimage_theta, nx, nz, ntheta, dtheta, eiproc->theta4, +dtheta*ltheta0, 0.0);  
    taperADCIG(eimage_theta, eimage_theta, nx, nz, ntheta, dtheta, eiproc->theta1, -dtheta*ltheta0, 0.0);

    float *eimage_taper_angle = cg_theta2lx_gpu(procGPU, eimage_theta, nx, nz, dx, dz, lx0, ntheta, dtheta, otheta, 7);
    free(eimage_theta);

return eimage_taper_angle;
}

void reciprocityAngle(float *eimage, EIProc *eiproc, int acqGeom, \
                        int nx, int nz, int nxb, int nzb, \
                        float dx, float dz, float ox, float oz, \
                        int ltheta0, float dtheta, int lx0, \
                        int iproc, int procGPU)
{
    if(!eiproc->makeReciprocalADCIG)    return;

    int nlx = 2*lx0+1;
    int   ntheta = 2*ltheta0 + 1;
    float otheta = -ltheta0*dtheta;

    float *eimage_theta        = CPU_zaloc1F(nx*nz*ntheta);

    GPU_WRAP_lx2theta(0, procGPU, 0, eimage, eimage_theta, nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);

    makeReciprocalADCIG(eimage_theta, eimage_theta, nx, nz, nxb, nzb, ntheta, acqGeom);

    float *eimage_reciprocal = cg_theta2lx_gpu(procGPU, eimage_theta, nx, nz, dx, dz, lx0, ntheta, dtheta, otheta, 7);
    free(eimage_theta);

    memset(eimage, 0, nx*nz*nlx*sizeof(float));
    memcpy(eimage, eimage_reciprocal, nx*nz*nlx*sizeof(float));
    free(eimage_reciprocal);

return;
}






void initializeFileNames(int iter)
{
    // Extended images file names
    sprintf(gradient_FileName                             , "gradient_Iter%d"                       , iter);
    sprintf(image_FileName                                , "image_Iter%d"                          , iter);
    sprintf(eimage_FileName                               , "eimage_Iter%d"                         , iter);
    sprintf(eimage_LambdaDip_FileName                     , "eimage_LambdaDip_Iter%d"               , iter);
    sprintf(eimage_LambdaDip_contracted_FileName          , "eimage_LambdaDip_contracted_Iter%d"    , iter);
    sprintf(eimage_LambdaDip_contracted_ref_FileName      , "eimage_LambdaDip_contracted_ref_Iter%d", iter);
    sprintf(eimage_taper_FileName                         , "eimage_taper_Iter%d"                   , iter);
    sprintf(adcigs_FileName                               , "adcigs_Iter%d"                         , iter);
    sprintf(adcigs_taper_angle_FileName                   , "adcigs_taper_angle_Iter%d"             , iter);
    sprintf(adcigs_filterECIG_FileName                    , "adcigs_filterECIG_Iter%d"              , iter);
    sprintf(Gocig_FileName                                , "Gocig_Iter%d"                          , iter);
    sprintf(Gocig_foc_FileName                            , "Gocig_contracted_Iter%d"               , iter);
    sprintf(Gocig_ref_FileName                            , "Gocig_contracted_ref_Iter%d"           , iter);
    sprintf(Docig_FileName                                , "Docig_Iter%d"                          , iter);
    sprintf(Docig_foc_FileName                            , "Docig_foc_Iter%d"                      , iter);
    sprintf(Docig_ref_FileName                            , "Docig_ref_Iter%d"                      , iter);
    sprintf(Hocig_FileName                                , "Hocig_Iter%d"                          , iter);
    sprintf(Hocig_foc_FileName                            , "Hocig_foc_Iter%d"                      , iter);
    sprintf(Hocig_ref_FileName                            , "Hocig_ref_Iter%d"                      , iter);
    sprintf(Vocig_FileName                                , "Vocig_Iter%d"                          , iter);
    sprintf(Vocig_foc_FileName                            , "Vocig_foc_Iter%d"                      , iter);
    sprintf(Vocig_ref_FileName                            , "Vocig_ref_Iter%d"                      , iter);
    sprintf(Hocig_dif_FileName                            , "Hocig_dif_Iter%d"                      , iter);
    sprintf(Hacig_dif_FileName                            , "Hacig_dif_Iter%d"                      , iter);
    sprintf(Vocig_dif_FileName                            , "Vocig_dif_Iter%d"                      , iter);
    sprintf(Vacig_dif_FileName                            , "Vacig_dif_Iter%d"                      , iter);
    sprintf(Hocig_proc_FileName                           , "Hocig_proc_Iter%d"                     , iter);
    sprintf(Vocig_proc_FileName                           , "Vocig_proc_Iter%d"                     , iter);
    sprintf(Hacig_FileName                                , "Hacig_Iter%d"                          , iter);
    sprintf(Hacig_foc_FileName                            , "Hacig_foc_Iter%d"                      , iter);
    sprintf(Hacig_ref_FileName                            , "Hacig_ref_Iter%d"                      , iter);
    sprintf(Vacig_FileName                                , "Vacig_Iter%d"                          , iter);
    sprintf(eimage_taper_angle_FileName                   , "eimage_taper_angle_Iter%d"             , iter);
    sprintf(eimage_filterECIG_FileName                    , "eimage_filterECIG_Iter%d"              , iter);
    sprintf(eimage_contracted_FileName                    , "eimage_contracted_Iter%d"              , iter);
    sprintf(adcigs_contracted_FileName                    , "adcigs_contracted_Iter%d"              , iter);
    sprintf(eimage_contracted_ref_FileName                , "eimage_contracted_ref_Iter%d"          , iter);
    sprintf(adcigs_contracted_ref_FileName                , "adcigs_contracted_ref_Iter%d"          , iter);
    sprintf(eimage_diff_FileName                          , "eimage_diff_Iter%d"                    , iter);
    sprintf(adcigs_diff_FileName                          , "adcigs_diff_Iter%d"                    , iter);

return;
}


