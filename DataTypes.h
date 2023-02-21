#ifndef DATA_TYPES_H
#define DATA_TYPES_H

// Structs
typedef struct TapeLX_ {
    float xpos;
    float lx1_i, lx2_i, lx3_i, lx4_i;
    float lx1_f, lx2_f, lx3_f, lx4_f;
    float z_i, z_f;    
} TapeLX;

typedef struct Shot_ {
    int ishot;
    int nshot;
    int nsou, nrec, nt;
    int ixMin, ixMax;
    int *ixs, *izs;
    float *xs, *zs;
    int *ixr, *izr;
    float *xr, *zr;
    float dt;
    float ot;
    float *source;
    float *seismogram;
} Shot;

typedef struct Dataset_ {
    int nshot;
    Shot *shot;
} Dataset;

typedef struct Parms_ 
{
    // GPU distribution
    int ngpu;
    int firstGPU;

    int RUN_MODE;

    int computeOnlyGrad_NoLineSearch;

    int performEBM;

    // Obj Func Explorer
    float OF_vel_ref_step;
    float OF_vel_mis_step;
    float OF_vel_ref_mag;
    float OF_vel_mis_mag;

    // Grid
    float LX;
    float LY;
    float LZ;
    float LT;
    float dx;
    float dy;
    float dz;
    float dt;
    float ox;
    float oy;
    float oz;

    // Imaging
    float maxSubOff_lx;
    float maxSubOff_lz;
    float maxSubOff_lt;
    int   jt_lt;
    int   jt;
    int   jt_EIC;

    // Inversion parameters
    int imageDomain;
    int fixedEImage;
    int secondTerm;
    int   n_iter;
    float alpha;
    float muFocus;
    float epsFocus;
    float lambdaNull;
    float eimageFocus_zMin;
    float eimageFocus_zMax;

    int saveReferenceData;
    int saveMisfitData;

    // shot modeling
    float f1, f2, f3, f4;
    float agc_time_window;
    float dxShot;
    float offMax;
    int   acquisitionGeom;
    float offPerc;


    // eimage filtering
    int applyBandass2EImage;
    float lambda1, lambda2, lambda3, lambda4;

    // Input model
    char executionDirectory[1024];
    char *inputExecPathName;

    int nxInput_vel;
    int nzInput_vel;
    float deltaXInput_vel;
    float deltaZInput_vel;
    float OrigXInput_vel;
    float OrigZInput_vel;
    float velScalar;    

    int   DA_removeDA;
    float DA_waterSurfaceVel;
    float DA_offMin, DA_offMax;

    char inputSeismicFileName[1024];
    char inputSeismicFileName_base[1024];
    char inputSeismicFileName_moni[1024];
    char filePathVpTrue[1024];
    char filePathVsTrue[1024];
    char filePathRhoTrue[1024];
    char filePathEpsTrue[1024];
    char filePathDltTrue[1024];
    char filePathDipTrue[1024];
    char filePathAzimTrue[1024];
    char filePathVpInit[1024];
    char filePathVsInit[1024];
    char filePathRhoInit[1024];
    char filePathEpsInit[1024];
    char filePathDltInit[1024];
    char filePathDipInit[1024];
    char filePathAzimInit[1024];
    char filePathVpDA[1024];
    char filePathRhoDA[1024];

    char outputSeismicFileName[1024];
    char filePathInputData[1024];
    char filePathInputSource[1024];
    char filePathVel[1024];
    char filePathDAVel[1024];

    // eimage angle tapering
    int   eimageTaper_applyTaperAngle;
    float eimageTaper_theta1;
    float eimageTaper_theta2;
    float eimageTaper_theta3;
    float eimageTaper_theta4;
    
    // eimage tapering on Z and X
    int   eimageTaper_applyZTaper;
    float eimageTaper_z0;
    float eimageTaper_z1;
    float eimageTaper_z2;
    float eimageTaper_z3;
    int   eimageTaper_applyXTaper;
    float eimageTaper_x1;
    float eimageTaper_x2;
    float eimageTaper_x3;
    float eimageTaper_x4;

    // eimage tapering on LX and TAU
    int   eimageTaper_applyLXTaper;
    int   eimageTaper_applyTauTaper;
    float eimageTaper_z_i;
    float eimageTaper_z_f;
    float eimageTaper_lx1_i;
    float eimageTaper_lx2_i;
    float eimageTaper_lx3_i;
    float eimageTaper_lx4_i;
    float eimageTaper_lx1_f;
    float eimageTaper_lx2_f;
    float eimageTaper_lx3_f;
    float eimageTaper_lx4_f;
    float eimageTaper_tau1_i;
    float eimageTaper_tau2_i;
    float eimageTaper_tau3_i;
    float eimageTaper_tau4_i;
    float eimageTaper_tau1_f;
    float eimageTaper_tau2_f;
    float eimageTaper_tau3_f;
    float eimageTaper_tau4_f;

    // gradient bandpass
    int gradProc_applyBandpass;
    float gradProc_LBDA1;
    float gradProc_LBDA2;
    float gradProc_LBDA3;
    float gradProc_LBDA4;

    // gradient tapering in Z
    int   gradProc_applyZTaper;
    float gradProc_z1;
    float gradProc_z2;
    float gradProc_z3;
    float gradProc_z4;
    
    // gradient tapering in X
    int   gradProc_applyXTaper;
    float gradProc_x1;
    float gradProc_x2;
    float gradProc_x3;
    float gradProc_x4;

    // gradient smoothing
    int gradProc_applySmoothing;
    float gradProc_halfLengthX;
    float gradProc_halfLengthZ;

    int   applyADFocusing;

    int makeReciprocalADCIG;

    TapeLX *tapelx;
    int ntape;

    int   encoding;
    int   codelength;
    float codestep;
    int   realizations;
    
    float EIMSH_azmMin;
    float EIMSH_azmMax;
    float EIMSH_dazm;

    float EIMSH_dipMin;
    float EIMSH_dipMax;
    float EIMSH_ddip;

    float EIMSH_dlambda;
    float EIMSH_lambdaMax;

    float EIMSH_dtau;
    float EIMSH_tauMax;
    float EIMSH_tauNull;

    int   EIMSH_CIGS_TAU;
    int   EIMSH_CIGS_ORTH;
    int   EIMSH_CIGS_VCIG;

    int   EIMSH_focusingQC;

    float EIMSH_dtheta;
    float EIMSH_thetaMax;

    int extendInitModel;
    float extendInitModelLeft_x;
    float extendInitModelRight_x;


    // ==============================================================
    // ======================  New parameters  ======================
    // ==============================================================
    
    int nxb, nyb, nzb;

    int normalizeInputDataAmplitudes;

    float AppertureExtension_x;

    // Parameters for application of top mute
    int   TopMute_apply;
    float TopMute_t0;
    float TopMute_t1;
    float TopMute_off0;
    float TopMute_off1;
    float TopMute_ramp;

    // Parameters for application of bottom mute
    int   BottomMute_apply;
    float BottomMute_t0;
    float BottomMute_t1;
    float BottomMute_off0;
    float BottomMute_off1;
    float BottomMute_ramp;

    float shotOX, shotOY, shotOZ;
    float shotDX, shotDY, shotDZ;
    int   shotNX, shotNY, shotNZ;

    float recOX, recOY, recOZ;
    float recDX, recDY, recDZ;
    int   recNX, recNY, recNZ;

    int recCoordRefFromShot;

    // Parameters for direct arrival (DA) attenuation
    int   DA_modelAndRemoveDA;

} Parms;


typedef struct EIProc_ {

    // focusing operator
    float muFocus, epsFocus, lambdaNull, tauNull;
    float zMin, zMax;

    // eimage filtering
    int applyBandass2EImage;
    float lambda1, lambda2, lambda3, lambda4;

    // eimage angle tapering
    float applyTaperAngle;
    float theta1;
    float theta2;
    float theta3;
    float theta4;

    // eimage tapering on Z
    int   applyZTaper;
    float z0;
    float z1;
    float z2;
    float z3;
    
    // eimage tapering on X
    int   applyXTaper;
    float x1;
    float x2;
    float x3;
    float x4;

    // eimage tapering on LX
    int   applyLXTaper;
    float z_i;
    float z_f;
    float lx1_i;
    float lx2_i;
    float lx3_i;
    float lx4_i;
    float lx1_f;
    float lx2_f;
    float lx3_f;
    float lx4_f;

    // eimage tapering on Tau
    int   applyTauTaper;
    float tau1_i;
    float tau2_i;
    float tau3_i;
    float tau4_i;
    float tau1_f;
    float tau2_f;
    float tau3_f;
    float tau4_f;

    int applyADFocusing;

    int makeReciprocalADCIG;

} EIProc;

typedef struct GRDProc_ {

    // gradient tapering in Z
    int   applyZTaper;
    float z1;
    float z2;
    float z3;
    float z4;
    
    // gradient tapering in X
    int   applyXTaper;
    float x1;
    float x2;
    float x3;
    float x4;

    // gradient smoothing
    int applySmoothing;
    float halfLengthX;
    float halfLengthZ;
    

} GRDProc;


typedef struct InvParms_ {

    int n_iter;
    float alpha;    

} InvParms;


typedef struct MSH_ {

    float LX, LY, LZ, LT;
    float ox, oy, oz, ot;
    float dx, dy, dz, dt;
    int   nx, ny, nz, nt;
    int   nxb, nyb, nzb;

    int   jt, jt_EIC, jt_lt;

} MSH;

typedef struct Frame_ {

    int ixMin, ixMax;
    int iyMin, iyMax;

    float x1, y1;
    float x2, y2;
    float x3, y3;
    float x4, y4;

} Frame;


typedef struct Eimage_MSH_ {

    float azmMin, azmMax, dazm;
    float dipMin, dipMax, ddip;
    float lambdaMax, dlambda;
    int   nlambda, ndip;
    int   nx_LambdaDip, nz_LambdaDip;
    float ox_LambdaDip, oz_LambdaDip;
    
    float dlx, dlz, dlv, dly, dlt;
    int   lx0, lz0, lv0, ly0, lt0;
    int   nlx, nlz, nlv, nly, nlt;

    int   tau0, ntau;
    float tauMax, dtau;

    int   ltheta0, ntheta;
    float otheta, dtheta;

} Eimage_MSH;


typedef struct AcquisitionGeom_
{
    int nshot, nrec;

    float shotOX, shotOY, shotOZ;
    float shotDX, shotDY, shotDZ;
    int   shotNX, shotNY, shotNZ;

    float recOX, recOY, recOZ;
    float recDX, recDY, recDZ;
    int   recNX, recNY, recNZ;

    int  recCoordRefFromShot;

    float f1, f2, f3, f4;    

} AcquisitionGeom;

typedef struct Shot3D_
{
    int ishot;
    int nshot;
    int nsou, nrec, nt;
    int ixMin, ixMax;
    int iyMin, iyMax;
    int  *ixs, *iys, *izs;
    float *xs,  *ys,  *zs;
    float *alpha_xs,  *alpha_ys,  *alpha_zs;
    float *coeff_sou;
    int  *ixr, *iyr, *izr;
    float *xr,  *yr,  *zr;
    float *alpha_xr,  *alpha_yr,  *alpha_zr;
    float *coeff_rec;
    float dt;
    float ot;
    float *source;
    float *seismogram;
    float maxFrequency;
    long long int nStationsPreceding;

    int *recMap_P2Z;
    int *recMap_P2PZ;

    // Parameters for direct arrival (DA) attenuation
    int   DA_removeDA;
    float DA_waterSurfaceVel;
    float DA_offMin;
    float DA_offMax;

    // Parameters for application of top mute
    int   TopMute_apply;
    float TopMute_t0;
    float TopMute_t1;
    float TopMute_off0;
    float TopMute_off1;
    float TopMute_ramp;

    // Parameters for application of bottom mute
    int   BottomMute_apply;
    float BottomMute_t0;
    float BottomMute_t1;
    float BottomMute_off0;
    float BottomMute_off1;
    float BottomMute_ramp;

} Shot3D;

typedef struct Dataset3D_ {
    FILE *fMet0;
    FILE *fMet1;
    FILE *fMet2;
    FILE *fData;
    FILE *oMet0;
    FILE *oMet1;
    FILE *oMet2;
    FILE *oData;
    int OdisseiaFileType;
    int nshot;
    int nshot_aux;
    Shot3D *shot;
    AcquisitionGeom acqGeom;

    int *shotMap_P2Z;
    int *shotMap_P2PZ;
} Dataset3D;

#endif