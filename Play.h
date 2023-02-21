char Vel_TWI_FileName[1024];
char gradient_FileName[1024];
char image_FileName[1024];
char eimage_FileName[1024];
char Gocig_FileName[1024];
char Gocig_foc_FileName[1024];
char Gocig_ref_FileName[1024];
char Docig_FileName[1024];
char Docig_foc_FileName[1024];
char Docig_ref_FileName[1024];
char Hocig_FileName[1024];
char Hocig_foc_FileName[1024];
char Hocig_ref_FileName[1024];
char Hacig_foc_FileName[1024];
char Hacig_ref_FileName[1024];
char Vocig_FileName[1024];
char Vocig_foc_FileName[1024];
char Vocig_ref_FileName[1024];
char Hocig_dif_FileName[1024];
char Hacig_dif_FileName[1024];
char Vocig_dif_FileName[1024];
char Vacig_dif_FileName[1024];
char Hocig_proc_FileName[1024];
char Vocig_proc_FileName[1024];
char Hacig_FileName[1024];
char Vacig_FileName[1024];
char eimage_LambdaDip_FileName[1024];
char eimage_LambdaDip_contracted_FileName[1024];
char eimage_LambdaDip_contracted_ref_FileName[1024];
char eimage_taper_FileName[1024];
char adcigs_FileName[1024];
char adcigs_filterECIG_FileName[1024];
char adcigs_taper_angle_FileName[1024];
char eimage_filterECIG_FileName[1024];
char eimage_taper_angle_FileName[1024];
char eimage_contracted_FileName[1024];
char eimage_contracted_ref_FileName[1024];
char adcigs_contracted_FileName[1024];
char adcigs_contracted_ref_FileName[1024];
char eimage_diff_FileName[1024];
char adcigs_diff_FileName[1024];

char refDataFileName[1024];
char misDataFileName[1024];
char resDataFileName[1024];


// TL-TWI file names
char Dvel_FileName[1024];

char Tcig_FileName[1024];
char Tcig_proc_FileName[1024];
char Tcig_ref_FileName[1024];
char Tcig_foc_FileName[1024];
char Tcig_dif_FileName[1024];

char image_base_FileName[1024];
char Hocig_base_FileName[1024];
char Hocig_base_proc_FileName[1024];
char Hocig_base_foc_FileName[1024];
char Vocig_base_FileName[1024];
char Vocig_base_proc_FileName[1024];
char Vocig_base_foc_FileName[1024];

char image_moni_FileName[1024];
char Hocig_moni_FileName[1024];
char Hocig_moni_proc_FileName[1024];
char Hocig_moni_foc_FileName[1024];
char Vocig_moni_FileName[1024];
char Vocig_moni_proc_FileName[1024];
char Vocig_moni_foc_FileName[1024];


#define MET1_BYTESPERSHOT    (  (long long int) ( 3*sizeof(int) + 3*sizeof(float) + 1*sizeof(long long int) )  )
#define MET2_BYTESPERSTATION (  (long long int) ( 3*sizeof(long long int) )  )

#define  ODISSEIA_FILE_TYPE_SEISMIC              1
#define  ODISSEIA_FILE_TYPE_P_VELOCITY           2
#define  ODISSEIA_FILE_TYPE_S_VELOCITY           3
#define  ODISSEIA_FILE_TYPE_DENSITY              4
#define  ODISSEIA_FILE_TYPE_SOURCE_SIGNATURE     5

#define  ODISSEIA_OPEN_MODE_CREATENEW            1
#define  ODISSEIA_OPEN_MODE_OPENEXISTING         2

// Parameters
typedef int parmType;
#define TYPE_INT 1
#define TYPE_FLOAT 2
#define TYPE_DOUBLE 3
#define TYPE_STRING 4
static char inputFileName[1024];
static char executionDirectory[1024];



#define FILENAME_MAXLENGTH 1024
#define FILETYPE_LENGTH 16

#define OL 5

#define id2(iz,ix)                 ((iz)+(ix)*nz)
#define id3(iz,ix,iy)              ((iz)+(ix)*nz+(iy)*nx*nz)

#define id3_eic_lx(lx,iz,ix)                  ( (lx)+lx0 + (iz)*nlx + (ix)*nz*nlx)
#define id3_eic_lx_lz(lx,lz,ix,iz)            ( (lx)+lx0 + ((lz)+lz0)*nlx + (iz)*nlx*nlz + (ix)*nz*nlx*nlz)
#define id3_eic_lx_lst(itheta,lx,iz,ix)       ( (itheta)+ltheta0 + ((lx)+lx0)*ntheta + (iz)*nlx*ntheta + (ix)*nz*nlx*ntheta)


#define id3_eic_th(itheta,iz,ix)   ( (itheta)+ltheta0 + (iz)*ntheta + (ix)*nz*ntheta )

#define min(a,b) ((a)<=(b)) ? (a) : (b)
#define max(a,b) ((a)>=(b)) ? (a) : (b)

// ------------------------------
// ------------------------------
// Finite Difference Coefficients
// ------------------------------
// ------------------------------


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

#define RUN_MODELING        0
#define RUN_TWI_REG         1
#define RUN_TWI_TL          2
#define RUN_TWI_OBJFUNCEXP 10

#define FILE_POINTERS_INPUT       0
#define FILE_POINTERS_OUTPUT      1


// Parameters
int readParms(Parms *parms, int iproc);
int getParm(char *parmName, void* parm, char *requirement, parmType TYPE);
int getStringParm(char *inputFileName, char *parmName, char* parm, char *requirement);
int getStringListParm(char *inputFileName, char *parmName, int *nitens, char parm[][1024], char *requirement);
int readCommandLineArg(void* parm, int argc, char** argv, char *argName, char *requirement, parmType TYPE);

// Mem alloc
void    zeroArray(float *array, int n);
void fftw_freeMem(fftw_complex *array, long long int n1);
void freeMem(void *array, size_t size);
fftw_complex* CPU_zalocFFTW1(long long int n1);
double* CPU_zaloc1D(long long int n1);
float*  CPU_zaloc1F(long long int n1);
float** CPU_zaloc2F(long long int n1, long long int n2);
int*    CPU_zaloc1I(long long int n1);
int**   CPU_zaloc2I(long long int n1, long long int n2);
Shot*   CPU_zaloc1Shot(int n1);
Shot3D* CPU_zaloc1Shot3D(int n1);

// Starting model and source
float* getVp_VofZ(int nx, int nxb, int  nz, int nzb, float dx, float dz, \
                  float vp0, float coef, float zref, float lref);
void   getVp_add1Reflector(float *m, int nx, int nxb, int  nz, int nzb, float dx, float dz, float dvp, float zref, float lref);
void   getVp_add1DippingReflector(float *m, int nx, int nxb, int  nz, int nzb, float dx, float dz, float ox, float oz, \
                                  float dvp, float xleft, float xright, float zleft, float zright, float width);
float* getVp_1Reflector(int nx, int nxb, int  nz, int nzb, float dx, float dz, float vp0, float dvp, float zref, float lref);
float* getVp_1Reflector_1GaussAnomaly(int nx, int nxb, int  nz, int nzb, float dx, float dz, \
                                      float vp0, float dvp, float zref, float lref, \
                                      float dvpGauss, float gaussRadius, float x0, float z0);
void getVp_addGaussAnomaly(float *vp, int nx, int nxb, int  nz, int nzb, float dx, float dz, \
                           float dvpGauss, float gaussRadius, float x0, float z0);
void   transformVel(float *vp, int nx, int nz, float dx, float dt);
void   untransformVel(float *vp, int nx, int nz, float dx, float dt);
void transformVelSlowness(float *m, long long int n);
void model_vel2SquaredVaga(float *m, int nx, int nz);
void model_squaredVaga2Vel(float *m, int nx, int nz);
float* getRicker(int nt, float dt, float fm, float t0);
void  getLame_1Reflector(float *lambda, float *mu, float *lambda2mu, float *b_c, float *b_d, \
                        int nx, int nxb, int  nz, int nzb, float dx, float dz, \
                        float vp0, float dvp, float vs0, float dvs, float rho0, float drho, float zref, float lref);
void getLame_1Reflector_1GaussAnomaly(float *lambda, float *mu, float *lambda2mu, float *bc, float *bd, \
                                      int nx, int nxb, int  nz, int nzb, float dx, float dz, \
                                      float vp0, float dvp, float vs0, float dvs, float rho, float drho, float zref, float lref, \
                                      float dvpGauss, float dvsGauss, float gaussRadius, float x0, float z0);

float* extendInputModel(float *inpModel, int nx, int nz, int nxInput, int nzInput, float scalar);
float* extendResampInputModel(float *inpModel, int nxInput, int nzInput, float dxInput, float dzInput, float oxInput, float ozInput, \
                              int nx, int nz, float dx, float dz, float ox, float oz, int nxb, int nzb, float scalar);

void fillBoundariesInUpdatedModel(float *model, int nxb, int nx, int nz, float x1, float x2, float dx);
void padBoundaries(float *model, int nx, int nz, int nxb, int nzb);
void fillGradientEdges(float *gradient, int nx, int nz, float dx, float x1, float x2);

// Dataset structure
void initDataset(Dataset *dataset, int nshot, float dx_rec, \
                 int nt, int nx, int nz, float dt, float dx, float dz, \
                 float offMax, float f1, float f2, float f3, float f4);

void finalizeDataset(Dataset *dataset);

// (Odisseia Shot) Read/Save Odisseia files
void saveDataset3D(char *fileName, Dataset3D *dataset);
void saveDataset3D_metaData(char fileName[], char fileNameBinary[], Dataset3D *dataset, int filePointers, int saveNewBinaryFile);
void saveDataset3D_binaryData(Dataset3D *dataset, int filePointers, int ishot);
void closeOdisseiaFiles(Dataset3D *dataset, int filePointers);
void initDataset3DfromFile(Dataset3D *dataset, Parms *parms);
void initTimeLapseDataset3DfromFile(Dataset3D *dataset_base, Dataset3D *dataset_moni, Parms *parms);
void initDataset3D(Dataset3D *dataset, MSH *msh, Parms *parms);
void initSeismogramAndSourceMemory(Dataset3D *dataset, int ishot);
void initSeismogramMemory(Dataset3D *dataset, int ishot);
void initSourceMemory(Dataset3D *dataset, int ishot);
void initShotFilteringParms(Dataset3D *dataset, Parms *parms);
int initShot3DfromFile(Dataset3D *dataset, int ishot);
void loadSeismogram2Dataset(Dataset3D *dataset, int ishot, int normalizeInputDataAmplitudes);
void initShot3D(Dataset3D *dataset, MSH *msh);
Shot3D* initCopyShot3D(Dataset3D *dataset_inp, int ishot);
void initCoordinatesAndWavelet(Dataset3D *dataset, MSH *msh, AcquisitionGeom *acqGeom);
void finalizeShot3D(Shot3D *shot);
int finalizeShot3DfromFile(Dataset3D *dataset, int ishot);
int finalizeSeismogramAndSource(Dataset3D *dataset, int ishot);
int finalizeShot3DMetaInfo(Shot3D *shot);
int finalizeShotSeismogramAndSource(Shot3D *shot);
void finalizeDataset3D(Dataset3D *dataset);
float *getBandpassPulse3D(int nt, float dt, float f1, float f2, float f3, float f4, float t0);

//Math
void modelBandpass(float *model, int nx, int nz, float dz, float k1, float k2, float k3, float k4);
void modelBandpass_Lambdas(float *model, int nx, int nz, float dz, float lambda1, float lambda2, float lambda3, float lambda4);

// Shot bounding box
int checkShotBoundaries3D(MSH *msh, Shot3D *shot);
int gridShot3D(MSH *msh, Shot3D *shot);
void resampleSourceAndRec(Shot3D *shot, float dt, int nt);



// (Odisseia IO) Initialize files in Odisseia format
// Open files in Odisseia format
void openOdisseiaFiles(FILE **fMet0, FILE **fMet1, FILE **fMet2, FILE **fData, char *fileNameMet0, int OdisseiaFileType, char mode[]);
void getFileNameOdisseia(char *fileNameMet, FILE *fMet0, int which);
void readOdisseiaMetaFile_get_nshot(FILE *fMet1, int *nshot);
void readOdisseiaMetaFile_get_shotInfo(FILE *fMet1, int ishot, int *nsou, int *nrec, int *nt, float *dt, float *ot, float *maxFreq, long long int *nStationsPreceding);
void readOdisseiaMetaFile_get_shotCoords(FILE *fMet1, FILE *fMet2, int ishot, \
                                         int nsou, int nrec, long long int nStationsPreceding, \
                                         float *xs, float *ys, float *zs, \
                                         float *xr, float *yr, float *zr);
void readSeismogramFromFile(FILE *fData, float *source, float *seismogram, int ishot, int nsou, int nrec, int nt, long long int nStationsPreceding);

void fillOdisseiaMetaFile_Seismic_ShotInfo(FILE *fMet1, Dataset3D *dataset);
void fillOdisseiaMetaFile_Seismic_ShotCoords(FILE *fMet2, Dataset3D *dataset);
void writeSeismogram2DataFile(FILE *fData, Dataset3D *dataset, int ishot2write);
void writeDataset2DataFile(FILE *fData, Dataset3D *dataset);
void writeSeismogram2File(FILE *fData, float *source, float *seismogram, int nsou, int nrec, int nt, long long int nStationsPreceding);
void createOdisseiaFilesNew(char *fileName, FILE **fMet0, FILE **fMet1, FILE **fMet2, FILE **fData, int OdisseiaFileType);
void createOdisseiaFilesNew_ExistingBinary(char *fileName, char *fileNameBinary, FILE **fMet0, FILE **fMet1, FILE **fMet2, FILE **fData, int OdisseiaFileType);
void initOdisseiaHeaderFile(FILE *fMet0, char *fileNameMet1, char *fileNameMet2, char *fileNameData, int OdisseiaFileType);
int appendOdisseiaMetaFile_Seismic_ShotInfo(FILE *fMet1, Dataset3D *dataset);
void appendOdisseiaMetaFile_Seismic_ShotCoords(FILE *fMet1, FILE *fMet2, Dataset3D *dataset);
void appendSeismogram2DataFile(FILE *fData, Dataset3D *dataset, int ishot2write);

// Shot structure
void initShot(Shot *shot, int ishot, int nshot, float dx_rec, int nt, \
              int nx, int nz, float dt, float dx, float dz, \
              float offMax, float f1, float f2, float f3, float f4);
void initSou(Shot *shot, int ishot, int nshot, int nt, int nx, int nz, float dt, float dx, float dz, \
             float offMax, float f1, float f2, float f3, float f4);
void initRec(Shot *shot, int ishot, int nshot, float dx_rec, int nx, int nz, float dx, float dz);
void finalizeShot(Shot *shot);
double addObjFunc(Shot *shot);
void addShots(Shot *shot_F, float coef_A, Shot *shot_A, float coef_B, Shot *shot_B);
void removeDA(Shot *shot, Shot *shot_da);
void removeDA_3D(Shot3D *shot, Shot3D *shot_da);
void copyData(Shot *shot_inp, Shot *shot_out);
void copyDataset(Dataset *dst, Dataset *src);
void applyHalfDerivative2Source(Shot *shot);
void applyHalfDerivative2Source3D(Shot3D *shot);
void applyHalfDerivative2Data(Shot *shot);
void applyHalfDerivative2Shot(Shot3D *shot);
void applyHalfDerivative2Dataset(Dataset3D *data);
void applyShotTaper(Shot *shot, int nx);
void applyRecTaper(Shot3D *shot, float perc, float offMax, int sensorSide);
void applyWaterBottomMute(Shot *shot, float waterBottomDepth, float taper);
float *getBandpassPulse(int nt, float dt, float f1, float f2, float f3, float f4, float t0);
void integrateRecTime(Shot *shot);
void applyBandpassToTrace(Shot *shot, float f1, float f2, float f3, float f4);
void multiplySeismogramByTime(Shot *shot, float power);
void multiplyTraceByTime(float *trace, int nt, float dt, float power);
void demultiplexSeismogram(float *seismogram, int nt, int nrec);
void normalizeDatasetAmplitudes(Dataset *dataset);
void AGC(Dataset *dataset, float time_window);
void makeEncodedDataset(int idev, Dataset *data_inp, Dataset *data_tmp, Dataset *data_enc, Dataset *data_res, int codelength, int codestep, int realizations);
int  shotEncoder(int idev, Dataset *data_in, Dataset *data_ou, int codelength, int codestep, int seed);
void correlateWavelets(Dataset *data);
void windowData(int idev, Dataset *dataset, float tmin, float tmax, float dt);
void writeSeisData2File(Dataset *dataset, int jshot, int iproc, char *fileNameSeismos, char *fileNameSources);

// Forward / Adjoint Acoustic
void stripBoundary(float *wholeArray, float *strippedArray, int nx, int nz, int mode);

void stripBoundary2(float *wholeArray, float *strippedArray, int nx, int nz, int mode);

void stripProp(float *wholeArray, float *strippedArray, int nx, int nz, int ixMin, int ixMax, int mode);

void modelShot(float *vp, Shot *shot, int jt, int nt, int nx, int nz, float dx, float dz, float dt, \
               float ox, float oz, float ot, int recordMovie, int ishotMovie);

void modelShot_GPU(int gpu_device, float *vp, Shot3D *shot, int jt, int nt, int nx, int nz, int nxb, int nzb, \
                   float dx, float dz, float dt, float ox, float oz, float ot, int recordMovie, int ishotMovie);

void modelShotBorn_GPU(int gpu_device, float *vp, float *ref, Shot *shot, int jt, \
                       int nt, int nx, int nz, float dx, float dz, float dt, \
                       float ox, float oz, float ot, int recordMovie, int ishotMovie);

void migrateShot(float *image, float *vp, Shot *shot, int jt, int nt, int nx, int nz, float dx, float dz, float dt, \
                 float ox, float oz, float ot, int recordMovie, int ishotMovie);

void emigrateShot_lx(float *image, float *eimage_shot, float *eimage_shot_theta, \
                     float *eimage, float *eimage_theta, float *vp, Shot *shot, \
                     int jt, int nt, int nx, int nz, float dx, float dz, float dt, \
                     float ox, float oz, float ot, int lx0, int ntheta, float dtheta, \
                     float otheta, int recordMovie, int ishotMovie);

void emigrateShot_lx_lz_lt_GPU(int gpu_device, float *image, float *eimage_shot, float *eimage_shot_theta, \
                         float *eimage, float *vp, Shot *shot, \
                         int jt, int nt, int nx, int nz, float dx, float dz, float dt, \
                         float ox, float oz, float ot, int lx0, int lz0, int lt0, int ntheta, float dtheta, float otheta, \
                         int recordMovie, int ishotMovie);

void emigrateShot_TLCIG_Mem_GPU(int gpu_device, float *image, float *ilumin, float *TLcig, float *vp, Shot *shot, 
                                long long int jt, long long int nt, long long int nx, long long int nz, float dx, float dz, float dt, \
                                float ox, float oz, float ot, long long int lt0, int recordMovie, int ishotMovie);

void emigrateShot_TAU_Mem_GPU(int gpu_device, float *image, float *ilumin, float *Tcig, float *vp, \
                              Shot3D *shot, MSH *msh, Eimage_MSH *eimsh, int recordMovie, int ishotMovie);

void emigrateShot_HOCIG_VOCIG_Mem_GPU(int gpu_device, float *image, float *ilumin, float *Hocig, float *Vocig, \
                                      float *vp, Shot3D *shot, MSH *msh, Eimage_MSH *eimsh, int recordMovie, int ishotMovie);

void emigrateShot_lx_lz_lt_Mem_GPU(int gpu_device, float *image, float *ilumin, float *eimage_shot, \
                                  float *eimage_shot_theta, float *eimage, float *vp, Shot *shot, int jt, \
                                  int nt, int nx, int nz, float dx, float dz, float dt, float ox, float oz, \
                                  float ot, int lx0, int lz0, int lt0, int jt_lt, int ntheta, float dtheta, \
                                  float otheta, int recordMovie, int ishotMovie);

void migrateExtBorn_lx(float *image1, float *image2, float *ilumin1, float *ilumin2, float *eref, float *vp, \
                       Shot *shot, int jt, int nt, int nx, int nz, float dx, float dz, \
                       float dt, float ox, float oz, float ot, int lx0, int ishotMovie, int recordMovie);

void migrateExtBorn_lx_lz_lt_GPU(float *imageSouSide, float *imageRecSide, \
                           float *imageSouSide_shot, float *imageRecSide_shot, \
                           float *iluminSouSide, float *iluminRecSide, \
                           float *eref, float *vp, Shot *shot, int jt, int nt, \
                           int nx, int nz, float dx, float dz, \
                           float dt, float ox, float oz, float ot, int lx0, int lz0, int lt0, \
                           int ishotMovie, int recordMovie);

void emigrateExplodingExtendedImages_HOCIG_VOCIG_Mem_GPU(int gpu_device, float *image,     float *ilumin, \
                                                         float *Hocig_background_original, float *Hocig_background_phaseShift, \
                                                         float *Hocig_residual_original,   float *Hocig_residual_phaseShift, \
                                                         float *Vocig_background_original, float *Vocig_background_phaseShift, \
                                                         float *Vocig_residual_original,   float *Vocig_residual_phaseShift, \
                                                         float *vp, long long int jt, int ntMax, \
                                                         long long int nt, long long int nx, long long int nz, \
                                                         float dx, float dz, float dt, float ox, float oz, float ot, \
                                                         long long int lx0, long long int lz0, int recordMovie);

void emigrateExplodingExtendedImages_DOCIG_Mem_GPU(int gpu_device, float *image, float *ilumin, \
                                                   float *Docig_background_original, float *Docig_background_phaseShift, \
                                                   float *Docig_residual_original,   float *Docig_residual_phaseShift, \
                                                   float *vp, long long int jt, int ntMax, \
                                                   long long int nt, long long int nx, long long int nz, \
                                                   float dx, float dz, float dt, float ox, float oz, float ot, \
                                                   long long int ndip, float dipMin, float ddip, long long int lambda0, int recordMovie);

void ExtendedBornModeling_TAU_Mem_GPU(int idev, float *Tcig, float *vp, Shot3D *shot, \
                                      MSH *msh, Eimage_MSH *eimsh, int ishotMovie, int recordMovie);

void migrateExtBorn_TAU_Mem_GPU(int idev, float *imageSouSide, float *imageRecSide, \
                                float *iluminSouSide, float *iluminRecSide, \
                                float *Tcig, float *vp, Shot3D *shot, \
                                MSH *msh, Eimage_MSH *eimsh, int ishotMovie, int recordMovie);

void migrateExtBorn_lx_lz_lt_Mem_GPU(int idev, float *imageSouSide, float *imageRecSide, \
                                     float *iluminSouSide, float *iluminRecSide, \
                                     float *eref, float *Heref, float *Veref, float *vp, Shot3D *shot, \
                                     MSH *msh, Eimage_MSH *eimsh, int ishotMovie, int recordMovie);


void migrateExtBorn_lx_lz_lt_Mem_GPU_OPTIM(int idev, float *imageSouSide, float *imageRecSide, \
                                           float *imageSouSide_shot, float *imageRecSide_shot, \
                                           float *iluminSouSide, float *iluminRecSide, \
                                           float *eref, float *Heref, float *Veref, float *vp, Shot *shot, \
                                           long long int jt, long long int jt_EIC, long long int nt, \
                                           long long int nx, long long int nz, float dx, float dz, \
                                           float dt, float ox, float oz, float ot, long long int lx0, long long int lz0, long long int lt0, \
                                           int ishotMovie, int recordMovie, int imageDomain);

void migrateExtBorn_lx_lz_lt_Mem_GPU_old(int idev, float *imageSouSide, float *imageRecSide, \
                                     float *imageSouSide_shot, float *imageRecSide_shot, \
                                     float *iluminSouSide, float *iluminRecSide, \
                                     float *eref, float *Heref, float *Veref, float *vp, Shot *shot, \
                                     long long int jt, long long int jt_EIC, long long int nt, \
                                     long long int nx, long long int nz, float dx, float dz, \
                                     float dt, float ox, float oz, float ot, long long int lx0, long long int lz0, long long int lt0, \
                                     int ishotMovie, int recordMovie, int imageDomain, int OPTIM);

// Forward / Adjoint Elastic
void modelShot_elastic(float *lambda, float *mu, float *lambda2mu, float *bc, float *bd, \
                       Shot *shot, int jt, int nt, int nx, int nz, \
                       float dx, float dz, float dt, float ox, float oz, float ot, int recordMovie, int ishotMovie);

//Born modeling
void modelShotBorn(float *ref, float *vp, Shot *shot, int jt, int nt, int nx, int nz, float dx, float dz, \
                   float dt, float ox, float oz, float ot);

void modelShotBorn_lx(float *eref, float *vp, Shot *shot, int jt, int nt, int nx, int nz, float dx, float dz, \
                   float dt, float ox, float oz, float ot, int lx0);

void modelShotBorn_lx_lz_lt_GPU(int idev, float *eref, float *Heref, float *Veref, float *vp, Shot *shot, int jt, int nt, int nx, int nz, float dx, float dz, \
                   float dt, float ox, float oz, float ot, int lx0, int lz0, int lt0, int recordMovie, int shotToRecord, int OPTIM);

void modelShotBorn_ampAngle(float *ref, float *vp, Shot *shot, int jt, int nt, int nx, int nz, float dx, float dz, \
                   float dt, float ox, float oz, float ot, int lx0, int ntheta);

// Acoustic Propagation
void propagate_Acoustic(float *wav1, float *wav2, float *lap, float *vp, int nx, int nz);

void propagate_Acoustic_abs(float *wav1, float *wav2, float *lap, float *vp, float *vp2dtsigma, int nx, int nz);

void propagate_Acoustic_abs2(float *wav1, float *wav2, float *lap, float *vp, float *sigma, float *sigmaInv, int nx, int nz);

void timeStep_Acoustic(float *u1, const float *u2, const float *lap, const float *vp2dt2dx2, int nx, int nz);

void timeStep_Acoustic_abs(float *u1, const float *u2, const float *lap, const float *vp2dt2dx2, const float *sigma, int nx, int nz);

void timeStep_Acoustic_abs2(float *u1, const float *u2, const float *lap, const float *vp2dt2dx2, \
                            const float *sigma, const float *sigmaInv, int nx, int nz);

void laplacian_10thOrder_2D_Acoustic_Iso(float *lap, const float *u2, int nx, int nz);

void injectSource_Acoustic(float *u, float *source, int it, int ix0, int iz0, int nx, int nz, float dx, float dz);

void injectSourceTaper_Acoustic(float *u, float *source, float taper, int it, int ix0, int iz0, int nx, int nz, float dx, float dz);

void recordData(Shot *shot, float *wav2, int nx, int nz, int it);

void injectData(Shot *shot, float *wav, int nx, int nz, float dx, float dz, int it);

void injectDataTaper(Shot *shot, float *wav, int nx, int nz, float dx, float dz, int it);

void imagingCondition(float *image, float *ws, float *wr, int nx, int nz);

void  eimagingCondition_lx(float *eimage, float *ws, float *wr, int nx, int nz, int lx0);

void eimagingCondition_lx_theta(float *eimage, float *eimage_shot, float *eimage_theta, float *eimage_shot_theta, \
                                float *ws, float *wr, int nx, float dx, \
                                int nz, float dz, int lx0, int ntheta, float dtheta, float otheta);

void lx2theta_forward(int add, float *eimage, float *eimage_theta, \
             int nx, float dx, int nz, float dz, \
             int lx0, int ntheta, float dtheta, float otheta);

void lx2theta_adjoint(int add, float *eimage, float *eimage_theta, \
                      int nx, float dx, int nz, float dz, \
                      int lx0, int ntheta, float dtheta, float otheta);

void lx2theta_adjoint_singleTheta(int add, float *eimage, float *eimage_theta, \
                      int nx, float dx, int nz, float dz, \
                      int lx0, int ntheta, float dtheta, float otheta, int itheta);

float* cg_theta2lx(float *img_th, int nx_, int nz_, float dx_, float dz_, \
                   int lx0_, int ntheta_, float dtheta_, float otheta_, int niter);

float* cg_theta2lx_gpu(int gpuDev, float *img_th, int nx_, int nz_, float dx_, float dz_, \
                       int lx0_, int ntheta_, float dtheta_, float otheta_, int niter);

void cg_VHocig2Docig(int procGPU, float *Docig, float *Hocig, float *Vocig, long long int nx, long long int nz, float dx, float dz, \
                     long long int lx0, long long int lv0, long long int ld0, long long int ndip, float ddip, float dipMin, int niter, int iproc);

void getAp_VHocig2Docig(int procGPU, float *Aph, float *Apv, float *ph, float *pv, long long int nx, long long int nz, float dx, float dz, \
                        long long int lx0, long long int lv0, long long int ld0, long long int ndip, float ddip, float dipMin, int iproc);

void cg_tanTheta2theta(float *tdcig, float *adcig, long long int nx, long long int nz, float dx, float dz, float xMin, float xMax, float zMin, float zMax, \
                       long long int ltheta0, float dtheta, long long int ltanth0, float dtanth, int niter, int iproc);
                       

// Elastic Propagation
void propagate_Elastic(float *vx, float *vz, float *Sxx, float *Szz, float *Sxz, \
                       float *lap, float *lambda, float *mu, float *b, int nx, int nz, float dx, float dt);

void propagate_Elastic_abs(float *vx, float *vz, float *Sxx, float *Szz, float *Sxz, \
                       float *lap, float *lambda, float *mu, float *lambda2mu, float *bc, \
                       float *bd, int nx, int nz, float dx, float dt, float *sigX, float *sigZ);

void laplacian_80thOrder_Vx_Elastic_Iso(float *lap, float *b, float *Sdx, float *Sdz, int nx, int nz);
void laplacian_80thOrder_Vz_Elastic_Iso(float *lap, float *b, float *Sdx, float *Sdz, int nx, int nz);
void laplacian_80thOrder_Szz_Elastic_Iso(float *lap, float *lambda, float *lambda2mu, float *vx, float *vz, int nx, int nz);
void laplacian_80thOrder_Sxx_Elastic_Iso(float *lap, float *lambda, float *lambda2mu, float *vx, float *vz, int nx, int nz);
void laplacian_80thOrder_Sxz_Elastic_Iso(float *lap, float *mu, float *vx, float *vz, int nx, int nz);
void laplacian_40thOrder_Vx_Elastic_Iso(float *lap, float *b, float *Sdx, float *Sdz, int nx, int nz);
void laplacian_40thOrder_Vz_Elastic_Iso(float *lap, float *b, float *Sdx, float *Sdz, int nx, int nz);
void laplacian_40thOrder_Szz_Elastic_Iso(float *lap, float *lambda, float *lambda2mu, float *vx, float *vz, int nx, int nz);
void laplacian_40thOrder_Sxx_Elastic_Iso(float *lap, float *lambda, float *lambda2mu, float *vx, float *vz, int nx, int nz);
void laplacian_40thOrder_Sxz_Elastic_Iso(float *lap, float *mu, float *vx, float *vz, int nx, int nz);
void laplacian_20thOrder_Vx_Elastic_Iso(float *lap, float *b, float *Sdx, float *Sdz, int nx, int nz);
void laplacian_20thOrder_Vz_Elastic_Iso(float *lap, float *b, float *Sdx, float *Sdz, int nx, int nz);
void laplacian_20thOrder_Szz_Elastic_Iso(float *lap, float *lambda, float *lambda2mu, float *vx, float *vz, int nx, int nz);
void laplacian_20thOrder_Sxx_Elastic_Iso(float *lap, float *lambda, float *lambda2mu, float *vx, float *vz, int nx, int nz);
void laplacian_20thOrder_Sxz_Elastic_Iso(float *lap, float *mu, float *vx, float *vz, int nx, int nz);
void laplacian_10thOrder_Vx_Elastic_Iso(float *lap, float *b, float *Sdx, float *Sdz, int nx, int nz);
void laplacian_10thOrder_Vz_Elastic_Iso(float *lap, float *b, float *Sdx, float *Sdz, int nx, int nz);
void laplacian_10thOrder_Szz_Elastic_Iso(float *lap, float *lambda, float *lambda2mu, float *vx, float *vz, int nx, int nz);
void laplacian_10thOrder_Sxx_Elastic_Iso(float *lap, float *lambda, float *lambda2mu, float *vx, float *vz, int nx, int nz);
void laplacian_10thOrder_Sxz_Elastic_Iso(float *lap, float *mu, float *vx, float *vz, int nx, int nz);
void timeStep_Elastic(float *u, float *lap, int nx, int nz, float dt, float dx);
void timeStep_Elastic_abs(float *u, float *lap, int nx, int nz, float dt, float dx, float *sigX, float *sigZ);
void injectSource_elastic_Szz(float *szz, float *source, int it, int ix0, int iz0, int nx, int nz, float dx, float dz);
void injectSource_elastic_Omni(float *sxx, float *szz, float *source, int it, int ix0, int iz0, int nx, int nz, float dx, float dz);
void recordData_Elastic_Omni(Shot *shot, float *Sxx, float *Szz, int nx, int nz, int it);
void recordData_Elastic_phi(Shot *shot, float *vx, float *vz, int nx, int nz, float dx, float dz, int it);
void recordData_Elastic_Vertical(Shot *shot, float *vz, int nx, int nz, int it);

// Wavefield operations
void wavefieldTimeSecondDerivative(float *secDer, float coef_2, float *inc_2, float coef_1, float *inc_1, int nx, int nz, float dt);

// Born propagation
void scatter_Acoustic(float *ref, float *sct, float coef_2, float *inc_2, float coef_1, float *inc_1, int nx, int nz, float dt);
void scatter_Acoustic_lx(float *ref, float *sct, float coef_2, float *inc_2, float coef_1, float *inc_1, int nx, int nz, float dt, int lx0);
void scatter_Acoustic_lx_opt(float *ref, float *sct, float *inc, int nx, int nz, float dt, int lx0);
void scatter_Acoustic_AmpAngle(float *ref, float *sct, float coef_2, float *inc_2, float coef_1, float *inc_1, int nx, int nz, float dt, int lx0);


// Absorbing boundary
void   getSigmas(float *sigmaX, float *sigmaZ, int nb, int operLen, int nx, int nz, float dx, float dz);
void getVelSigma(float **sigma, float **sigmaInv, float *sigX, float *sigZ, int nx, int nz, float dx, float dt);

// Outputs
void   initExecPathName(char *inputExecPathName);
void   initFiles1d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1);
void   initFiles2d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2);
void   initFiles3d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2, int n3, float d3, float o3);
void   initFiles4d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2, int n3, float d3, float o3, \
                                                         int n4, float d4, float o4);
void   initFiles5d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2, int n3, float d3, float o3, \
                                                         int n4, float d4, float o4, int n5, float d5, float o5);

void   initFilesSmart1d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1);
void   initFilesSmart2d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2);
void   initFilesSmart3d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2, int n3, float d3, float o3);
void   initFilesSmart4d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2, int n3, float d3, float o3, \
                                                         int n4, float d4, float o4);
void   initFilesSmart5d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2, int n3, float d3, float o3, \
                                                         int n4, float d4, float o4, int n5, float d5, float o5);

void outputModel2d(char *name, float *inp, long long int n1, float d1, float o1, long long int n2, float d2, float o2, int n1b, int n2b);

void   initFilesSmart4d_J(char *name, FILE **fhdr, FILE **fbin, \
                          int n1, float d1, float o1, int j1, \
                          int n2, float d2, float o2, int j2, \
                          int n3, float d3, float o3, int j3, \
                          int n4, float d4, float o4, int j4);

void   createFiles(char *fileName, FILE **fBin, FILE **fHdr, int ndim, int *nsamp, float *dintv);
void   createFilesNew(char *fileName, FILE **fBin, FILE **fHdr, int ndim, int *nsamp, float *dintv, float *orig);
void   createFilesSmart(char *fileName, FILE **fBin, FILE **fHdr, int ndim, int *nsamp, float *dintv, float *orig);
FILE*  openFile(int *err, char *fileName, char *mode);
void   fillHeaderFile(FILE *fHdr, char *fileNameBin, int ndim, int *nsamp, float *dintv);
void   fillHeaderFileNew(FILE *fHdr, char *fileNameBin, int ndim, int *nsamp, float *dintv, float *orig);
void   writeSamples(const float *array, FILE *fBin, size_t ns);
int    fileStrCat(char *sout, const char *s1, const char *s2);
void   output1d(char *name, float *array, int n1, float d1, float o1);
void   output2d(char *name, float *array, int n1, float d1, float o1, int n2, float d2, float o2);
void   output3d(char *name, float *array, int n1, float d1, float o1, int n2, float d2, float o2, int n3, float d3, float o3);
void   outputSmart1d(char *name, float *array, long long int n1, float d1, float o1);
void   outputSmart2d(char *name, float *array, long long int n1, float d1, float o1, long long int n2, float d2, float o2);
void   outputSmart3d(char *name, float *array, long long int n1, float d1, float o1, long long int n2, float d2, float o2, long long int n3, float d3, float o3);
void   outputSmart4d(char *name, float *array, long long int n1, float d1, float o1, long long int n2, float d2, float o2, \
                     long long int n3, float d3, float o3, long long int n4, float d4, float o4);
void   outputSmart5d(char *name, float *array, long long int n1, float d1, float o1, long long int n2, float d2, float o2, \
                     long long int n3, float d3, float o3, long long int n4, float d4, float o4, long long int n5, float d5, float o5);
void  outputSmart4d_J(char *name, float *array, \
                      long long int n1, float d1, float o1, long long int j1, \
                      long long int n2, float d2, float o2, long long int j2, \
                      long long int n3, float d3, float o3, long long int j3, \
                      long long int n4, float d4, float o4, long long int j4);

// Inputs
void readImageFile(char *filename, float *image, int n);
void readShot(char *seismogramFileName, char *sourceFileName, Shot *shot);
float* readModel(char *arqName, int nx, int nz);

// SEGY reading
void initDatasetFromSEGY(FILE *fp, char *filename, Dataset *dataset);
int indexShotGathersFromSEGY(FILE *fp, char *filename, int **shotIndexes, int **tracesPerShot, int *nt, float *dt, int *ntraces);

// Inversion
float* cg(float *image, Dataset *dataset, float *vp_bg, int jt, int nt, int nx, int nz, float dx, \
          float dz, float dt, float ox, float oz, float ot, int lx0, int niter);
float* cg_psf(float *image, float *image_psf, int **positionPSF, int nx_, int nz_, float dx_, \
                 float dz_, float dt_, float ox_, float oz_, int niter, float epsilon);
int** definePositionPSF(int nx, int nz);
void createPointScatterersModel(float *m, float dv, int **positionPSF, int nx, int nz);
void addWhiteNoise2EImagePSF(float *eimage_psf, int nx, int nz, int lx0, float epsilon);

// Slicing offset/angle planes
float *pullOffsetPlane(float *eimage, int nx, int nz, int lx0, int lx);
void   pushOffsetPlane(float *eimage, float *plane, int nx, int nz, int lx0, int lx);
float *pullThetaPlane(float *eimage_theta, int nx, int nz, int ltheta0, int itheta);
void   pushThetaPlane(float *eimage_theta, float *thetaPlane, int nx, int nz, int ltheta0, int itheta);

// Operations
void accumulate(float *out, const float *inp, size_t n);

// Analytic Trace
void hilbert(float *trace, int nt, float dt);

void instantaneousPhase(float *trace, int nt, float dt);

float* transform2InstantaneousPhase(float *traces, int nt, float dt, int ntrace);

void envelope(float *trace, int nt, float dt);

float* transform2Envelope(float *traces, int nt, float dt, int ntrace);

// Math
float sinc(float x);

float mBFZ(float x);

float kaiser(float x, float r, float b);

void applyHalfDerivative(float *inp, int nt, float dt, float sign);

void bandpass(float *inp, int nt, float dt, float f1, float f2, float f3, float f4);

void modelSmooth_Target(float *out, float *inp, int nx, int nz, int rx, int rz, int min_z, int max_z);

void modelSmooth(float *out, float *inp, int nx, int nz, int rx, int rz);

double getArrayL2Norm(float *array, long long int n);

void normalizeArray(float *array, long long int n);

void normalizeArraysJointly(float *array1, long long int n1, float *array2, long long int n2);

void addArrays(float *out, float alpha1, float *inp1, float alpha2, float *inp2, long long int n);

void divideArraysPow(float *out, float pow1, float *inp1, float pow2, float *inp2, long long int n);

void multiplyArrays(float *out, float *inp1, float *inp2, long long int n);

void takeSquareRootArrays(float *out, float *inp, long long int n);

void takeAbsluteValueArrays(float *out, float *inp, long long int n);

void multiplyArrayByScalar(float alpha, float *out, long long int n);

void divideArrays(float *out, float *inp1, float *inp2, long long int n);

void divideArraysWhiteNoise(float *out, float *num, float *den, long long int n, float noiseLevel);

double computeL2NormArray(float *array, long long int n);

float getVolumeRMS(float *array, int n);

void applyLinearWeight_TPCIG(float *tpcig_dso, float *tpcig, float dp, int nx, int ntau, int lp0);

void applyLinearWeight_ODCIG(float *eimage_dso, float *eimage, float dx, int nx, int nz, int lx0);

void applyLinearDSO_ODCIG(float *eimage_dso, float *eimage, float dx, int nx, int nz, int lx0, float linear_coefficient);

void applyAngleDerivative(float *out, float *eimage_theta_phase, int nx, int nz, int ntheta, float dtheta);

void getThetaTrace(float *eimage_theta, float *thetaTrace, int nx, int nz, int ntheta, int itheta, int ix);

int  interpTraceLinear(float *inp, float *out, int nin, float din, int nou, float dou);

void smoothTrace(float *gather, int ntrc, int nsmp, int nsmt);

void transp(float *ou, float *in, int n1, int n2);

void getFlatADCIG(float *out, float *eimage_theta, int nx, int nz, int ntheta, float dz, float dtheta, float oz, float theta_max);

void taperADCIG(float *out, float *eimage_theta, int nx, int nz, int ntheta, float dtheta, float thetaIni, float thetaFin, float power);

void makeReciprocalADCIG(float *out, float *inp, int nx, int nz, int nxb, int nzb, int ntheta, int acqGeom);

void applyCosADCIG(float *out, float *eimage_theta, int nx, int nz, int ntheta, float dtheta);

void flatSpikeADCIG(float *out, int nx, int nz, float dz, int ntheta, float dtheta, float z);

float* getADCIGDipGrad(float *eimage_theta, int nx, int nz, int ntheta, float dx, float dtheta);

void ODCIG_localSlantStack(int adj, float *eimage_localSlantStack, float *eimage, \
                           int nx, float dx, int nz, float dz, int lx0, int ntheta, float dtheta);

void ODCIG_MoveIn_LocalSlantStack(float *eimage_SlantStack, float *eimage_SlantStack_deMO, \
                                int nx, float dx, int nz, float dz, int lx0, \
                                float contraction_factor, float lx_null, int ntheta, float dtheta);

void lxMoveIn(int add, float *eimage, float *eimage_movein, \
              int nx, float dx, int nz, float dz, int lx0, \
              int ntheta, float dtheta, \
              float contraction_factor, float lx_null);    

void lx2DeMoveout_forward(int add, float *eimage, float *eimage_deMO, \
                          int nx, float dx, int nz, float dz, int lx0, \
                          float xMin, float xMax, float zMin, float zMax,
                          float contraction_factor, float lx_null);

void lx2DeMoveout_forward2(int add, float *eimage, float *eimage_deMO, \
                          int nx, float dx, int nz, float dz, int lx0, \
                          float contraction_factor, float lx_null);

void lx2DeMoveout_forward3(int add, float *eimage, float *eimage_deMO, \
                          int nx, float dx, int nz, float dz, int lx0, \
                          float contraction_factor, float lx1, float lx2, float lx3, float lx4);

void lx2Contraction_forward(int add, float *eimage, float *eimage_contracted, \
                            int nx, float dx, int nz, int lx0, float contraction_factor);

// void modelTaperingX(float *model, float *model_tapered, \
//                     int nx, int nz, float dx, float x0, float x1);

void modelTaperingZ(float *model, float *model_tapered, \
                    int nx, int nz, float dz, float z0,
                    float z1);

void modelTaperingX(float *model, float *model_tapered, \
                    int nx, int nz, float dx, float x1, \
                    float x2, float x3, float x4);

void lxTapering(float *eimage, float *eimage_tapered, int nx, int nz, float dz, int lx0, int lz0, float z0, float z1);

void lxTapering_x(float *eimage, float *eimage_tapered, \
                         int nx, int nz, float dx, int lx0, \
                         float x1, float x2, float x3, float x4);

void applyIlum2EImage(float *eimage, float *eimage_ilumin, \
                      float *image, float *image_ilumin, float *ilumin, \
                      int nx, int nz, float dx, int lx0, float epsilon);
void lxTapering_lx(float *eimage, float *eimage_tapered, \
                   int nx, int nz, float dx, int lx0, \
                   float lx1, float lx2, float lx3, float lx4);

void lxNullTapering(float *eimage, float *eimage_tapered, \
                    int nx, int nz, float dx, int lx0, float lxNull);

void lxTapering_lx_varXZ(float *eimage, float *eimage_tapered, \
                         int nx, int nz, float dx, int lx0, \
                         TapeLX *tapelx, int ntape);               

void lxTapering_lx_varDepth(float *eimage, float *eimage_tapered, \
                            int nx, int nz, float dx, int lx0, \
                            float lx1_i, float lx2_i, float lx3_i, float lx4_i, \
                            float lx1_f, float lx2_f, float lx3_f, float lx4_f, \
                            float z_i, float z_f);

void lxTapering_MultiDim(float *eimage, float *eimage_tapered, \
                         int nx, int nz, float dx, float dz, int lx0, \
                         float z0, float z1, \
                         float x0, float x1);

void lxApplyBandpass(float *eimage, int nx, int nz, float dz, int lx0, float lambda1, float lambda2, float lambda3, float lambda4);

void stacking_eimage_threads(int nthreads, float *eimage_input, float *eimage_reduced, int nx, int nz, int lx0);

void cpyAngle(float *eimage_theta_temp, float *eimage_theta, int nx, int nz, int ntheta, int itheta);
void cpyAngleODCIG(float *eimage_temp, float *eimage_ampAngle, int nx, int nz, int lx0, int itheta);

float  fsignf(float x);
double fsign(double x);
int     sign(int x);

void transposeArray(float *array, int n1, int n2);

float findParabolaMin(float x0, float x1, float x2, float y0, float y1, float y2, float *ymin);
void getParabolaParams(float *a, float *b, float *c, double x0, double x1, double x2, double y0, double y1, double y2);

float* convolution(float *a, float *b, int na, int nb);
float* correlation(float *a, float *b, int na, int nb);
float* codeGen(int n1, int j1, int ncodes, int seed);

void resamp(float *trace, float *trace_resamp, int nt, float dt, int nt_resamp, float dt_resamp);

// Functions to transform to geologic dip
float* VHocigs2Gocig_old(float *Hocig, float *Vocig, int nx, int nz, float dx, float dz, \
                     int lx0, int lv0, float dipMin, float ddip, float dipMax, float lambdaMax, \
                     int *lg0_out, float *dlg_out, int *nDip_out);
void VHocig2Docig(float *Docig, float *Hocig, float *Vocig, \
                  long long int nx, long long int nz, long long int ndip, \
                  float dx, float dz, float ddip, float dipMin, \
                  long long int lx0, long long int lv0, long long int ld0, int iproc);
void Docig2VHocig(float *Docig, float *Hocig, float *Vocig, \
                  long long int nx, long long int nz, long long int ndip, \
                  float dx, float dz, float ddip, float dipMin, \
                  long long int lx0, long long int lv0, long long int ld0, int iproc);                  
void VHocig2Docig_old(float *Docig, float *Hocig, float *Vocig, float dip, float dipRange, \
                    long long int nx, long long int nz, float dx, float dz, \
                    long long int lx0, long long int lv0, long long int ld0, int iproc);
void Docig2VHocig_old(float *Docig, float *Hocig, float *Vocig, float dip, float dipRange, \
                      long long int nx, long long int nz, float dx, float dz, \
                      long long int lx0, long long int lv0, long long int ld0, int iproc);
float* VHocig2Gocig(float *Hocig, float *Vocig, long long int nx, long long int nz, float dx, float dz, \
                    long long int lx0, long long int lv0, long long int lz0, int iproc);
void Gocig2VHocig(float *Gocig, float *Hocig, float *Vocig, long long int nx, long long int nz, float dx, float dz, \
                  long long int lx0, long long int lv0, long long int lz0, int iproc);
float* VHocig2Gocig_simple(float **Hocig_new, float *Hocig, float *Vocig, int nx, int nz, float dx, float dz, int lx0, int lv0, int lz0, int iproc);
void eimage_NonGeologicDipFiltering(float *eimage, int nx, int nz, float dx, float dz, int lx0, int lz0);
void VHocig_dispersalCorrection(float *Hocig, float *Vocig, int nx, int nz, float dx, float dz, int lx0, int lv0);

// Functions to perform matrix transposition
float* transp_dim3(float *m_input, long long int n1, long long int n2, long long int n3, int plane);
float* transp_dim4(float *m_input, long long int n1, long long int n2, long long int n3, long long int n4, int plane);

// Functions to apply FFTW2D and FFTW3D
void applyFFTW2D_forwardR2I(fftw_complex *out, float *inp, int n1, int n2);
void applyFFTW2D_backwardI2R(float *out, fftw_complex *inp, long long int n1, long long int n2);
fftw_complex *applyFFTW3D_forwardR2I(float *inp, int n1, int n2, int n3);
float *applyFFTW3D_backwardI2R(fftw_complex *inp, int n1, int n2, int n3);

// Taup Transform
void dsoTauP(float *adcigtp, float *adcigtp_dso, int nx, int np, int ntau, float dp, float dtau);

void ADCIG2TauP_forward(float *adcigtp, float *eimage_theta, int nx, int nz, int ntheta, int np, int ntau, float dx, float dtheta, float dp, float dtau);

void ADCIG2TauP_adjoint(float *adcigtp, float *eimage_theta, int nx, int nz, int ntheta, int np, int ntau, float dx, float dtheta, float dp, float dtau);

float* cg_taup2adcig(float *taup, int nx_, int nz_, float dx_, float dz_, \
                     int np, int ntau, int ntheta_, float dp, float dtau, \
                     float dtheta_, int niter);

// Focusing operator
void FocusingOperator_lx_lz_lt(int adj, float *eimage_diff, float *eimage_contracted, float *eimage_contracted_ref, float *eimage, \
                               float zMin, float zMax, \
                               float maxDepth, long long int nx, long long int nz, float dx, float dz, float ox, float oz, \
                               int lx0, int lz0, int lt0, int ltheta0, float dtheta, float otheta, \
                               float dipMin, float dipMax, float ddip, float lambdaMax, float dlambda, int applyADFocusing, \
                               int iter, float mu, float epsFocus, float lambdaNull, \
                               int outputs, int iproc, int nproc);

void FocusingOperator_AD(int adj, float *eimage_diff, float *eimage_contracted, float *eimage_contracted_ref, float *eimage, \
                         float zMin, float zMax, float maxDepth, \
                         long long int nx, long long int nz, float dx, float dz, float ox, float oz, \
                         int lx0, int ltheta0, float dtheta, float otheta, \
                         int iter, float mu, float epsFocus, float lambdaNull, int outputs, int iproc, int nproc);

void FocusingOperator_Gocig(int adj, float *eimage_diff, float *eimage_contracted, float *eimage_contracted_ref, float *eimage, \
                            float zMin, float zMax, float maxDepth, \
                            long long int nx, long long int nz, float dx, float dz, float ox, float oz, long long int nxb, long long int nzb, \
                            long long int lx0, long long int lz0, long long int ltheta0, float dtheta, float otheta, \
                            float lambdaMax, float dipMin, float dipMax, float ddip, \
                            int iter, float muFocus, float epsFocus, float lambdaNull, \
                            int outputs, int iproc, int nproc);

void FocusingOperator_Tcig(int adj, float *Tcig_foc, float *Tcig_ref, float *Tcig, float zMin, float zMax, \
                           long long int nx, long long int nz, float dx, float dz, float ox, float oz, \
                           long long int nxb, long long int nzb, long long int tau0, \
                           float muFocus, float tauNull, int outputs, int iproc, int nproc);

void FocusingOperator_Docig(int adj, float *eimage_diff, float *eimage_contracted, float *eimage_contracted_ref, float *eimage, \
                            float zMin, float zMax, long long int nx, long long int nz, float dx, float dz, float ox, float oz, \
                            long long int nxb, long long int nzb, long long int lx0, long long int ndip, float dipMin, float ddip, \
                            float muFocus, float epsFocus, float lambdaNull, int outputs, int iproc, int nproc);

void FocusingOperator_TL_Docig(int adj, float *eimage_contracted, float *eimage, \
                               float zMin, float zMax, long long int nx, long long int nz, float dx, float dz, float ox, float oz, \
                               long long int nxb, long long int nzb, long long int lx0, long long int ndip, float dipMin, float ddip, \
                               float muFocus, float epsFocus, float lambdaNull, int outputs, int iproc, int nproc);

void FocusingOperator_lambda(int adj, float *eimage_diff, float *eimage_contracted, float *eimage_contracted_ref, float *eimage, \
                             float zMin, float zMax, float maxDepth, long long int nx, long long int nz, float dx, float dz, float ox, float oz, \
                             int nxb, int nzb, int lx0, int ltheta0, float dtheta, float otheta, \
                             int iter, float muFocus, float epsFocus, float lambdaNull, \
                             int outputs, int iproc, int nproc);

void FocusingOperator_simple_lx(int adj, float *eimage_contracted, float *eimage_input, float zMin, float zMax, \
                                 float maxDepth, long long int nx, long long int nz, float dx, float dz, float ox, float oz, \
                                 int lx0, int ltheta0, float dtheta, float otheta, int iter, float mu, float eps, int iproc, int nproc);


// ModMig
void BornDataModeling(Dataset *dataset, float *vp, float *ref, MSH *msh, \
                      int acquisitionGeom, float offMax, float perc, \
                      float agc_time_window, int iproc, int nproc);

void AcousticDataModeling(Dataset3D *dataset, Dataset3D *dataset_da, float *vp, float *vp_da, \
                          MSH *msh, Parms *parms, int iproc, int nproc, int readInputData);

void migration(Dataset *dataset, float *image, float *eimage, float *vp_bg,\
               int nx, int nz, int lx0, int lz0, int lt0, int jt_lt, int ltheta0, \
               float dx, float dz, float dtheta, \
               float ox, float oz, float otheta, \
               int jt, int ism, int recordMovie, int iproc, int nproc);

void migration_TLCIG(Dataset *dataset, float *image, float *TLcig, float *vp_bg,\
                     long long int nx, long long int nz, \
                     float dx, float dz, float ox, float oz, long long int lt0, \
                     long long int jt, int ism, int recordMovie, int iproc, int nproc);

void migration_TAU(Dataset3D *dataset, float *image, float *ilumin, float *Tcig, float *vp_bg,\
                   MSH *msh, Eimage_MSH *eimsh, Parms *parms, int ism, int recordMovie, int iproc, int nproc);

void migration_HOCIG_VOCIG(Dataset3D *dataset, float *image, float *Hocig, float *Vocig, float *vp_bg,\
                           MSH *msh, Eimage_MSH *eimsh, Parms *parms, int ism, int recordMovie, int iproc, int nproc);

double extendedBornModeling(Dataset *dataset, float *eimage, float *Hocig, float *Vocig, float *vp_bg, \
                            long long int nx, long long int nz, long long int lx0, long long int lz0, long long int lt0, \
                            float dx, float dz, float ox, float oz, \
                            long long int jt, float perc, float offMax, int acquisitionGeom, \
                            char *data_FileName, int iproc, int nproc);

void ExtendedBornModeling_TAU(Dataset3D *dataset, float *vp_bg, float *Tcig, MSH *msh, Eimage_MSH *eimsh, \
                              Parms *parms, int ism, int recordMovie, int iproc, int nproc);

void BornGradient_DataDomain_TAU(Dataset3D *dataset, float *gradient, float *ilumin, float *vp_bg, \
                                 float *Tcig, GRDProc *grdProc, MSH *msh, Eimage_MSH *eimsh, \
                                 Parms *parms, int ism, int recordMovie, int iter, int iproc, int nproc);

void BornGradient_DataDomain(Dataset3D *dataset, float *gradient, float *vp_bg, \
                             float *eimage, float *Heref, float *Veref, GRDProc *grdProc, MSH *msh, Eimage_MSH *eimsh, \
                             Parms *parms, int ism, int recordMovie, int iter, int iproc, int nproc);

// void BornGradient_DataDomain(Dataset3D *dataset, float *gradient, float *vp_bg, \
                             float *eimage, float *Heref, float *Veref, GRDProc *grdProc, \
                             long long int nx, long long int nz, long long int lx0, long long int lz0, long long int lt0, \
                             float dx, float dz, float ox, float oz, \
                             long long int jt, long long int jt_EIC, int iter, int imageDomain, \
                             int iproc, int nproc);

void BornGradient_ExtExpRef(float *gradient, float *vp_bg, GRDProc *grdProc, \
                            float *Hocig_background_phaseShift, float *Hocig_background_original, \
                            float *Hocig_residual_phaseShift, float *Hocig_residual_original, \
                            float *Vocig_background_phaseShift, float *Vocig_background_original, \
                            float *Vocig_residual_phaseShift, float *Vocig_residual_original, \
                            long long int lx0, long long int lz0, \
                            long long int nx, long long int nz, long long int nt, \
                            float dx, float dz, float dt, float ox, float oz, float ot, \
                            long long int jt, int iter, int iproc, int nproc);

void BornGradient_ExtExpRef_Docig(float *gradient, float *vp_bg, GRDProc *grdProc, \
                                  float *Docig_background_phaseShift, float *Docig_background_original, \
                                  float *Docig_residual_phaseShift, float *Docig_residual_original, \
                                  long long int lambda0, long long int ndip, float dipMin, float ddip, \
                                  long long int nx, long long int nz, long long int nt, \
                                  float dx, float dz, float dt, float ox, float oz, float ot, \
                                  long long int jt, int iter, int iproc, int nproc);

// Wrap
float computeResiduals(Dataset *dataset_input, Dataset *dataset_residuals, Dataset *dataset_reference, Dataset *dataset_misfit, \
                       float *eimage_contracted, float *eimage_contracted_ref, float *eimage_diff_diff, float *eimage_residuals_contracted, 
                       MSH *msh, float *vp_bg, EIProc *eiproc, Eimage_MSH *eimage_msh, float offMax, int acquisitionGeom, float perc, \
                       float zMin, float zMax, float maxDepth, float muFocus, float epsFocus, float lambdaNull, int iter, int ism, int outputs, \
                       int imageDomain, int iproc, int nproc);

float computeResiduals_AD(Dataset *dataset_input, Dataset *dataset_residuals, Dataset *dataset_reference, Dataset *dataset_misfit, \
                          float *Hocig_contracted, float *Hocig_contracted_ref, float *Vocig_contracted, float *Vocig_contracted_ref, \
                          float *Hocig_diff, float *Vocig_diff, \
                          MSH *msh, float *vp_bg, EIProc *eiproc, Eimage_MSH *eimage_msh, float offMax, int acquisitionGeom, float perc, \
                          float zMin, float zMax, float maxDepth, float muFocus, float epsFocus, float lambdaNull, int iter, int ism, int outputs, \
                          int imageDomain, int computingGradient, int iproc, int nproc);

float computeResiduals_ORTH(Dataset *dataset_input, Dataset *dataset_residuals, Dataset *dataset_reference, Dataset *dataset_misfit, \
                            float *Hocig_contracted, float *Hocig_contracted_ref, float *Hocig_normalizer, \
                            float *Vocig_contracted, float *Vocig_contracted_ref, float *Vocig_normalizer, \
                            float *Hocig_diff, float *Vocig_diff, \
                            MSH *msh, float *vp_bg, EIProc *eiproc, Eimage_MSH *eimage_msh, float offMax, int acquisitionGeom, float perc, \
                            float zMin, float zMax, float maxDepth, float muFocus, float epsFocus, float lambdaNull, int iter, int ism, int outputs, \
                            int imageDomain, int computingGradient, int iproc, int nproc);

float computeResiduals_ORTH_PopCIG(Dataset *dataset_input, \
                                   float *Hocig_background_phaseShift, float *Hocig_background_original, \
                                   float *Hocig_residual_phaseShift, float *Hocig_residual_original, \
                                   float *Vocig_background_phaseShift, float *Vocig_background_original, \
                                   float *Vocig_residual_phaseShift, float *Vocig_residual_original, \
                                   MSH *msh, float *vp_bg, EIProc *eiproc, Eimage_MSH *eimage_msh, \
                                   float offMax, int acquisitionGeom, float perc, float zMin, float zMax, float maxDepth, \
                                   float muFocus, float epsFocus, float lambdaNull, int iter, int ism, int outputs, \
                                   int computingGradient, int focusingQC, int iproc, int nproc);

float computeResiduals_ORTH_PopDocig(Dataset *dataset_input, \
                                     float *Docig_background_phaseShift, float *Docig_background_original, \
                                     float *Docig_residual_phaseShift, float *Docig_residual_original, \
                                     MSH *msh, float *vp_bg, EIProc *eiproc, Eimage_MSH *eimage_msh, \
                                     float offMax, int acquisitionGeom, float perc, float zMin, float zMax, float maxDepth, \
                                     float muFocus, float epsFocus, float lambdaNull, int iter, int ism, int outputs, \
                                     int computingGradient, int focusingQC, int iproc, int nproc);

float computeObjFunc(Dataset *dataset_misfit, Dataset *dataset_reference, float *eimage_contracted, float *vp_bg, \
                     int nx, int nz, int nt, int lx0, int lz0, int lt0, int jt_lt, float dx, float dz, float dt, float ox, float oz, \
                     float ot, float offMax, int acquisitionGeom, float perc, int jt, int imageDomain, int iproc, int nproc, int iter);

float computeODCIGs(Dataset *dataset_input, Dataset *dataset_reference, float *eimage_contracted, float *eimage_contracted_ref, \
                    float *eimage_diff, float *vp_bg, \
                    int nx, int nz, int nt, int lx0, int lz0, int lt0, int jt_lt, int ltheta0, \
                    EIProc *eiproc, Eimage_MSH *eimage_msh, float dx, float dz, float dt, float dtheta, float ox, float oz, \
                    float ot, float otheta, float offMax, int acquisitionGeom, float perc, \
                    float zMin, float zMax, float maxDepth, float muFocus, float eps, float lambdaNull, int iter, int ism, int jt, int outputs, \
                    int iproc, int nproc);


// Filtering ECIG
void normalizeHocigByEnvelope(float *Hocig, float *Hocig_normalizer, long long int lx0, long long int nx, long long int nz, float dx);

void normalizeVocigByEnvelope(float *Vocig, float *Vocig_normalizer, long long int lv0, long long int nx, long long int nz, float dz);

void filterTCIG(float *eimage, EIProc *eiproc, Eimage_MSH *emsh, MSH *msh);

void filterECIG(float *eimage, EIProc *eiproc, Eimage_MSH *emsh, MSH *msh);

float* filterECIG_Angle(float *eimage, EIProc *eiproc, \
                        int nx, int nz, int nxb, int nzb, float dx, float dz, \
                        float ox, float oz, int ltheta0, float dtheta, int lx0, \
                        int iproc, int procGPU);

void reciprocityAngle(float *eimage, EIProc *eiproc, int acqGeom, \
                        int nx, int nz, int nxb, int nzb, \
                        float dx, float dz, float ox, float oz, \
                        int ltheta0, float dtheta, int lx0, \
                        int iproc, int procGPU);

// Log
void createFileLog();
void openFileLog(FILE **fp);
void log_LineSearch(float alpha0, float alpha1, float alpha2, float objFunc0,  float objFunc1,  float objFunc2, float alphaMin, int iter);
void log_InitIter(int iter, time_t time_now, time_t time_begin);
void log_endIter(int iter, time_t time_now, time_t time_begin, time_t time_begin_iter);



// CIG processing
void get_HVocig_ExtExpRef_fromDocig(float *Docig_background, float *Docig_residual, \
                                    float *Docig_background_phaseShift, float *Docig_background_original, \
                                    float *Docig_residual_phaseShift, float *Docig_residual_original, \
                                    long long int lx0, long long int ndip, float dipMin, float ddip, \
                                    long long int nx, long long int nz, long long int nt, \
                                    float dx, float dz, float dt, float ox, float oz, float ot, int iproc, int nproc, float *vp_bg);

void get_HVocig_ExtExpRef(float *Hocig_background, float *Vocig_background, \
                          float *Hocig_residual,   float *Vocig_residual, \
                          float *Hocig_background_phaseShift, float *Hocig_background_original, \
                          float *Hocig_residual_phaseShift, float *Hocig_residual_original, \
                          float *Vocig_background_phaseShift, float *Vocig_background_original, \
                          float *Vocig_residual_phaseShift, float *Vocig_residual_original, \
                          long long int lx0, long long int lz0, long long int nx, long long int nz, long long int nt, \
                          float dx, float dz, float dt, float ox, float oz, float ot, int iproc, int nproc, float *vp_bg);


void initializeFileNames(int iter);

void initializeTLFileNames(int iter);



