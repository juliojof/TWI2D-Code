#ifndef SHOT_3D_H
#define SHOT_3D_H

// #include "Sensor.hpp"
#include "DataTypes.h"
// #include "parmHandling.h"

#define SOURCETYPE_ZEROPHASE          0
#define SOURCETYPE_MINIMUMPHASE       1
#define SOURCETYPE_FROM_DIRECTARRIVAL 2

#define FILE_POINTERS_INPUT       0
#define FILE_POINTERS_OUTPUT      1

#define MODE_CONVOLUTION_DECONVOLVE 0
#define MODE_CONVOLUTION_CONVOLVE 1

// Data types

typedef struct AcquisitionGeom_
{
    int nshot, nrec;

    float shotOX, shotOY, shotOZ;
    float shotDX, shotDY, shotDZ;
    int   shotNX, shotNY, shotNZ;

    float recOX, recOY, recOZ;
    float recDX, recDY, recDZ;
    int   recNX, recNY, recNZ;

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






// Write/Read dataset to/from file
void saveDataset3D(char *fileName, Dataset3D *dataset);
void saveDataset3D_metaData(char fileName[], char fileNameBinary[], Dataset3D *dataset, int filePointers, int saveNewBinaryFile);
void saveDataset3D_binaryData(Dataset3D *dataset, int filePointers, int ishot);
void closeOdisseiaFiles(Dataset3D *dataset, int filePointers);
void initDataset3DfromFile(Dataset3D *dataset, Parms *parms);
void initDataset3DfromSEGY_original(Dataset3D *dataset, int nshot, int nshot_out, int shot_step, \
                           int *nrec, int *nsou, int *gatherIni, int *gatherFin, \
                           int nt, float dt, float ot, float maxFreq, \
                           float *xs, float *ys, float *zs, \
                           float *xr, float *yr, float *zr);
void initDataset3DfromSEGY(Dataset3D *dataset, int nshot, int nshot_out, int shot_step, \
                           int *nrec, int *nsou, int **gatherTraces, \
                           int nt, float dt, float ot, float maxFreq, \
                           float *xs, float *ys, float *zs, \
                           float *xr, float *yr, float *zr);

void initDataset3D_PZSum(Dataset3D *dataset_P, Dataset3D *dataset_Z, Dataset3D *dataset_PZ);
void PZSum(Dataset3D *dataset_P, Dataset3D *dataset_Z, Dataset3D *dataset_PZ);
void PZSum_shot(Shot3D *shot_P, Shot3D *shot_Z, Shot3D *shot_PZ, float factor_P, float factor_Z, float *source);
void PZSum_sensor(Shot3D *shot_P, Shot3D *shot_Z, Shot3D *shot_PZ, float *amplitude_P, float *amplitude_Z, \
                  long long int *trace2Station, long long int *traceNumber, float *source);

// Sensor* mapSensors(Dataset3D *dataset, long long int *nsensors, long long int **trace2Station);
void getShotRecNum(Dataset3D *dataset, long long int traceNum, long long int *ishot_, long long int *irec_);
// float*  get_SignatureFromDA_Sensor(Dataset3D *dataset, Sensor *sensor, long long int nsensors, long long int *trace2Station, Parms *parms, int applyCosineFactor);

float*  get_SignatureFromDA(Dataset3D *dataset, Parms *parms);
void Take1stDerivative_singleTrace_TimeDomain(float *trace_time, int nt, float dt, float ot);
void ApplyTaperToTrace(float *trace_time, int nt, float dt, float ot, float t_ini, float t_fin, float taper_length);
float getMaximumAmplitude_singleTrace(float *trace_time, int nt);
void ZeroOutNoisyTrace(Shot3D *shot);
void Take1stDerivative_shot_TimeDomain(Shot3D *shot);
float Take1stDerivative_singleTrace(float *trace_time, int nt, float dt, float ot);
void Take1stDerivative_shot(Shot3D *shot, float multiplyingFactor);
void Deconvolve_shot(Shot3D *shot, float *trace_time, int mode_convolution, float white_noise_perc);
void Deconvolve_Data(Dataset3D *dataset, float *trace_time, int mode_convolution, float white_noise_perc);
long long int buildSourceFromDirectArrival(float *source, float *seismogram, int nt, float dt, float ot, int ntrc, \
                                           float *xs, float *ys, float *zs, float *xr, float *yr, float *zr, float offMin, float offMax, \
                                           float srcWin_ti, float srcWin_ti_taper, float srcWin_tf, float srcWin_tf_taper, float velWater);

float* timeShiftDA(float *trace_inp, int nt, float dt, float ot, float xs, float ys, float zs, \
                   float xr, float yr, float zr, float offMin, float offMax, \
                   float srcWin_ti, float srcWin_ti_taper, float srcWin_tf, float srcWin_tf_taper, float velWater, int applyCosineFactor);

void finalizeDataset3DfromSEGY(Dataset3D *dataset, int shot_step);
int initShot3DfromFile(Dataset3D *dataset, int ishot);
void initSeismogramAndSourceMemory(Dataset3D *dataset, int ishot);
void initSeismogramMemory(Dataset3D *dataset, int ishot);
void initSourceMemory(Dataset3D *dataset, int ishot);
void loadSeismogram2Dataset(Dataset3D *dataset, int ishot, int normalizeInputDataAmplitudes);
int finalizeShot3DfromFile(Dataset3D *dataset, int ishot);
int finalizeSeismogramAndSource(Dataset3D *dataset, int ishot);
int finalizeShot3DMetaInfo(Shot3D *shot);
int finalizeShotSeismogramAndSource(Shot3D *shot);

int getShotWithRestrictedAzimuths(Shot3D *shot_out, Shot3D *shot_inp, double maxAzimDegrees);
double findAzimuth(float xs, float ys, float xr, float yr);

void ZeroUpOrDown(Shot3D *shot, int option);
void UpOrDown_Coefficients(Shot3D *shot, float coeff_up, float coeff_down);
void pickReceiver(Shot3D *shot, int irec_pick);

void initDataset3D(Dataset3D *dataset, MSH *msh, Parms *parms);

void initShotFilteringParms(Dataset3D *dataset, Parms *parms);

void initGeomParms(AcquisitionGeom *acqGeom, Parms *parms);

Shot3D* initCopyShot3D(Dataset3D *dataset_inp, int ishot);

void initShot3D(Dataset3D *dataset, MSH *msh);

void initCoordinatesAndWavelet(Dataset3D *dataset, MSH *msh, AcquisitionGeom *acqGeom);

void finalizeDataset3D(Dataset3D *dataset);

void integrateRecTime3D(Shot3D *shot);

void AGC3D(Dataset3D *dataset, float time_window);

void applyRecTaper3D(Shot3D *shot, float perc, float offMax, int sensorSide);

void finalizeShot3D(Shot3D *shot);

double addObjFunc3D(Shot3D *shot);

void removeDA3D(Shot3D *shot, Shot3D *shot_da);

void copyData3D(Shot3D *shot_inp, Shot3D *shot_out);

void applyHalfDerivative2Source3D(Shot3D *shot);

void applyHalfDerivative2Data3D(Shot3D *shot);

float *getBandpassPulse3D(int nt, float dt, float f1, float f2, float f3, float f4, float t0);

void applyBandpassToTrace3D(Shot3D *shot, float f1, float f2, float f3, float f4);
void applyBandpassToSourceAndSeismogram(Shot3D *shot, float f1, float f2, float f3, float f4);

void multiplySeismogramByTime3D(Shot3D *shot, float power);

void multiplyTraceByTime3D(float *trace, int nt, float dt, float power);

void demultiplexSeismogram3D(float *seismogram, int nt, int nrec);

void makeEncodedDataset3D(int idev, Dataset3D *data_inp, Dataset3D *data_tmp, Dataset3D *data_enc, Dataset3D *data_res, int codelength, int codestep, int realizations);

void windowData3D(int idev, Dataset3D *dataset, float tmin, float tmax, float dt);

void copyDataset3D(Dataset3D *dst, Dataset3D *src);

int shotEncoder3D(int idev, Dataset3D *data_in, Dataset3D *data_ou, int codelength, int codestep, int seed);

void correlateWavelets3D(Dataset3D *data);

void writeSeisData2File3D(Dataset3D *dataset, int jshot, int iproc, char *fileNameSeismos, char *fileNameSources);

int dumpShotCoordinates3D(MSH *msh, Shot3D *shot);

int checkShotBoundaries3D(MSH *msh, Shot3D *shot);

int gridShot3D(MSH *msh, Shot3D *shot);

void computeDataBoundingBox(Dataset3D *dataset, float p1[], float p2[], float p3[], float p4[], \
                            float skirtLeft, float skirtRight, float skirtUp, float skirtDown, \
                            int applyCoordsShift);

void findDataGeometry_getShotLineDir(Dataset3D *dataset, float Vec_shotLine[]);
void findDataGeometry_getRecvLineDir(Dataset3D *dataset, float Vec_recvLine[]);
void changeStationHorizontalCoordinates(Dataset3D *dataset, double P0[], double angle);

void changeStationHorizontalCoordinates2NewFrame(Dataset3D *dataset, double p1[], double p2[], double p3[], double p4[], double Pref_x, double Pref_y);

void resampleSourceAndRec(Shot3D *shot, float dt, int nt);

void applyTopMute_3D(Shot3D *shot);
void attenuateDA_3D(Shot3D *shot);
void modelDA_3D(Shot3D *shot, float *seisDA, int ghostSign);
int  removeDA_3D(Shot3D *shot, float *seisDA);

void removeNoisyTraces(Shot3D *shot, float t0, float t1);

#endif