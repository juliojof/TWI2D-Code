

#ifndef IO_H
#define IO_H

#include <fftw3.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>   // stat
#include <stdbool.h>    // bool type
// #include "shot.h"
// #include "shot_3D.h"
// #include "Math_Functions.h"

#define FILENAME_MAXLENGTH 1024
#define FILETYPE_LENGTH 16

#define MET1_BYTESPERSHOT    (  (long long int) ( 3*sizeof(int) + 3*sizeof(float) + 1*sizeof(long long int) )  )
#define MET2_BYTESPERSTATION (  (long long int) ( 3*sizeof(long long int) )  )

#define  ODISSEIA_FILE_TYPE_SEISMIC              1
#define  ODISSEIA_FILE_TYPE_P_VELOCITY           2
#define  ODISSEIA_FILE_TYPE_S_VELOCITY           3
#define  ODISSEIA_FILE_TYPE_DENSITY              4
#define  ODISSEIA_FILE_TYPE_SOURCE_SIGNATURE     5

#define  ODISSEIA_OPEN_MODE_CREATENEW            1
#define  ODISSEIA_OPEN_MODE_OPENEXISTING         2


void initExecPathName(char *inputExecPathName);

void readSmart(char name[], float *array, size_t nitems);

void output1d(char *name, float *array, int n1, float d1, float o1);
void output2d(char *name, float *array, int n1, float d1, float o1, int n2, float d2, float o2);
void output3d(char *name, float *array, int n1, float d1, float o1, int n2, float d2, float o2, int n3, float d3, float o3);

void outputSmart1d(char *name, float *array, long long int n1, float d1, float o1);
void outputSmart2d(char *name, float *array, long long int n1, float d1, float o1, long long int n2, float d2, float o2);
void outputSmart3d(char *name, float *array, long long int n1, float d1, float o1, long long int n2, float d2, float o2, long long int n3, float d3, float o3);
void outputSmart4d(char *name, float *array, long long int n1, float d1, float o1, long long int n2, float d2, float o2, \
                   long long int n3, float d3, float o3, long long int n4, float d4, float o4);


void outputSmart4d_J(char *name, float *array, \
                     long long int n1, float d1, float o1, long long int j1, \
                     long long int n2, float d2, float o2, long long int j2, \
                     long long int n3, float d3, float o3, long long int j3, \
                     long long int n4, float d4, float o4, long long int j4);

void outputSmart5d(char *name, float *array, long long int n1, float d1, float o1, long long int n2, float d2, float o2, \
                   long long int n3, float d3, float o3, long long int n4, float d4, float o4, long long int n5, float d5, float o5);

void initFiles1d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1);
void initFiles2d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2);
void initFiles3d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2, int n3, float d3, float o3);
void initFiles4d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2, int n3, float d3, float o3, int n4, float d4, float o4);
void initFiles5d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2, int n3, float d3, float o3, \
                 int n4, float d4, float o4, int n5, float d5, float o5);

void initFilesSmart1d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1);
void initFilesSmart2d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2);
void initFilesSmart3d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2, int n3, float d3, float o3);
void initFilesSmart4d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2, int n3, float d3, float o3, \
                                                            int n4, float d4, float o4);

void   initFilesSmart4d_J(char *name, FILE **fhdr, FILE **fbin, \
                          int n1, float d1, float o1, int j1, \
                          int n2, float d2, float o2, int j2, \
                          int n3, float d3, float o3, int j3, \
                          int n4, float d4, float o4, int j4);

void   initFilesSmart5d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2, int n3, float d3, float o3, \
                                                         int n4, float d4, float o4, int n5, float d5, float o5);


void outputModel2d(char *name, float *inp, long long int n1, float d1, float o1, long long int n2, float d2, float o2, int n1b, int n2b);
void outputModel3d(char *name, float *inp, long long int n1, float d1, float o1, long long int n2, float d2, float o2, long long int n3, float d3, float o3, long long int n1b, long long int n2b, long long int n3b);

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

// Initialize files in Odisseia format
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


void createFilesNew(char *fileName, FILE **fBin, FILE **fHdr, int ndim, int *nsamp, float *dintv, float *orig);

void createFilesSmart(char *fileName, FILE **fBin, FILE **fHdr, int ndim, int *nsamp, float *dintv, float *orig);

void createFiles(char *fileName, FILE **fBin, FILE **fHdr, int ndim, int *nsamp, float *dintv);

void createFileLog();

void openFileLog(FILE **fp);

void log_LineSearch(float alpha0, float alpha1, float alpha2, float objFunc0,  float objFunc1,  float objFunc2, float alphaMin, int iter);

void log_InitIter(int iter, time_t time_now, time_t time_begin);

void log_endIter(int iter, time_t time_now, time_t time_begin, time_t time_begin_iter);

FILE* openFile(int *err, char *fileName, char *mode);
bool file_exists (char *filename);
void chmodFile(char mode[], char fileName[]);
void output_permissions(char *fileName);

void fillHeaderFile(FILE *fHdr, char *fileNameBin, int ndim, int *nsamp, float *dintv);

void fillHeaderFileNew(FILE *fHdr, char *fileNameBin, int ndim, int *nsamp, float *dintv, float *orig);

void writeSamples(const float *array, FILE *fBin, size_t ns);

int fileStrCat(char *sout, const char *s1, const char *s2);

void fftw_freeMem(fftw_complex *array, long long int n1);

void freeMem(void *array, size_t size);

fftw_complex* CPU_zalocFFTW1(long long int n1);

double* CPU_zaloc1D(long long int n1);

float* CPU_zaloc1F(long long int n1);

float** CPU_zaloc2F(long long int n2, long long int n1);

int* CPU_zaloc1I(long long int n1);

int** CPU_zaloc2I(long long int n2, long long int n1);

Shot* CPU_zaloc1Shot(int n1);

Shot3D* CPU_zaloc1Shot3D(int n1);

void zeroArray(float *array, int n);

void readImageFile(char *filename, float *image, int n);

void readShot(char *seismogramFileName, char *sourceFileName, Shot *shot);
void readShot3D(char *seismogramFileName, char *sourceFileName, Shot3D *shot);

float* readModel(char *arqName, int nx, int nz);
float* readModel3D(char *arqName, int nx, int ny, int nz);

int indexShotGathersFromSEGY(FILE *fp, char *filename, int **shotIndexes, int **tracesPerShot, int *nt, float *dt, int *ntraces);

void backSpace(int ntimes);


#endif