char execPathName[FILENAME_MAXLENGTH];

FILE *fLOG;

void initExecPathName(char *inputExecPathName) {
    strncpy(execPathName, inputExecPathName, FILENAME_MAXLENGTH);
    printf("\n\n\n  Printing Exec Path Name: \n");
    printf("    execPathName:     %s  \n", execPathName);
    printf("\n\n\n");
}

// Save seismogram in file
void output1d(char *name, float *array, int n1, float d1, float o1) {
    FILE *fhdr, *fbin;
    initFiles1d(name, &fhdr, &fbin, n1, d1, o1);
    size_t ns = n1;
    writeSamples(array, fbin, ns);
    fclose(fhdr);
    fclose(fbin);
}
void output2d(char *name, float *array, int n1, float d1, float o1, int n2, float d2, float o2) {
    FILE *fhdr, *fbin;
    initFiles2d(name, &fhdr, &fbin, n1, d1, o1, n2, d2, o2);
    size_t ns = n1*n2;
    writeSamples(array, fbin, ns);
    fclose(fhdr);
    fclose(fbin);
}
void output3d(char *name, float *array, int n1, float d1, float o1, int n2, float d2, float o2, int n3, float d3, float o3) {
    FILE *fhdr, *fbin;
    initFiles3d(name, &fhdr, &fbin, n1, d1, o1, n2, d2, o2, n3, d3, o3);
    size_t ns = n1*n2*n3;
    writeSamples(array, fbin, ns);
    fclose(fhdr);
    fclose(fbin);
}
void outputSmart1d(char *name, float *array, long long int n1, float d1, float o1)
{
    FILE *fhdr, *fbin;
    initFilesSmart1d(name, &fhdr, &fbin, n1, d1, o1);
    size_t ns = n1;
    writeSamples(array, fbin, ns);
    fclose(fhdr);
    fclose(fbin);
}
void outputSmart2d(char *name, float *array, long long int n1, float d1, float o1, long long int n2, float d2, float o2) {
    FILE *fhdr, *fbin;
    initFilesSmart2d(name, &fhdr, &fbin, n1, d1, o1, n2, d2, o2);
    size_t ns = n1*n2;
    writeSamples(array, fbin, ns);
    fclose(fhdr);
    fclose(fbin);
}
void outputSmart3d(char *name, float *array, long long int n1, float d1, float o1, long long int n2, float d2, float o2, long long int n3, float d3, float o3) {
    FILE *fhdr, *fbin;
    initFilesSmart3d(name, &fhdr, &fbin, n1, d1, o1, n2, d2, o2, n3, d3, o3);
    size_t ns = n1*n2*n3;
    writeSamples(array, fbin, ns);
    fclose(fhdr);
    fclose(fbin);
}
void outputSmart4d(char *name, float *array, long long int n1, float d1, float o1, long long int n2, float d2, float o2, \
                   long long int n3, float d3, float o3, long long int n4, float d4, float o4) {
    FILE *fhdr, *fbin;
    initFilesSmart4d(name, &fhdr, &fbin, n1, d1, o1, n2, d2, o2, n3, d3, o3, n4, d4, o4);
    size_t ns = n1*n2*n3*n4;
    writeSamples(array, fbin, ns);
    fclose(fhdr);
    fclose(fbin);
}
void outputSmart4d_J(char *name, float *array, \
                     long long int n1, float d1, float o1, long long int j1, \
                     long long int n2, float d2, float o2, long long int j2, \
                     long long int n3, float d3, float o3, long long int j3, \
                     long long int n4, float d4, float o4, long long int j4) {
    FILE *fhdr, *fbin;
    initFilesSmart4d_J(name, &fhdr, &fbin, n1, d1, o1, j1, n2, d2, o2, j2, n3, d3, o3, j3, n4, d4, o4, j4);

    int n1_ = n1/j1;
    int n2_ = n2/j2;
    int n3_ = n3/j3;
    int n4_ = n4/j4;
    float *tmp = CPU_zaloc1F(n1_*n2_*n3_*n4_);

    int i1, i2, i3, i4;
    for(i4=0;i4<n4_;i4++)
      for(i3=0;i3<n3_;i3++)
        for(i2=0;i2<n2_;i2++)
          for(i1=0;i1<n1_;i1++)
          {
              int i_ =  i1     + i2*n1_     + i3*n2_*n1_    + i4*n3_*n2_*n1_;
              int i  = (i1*j1) + (i2*j2)*n1 + (i3*j3)*n2*n1 + (i4*j4)*n3*n2*n1;
              if( i1*j1>n1  ||  i2*j2>n2  ||  i3*j3>n3  ||  i4*j4>n4)    continue;
              tmp[i_] = array[i];
          }

    size_t ns = n1_ * n2_ * n3_ * n4_;
    writeSamples(tmp, fbin, ns);
    fclose(fhdr);
    fclose(fbin);
    free(tmp);
}
void outputSmart5d(char *name, float *array, long long int n1, float d1, float o1, long long int n2, float d2, float o2, \
                   long long int n3, float d3, float o3, long long int n4, float d4, float o4, long long int n5, float d5, float o5) {
    FILE *fhdr, *fbin;
    initFilesSmart5d(name, &fhdr, &fbin, n1, d1, o1, n2, d2, o2, n3, d3, o3, n4, d4, o4, n5, d5, o5);
    size_t ns = n1*n2*n3*n4*n5;
    writeSamples(array, fbin, ns);
    fclose(fhdr);
    fclose(fbin);
}

void initFiles1d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1) {
    int   ndim = 1;
    int   nsamp[1];
    nsamp[0] = n1;
    float dintv[1];
    dintv[0] = d1;
    float  orig[1];
    orig[0] = o1;
    createFilesNew(name, fbin, fhdr, ndim, nsamp, dintv, orig);
return;
}
void initFiles2d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2) {
    int   ndim = 2;
    int   nsamp[2];
    nsamp[0] = n1;    nsamp[1] = n2;
    float dintv[2];
    dintv[0] = d1;    dintv[1] = d2;
    float  orig[2];
    orig[0] = o1;     orig[1] = o2;
    createFilesNew(name, fbin, fhdr, ndim, nsamp, dintv, orig);
return;
}
void   initFiles3d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2, int n3, float d3, float o3) {
    int   ndim = 3;
    int   nsamp[3];
    nsamp[0] = n1;    nsamp[1] = n2;    nsamp[2] = n3;
    float dintv[3];
    dintv[0] = d1;    dintv[1] = d2;    dintv[2] = d3;
    float orig[3];
    orig[0] = o1;     orig[1] = o2;     orig[2] = o3;
    createFilesNew(name, fbin, fhdr, ndim, nsamp, dintv, orig);
return;
}

void   initFiles4d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2, int n3, float d3, float o3, \
                                                         int n4, float d4, float o4) {
    int   ndim = 4;
    int   nsamp[4];
    nsamp[0] = n1;    nsamp[1] = n2;    nsamp[2] = n3;    nsamp[3] = n4;
    float dintv[4];
    dintv[0] = d1;    dintv[1] = d2;    dintv[2] = d3;    dintv[3] = d4;
    float orig[4];
    orig[0] = o1;     orig[1] = o2;     orig[2] = o3;     orig[3] = o4;
    createFilesNew(name, fbin, fhdr, ndim, nsamp, dintv, orig);
return;
}
void   initFiles5d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2, int n3, float d3, float o3, \
                                                         int n4, float d4, float o4, int n5, float d5, float o5) {
    int   ndim = 5;
    int   nsamp[5];
    nsamp[0] = n1;    nsamp[1] = n2;    nsamp[2] = n3;    nsamp[3] = n4;    nsamp[4] = n5;
    float dintv[5];
    dintv[0] = d1;    dintv[1] = d2;    dintv[2] = d3;    dintv[3] = d4;    dintv[4] = d5;
    float orig[5];
    orig[0] = o1;     orig[1] = o2;     orig[2] = o3;     orig[3] = o4;     orig[4] = o5;
    createFilesNew(name, fbin, fhdr, ndim, nsamp, dintv, orig);
return;
}
void initFilesSmart1d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1) {
    int   ndim = 1;
    int   nsamp[1];
    nsamp[0] = n1;
    float dintv[1];
    dintv[0] = d1;
    float  orig[1];
    orig[0] = o1;
    createFilesSmart(name, fbin, fhdr, ndim, nsamp, dintv, orig);
return;
}
void initFilesSmart2d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2) {
    int   ndim = 2;
    int   nsamp[2];
    nsamp[0] = n1;    nsamp[1] = n2;
    float dintv[2];
    dintv[0] = d1;    dintv[1] = d2;
    float  orig[2];
    orig[0] = o1;     orig[1] = o2;
    createFilesSmart(name, fbin, fhdr, ndim, nsamp, dintv, orig);
return;
}
void   initFilesSmart3d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2, int n3, float d3, float o3) {
    int   ndim = 3;
    int   nsamp[3];
    nsamp[0] = n1;    nsamp[1] = n2;    nsamp[2] = n3;
    float dintv[3];
    dintv[0] = d1;    dintv[1] = d2;    dintv[2] = d3;
    float orig[3];
    orig[0] = o1;     orig[1] = o2;     orig[2] = o3;
    createFilesSmart(name, fbin, fhdr, ndim, nsamp, dintv, orig);
return;
}

void   initFilesSmart4d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2, int n3, float d3, float o3, \
                                                         int n4, float d4, float o4) {
    int   ndim = 4;
    int   nsamp[4];
    nsamp[0] = n1;    nsamp[1] = n2;    nsamp[2] = n3;    nsamp[3] = n4;
    float dintv[4];
    dintv[0] = d1;    dintv[1] = d2;    dintv[2] = d3;    dintv[3] = d4;
    float orig[4];
    orig[0] = o1;     orig[1] = o2;     orig[2] = o3;     orig[3] = o4;
    createFilesSmart(name, fbin, fhdr, ndim, nsamp, dintv, orig);
return;
}
void   initFilesSmart4d_J(char *name, FILE **fhdr, FILE **fbin, \
                          int n1, float d1, float o1, int j1, \
                          int n2, float d2, float o2, int j2, \
                          int n3, float d3, float o3, int j3, \
                          int n4, float d4, float o4, int j4) {
    int   ndim = 4;
    int   nsamp[4];
    nsamp[0] = n1/j1;    nsamp[1] = n2/j2;    nsamp[2] = n3/j3;    nsamp[3] = n4/j3;
    float dintv[4];
    dintv[0] = d1*j1;    dintv[1] = d2*j2;    dintv[2] = d3*j3;    dintv[3] = d4*j4;
    float orig[4];
    orig[0] = o1;     orig[1] = o2;     orig[2] = o3;     orig[3] = o4;
    createFilesSmart(name, fbin, fhdr, ndim, nsamp, dintv, orig);
return;
}
void   initFilesSmart5d(char *name, FILE **fhdr, FILE **fbin, int n1, float d1, float o1, int n2, float d2, float o2, int n3, float d3, float o3, \
                                                         int n4, float d4, float o4, int n5, float d5, float o5) {
    int   ndim = 5;
    int   nsamp[5];
    nsamp[0] = n1;    nsamp[1] = n2;    nsamp[2] = n3;    nsamp[3] = n4;    nsamp[4] = n5;
    float dintv[5];
    dintv[0] = d1;    dintv[1] = d2;    dintv[2] = d3;    dintv[3] = d4;    dintv[4] = d5;
    float orig[5];
    orig[0] = o1;     orig[1] = o2;     orig[2] = o3;     orig[3] = o4;     orig[4] = o5;
    createFilesSmart(name, fbin, fhdr, ndim, nsamp, dintv, orig);
return;
}

void outputModel2d(char *name, float *inp, long long int n1, float d1, float o1, long long int n2, float d2, float o2, int n1b, int n2b) {
    int n1_ = n1-2*n1b;
    int n2_ = n2-2*n2b;
    float o1_ = o1 + n1b*d1;
    float o2_ = o2 + n2b*d2;
    float *out = CPU_zaloc1F(n1_*n2_);
    int i1, i2;
    for(i2=0; i2<n2_; i2++)
    {
        int shiftInp = (i2+n2b)*n1 + n1b;
        int shiftOut = i2*n1_;
        memcpy(out+shiftOut, inp+shiftInp, n1_*sizeof(float));
    }
    FILE *fhdr, *fbin;
    initFilesSmart2d(name, &fhdr, &fbin, n1_, d1, o1_, n2_, d2, o2_);
    size_t ns = n1_*n2_;
    writeSamples(out, fbin, ns);
    fclose(fhdr);
    fclose(fbin);
}


void createFilesNew(char *fileName, FILE **fBin, FILE **fHdr, int ndim, int *nsamp, float *dintv, float *orig) {
    int err;
    char tmpHdr[FILETYPE_LENGTH];
    char tmpBin[FILETYPE_LENGTH];
    char fileNameBin[FILENAME_MAXLENGTH];
    char fileNameHdr[FILENAME_MAXLENGTH];
    //strncpy(fileNameHdr, fileName, FILENAME_MAXLENGTH);
    //strncpy(fileNameBin, fileName, FILENAME_MAXLENGTH);
    strncpy(tmpHdr, ".hdr", FILETYPE_LENGTH);
    strncpy(tmpBin, ".bin", FILETYPE_LENGTH);
    //strncat( fileNameHdr, tmpHdr, FILETYPE_LENGTH);
    //strncat( fileNameBin, tmpBin, FILETYPE_LENGTH);
    //*
    err = fileStrCat(fileNameHdr, fileName, tmpHdr);
    if(err!=0)    printf("\n Error: Could not concatenate fileNameHdr for file: %s\n", fileName);
    err = fileStrCat(fileNameBin, fileName, tmpBin);
    if(err!=0)    printf("\n Error: Could not concatenate fileNameBin for file: %s\n", fileName);
    //*/

    *fHdr = openFile(&err, fileNameHdr, "w+");
    *fBin = openFile(&err, fileNameBin, "w+");

    fillHeaderFileNew(*fHdr, fileNameBin, ndim, nsamp, dintv, orig);

return;
}

void createFilesSmart(char *fileName, FILE **fBin, FILE **fHdr, int ndim, int *nsamp, float *dintv, float *orig)
{    
    int err;
    char tmpHdr[FILETYPE_LENGTH];
    char tmpBin[FILETYPE_LENGTH];
    char fileNameBinNoDir[FILENAME_MAXLENGTH];
    char fileNameBin[FILENAME_MAXLENGTH];
    char fileNameHdr[FILENAME_MAXLENGTH];
    char fileNameTmp[FILENAME_MAXLENGTH];
    char fileNameNoDirTmp[FILENAME_MAXLENGTH];
    strncpy(fileNameHdr, execPathName, FILENAME_MAXLENGTH);
    strncpy(fileNameBin, execPathName, FILENAME_MAXLENGTH);
    err = fileStrCat(fileNameTmp, execPathName, fileName);
    if(err!=0)    printf("\n Error: Could not concatenate fileNameTmp for file: %s\n", fileName);
    err = fileStrCat(fileNameNoDirTmp, "./", fileName);
    if(err!=0)    printf("\n Error: Could not concatenate fileNameTmp for file: %s\n", fileName);

    strncpy(tmpHdr, ".hdr", FILETYPE_LENGTH);
    strncpy(tmpBin, ".bin", FILETYPE_LENGTH);
    
    err = fileStrCat(fileNameHdr,      fileNameTmp, tmpHdr);
    if(err!=0)    printf("\n Error: Could not concatenate fileNameHdr for file: %s\n", fileName);
    err = fileStrCat(fileNameBin,      fileNameTmp, tmpBin);
    if(err!=0)    printf("\n Error: Could not concatenate fileNameBin for file: %s\n", fileName);
    err = fileStrCat(fileNameBinNoDir, fileNameNoDirTmp,    tmpBin);
    if(err!=0)    printf("\n Error: Could not concatenate fileNameBinNoDir for file: %s\n", fileNameNoDirTmp);

    *fHdr = openFile(&err, fileNameHdr, "w+");
    *fBin = openFile(&err, fileNameBin, "w+");

    // printf("\n\n\n  Printing smart paths: \n");
    // printf("    fileName:         %s  \n", fileName);
    // printf("    fileNameBinNoDir: %s  \n", fileNameBinNoDir);
    // printf("    execPathName:     %s  \n", execPathName);
    // printf("    fileNameHdr:      %s  \n", fileNameHdr);
    // printf("    fileNameBin:      %s  \n", fileNameBin);
    // printf("\n\n\n");

    fillHeaderFileNew(*fHdr, fileNameBinNoDir, ndim, nsamp, dintv, orig);

return;
}

void createFiles(char *fileName, FILE **fBin, FILE **fHdr, int ndim, int *nsamp, float *dintv) {
    int err;
    char tmpHdr[FILETYPE_LENGTH];
    char tmpBin[FILETYPE_LENGTH];
    char fileNameBin[FILENAME_MAXLENGTH];
    char fileNameHdr[FILENAME_MAXLENGTH];
    //strncpy(fileNameHdr, fileName, FILENAME_MAXLENGTH);
    //strncpy(fileNameBin, fileName, FILENAME_MAXLENGTH);
    strncpy(tmpHdr, ".hdr", FILETYPE_LENGTH);
    strncpy(tmpBin, ".bin", FILETYPE_LENGTH);
    //strncat( fileNameHdr, tmpHdr, FILETYPE_LENGTH);
    //strncat( fileNameBin, tmpBin, FILETYPE_LENGTH);
    //*
    err = fileStrCat(fileNameHdr, fileName, tmpHdr);
    if(err!=0)    printf("\n Error: Could not concatenate fileNameHdr for file: %s\n", fileName);
    err = fileStrCat(fileNameBin, fileName, tmpBin);
    if(err!=0)    printf("\n Error: Could not concatenate fileNameBin for file: %s\n", fileName);
    //*/

    *fHdr = openFile(&err, fileNameHdr, "w");
    *fBin = openFile(&err, fileNameBin, "w");

    fillHeaderFile(*fHdr, fileNameBin, ndim, nsamp, dintv);

return;
}

void createFileLog()
{
    openFileLog(&fLOG);
}

void openFileLog(FILE **fp)
{    
    int err;
    char tmpLog[FILETYPE_LENGTH];
    char fileNameLog[FILENAME_MAXLENGTH];
    char fileNameTmp[FILENAME_MAXLENGTH];
    strncpy(fileNameLog, execPathName, FILENAME_MAXLENGTH);

    err = fileStrCat(fileNameTmp, execPathName, "LOG");
    if(err!=0)    printf("\n Error: Could not concatenate fileNameTmp for LOG file\n");

    strncpy(tmpLog, ".txt", FILETYPE_LENGTH);
    
    err = fileStrCat(fileNameLog, fileNameTmp, tmpLog);
    if(err!=0)    printf("\n Error: Could not concatenate fileNameLog for LOG file\n");
    

    *fp = openFile(&err, fileNameLog, "w+");

    printf("\n\n\n  Printing log paths: \n");
    printf("    execPathName:     %s  \n", execPathName);
    printf("    fileNameLog:      %s  \n", fileNameLog);
    printf("    fileNameTmp:      %s  \n", fileNameTmp);
    printf("\n\n\n");

return;
}

void log_LineSearch(float alpha0, float alpha1, float alpha2, float objFunc0,  float objFunc1,  float objFunc2, float alphaMin, int iter)
{
    fprintf(fLOG, "\n\n Result of line search in iteration %d:", iter);
    fprintf(fLOG, "  \n    alpha0   = %f    alpha1   = %f    alpha2   = %f",   alpha0,   alpha1,   alpha2);
    fprintf(fLOG, "  \n    objFunc0 = %g    objFunc1 = %g    objFunc2 = %g", objFunc0, objFunc1, objFunc2);
    fprintf(fLOG, "  \n    objFunc0 = %g    objFunc1 = %g    objFunc2 = %g", 1.0, objFunc1/objFunc0, objFunc2/objFunc0);    printf("\n\n");
    fprintf(fLOG, "\n  >>> Optimum step-length: alphaMin = %f <<< \n\n", alphaMin);
}

void log_InitIter(int iter, time_t time_now, time_t time_begin)
{
    time_t time_elapsed = time_now-time_begin;
    fprintf(fLOG, "\n\n   ********************************* \n");
    fprintf(fLOG,     "   ***** Starting iteration %d ***** ", iter);
    fprintf(fLOG,   "\n   ********************************* \n\n");
    fprintf(fLOG, "\n Total elapsed time so far: %d seconds. \n\n", ((int) time_elapsed));
}
void log_endIter(int iter, time_t time_now, time_t time_begin, time_t time_begin_iter)
{
    time_t time_elapsed           = time_now-time_begin;
    time_t time_elapsed_iteration = time_now-time_begin_iter;
    fprintf(fLOG, "\n Total elapsed time so far: %d seconds. \n", ((int) time_elapsed));
    fprintf(fLOG, "\n Total elapsed time in this iteration: %d seconds. \n", ((int) time_elapsed_iteration));
    fprintf(fLOG," \n\n   =) =) =)    End of iteration %d    =) =) =)   \n\n", iter);
}

FILE* openFile(int *err, char *fileName, char *mode) {
    *err = 0;
    FILE *fp = NULL;
    fp = fopen(fileName, mode);
    if(fp==NULL) {
        printf("\n Error: Could not open file %s in mode %s \n ", fileName, mode);
        *err = -1;
    }
    //else {
        //printf("\n Log: File %s was successfully opened in mode %s\n", fileName, mode);
    //}
return fp;
}

void fillHeaderFile(FILE *fHdr, char *fileNameBin, int ndim, int *nsamp, float *dintv) {
    int idim;
    for(idim=0; idim<ndim; idim++)
        fprintf(fHdr, " n%d=%d \n d%d=%f \n o%d=0 \n x%d=0 \n", idim+1, nsamp[idim], idim+1, dintv[idim], idim+1, idim+1);
    fprintf(fHdr, " esize=%d \n in=\"%s\" \n data_format=\"native_float\" \n", ((int) sizeof(float)), fileNameBin);
}
void fillHeaderFileNew(FILE *fHdr, char *fileNameBin, int ndim, int *nsamp, float *dintv, float *orig) {
    int idim;
    for(idim=0; idim<ndim; idim++)
        fprintf(fHdr, " n%d=%d \n d%d=%f \n o%d=%f \n x%d=0 \n", idim+1, nsamp[idim], idim+1, dintv[idim], idim+1, orig[idim], idim+1);
    fprintf(fHdr, " esize=%d \n in=\"%s\" \n data_format=\"native_float\" \n", ((int) sizeof(float)), fileNameBin);
}
void writeSamples(const float *array, FILE *fBin, size_t ns) {
    fwrite(array, sizeof(float), ns, fBin);
return;
}
int fileStrCat(char *sout, const char *s1, const char *s2) {
    
    int  err=0;
    char stmp[FILENAME_MAXLENGTH];

    if ( (strlen(s1)+1) > sizeof(stmp) ) {
        printf("\n Error: Could not copy strings s1:%s to stmp:%s. Size of stmp does not fit size of s1.\n", s1, stmp);
        err -= 1;
    }
    else {
        strncpy(stmp, s1, sizeof(stmp));
    }
    if ( (strlen(s2)+1) > (sizeof(stmp)-strlen(stmp)) ) {
        printf("\n Error: Could not concatenate strings s1:%s and s2:%s. Size of s1 does not fit size of s2.\n", s1, s2);
        err -= 1;
    }    
    else {
        strncat(stmp, s2, strlen(s2));
    }
    strncpy(sout, stmp, sizeof(stmp));
return err;
}

void fftw_freeMem(fftw_complex *array, long long int n1) {
    fftw_free(array);
    AllocatedMem_CPU -= 2*n1*sizeof(float);
}
void freeMem(void *array, size_t size) {
    free(array);
    AllocatedMem_CPU -= size;
}
fftw_complex* CPU_zalocFFTW1(long long int n1)
{
    fftw_complex *array = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n1);
    AllocatedMem_CPU += 2*n1*sizeof(float);
    long long int i;
    for(i=0; i<n1; i++) {
        array[i][0] = 0.0f;
        array[i][1] = 0.0f;
    }
return array;
}
double* CPU_zaloc1D(long long int n1) {
    double *m;
    m = (double *) calloc(n1, sizeof(double));
return m;
}
float* CPU_zaloc1F(long long int n1) {
    float *m;
    m = (float *) calloc(n1, sizeof(float));
    AllocatedMem_CPU += n1*sizeof(float);
return m;
}
float** CPU_zaloc2F(long long int n2, long long int n1) {
    float **m;
    long long int i2;
    m = (float **) malloc(n2*sizeof(float*));
    for(i2=0; i2<n2; i2++) {
        m[i2] = (float *) calloc(n1, sizeof(float));
    }
return m;
}
int* CPU_zaloc1I(long long int n1) {
    int *m;
    m = (int *) calloc(n1, sizeof(int));
return m;
}
int** CPU_zaloc2I(long long int n2, long long int n1) {
    int **m;
    long long int i2;
    m = (int **) malloc(n2*sizeof(int*));
    for(i2=0; i2<n2; i2++) {
        m[i2] = (int *) calloc(n1, sizeof(int));
    }
return m;
}
Shot* CPU_zaloc1Shot(int n1) {
    Shot *m;
    m = (Shot *) calloc(n1, sizeof(Shot));
return m;
}
Shot3D* CPU_zaloc1Shot3D(int n1) {
    Shot3D *m;
    m = (Shot3D *) calloc(n1, sizeof(Shot3D));
return m;
}


void zeroArray(float *array, int n) {
    int i;
    for(i=0; i<n; i++)
        array[i] = 0.0f;
}


void readImageFile(char *filename, float *image, int n) {
    FILE *fp = NULL;
    fp = fopen(filename, "r");
    if(fp==NULL)    printf("Could not open image file");
    int nread = fread(image, sizeof(float), n, fp);
    printf("\n\n  Read %d samples from image file. n=%d \n\n", nread, n);
    fclose(fp);
}


void readShot(char *seismogramFileName, char *sourceFileName, Shot *shot) {

    int ishot = shot->ishot;
    int nrec  = shot->nrec;
    int nt    = shot->nt;

    FILE *fp1 = NULL;
    FILE *fp2 = NULL;
    fp1 = fopen(seismogramFileName, "r");  if(fp1==NULL)  printf("\n\n Cannot open seismogram file\n\n");
    fp2 = fopen(    sourceFileName, "r");  if(fp2==NULL)  printf("\n\n Cannot open source file\n\n");

    fseek(fp1, ishot*nrec*nt*sizeof(float), 0);

    fread((shot->seismogram), sizeof(float), nrec*nt, fp1);
    fread((shot->source)    , sizeof(float),      nt, fp2);

    fclose(fp1);
    fclose(fp2);

return;    
}

float* readModel(char *arqName, int nx, int nz) {
    
    if(access(arqName, F_OK) != 0)     return NULL;

    float *model = CPU_zaloc1F(nx*nz);
    FILE *fp = fopen(arqName,"r");
    fread(model, sizeof(float), nx*nz, fp);

return model;
}


/*
float* readModelInterp(char *arqName, int nx, int nz, float dx, float dz, int nx_, int nz_, float dx_, float dz_) {
    float *model_input  = CPU_zaloc1F(nx*nz);
    float *model_interp = CPU_zaloc1F(nx_*nz_);
    FILE *fp = fopen(arqName,"r");
    fread(model_input, sizeof(float), nx*nz, fp);

    int ix, iz;
    int ix_, iz_;
    int i_, i_x0z0, i_x0z1, i_x1z0, ix1z1; 
    float x, z, x0, z0, x1, z1;
    float px0, px1, pz0, pz1, ptot;

    for(ix_=0; ix_<nx_; ix_++)
    {
        for(iz_=0; iz_<nz_; iz_++)
        {
            x = ix_*dx_;
            z = iz_*dz_;
            ix = floorf(x/dx);
            iz = floorf(z/dz);

            x0 = ix*dx;
            z0 = iz*dz;
            x1 = (ix+1)*dx;
            z1 = (iz+1)*dz;

            i_ = iz_ + ix_*nz_;

            i_x0z0  = iz   +  ix    *nz;
            i_x0z1  = iz+1 +  ix    *nz;
            i_x1z0  = iz   + (ix+1) *nz;
            i_x1z1  = iz+1 + (ix+1) *nz;

            px0 = (x1-x)/dx;
            px1 = (x-x0)/dx;
            pz0 = (z1-z)/dz;
            pz1 = (z-z0)/dz;
            ptot = px0 + px1 + pz0 + pz1;
            px0 /= ptot;
            px1 /= ptot;
            pz0 /= ptot;
            pz1 /= ptot;

            model_interp[i_] = px0 * (pz0*model[i_x0z0] + pz1*model[i_x0z1]) + px1 * (pz0*model[i_x1z0] + pz1*model[i_x1z1]);

        }
    }

    free(model_input);

return model_interp;
}
*/

//*
int indexShotGathersFromSEGY(FILE *fp, char *filename, int **shotIndexes, int **tracesPerShot, int *nt, float *dt, int *ntraces)
{
    // Opening seismic data file
    fp = fopen(filename, "r");
    if(fp==NULL)
    { 
        printf("\n\n  ERROR: could not open seismic file <%s>. Exiting program.  \n\n", filename);
        exit(-1);
    }

    // reading number of samples and sample interval in data trace
    short int ns;
    short int ds;
    short int traceSorting;
    fseek(fp, 3217, SEEK_SET);
    fread(&ds, 2, 1, fp);
    fseek(fp, 3221, SEEK_SET);
    fread(&ns, 2, 1, fp);
    fseek(fp, 3229, SEEK_SET);
    fread(&traceSorting, 2, 1, fp);

    // finding the size of the entire segy file
    // fseek(fp, 0L, SEEK_END);
    fseek(fp, 0, SEEK_END);
    size_t fileSize = ftell(fp);
    rewind(fp);

    size_t bytesPerTrace  = 240 + ns*4;
    size_t numberOfTraces = (fileSize-3600) / bytesPerTrace;
    if( ((fileSize-3600)%bytesPerTrace) != 0 ) {
        printf(" \n\n ERROR: (fileSize-3600) is not multiple of bytesPerTrace. Exiting program. \n\n ");
        exit(-1);
    }

    int nshot_temp = 10000000;
    int *shotIndexes_temp   = CPU_zaloc1I(nshot_temp);
    int *tracesPerShot_temp = CPU_zaloc1I(nshot_temp);

    int itr;
    int shot_curr=-9999, shot_trc=-9999, shot_num=0;
    int tracesInShot = 0;
    for(itr=0; itr<numberOfTraces; itr++) 
    {
        tracesInShot++;
        size_t pos = 3600 + bytesPerTrace * itr + 9;
        fseek(fp, pos, SEEK_SET);
        fread(&shot_trc, 2, 1, fp);
        if( (shot_curr!=shot_trc)  &&  itr>0)
        {
            tracesPerShot_temp[shot_num] = tracesInShot;
            shotIndexes_temp[shot_num]   = itr;
            shot_num++;
            tracesInShot = 0;
        }
        shot_curr = shot_trc;
    }

    int nshot = shot_num + 1;
    *shotIndexes   = CPU_zaloc1I(nshot);
    *tracesPerShot = CPU_zaloc1I(nshot);

    int ishot;
    for(ishot=0; ishot<nshot; ishot++)
    {
        (*shotIndexes)[ishot]   = shotIndexes_temp[ishot];
        (*tracesPerShot)[ishot] = tracesPerShot_temp[ishot];
    }
    free(shotIndexes_temp);
    free(tracesPerShot_temp);

    *nt = ns;
    *dt = ds;
    *ntraces = numberOfTraces;

return nshot;
}

//*/

////////// Odisseia functions
void openOdisseiaFiles(FILE **fMet0, FILE **fMet1, FILE **fMet2, FILE **fData, char *fileNameMet0_inp, int OdisseiaFileType, char mode[]) {
    
    int err;
    char fileNameMet0[FILENAME_MAXLENGTH];
    char fileNameMet1[FILENAME_MAXLENGTH];
    char fileNameMet2[FILENAME_MAXLENGTH];
    char fileNameData[FILENAME_MAXLENGTH];

    // printf("\n  Opening Odisseia file <%s> in mode  %s  \n", fileNameMet0_inp, mode);

    fileStrCat(fileNameMet0, fileNameMet0_inp, ".odmet0");

    // output_permissions(fileNameMet0);

    // printf("\n  Opening Odisseia Met0 file <%s>. \n", fileNameMet0);

    // chmodFile("0777", fileNameMet0);
    // output_permissions(fileNameMet0);
    *fMet0 = openFile(&err, fileNameMet0, "r");
    // output_permissions(fileNameMet0);

    getFileNameOdisseia(fileNameMet1, *fMet0, 1); //backhere
    getFileNameOdisseia(fileNameMet2, *fMet0, 2);
    getFileNameOdisseia(fileNameData, *fMet0, 3);
    // output_permissions(fileNameMet1);
    // output_permissions(fileNameMet2);
    // output_permissions(fileNameData);
    *fMet1 = openFile(&err, fileNameMet1, mode);  
    int nshot_test;
    fread(&nshot_test,  sizeof(int), 1, *fMet1);
    rewind(*fMet1);
    *fMet2 = openFile(&err, fileNameMet2, mode);
    *fData = openFile(&err, fileNameData, mode);
    // *fMet1 = openFile(&err, fileNameMet1, "r");
    // *fMet2 = openFile(&err, fileNameMet2, "r");
    // *fData = openFile(&err, fileNameData, "r");

    // printf("\n\n ()()()()()()()()  Opening input files for the first time: ");
    // printf("\n fileNameMet0 = %s     fMet0 = %p", fileNameMet0, *fMet0);
    // printf("\n fileNameMet1 = %s     fMet1 = %p", fileNameMet1, *fMet1);
    // printf("\n fileNameMet2 = %s     fMet2 = %p", fileNameMet2, *fMet2);
    // printf("\n fileNameData = %s     fData = %p", fileNameData, *fData);
    // printf("\n ()()()()()()()()()()()()()()()()()()()()()()()()()()()()()() \n\n");

return;
}
void getFileNameOdisseia(char *fileName, FILE *fMet0, int which)
{
    char inputParm[1024];
    char parmValue[1024];
    
    char *parmName;
    if(which==1)    parmName = "MetaDataFile1";
    if(which==2)    parmName = "MetaDataFile2";
    if(which==3)    parmName = "DataFile";
    rewind(fMet0);

    while(!feof(fMet0))
    {
        fscanf(fMet0, "%s", inputParm);
        if( (strncmp(inputParm, parmName, strlen(parmName))) == 0) {
            fscanf(fMet0, "%s", parmValue);
        }
    }
    strncpy(fileName, parmValue, 1024);

return;    
}
void readOdisseiaMetaFile_get_nshot(FILE *fMet1, int *nshot)
{
    fread(nshot, sizeof(int), 1, fMet1);
return;
}
void readOdisseiaMetaFile_get_shotInfo(FILE *fMet1, int ishot, int *nsou, int *nrec, int *nt, float *dt, float *ot, float *maxFreq, long long int *nStationsPreceding)
{
    long long int position = (3*sizeof(int) + 3*sizeof(float) + 1*sizeof(long long int)) * ishot + 1*sizeof(int);
    fseek(fMet1, position, SEEK_SET);
    fread(nsou              ,   sizeof(int)          , 1, fMet1);
    fread(nrec              ,   sizeof(int)          , 1, fMet1);
    fread(nt                ,   sizeof(int)          , 1, fMet1);
    fread(dt                ,   sizeof(float)        , 1, fMet1);
    fread(ot                ,   sizeof(float)        , 1, fMet1);
    fread(maxFreq           ,   sizeof(float)        , 1, fMet1);
    fread(nStationsPreceding,   sizeof(long long int), 1, fMet1);
return;
}
void readOdisseiaMetaFile_get_shotCoords(FILE *fMet1, FILE *fMet2, int ishot, \
                                         int nsou, int nrec, long long int nStationsPreceding, \
                                         float *xs, float *ys, float *zs, \
                                         float *xr, float *yr, float *zr)
{
    long long int position = 3*sizeof(long long int) * nStationsPreceding;
    fseek(fMet2, position, SEEK_SET);

    int isou;
    for(isou=0; isou<nsou; isou++)
    {
        long long int xs_lli, ys_lli, zs_lli;
        double        xs_dbl, ys_dbl, zs_dbl;
        fread(&xs_lli, sizeof(long long int), 1, fMet2);
        fread(&ys_lli, sizeof(long long int), 1, fMet2);
        fread(&zs_lli, sizeof(long long int), 1, fMet2);

        xs_dbl = (double) xs_lli;
        ys_dbl = (double) ys_lli;
        zs_dbl = (double) zs_lli;

        xs[isou] = xs_dbl/100.0;
        ys[isou] = ys_dbl/100.0;
        zs[isou] = zs_dbl/100.0;
        // printf("\n\n\n  >>>>>>> <<<<<< ");
        // printf("\n nStationsPreceding=%lld     xs_lli=%lld  ys_lli=%lld  zs_lli=%lld \n", nStationsPreceding, xs_lli, ys_lli, zs_lli);
        // printf("\n  >>>>>>> Reading from file: ishot:%d      xs=%f  ys=%f  zs=%f <<<<<<<< \n", ishot, xs[isou], ys[isou], zs[isou]);
        // printf("\n fMet1=%p  fMet2=%p\n", fMet1, fMet2);
        // printf("  >>>>>>> <<<<<< \n\n\n");
    }
    int irec;
    for(irec=0; irec<nrec; irec++)
    {
        long long int xr_lli, yr_lli, zr_lli;
        double        xr_dbl, yr_dbl, zr_dbl;
        fread(&xr_lli, sizeof(long long int), 1, fMet2);
        fread(&yr_lli, sizeof(long long int), 1, fMet2);
        fread(&zr_lli, sizeof(long long int), 1, fMet2);

        xr_dbl = (double) xr_lli;
        yr_dbl = (double) yr_lli;
        zr_dbl = (double) zr_lli;

        xr[irec] = xr_dbl/100.0;
        yr[irec] = yr_dbl/100.0;
        zr[irec] = zr_dbl/100.0;
        // printf("\n  Reading from file:  xr[irec]=%f  yr[irec]=%f  zr[irec]=%f  \n", xr[irec], yr[irec], zr[irec]);
    }

return;
}
void fillOdisseiaMetaFile_Seismic_ShotInfo(FILE *fMet1, Dataset3D *dataset)
{
    int nshot = dataset->nshot;
    fwrite(&nshot, sizeof(int), 1, fMet1);

    long long int nStationsPreceding = 0;
    int ishot;

    Shot3D *shot;
    
    shot = &((dataset->shot)[0]);
    // printf("\n Output dimensions: nshot=%d  nsou=%d  nrec=%d  nt=%d  dt=%f  ot=%f", nshot, shot->nsou, shot->nrec, shot->nt, shot->dt, shot->ot);

    for(ishot=0; ishot<nshot; ishot++)
    {        
        shot = &((dataset->shot)[ishot]);

        int   nsou    = shot->nsou;
        int   nrec    = shot->nrec;
        int   nt      = shot->nt;
        float dt      = shot->dt;
        float ot      = shot->ot;
        float maxFreq = shot->maxFrequency;

        fwrite(&nsou              , sizeof(int)          , 1, fMet1);
        fwrite(&nrec              , sizeof(int)          , 1, fMet1);
        fwrite(&nt                , sizeof(int)          , 1, fMet1);
        fwrite(&dt                , sizeof(float)        , 1, fMet1);
        fwrite(&ot                , sizeof(float)        , 1, fMet1);
        fwrite(&maxFreq           , sizeof(float)        , 1, fMet1);
        fwrite(&nStationsPreceding, sizeof(long long int), 1, fMet1);

        nStationsPreceding += (nsou+nrec);
    }

}
void fillOdisseiaMetaFile_Seismic_ShotCoords(FILE *fMet2, Dataset3D *dataset)
{
    int nshot = dataset->nshot;

    int ishot;
    for(ishot=0; ishot<nshot; ishot++)
    {
        Shot3D *shot = &((dataset->shot)[ishot]);
        
        int   nsou    = shot->nsou;
        int   nrec    = shot->nrec;

        int isou;
        for(isou=0; isou<nsou; isou++)
        {
            long long int xs = 100 * ((long long int) shot->xs[isou]);
            long long int ys = 100 * ((long long int) shot->ys[isou]);
            long long int zs = 100 * ((long long int) shot->zs[isou]);
            fwrite(&xs, sizeof(long long int), 1, fMet2);
            fwrite(&ys, sizeof(long long int), 1, fMet2);
            fwrite(&zs, sizeof(long long int), 1, fMet2);
        }
        int irec;
        for(irec=0; irec<nrec; irec++)
        {
            long long int xr = 100 * ((long long int) shot->xr[irec]);
            long long int yr = 100 * ((long long int) shot->yr[irec]);
            long long int zr = 100 * ((long long int) shot->zr[irec]);
            fwrite(&xr, sizeof(long long int), 1, fMet2);
            fwrite(&yr, sizeof(long long int), 1, fMet2);
            fwrite(&zr, sizeof(long long int), 1, fMet2);
        }
    }

}

void writeSeismogram2DataFile(FILE *fData, Dataset3D *dataset, int ishot2write)
{

    // Check if ishot2write is withing allowed range
    int nshot = dataset->nshot;
    if(ishot2write<0 || ishot2write>=nshot)    return;

    // Compute the offset
    long long int offset = 0;
    int ishot;
    for(ishot=0; ishot<ishot2write; ishot++)
    {
        Shot3D *shot = &((dataset->shot)[ishot]);
        int nsou = shot->nsou;
        int nrec = shot->nrec;
        int nt   = shot->nt;
        offset += (nsou+nrec) * nt * sizeof(float);
    }

    // Position file pointer at the correct offset
    fseek(fData, offset, 0);

    // Read information about ishot2write
    Shot3D *shot = &((dataset->shot)[ishot2write]);        
    int nsou = shot->nsou;
    int nrec = shot->nrec;
    int nt   = shot->nt;

    // Write source and seismogram of the ishot2write
    fwrite((shot->source)    , sizeof(float), nsou*nt, fData);
    fwrite((shot->seismogram), sizeof(float), nrec*nt, fData);

return;
}



int appendOdisseiaMetaFile_Seismic_ShotInfo(FILE *fMet1, Dataset3D *dataset)
{
    int ishot;
    Shot3D *shot;    
    shot = &((dataset->shot)[0]);

    // Stores original value of nshot
    int nshot_orig;
    rewind(fMet1); //printf("\n  WILL APPEND:  fMet1=%p  \n", fMet1);
    fread(&nshot_orig,  sizeof(int), 1, fMet1);
    
    int nshot     = dataset->nshot;
    int nshot_new = nshot_orig + nshot;
    rewind(fMet1);
    fwrite(&nshot_new,  sizeof(int), 1, fMet1);
    
    dataset->nshot_aux = nshot_orig;


    // long long int singleShotBytes_Met1 = (3*sizeof(int) + 3*sizeof(float) + 1*sizeof(long long int));
    long long int position_lastShot    = MET1_BYTESPERSHOT * (nshot_orig-1) + 1*sizeof(int);
    long long int position_newShot     = MET1_BYTESPERSHOT *  nshot_orig    + 1*sizeof(int);

    // printf("\n\n  nshot_orig=%d   MET1_BYTESPERSHOT=%lld  position_lastShot=%lld  position_newShot=%lld \n\n", \
    nshot_orig, MET1_BYTESPERSHOT, position_lastShot, position_newShot);

    int   nsou, nrec, nt;
    float dt, ot, maxFreq;
    long long int nStationsPreceding;

    fseek(fMet1, position_lastShot, SEEK_SET);
    fread(&nsou              , sizeof(int)          , 1, fMet1);
    fread(&nrec              , sizeof(int)          , 1, fMet1);
    fread(&nt                , sizeof(int)          , 1, fMet1);
    fread(&dt                , sizeof(float)        , 1, fMet1);
    fread(&ot                , sizeof(float)        , 1, fMet1);
    fread(&maxFreq           , sizeof(float)        , 1, fMet1);
    fread(&nStationsPreceding, sizeof(long long int), 1, fMet1);

    nStationsPreceding += nsou + nrec;
    
    if(dt!=shot->dt  ||  nt!=shot->nt)
    {
        printf("\n\n\n  Problems to append new data to original data. \n  Number of samples or sampling ratio are not compatible:  \n");
        printf("  Original data      :  nt=%d  dt=%f  \n", nt, dt);
        printf("  Data to be appended:  nt=%d  dt=%f  \n\n", shot->nt, shot->dt);
        printf(" ATTENTION:  This dataset will not be appended.  \n\n\n");
        return -1;
    }
    
    fseek(fMet1, position_newShot, SEEK_SET);
    for(ishot=0; ishot<nshot; ishot++)
    {
        shot = &((dataset->shot)[ishot]);

        shot->nStationsPreceding = nStationsPreceding;

        nsou    = shot->nsou;
        nrec    = shot->nrec;
        nt      = shot->nt;
        dt      = shot->dt;
        ot      = shot->ot;
        maxFreq = shot->maxFrequency;

        fwrite(&nsou              , sizeof(int)          , 1, fMet1);
        fwrite(&nrec              , sizeof(int)          , 1, fMet1);
        fwrite(&nt                , sizeof(int)          , 1, fMet1);
        fwrite(&dt                , sizeof(float)        , 1, fMet1);
        fwrite(&ot                , sizeof(float)        , 1, fMet1);
        fwrite(&maxFreq           , sizeof(float)        , 1, fMet1);
        fwrite(&nStationsPreceding, sizeof(long long int), 1, fMet1);
        
        nStationsPreceding += (nsou+nrec);
    }

return 0;
}
void appendOdisseiaMetaFile_Seismic_ShotCoords(FILE *fMet1, FILE *fMet2, Dataset3D *dataset)
{    
    int nshot      = dataset->nshot;
    int nshot_orig = dataset->nshot_aux;

    // long long int nStationsPreceding;
    // long long int position = MET1_BYTESPERSHOT * nshot_orig - sizeof(long long int);
    // fseek(fMet1, position, SEEK_SET);
    // fread(&nStationsPreceding, sizeof(long long int), 1, fMet1);

    Shot3D *shot = &((dataset->shot)[0]);
    long long int position = MET2_BYTESPERSTATION * shot->nStationsPreceding;
    fseek(fMet2, position, SEEK_SET);

    int ishot;
    for(ishot=0; ishot<nshot; ishot++)
    {
        shot = &((dataset->shot)[ishot]);
        
        int   nsou    = shot->nsou;
        int   nrec    = shot->nrec;

        int isou;
        for(isou=0; isou<nsou; isou++)
        {
            long long int xs = 100 * ((long long int) shot->xs[isou]);    //printf("\n  ishot=%d  xs=%d  \n", ishot, xs);
            long long int ys = 100 * ((long long int) shot->ys[isou]);    //printf("\n  ishot=%d  ys=%d  \n", ishot, ys);
            long long int zs = 100 * ((long long int) shot->zs[isou]);
            fwrite(&xs, sizeof(long long int), 1, fMet2);
            fwrite(&ys, sizeof(long long int), 1, fMet2);
            fwrite(&zs, sizeof(long long int), 1, fMet2);
        }
        int irec;
        for(irec=0; irec<nrec; irec++)
        {
            long long int xr = 100 * ((long long int) shot->xr[irec]);  //if(irec==0)  printf("\n  ishot=%d  xr=%d  \n", ishot, xr);
            long long int yr = 100 * ((long long int) shot->yr[irec]);
            long long int zr = 100 * ((long long int) shot->zr[irec]);
            fwrite(&xr, sizeof(long long int), 1, fMet2);
            fwrite(&yr, sizeof(long long int), 1, fMet2);
            fwrite(&zr, sizeof(long long int), 1, fMet2);
        }
    }
}
void appendSeismogram2DataFile(FILE *fData, Dataset3D *dataset, int ishot2write)
{
    // Check if ishot2write is withing allowed range
    int ishot;
    int nshot      = dataset->nshot;
    int nshot_orig = dataset->nshot_aux;

    if(ishot2write<0 || ishot2write>=nshot)    return;

    // Finding nStationsPreceding
    // long long int nStationsPreceding;
    // long long int position = MET1_BYTESPERSHOT * (nshot_orig + ishot2write) - sizeof(long long int);
    // fseek(fMet1, position, SEEK_SET);
    // fread(&nStationsPreceding, sizeof(long long int), 1, fMet1);

    // Position file pointer at the correct byte
    Shot3D *shot = &((dataset->shot)[ishot2write]);
    int nsou = shot->nsou;
    int nrec = shot->nrec;
    int nt   = shot->nt;
    long long int position = shot->nStationsPreceding * nt * sizeof(float);
    fseek(fData, position, SEEK_SET);

    // Write source and seismogram of the ishot2write
    fwrite((shot->source)    , sizeof(float), nsou*nt, fData);
    fwrite((shot->seismogram), sizeof(float), nrec*nt, fData);

return;
}



void writeDataset2DataFile(FILE *fData, Dataset3D *dataset)
{
    // Check if ishot2write is withing allowed range
    int nshot = dataset->nshot;    
    int ishot;
    for(ishot=0; ishot<nshot; ishot++)
    {
        // Read information about ishot2write
        Shot3D *shot = &((dataset->shot)[ishot]);        
        int nsou = shot->nsou;
        int nrec = shot->nrec;
        int nt   = shot->nt;
        
        // Write source and seismogram of the ishot2write
        fwrite((shot->source)    , sizeof(float), nsou*nt, fData);
        fwrite((shot->seismogram), sizeof(float), nrec*nt, fData);
    }

return;
}
void writeSeismogram2File(FILE *fData, float *source, float *seismogram, int nsou, int nrec, int nt, long long int nStationsPreceding)
{
    // Find position from where to read
    long long int position = nt * sizeof(float) * nStationsPreceding;
    fseek(fData, position, 0);
        
    // Read source and seismogram of the ishot2write
    fwrite(source    , sizeof(float), nsou*nt, fData);
    fwrite(seismogram, sizeof(float), nrec*nt, fData);
return;
}
void readSeismogramFromFile(FILE *fData, float *source, float *seismogram, int ishot, int nsou, int nrec, int nt, long long int nStationsPreceding)
{
    // Find position from where to read
    long long int position = nt * sizeof(float) * nStationsPreceding;
    fseek(fData, position, 0);
        
    // Read source and seismogram of the ishot2write
    fread(source    , sizeof(float), nsou*nt, fData);
    fread(seismogram, sizeof(float), nrec*nt, fData);
return;
}


void readSmart(char name[], float *array, size_t nitems)
{
    char fileName[1024];
    sprintf( (char *)  fileName,  "%s%s.bin", execPathName, name);
    
    FILE *fp = fopen(fileName, "r");
    fread(array, sizeof(float), nitems, fp);
    fclose(fp);
}

void createOdisseiaFilesNew_ExistingBinary(char *fileName, char *fileNameBinary, FILE **fMet0, FILE **fMet1, FILE **fMet2, FILE **fData, int OdisseiaFileType) {
    int err;
    char tmpMet0[FILETYPE_LENGTH];
    char tmpMet1[FILETYPE_LENGTH];
    char tmpMet2[FILETYPE_LENGTH];
    char tmpData[FILETYPE_LENGTH];
    char fileNameMet0[FILENAME_MAXLENGTH];
    char fileNameMet1[FILENAME_MAXLENGTH];
    char fileNameMet2[FILENAME_MAXLENGTH];
    char fileNameData[FILENAME_MAXLENGTH];
    
    strncpy(tmpMet0, ".odmet0", FILETYPE_LENGTH);
    strncpy(tmpMet1, ".odmet1", FILETYPE_LENGTH);
    strncpy(tmpMet2, ".odmet2", FILETYPE_LENGTH);
    if(OdisseiaFileType == ODISSEIA_FILE_TYPE_SEISMIC)
        strncpy(tmpData, ".odseis", FILETYPE_LENGTH);
    if(OdisseiaFileType == ODISSEIA_FILE_TYPE_P_VELOCITY)
        strncpy(tmpData, ".odvelp", FILETYPE_LENGTH);
    if(OdisseiaFileType == ODISSEIA_FILE_TYPE_S_VELOCITY)
        strncpy(tmpData, ".odvels", FILETYPE_LENGTH);
    if(OdisseiaFileType == ODISSEIA_FILE_TYPE_DENSITY)
        strncpy(tmpData, ".oddens", FILETYPE_LENGTH);
    if(OdisseiaFileType == ODISSEIA_FILE_TYPE_SOURCE_SIGNATURE)
        strncpy(tmpData, ".odsign", FILETYPE_LENGTH);

    err = fileStrCat(fileNameMet0, fileName, tmpMet0);
    if(err!=0)    printf("\n Error: Could not concatenate fileNameMet0 for file: %s\n", fileName);
    err = fileStrCat(fileNameMet1, fileName, tmpMet1);
    if(err!=0)    printf("\n Error: Could not concatenate fileNameMet1 for file: %s\n", fileName);
    err = fileStrCat(fileNameMet2, fileName, tmpMet2);
    if(err!=0)    printf("\n Error: Could not concatenate fileNameMet2 for file: %s\n", fileName);
    err = fileStrCat(fileNameData, fileNameBinary, tmpData);
    if(err!=0)    printf("\n Error: Could not concatenate fileNameDataBinary for file: %s\n", fileName);

    *fMet0 = openFile(&err, fileNameMet0, "w+");
    *fMet1 = openFile(&err, fileNameMet1, "w+");
    *fMet2 = openFile(&err, fileNameMet2, "w+");
    *fData = openFile(&err, fileNameData, "r+");

    initOdisseiaHeaderFile(*fMet0, fileNameMet1, fileNameMet2, fileNameData, OdisseiaFileType);

return;
}
void createOdisseiaFilesNew(char *fileName, FILE **fMet0, FILE **fMet1, FILE **fMet2, FILE **fData, int OdisseiaFileType) {
    int err;
    char tmpMet0[FILETYPE_LENGTH];
    char tmpMet1[FILETYPE_LENGTH];
    char tmpMet2[FILETYPE_LENGTH];
    char tmpData[FILETYPE_LENGTH];
    char fileNameMet0[FILENAME_MAXLENGTH];
    char fileNameMet1[FILENAME_MAXLENGTH];
    char fileNameMet2[FILENAME_MAXLENGTH];
    char fileNameData[FILENAME_MAXLENGTH];
    
    strncpy(tmpMet0, ".odmet0", FILETYPE_LENGTH);
    strncpy(tmpMet1, ".odmet1", FILETYPE_LENGTH);
    strncpy(tmpMet2, ".odmet2", FILETYPE_LENGTH);
    if(OdisseiaFileType == ODISSEIA_FILE_TYPE_SEISMIC)
        strncpy(tmpData, ".odseis", FILETYPE_LENGTH);
    if(OdisseiaFileType == ODISSEIA_FILE_TYPE_P_VELOCITY)
        strncpy(tmpData, ".odvelp", FILETYPE_LENGTH);
    if(OdisseiaFileType == ODISSEIA_FILE_TYPE_S_VELOCITY)
        strncpy(tmpData, ".odvels", FILETYPE_LENGTH);
    if(OdisseiaFileType == ODISSEIA_FILE_TYPE_DENSITY)
        strncpy(tmpData, ".oddens", FILETYPE_LENGTH);
    if(OdisseiaFileType == ODISSEIA_FILE_TYPE_SOURCE_SIGNATURE)
        strncpy(tmpData, ".odsign", FILETYPE_LENGTH);

    err = fileStrCat(fileNameMet0, fileName, tmpMet0);
    if(err!=0)    printf("\n Error: Could not concatenate fileNameMet0 for file: %s\n", fileName);
    err = fileStrCat(fileNameMet1, fileName, tmpMet1);
    if(err!=0)    printf("\n Error: Could not concatenate fileNameMet1 for file: %s\n", fileName);
    err = fileStrCat(fileNameMet2, fileName, tmpMet2);
    if(err!=0)    printf("\n Error: Could not concatenate fileNameMet2 for file: %s\n", fileName);
    err = fileStrCat(fileNameData, fileName, tmpData);
    if(err!=0)    printf("\n Error: Could not concatenate fileNameData for file: %s\n", fileName);

    *fMet0 = openFile(&err, fileNameMet0, "w+");
    *fMet1 = openFile(&err, fileNameMet1, "w+");
    *fMet2 = openFile(&err, fileNameMet2, "w+");
    *fData = openFile(&err, fileNameData, "w+");

    initOdisseiaHeaderFile(*fMet0, fileNameMet1, fileNameMet2, fileNameData, OdisseiaFileType);

return;
}
void initOdisseiaHeaderFile(FILE *fMet0, char *fileNameMet1, char *fileNameMet2, char *fileNameData, int OdisseiaFileType) {
    // printf("\n  Writting this on fMet0: \n");
    // printf(" OdisseiaFileType %d \n MetaDataFile1 %s \n MetaDataFile2 %s \n DataFile %s \n data_format native_float \n", \
    //         OdisseiaFileType, fileNameMet1, fileNameMet2, fileNameData);
    fprintf(fMet0, " OdisseiaFileType %d \n MetaDataFile1 %s \n MetaDataFile2 %s \n DataFile %s \n data_format native_float \n", \
            OdisseiaFileType, fileNameMet1, fileNameMet2, fileNameData);
return;
}