#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
// #include "Odisseia_defines.h"
#include "shot_3D.h"
// #include "IO.h"
// #include "Math_Functions.h"
#include "DataTypes.h"
// #include "parmHandling.h"
// #include <unistd.h>
// #include "PropagationFunctions.h"

void saveDataset3D(char *fileName, Dataset3D *dataset)
{
    createOdisseiaFilesNew(fileName, &(dataset->fMet0), &(dataset->fMet1), &(dataset->fMet2), &(dataset->fData), ODISSEIA_FILE_TYPE_SEISMIC);
    fillOdisseiaMetaFile_Seismic_ShotInfo(dataset->fMet1, dataset);
    fillOdisseiaMetaFile_Seismic_ShotCoords(dataset->fMet2, dataset);
    writeDataset2DataFile(dataset->fData, dataset);
    fclose(dataset->fMet0);
    fclose(dataset->fMet1);
    fclose(dataset->fMet2);
    fclose(dataset->fData);
return;
}

void saveDataset3D_metaData(char fileName[], char fileNameBinary[], Dataset3D *dataset, int filePointers, int saveNewBinaryFile)
{
    FILE **pMet0, **pMet1, **pMet2, **pData;
    if(filePointers==FILE_POINTERS_INPUT)
    {
        pMet0 =  &(dataset->fMet0);
        pMet1 =  &(dataset->fMet1);
        pMet2 =  &(dataset->fMet2);
        pData =  &(dataset->fData);
    }
    if(filePointers==FILE_POINTERS_OUTPUT)
    {
        pMet0 =  &(dataset->oMet0);
        pMet1 =  &(dataset->oMet1);
        pMet2 =  &(dataset->oMet2);
        pData =  &(dataset->oData);
    }
    if(saveNewBinaryFile)
        createOdisseiaFilesNew(fileName, pMet0, pMet1, pMet2, pData, ODISSEIA_FILE_TYPE_SEISMIC);
    else
        createOdisseiaFilesNew_ExistingBinary(fileName, fileNameBinary, pMet0, pMet1, pMet2, pData, ODISSEIA_FILE_TYPE_SEISMIC);

    fillOdisseiaMetaFile_Seismic_ShotInfo(*pMet1, dataset);
    fillOdisseiaMetaFile_Seismic_ShotCoords(*pMet2, dataset);
    
return;
}
void saveDataset3D_binaryData(Dataset3D *dataset, int filePointers, int ishot)
{
    FILE *pMet0, *pMet1, *pMet2, *pData;
    if(filePointers==FILE_POINTERS_INPUT)     pData =  dataset->fData;
    if(filePointers==FILE_POINTERS_OUTPUT)    pData =  dataset->oData;

    Shot3D *shot = &((dataset->shot)[ishot]);
    writeSeismogram2File(pData, shot->source, shot->seismogram, shot->nsou, shot->nrec, shot->nt, shot->nStationsPreceding);
return;
}
void closeOdisseiaFiles(Dataset3D *dataset, int filePointers)
{
    FILE *pMet0, *pMet1, *pMet2, *pData;
    if(filePointers==FILE_POINTERS_INPUT) 
    {
        pMet0 =  dataset->fMet0;
        pMet1 =  dataset->fMet1;
        pMet2 =  dataset->fMet2;
        pData =  dataset->fData;
    }
    if(filePointers==FILE_POINTERS_OUTPUT)
    {
        pMet0 =  dataset->oMet0;
        pMet1 =  dataset->oMet1;
        pMet2 =  dataset->oMet2;
        pData =  dataset->oData;
    }
    fclose(pMet0);
    fclose(pMet1);
    fclose(pMet2);
    fclose(pData);
}

void initDataset3DfromFile(Dataset3D *dataset, Parms *parms)
{    
    // int OdisseiaFileType;
    printf("\n Input seismic file: <%s> \n", parms->inputSeismicFileName);
    openOdisseiaFiles(&(dataset->fMet0), &(dataset->fMet1), &(dataset->fMet2), \
                      &(dataset->fData), parms->inputSeismicFileName, \
                      (dataset->OdisseiaFileType), "r");
    int nshot;
    readOdisseiaMetaFile_get_nshot(dataset->fMet1, &nshot);
    dataset->nshot = nshot;
    dataset->shot  = CPU_zaloc1Shot3D(nshot);

    initShotFilteringParms(dataset, parms);

return;
}

void initDataset3DfromSEGY(Dataset3D *dataset, int nshot, int nshot_out, int shot_step, \
                           int *nrec, int *nsou, int **gatherTraces, \
                           int nt, float dt, float ot, float maxFreq, \
                           float *xs, float *ys, float *zs, \
                           float *xr, float *yr, float *zr)
{   
    // printf("\n [initDataset3DfromSEGY] nshot: %d \n", nshot);
    // fflush(stdout);

    dataset->nshot = nshot_out;
    dataset->shot  = CPU_zaloc1Shot3D(nshot_out);

    int nStationsPreceding = 0;
    int ishot;
    for(ishot=0; ishot<nshot; ishot+=shot_step)
    {
        int ishot_out = (ishot/shot_step);
        Shot3D *shot = &((dataset->shot)[ishot_out]);

        int ns = nsou[ishot];
        int nr = nrec[ishot];

        shot->ishot        = ishot_out;
        shot->nshot        = nshot_out;
        shot->nsou         = ns;
        shot->nrec         = nr;
        shot->nt           = nt;
        shot->dt           = dt;
        shot->ot           = ot;
        shot->maxFrequency = maxFreq;

        shot->xs  = CPU_zaloc1F(ns);
        shot->ys  = CPU_zaloc1F(ns);
        shot->zs  = CPU_zaloc1F(ns);
        shot->xr  = CPU_zaloc1F(nr);
        shot->yr  = CPU_zaloc1F(nr);
        shot->zr  = CPU_zaloc1F(nr);

        shot->nStationsPreceding = nStationsPreceding;

        int trc_i = gatherTraces[ishot][0];
        shot->xs[0] = xs[trc_i];
        shot->ys[0] = ys[trc_i];
        shot->zs[0] = zs[trc_i];
        int ir;
        for(ir=0; ir<nr; ir++)
        {
            int trc = gatherTraces[ishot][ir];
            shot->xr[ir] = xr[trc];
            shot->yr[ir] = yr[trc];
            shot->zr[ir] = zr[trc];
        }
        nStationsPreceding += ns + nr;
    }

return;
}
void initDataset3DfromSEGY_original(Dataset3D *dataset, int nshot, int nshot_out, int shot_step, \
                           int *nrec, int *nsou, int *gatherIni, int *gatherFin, \
                           int nt, float dt, float ot, float maxFreq, \
                           float *xs, float *ys, float *zs, \
                           float *xr, float *yr, float *zr)
{   
    printf("\n [initDataset3DfromSEGY] nshot: %d \n", nshot);
    fflush(stdout);

    dataset->nshot = nshot_out;
    dataset->shot  = CPU_zaloc1Shot3D(nshot_out);

    int nStationsPreceding = 0;
    int ishot;
    for(ishot=0; ishot<nshot; ishot+=shot_step)
    {
        int ishot_out = (ishot/shot_step);
        Shot3D *shot = &((dataset->shot)[ishot_out]);

        int ns = nsou[ishot];
        int nr = nrec[ishot];

        printf("\n [initDataset3DfromSEGY] ns: %d    nr: %d\n", ns, nr);

        shot->ishot        = ishot_out;
        shot->nshot        = nshot_out;
        shot->nsou         = ns;
        shot->nrec         = nr;
        shot->nt           = nt;
        shot->dt           = dt;
        shot->ot           = ot;
        shot->maxFrequency = maxFreq;

        shot->xs  = CPU_zaloc1F(ns);  printf("\n  shot->xs=%p  ns=%d\n", shot->xs, ns);
        shot->ys  = CPU_zaloc1F(ns);
        shot->zs  = CPU_zaloc1F(ns);
        shot->xr  = CPU_zaloc1F(nr);
        shot->yr  = CPU_zaloc1F(nr);
        shot->zr  = CPU_zaloc1F(nr);

        shot->nStationsPreceding = nStationsPreceding;

        int trc_i = gatherIni[ishot];
        int trc_f = gatherFin[ishot];
        printf("\n [initDataset3DfromSEGY] trc_i: %d    trc_f: %d\n", trc_i, trc_f);
        shot->xs[0] = xs[trc_i];
        shot->ys[0] = ys[trc_i];
        shot->zs[0] = zs[trc_i];
        int ir;
        for(ir=0; ir<nr; ir++)
        {
            shot->xr[ir] = xr[trc_i + ir];
            shot->yr[ir] = yr[trc_i + ir];
            shot->zr[ir] = zr[trc_i + ir];
        }
        printf("\n Checking: nStationsPreceding=%d  and shot->nStationsPreceding=%lld    difference=%d \n", \
                   nStationsPreceding, shot->nStationsPreceding, (int) (nStationsPreceding - shot->nStationsPreceding) );
        nStationsPreceding += ns + nr;
    }

return;
}

//* imhere
void initDataset3D_PZSum(Dataset3D *dataset_P, Dataset3D *dataset_Z, Dataset3D *dataset_PZ)
{
    int nshot_P  = dataset_P->nshot;
    int nshot_Z  = dataset_Z->nshot;

    // Find number of matching shots between P and Z
    dataset_P->shotMap_P2Z  = CPU_zaloc1I(nshot_P);
    dataset_P->shotMap_P2PZ = CPU_zaloc1I(nshot_P);

    int ishot, ishot_P, ishot_Z, ishot_PZ;
    int nshot_PZ = 0;
    for(ishot_P=0; ishot_P<nshot_P; ishot_P++)
    {
        Shot3D *shot_P = &((dataset_P->shot)[ishot_P]);
        
        // Find corresponding shot in Z
        ishot_Z  = -1;
        ishot_PZ = -1;
        for(ishot=0; ishot<nshot_Z; ishot++)
        {
            Shot3D *shot_Z = &((dataset_Z->shot)[ishot]);
            if( shot_P->xs[0] == shot_Z->xs[0]  &&  shot_P->ys[0] == shot_Z->ys[0]  &&  shot_P->zs[0] == shot_Z->zs[0])
            {
                ishot_Z = ishot;
                ishot_PZ = nshot_PZ;
                nshot_PZ++;
                break;
            }
        }
        // Check if the same ishot_Z has been selected previously
        for(ishot=0; ishot<ishot_P; ishot++)
        {
            if( (ishot_Z==dataset_P->shotMap_P2Z[ishot])  &&  (ishot_Z!=-1) ) 
            {
                ishot_Z  = -1;
                ishot_PZ = -1;
                nshot_PZ--;
                break;
            }
        }
        dataset_P->shotMap_P2Z [ishot_P] = ishot_Z;
        dataset_P->shotMap_P2PZ[ishot_P] = ishot_PZ;
        // printf("\n Building shot map - ishot_P=%d  ishot_Z=%d  ishot_PZ=%d  nshot_Z=%d\n", ishot_P, ishot_Z, ishot_PZ, nshot_Z);  fflush(stdout);
    }
    // End of shot mapping


    // Initialize PZ dataset
    dataset_PZ->nshot = nshot_PZ;
    dataset_PZ->shot  = CPU_zaloc1Shot3D(nshot_PZ);

    long long int nStationsPreceding = 0;
    for(ishot_P=0; ishot_P<nshot_P; ishot_P++)
    {
        if( (dataset_P->shotMap_P2Z)[ishot_P] == -1 )    continue;

        ishot_Z  = (dataset_P->shotMap_P2Z )[ishot_P];
        ishot_PZ = (dataset_P->shotMap_P2PZ)[ishot_P];

        Shot3D *shot_P  = &((dataset_P->shot) [ishot_P ]);
        Shot3D *shot_Z  = &((dataset_Z->shot) [ishot_Z ]);
        Shot3D *shot_PZ = &((dataset_PZ->shot)[ishot_PZ]);

        int nrec_P  = shot_P->nrec;
        int nrec_Z  = shot_Z->nrec;
        int nrec_PZ = 0;
        int irec, irec_P, irec_Z, irec_PZ;
        // Find matching receivers between P and Z
        shot_P->recMap_P2Z  = CPU_zaloc1I(nrec_P);
        shot_P->recMap_P2PZ = CPU_zaloc1I(nrec_P);
        for(irec_P=0; irec_P<nrec_P; irec_P++)
        {
            irec_Z  = -1;
            irec_PZ = -1;
            for(irec=0; irec<nrec_Z; irec++)
            {
                if( shot_P->xr[irec_P] == shot_Z->xr[irec]  &&  shot_P->yr[irec_P] == shot_Z->yr[irec]  &&  shot_P->zr[irec_P] == shot_Z->zr[irec])
                {
                    irec_Z  = irec;
                    irec_PZ = nrec_PZ;
                    nrec_PZ++;
                    break;
                }  
            }
            shot_P->recMap_P2Z [irec_P] = irec_Z;
            shot_P->recMap_P2PZ[irec_P] = irec_PZ;            
        }
        shot_PZ->nrec               = 2*nrec_PZ; // Multiplied by two (one for upgoing and one for downgoing wavefields)
        shot_PZ->nsou               = shot_P->nsou;
        shot_PZ->ishot              = ishot_PZ;
        shot_PZ->nshot              = nshot_PZ;
        shot_PZ->nt                 = shot_P->nt;
        shot_PZ->dt                 = shot_P->dt;
        shot_PZ->ot                 = shot_P->ot;
        shot_PZ->maxFrequency       = shot_P->maxFrequency;
        shot_PZ->nStationsPreceding = nStationsPreceding;
        nStationsPreceding += shot_PZ->nrec + shot_PZ->nsou;

        shot_PZ->xs  = CPU_zaloc1F(shot_PZ->nsou);
        shot_PZ->ys  = CPU_zaloc1F(shot_PZ->nsou);
        shot_PZ->zs  = CPU_zaloc1F(shot_PZ->nsou);
        shot_PZ->xr  = CPU_zaloc1F(shot_PZ->nrec);
        shot_PZ->yr  = CPU_zaloc1F(shot_PZ->nrec);
        shot_PZ->zr  = CPU_zaloc1F(shot_PZ->nrec);
        
        // This cannot stay here, otherwise will allocate memory for all shots 
        // "simultaneously", which is too much memory.
        // shot_PZ->seismogram = CPU_zaloc1F(shot_PZ->nrec * shot_PZ->nt);
        // shot_PZ->source     = CPU_zaloc1F(shot_PZ->nsou * shot_PZ->nt);

        int is;
        for(is=0; is<shot_PZ->nsou; is++)
        {
            shot_PZ->xs[is] = shot_P->xs[is];
            shot_PZ->ys[is] = shot_P->ys[is];
            shot_PZ->zs[is] = shot_P->zs[is];
        }
        for(irec_P=0; irec_P<nrec_P; irec_P++)
        {
            irec_PZ = shot_P->recMap_P2PZ[irec_P];
            int operSign;
            int i_PZ;
            if(irec_PZ != -1)
            {
                operSign = +1;
                i_PZ = irec_PZ + 0 * nrec_PZ;
                shot_PZ->xr[i_PZ] = shot_P->xr[irec_P];
                shot_PZ->yr[i_PZ] = shot_P->yr[irec_P];
                shot_PZ->zr[i_PZ] = shot_P->zr[irec_P] * operSign;

                i_PZ = irec_PZ + 1 * nrec_PZ;
                operSign = -1;
                shot_PZ->xr[i_PZ] = shot_P->xr[irec_P];
                shot_PZ->yr[i_PZ] = shot_P->yr[irec_P];
                shot_PZ->zr[i_PZ] = shot_P->zr[irec_P] * operSign;
            }
        }
    }

return;
}


void PZSum_sensor(Shot3D *shot_P, Shot3D *shot_Z, Shot3D *shot_PZ, float *amplitude_P, float *amplitude_Z, \
                  long long int *trace2Station, long long int *traceNumber, float *source)
{

    int nrec_P  = shot_P->nrec;
    int nrec_Z  = shot_Z->nrec;
    int nrec_PZ = shot_PZ->nrec/2;
    int nt = shot_P->nt;
    int it;
    int irec, irec_P, irec_Z, irec_PZ;
    int operSign;

    // Allocating PZ seismogram and source arrays
    shot_PZ->seismogram = CPU_zaloc1F(shot_PZ->nrec * shot_PZ->nt);
    shot_PZ->source     = CPU_zaloc1F(shot_PZ->nsou * shot_PZ->nt);
        
    // Filling seismogram array
    for(irec_P=0; irec_P<nrec_P; irec_P++)
    {
        long long int stationNumber = trace2Station[*traceNumber];
        float amplitudeFactor_Z, amplitudeFactor_P;
        if(amplitude_Z[stationNumber]>0.0) {
            amplitudeFactor_P = 1.0;
            amplitudeFactor_Z = amplitude_P[stationNumber]/amplitude_Z[stationNumber];
        }   
        else {
            amplitudeFactor_P = 0.0;
            amplitudeFactor_Z = 0.0;
        }
        irec_Z  = shot_P->recMap_P2Z [irec_P];
        irec_PZ = shot_P->recMap_P2PZ[irec_P];
        if(irec_Z != -1)
        {
            int i_P, i_Z, i_PZ;
            for(it=0; it<nt; it++)
            {
                i_P  = it + irec_P*nt;
                i_Z  = it + irec_Z*nt;

                i_PZ = it + (irec_PZ+0*nrec_PZ)*nt;
                operSign = +1;
                shot_PZ->seismogram[i_PZ] = amplitudeFactor_P * shot_P->seismogram[i_P] \
                                          + operSign * amplitudeFactor_Z * shot_Z->seismogram[i_Z];

                i_PZ = it + (irec_PZ+1*nrec_PZ)*nt;
                operSign = -1;
                shot_PZ->seismogram[i_PZ] = amplitudeFactor_P * shot_P->seismogram[i_P] \
                                          + operSign * amplitudeFactor_Z * shot_Z->seismogram[i_Z];
            }
        }

        *traceNumber += 1;
    }

    // Filling source array
    int is;
    for(is=0; is<shot_PZ->nsou; is++)
        for(it=0; it<nt; it++) 
        {
            int i = it + is*nt;
            if(source==NULL)
                shot_PZ->source[i] = shot_P->source[i];
            else
                shot_PZ->source[i] = source[it];
        }


return;
}

void PZSum_shot(Shot3D *shot_P, Shot3D *shot_Z, Shot3D *shot_PZ, float factor_P, float factor_Z, float *source)
{

    int nrec_P  = shot_P->nrec;
    int nrec_Z  = shot_Z->nrec;
    int nrec_PZ = shot_PZ->nrec/2;
    int nt = shot_P->nt;
    int it;
    int irec, irec_P, irec_Z, irec_PZ;
    int operSign;

    // Allocating PZ seismogram and source arrays
    shot_PZ->seismogram = CPU_zaloc1F(shot_PZ->nrec * shot_PZ->nt);
    shot_PZ->source     = CPU_zaloc1F(shot_PZ->nsou * shot_PZ->nt);
        
    // Filling seismogram array
    for(irec_P=0; irec_P<nrec_P; irec_P++)
    {
        // printf("\n  mapping rec indexes \n");  fflush(stdout);
        irec_Z  = shot_P->recMap_P2Z [irec_P];
        irec_PZ = shot_P->recMap_P2PZ[irec_P];
        // printf("\n irec_P=%d  irec_Z=%d  irec_PZ=%d  nrec_P=%d  nrec_Z=%d  nrec_PZ=%d", irec_P, irec_Z, irec_PZ, nrec_P, nrec_Z, nrec_PZ);  fflush(stdout);
        if(irec_Z != -1)
        {
            int i_P, i_Z, i_PZ;
            for(it=0; it<nt; it++)
            {
                i_P  = it + irec_P*nt;
                i_Z  = it + irec_Z*nt;

                i_PZ = it + (irec_PZ+0*nrec_PZ)*nt;
                operSign = +1;
                shot_PZ->seismogram[i_PZ] = factor_P * shot_P->seismogram[i_P] + operSign * factor_Z * shot_Z->seismogram[i_Z];

                i_PZ = it + (irec_PZ+1*nrec_PZ)*nt;
                operSign = -1;
                shot_PZ->seismogram[i_PZ] = 1.0 * (factor_P * shot_P->seismogram[i_P] + operSign * factor_Z * shot_Z->seismogram[i_Z]);
            }
        }
    }

    // printf("\n  Finished sum on receivers \n");  fflush(stdout);
    // Filling source array
    int is;
    for(is=0; is<shot_PZ->nsou; is++)
        for(it=0; it<nt; it++) 
        {
            int i = it + is*nt;
            if(source==NULL)
                shot_PZ->source[i] = shot_P->source[i];
            else
                shot_PZ->source[i] = source[it];
        }


return;
}

void PZSum(Dataset3D *dataset_P, Dataset3D *dataset_Z, Dataset3D *dataset_PZ)
{

    int nshot_P  = dataset_P->nshot;
    int nshot_Z  = dataset_Z->nshot;
    int nshot_PZ = dataset_PZ->nshot;

    int ishot, ishot_P, ishot_Z, ishot_PZ;
    nshot_PZ = 0;
    long long int nStationsPreceding;

    for(ishot_P=0; ishot_P<nshot_P; ishot_P++)
    {
        if( (dataset_P->shotMap_P2Z)[ishot_P] == -1 )    continue;

        ishot_Z  = (dataset_P->shotMap_P2Z )[ishot_P];
        ishot_PZ = (dataset_P->shotMap_P2PZ)[ishot_P];

        Shot3D *shot_P  = &((dataset_P->shot) [ishot_P ]);
        Shot3D *shot_Z  = &((dataset_Z->shot) [ishot_Z ]);
        Shot3D *shot_PZ = &((dataset_PZ->shot)[ishot_PZ]);

        int nrec_P  = shot_P->nrec;
        int nrec_Z  = shot_Z->nrec;
        int nrec_PZ = shot_PZ->nrec;
        int nt = shot_P->nt;
        int it;
        int irec, irec_P, irec_Z, irec_PZ;
        int operSign;
        
        // Filling seismogram array
        for(irec_P=0; irec_P<nrec_P; irec_P++)
        {
            irec_Z  = shot_P->recMap_P2Z [irec_P];
            irec_PZ = shot_P->recMap_P2PZ[irec_P];
            if(irec_Z != -1)
            {
                int i_P, i_Z, i_PZ;
                for(it=0; it<nt; it++)
                {
                    i_P = it + irec_P*nt;
                    i_Z = it + irec_Z*nt;
                
                    i_PZ = it + (irec_PZ+0*nrec_PZ)*nt;
                    operSign = +1;
                    shot_PZ->seismogram[i_PZ] = shot_P->seismogram[i_P] + operSign*shot_Z->seismogram[i_Z];

                    i_PZ = it + (irec_PZ+0*nrec_PZ)*nt;
                    operSign = -1;
                    shot_PZ->seismogram[i_PZ] = shot_P->seismogram[i_P] + operSign*shot_Z->seismogram[i_Z];

                }
            }
        }
        // Filling source array
        int is;
        for(is=0; is<shot_PZ->nsou; is++)
            for(it=0; it<nt; it++) 
            {
                int i = it + is*nt;
                shot_PZ->source[i] = shot_P->source[i];
            }


    }

return;
}


/*
Sensor* mapSensors(Dataset3D *dataset, long long int *nsensors, long long int **trace2Station)
{
    int nshot  = dataset->nshot;
    int ishot, nrec, irec;

    long long int totalNumberOfTraces = 0;
    for(ishot=0; ishot<nshot; ishot++)
        totalNumberOfTraces += (dataset->shot)[ishot].nrec;
    
    *trace2Station = (long long int *) calloc( totalNumberOfTraces, sizeof(long long int) );

    long long int nstations = 100*nshot; // this number was arbitrarily chosen, but should be large enough for most OBS datasets
    Sensor *sensor = new Sensor[nstations];

    printf("\n DEBUGGING: totalNumberOfTraces=%lld  nshot=%d  nstations=%lld \n", totalNumberOfTraces, nshot, nstations);  fflush(stdout);

    long long int ist;
    long long int nst = 0;
    long long int traceNumber = 0;
    printf("\n Extracting station from shot:        "); fflush(stdout);
    for(ishot=0; ishot<nshot; ishot++)
    {
        if(ishot%100==0)  { 
            backSpace(7); fflush(stdout);
            printf("%7d", ishot); fflush(stdout); 
        }
        Shot3D *shot = &((dataset->shot)[ishot]);
        nrec = shot->nrec;
        for(irec=0; irec<nrec; irec++)
        {
            long long int stationNumber = -1;
            double x = shot->xr[irec];
            double y = shot->yr[irec];
            double z = shot->zr[irec];
            
            for(ist=0; ist<nst; ist++)
            {
                stationNumber = sensor[ist].addTrace(x, y, z, traceNumber);
                if(stationNumber==ist) {
                    (*trace2Station)[traceNumber] = stationNumber;
                    break;
                }
            }
            if(stationNumber==-1)
            {
                sensor[nst] = Sensor(x, y, z, nst, traceNumber);
                stationNumber = sensor[nst].stationNumber;
                (*trace2Station)[traceNumber] = stationNumber;
                nst++;
            }
            // printf("\n  (*trace2Station)[%lld] = %lld  \n", traceNumber, (*trace2Station)[traceNumber]);  fflush(stdout);
            traceNumber++;
        } // irec loop
    } // ishot loop

    *nsensors = nst;

    printf("\n DEBUGGING 2: nst=%lld \n", nst);  fflush(stdout);

return sensor;
}
float*  get_SignatureFromDA_Sensor(Dataset3D *dataset, Sensor *sensor, long long int nsensors, long long int *trace2Station, Parms *parms, int applyCosineFactor)
{
 
    Shot3D *shot  = &((dataset->shot)[0]);
    int   nt = shot->nt;
    float dt = shot->dt;
    float ot = shot->ot;

    float *trace_time = CPU_zaloc1F(nt*nsensors);
    
    int nshot  = dataset->nshot;
    int ishot, it, nrec, irec;
    float velWater = 1490.0;
    float offMin=parms->offMin, offMax=parms->offMax;
    float srcWin_ti       = parms->srcWin_ti;
    float srcWin_ti_taper = parms->srcWin_ti_taper;
    float srcWin_tf       = parms->srcWin_tf;
    float srcWin_tf_taper = parms->srcWin_tf_taper;

    long long int traceNum, isensor;

    long long int *n_contributions = (long long int *)  calloc(nsensors, sizeof(long long int));
    traceNum = 0;
    printf("\n Extracting signature from shot:        "); fflush(stdout);
    for(ishot=0; ishot<nshot; ishot++)
    {
        if(ishot%100==0)  { 
            backSpace(7); fflush(stdout);
            printf("%7d", ishot); fflush(stdout); 
        }
        loadSeismogram2Dataset(dataset, ishot, 0);
        Shot3D *shot  = &((dataset->shot)[ishot]);
        nrec = shot->nrec;
        for(irec=0; irec<nrec; irec++)  
        {
            long long int shiftMem = irec * nt;
            float *trace_TimeShift = timeShiftDA(shot->seismogram+shiftMem, nt, dt, ot, \
                                                 shot->xs[0], shot->ys[0], shot->zs[0], \
                                                 shot->xr[irec], shot->yr[irec], shot->zr[irec], \
                                                 offMin, offMax, srcWin_ti, srcWin_ti_taper, \
                                                 srcWin_tf, srcWin_tf_taper, velWater, applyCosineFactor);
            if(trace_TimeShift != NULL)
            {
                isensor = trace2Station[traceNum];
                for(it=0; it<nt; it++)
                    trace_time[it + isensor*nt] += trace_TimeShift[it];
                free(trace_TimeShift);
                n_contributions[isensor] += 1;
            }
            traceNum++;
        }
        finalizeSeismogramAndSource(dataset, ishot);
    }
    for(isensor=0; isensor<nsensors; isensor++)
    {
        if(n_contributions[isensor]>0)
            for(it=0; it<nt; it++)    trace_time[it + isensor*nt] /= n_contributions[isensor];
    }
    free(n_contributions);

return trace_time;
}
//*/
void getShotRecNum(Dataset3D *dataset, long long int traceNum, long long int *ishot_, long long int *irec_)
{
    int nshot  = dataset->nshot;
    int ishot, nrec, irec;

    long long int tracesRead = 0;
    for(ishot=0; ishot<nshot; ishot++)
    {
        nrec = (dataset->shot)[ishot].nrec;
        long long int remainder = traceNum-tracesRead;
        if(remainder<nrec) {
            *ishot_ = ishot;
            *irec_  = remainder;
        }   
        else tracesRead += nrec;
    }
return;
}
float*  get_SignatureFromDA(Dataset3D *dataset, Parms *parms)
{
 
    Shot3D *shot  = &((dataset->shot)[0]);
    int   nt = shot->nt;
    float dt = shot->dt;
    float ot = shot->ot;

    float *trace_time = CPU_zaloc1F(nt);
    
    int nshot  = dataset->nshot;
    int ishot, it, nrec, irec, ifreq, nfreq;
    float dfreq;
    float velWater = 1490.0;
    float offMin=parms->offMin, offMax=parms->offMax;
    float srcWin_ti       = parms->srcWin_ti;
    float srcWin_ti_taper = parms->srcWin_ti_taper;
    float srcWin_tf       = parms->srcWin_tf;
    float srcWin_tf_taper = parms->srcWin_tf_taper;

    long long int n_contributions = 0;

    // Forming trace representative of the dataset
    printf("\n Extracting signature from shot:        "); fflush(stdout);
    for(ishot=0; ishot<nshot; ishot++)
    {
        if(ishot%100==0)  { 
            backSpace(7); fflush(stdout);
            printf("%7d", ishot); fflush(stdout); 
        }
        loadSeismogram2Dataset(dataset, ishot, 0);
        Shot3D *shot  = &((dataset->shot)[ishot]);
        nrec = shot->nrec;
        n_contributions += buildSourceFromDirectArrival(trace_time, shot->seismogram, nt, dt, ot, nrec, \
                                                        shot->xs, shot->ys, shot->zs, shot->xr, shot->yr, shot->zr, \
                                                        offMin, offMax, srcWin_ti, \
                                                        srcWin_ti_taper, srcWin_tf, srcWin_tf_taper, velWater);
        finalizeSeismogramAndSource(dataset, ishot);
    }

    if(n_contributions>0)
        for(it=0; it<nt; it++)    trace_time[it] /= n_contributions;

return trace_time;
}

void Take1stDerivative_singleTrace_TimeDomain(float *trace_time, int nt, float dt, float ot)
{
    
    float coef[6];
    coef[0] = 0.0;
    coef[1] = (+1.66666667);
    coef[2] = (-0.47619048);
    coef[3] = (+0.11904762);
    coef[4] = (-0.01984127);
    coef[5] = (+0.00158730);

    int it, ifreq, i;
    float *trace_tmp = CPU_zaloc1F(nt);
    int DL = 5;
    for(it=DL; it<nt-DL; it++) {
        for(i=1; i<=DL; i++)
            trace_tmp[it] +=  coef[i] * (trace_time[it+i] - trace_time[it-i]) / dt;
    }
    for(it=1; it<nt; it++) {
        trace_time[it] = trace_tmp[it-1];
    }
    
return;
}

void ApplyTaperToTrace(float *trace_time, int nt, float dt, float ot, float t_ini, float t_fin, float taper_length)
{
    
    float taper, t;

    if(t_ini == -999999.0)    t_ini = ot;              // trace start
    if(t_fin == -999999.0)    t_fin = ot + (nt-1)*dt;  // trace end

    int it;
    for(it=0; it<nt; it++) {
        t = ot + it*dt;
        if(t<t_ini || t>t_fin)               taper = 0.0;
        else if( (t-t_ini)<taper_length ) {
            float arg = 0.5*M_PI * (t-t_ini)/taper_length;
            taper = sinf(arg) * sinf(arg);
        }
        else if( (t_fin-t)<taper_length ) {
            float arg = 0.5*M_PI * (t_fin-t)/taper_length;
            taper = sinf(arg) * sinf(arg);
        }
        
        trace_time[it] *= taper;
    }
    
return;
}

float getMaximumAmplitude_singleTrace(float *trace_time, int nt)
{
    int it;
    float maxAmp = 0.0;
    for(it=0; it<nt; it++) {
        if( maxAmp<fabsf(trace_time[it]) )    maxAmp = fabsf(trace_time[it]);
    }
return maxAmp;
}

void ZeroOutNoisyTrace(Shot3D *shot)
{
    int   nt = shot->nt;
    float dt = shot->dt;
    float ot = shot->ot;
    int   irec, it;
    int   nrec  = shot->nrec;

    for(irec=0; irec<nrec; irec++)
    {
        // Find average value of first 0.2 seconds of the trace
        double average_top    = 0.0;
        double average_bottom = 0.0;
        
        int nt0 = (int) (0.5 + 0.2/dt);
        for(it=0; it<nt0; it++)
            average_top += fabs( (shot->seismogram)[it + irec*nt] ) / nt0;

        int it1 = nt/2;
        int nt1 = (nt-nt/2);
        for(it=it1; it<nt; it++)
            average_bottom += fabs( (shot->seismogram)[it + irec*nt] ) / nt1;
        
        if(average_top >= 0.8*average_bottom)
        {
            for(it=0; it<nt; it++)
                (shot->seismogram)[it + irec*nt] = 0.0;
        }
    }
        
return;
}
void Take1stDerivative_shot_TimeDomain(Shot3D *shot)
{
    int   nt = shot->nt;
    float dt = shot->dt;
    float ot = shot->ot;
    int irec;
    int nrec  = shot->nrec;
    // Taking derivative of seismogram
    for(irec=0; irec<nrec; irec++)
        Take1stDerivative_singleTrace_TimeDomain((shot->seismogram)+irec*nt, nt, dt, ot);
return;
}
/*
void Take1stDerivative_shot(Shot3D *shot, float multiplyingFactor)
{
    int nshot  = shot->nshot;
    int ishot  = shot->ishot;

    int   nt = shot->nt;
    float dt = shot->dt;
    float ot = shot->ot;

    int it, irec, ifreq;
    int nrec  = shot->nrec;
    int nfreq = nt;
    float dfreq = 1.0f/(nt*dt);

    fftw_complex *seismo_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt * nrec);
    applyFFTW1D_forwardR2I(seismo_freq, shot->seismogram, nt, nrec);
    
    // Taking derivative of seismogram
    float RE, IM;
    for(irec=0; irec<nrec; irec++)
    {
        float freq = ifreq * dfreq;
        for(ifreq=0; ifreq<=nfreq/2; ifreq++)
        {
            int i = ifreq + irec*nfreq;
            seismo_freq[i][0] *= (+freq) * multiplyingFactor;
            seismo_freq[i][1] *= (+freq) * multiplyingFactor;

            RE = seismo_freq[i][0];
            IM = seismo_freq[i][1];
            seismo_freq[i][0] *= (+freq) * (-IM) * multiplyingFactor;
            seismo_freq[i][1] *= (+freq) * (+RE) * multiplyingFactor;

            if(ifreq>0 && ifreq<nfreq/2) {
                int i = nfreq-ifreq + irec*nfreq;
                RE = seismo_freq[i][0];
                IM = seismo_freq[i][1];
                seismo_freq[i][0] *= (-freq) * (-IM) * multiplyingFactor;
                seismo_freq[i][1] *= (-freq) * (+RE) * multiplyingFactor;
            }
        }
    }

    applyFFTW1D_backwardI2R(shot->seismogram, seismo_freq, nt, nrec);
    fftw_free(seismo_freq);

return;
}
*/

void Deconvolve_shot(Shot3D *shot, float *trace_time, int mode_convolution, float white_noise_perc)
{
    int nshot  = shot->nshot;
    int ishot  = shot->ishot;

    int   nt = shot->nt;
    float dt = shot->dt;
    float ot = shot->ot;

    int it, irec, ifreq;
    int nrec  = shot->nrec;
    int nfreq = nt;
    float dfreq = 1.0f/(nt*dt);

    printf("\n WWW 01  nt=%d dt=%f ot=%f nrec=%d dfreq=%f \n", nt, dt, ot, nrec, dfreq); fflush(stdout);

    // Transforming to frequency-domain
    fftw_complex *trace_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt);
    applyFFTW1D_forwardR2I(trace_freq, trace_time, nt, 1);
    float *amp_freq = CPU_zaloc1F(nfreq);
    float *phi_freq = CPU_zaloc1F(nfreq);
        
    float maxAmp = 0.0;
    float RE, IM;
    for(ifreq=0; ifreq<nfreq; ifreq++)
    {
        RE = trace_freq[ifreq][0];
        IM = trace_freq[ifreq][1];
        amp_freq[ifreq] = sqrtf(RE*RE + IM*IM);
        phi_freq[ifreq] = atan2f(IM,RE);
        if(maxAmp < amp_freq[ifreq])    maxAmp = amp_freq[ifreq];
    }

    fftw_complex *seismo_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt * nrec);
    applyFFTW1D_forwardR2I(seismo_freq, shot->seismogram, nt, nrec);
        
    // Deconvolving seismogram
    for(irec=0; irec<nrec; irec++)
    {
        float *amp = CPU_zaloc1F(nfreq);
        float *phi = CPU_zaloc1F(nfreq);
        for(ifreq=0; ifreq<nfreq; ifreq++)
        {
            int i = ifreq + irec*nfreq;
            float RE = seismo_freq[i][0];
            float IM = seismo_freq[i][1];
                
            if( mode_convolution == MODE_CONVOLUTION_DECONVOLVE) 
            {
                float WN = white_noise_perc * maxAmp;
                amp[ifreq] = sqrtf(RE*RE + IM*IM) / (amp_freq[ifreq] + WN);
                phi[ifreq] = atan2f(IM,RE) - phi_freq[ifreq];
            }
            else if( mode_convolution == MODE_CONVOLUTION_CONVOLVE) 
            {
                amp[ifreq] = sqrtf(RE*RE + IM*IM) * amp_freq[ifreq];
                phi[ifreq] = atan2f(IM,RE) + phi_freq[ifreq];
            }

            seismo_freq[i][0] = amp[ifreq] * cosf(phi[ifreq]);
            seismo_freq[i][1] = amp[ifreq] * sinf(phi[ifreq]);
        }
        free(amp);
        free(phi);
    }

    applyFFTW1D_backwardI2R(shot->seismogram, seismo_freq, nt, nrec);

    fftw_free(seismo_freq);    
    fftw_free(trace_freq);
    free(amp_freq);
    free(phi_freq);

return;
}

void Deconvolve_Data(Dataset3D *dataset, float *trace_time, int mode_convolution, float white_noise_perc)
{
 
    Shot3D *shot  = &((dataset->shot)[0]);
    int   nt = shot->nt;
    float dt = shot->dt;
    float ot = shot->ot;

    int nshot  = dataset->nshot;
    int ishot, it, nrec, irec, ifreq, nfreq;
    float dfreq;

    // Transforming to frequency-domain
    fftw_complex *trace_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt);
    applyFFTW1D_forwardR2I(trace_freq, trace_time, nt, 1);
    float *amp_freq = CPU_zaloc1F(nfreq);
    float *phi_freq = CPU_zaloc1F(nfreq);
        
    float maxAmp = 0.0;
    float RE, IM;
    for(ifreq=0; ifreq<nfreq; ifreq++)
    {
        RE = trace_freq[ifreq][0];
        IM = trace_freq[ifreq][1];
        amp_freq[ifreq] = sqrtf(RE*RE + IM*IM);
        phi_freq[ifreq] = atan2f(IM,RE);
        if(maxAmp < amp_freq[ifreq])    maxAmp = amp_freq[ifreq];
    }

    // Deconvolution loop
    for(ishot=0; ishot<nshot; ishot++)
    {
        Shot3D *shot  = &((dataset->shot)[ishot]);

        nrec  = shot->nrec;
        nfreq = nt;
        dfreq = 1.0f/(nt*dt);

        fftw_complex *seismo_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt * nrec);

        applyFFTW1D_forwardR2I(seismo_freq, shot->seismogram, nt, nrec);
        
        // Deconvolving seismogram
        for(irec=0; irec<nrec; irec++)
        {
            float *amp = CPU_zaloc1F(nfreq);
            float *phi = CPU_zaloc1F(nfreq);
            for(ifreq=0; ifreq<nfreq; ifreq++)
            {
                int i = ifreq + irec*nfreq;
                float RE = seismo_freq[i][0];
                float IM = seismo_freq[i][1];
                
                if( mode_convolution == MODE_CONVOLUTION_DECONVOLVE) 
                {
                    float WN = white_noise_perc * maxAmp;
                    amp[ifreq] = sqrtf(RE*RE + IM*IM) / (amp_freq[ifreq] + WN);
                    phi[ifreq] = atan2f(IM,RE) - phi_freq[ifreq];
                }
                else if( mode_convolution == MODE_CONVOLUTION_CONVOLVE) 
                {
                    amp[ifreq] = sqrtf(RE*RE + IM*IM) * amp_freq[ifreq];
                    phi[ifreq] = atan2f(IM,RE) + phi_freq[ifreq];
                }

                seismo_freq[i][0] = amp[ifreq] * cosf(phi[ifreq]);
                seismo_freq[i][1] = amp[ifreq] * sinf(phi[ifreq]);
            }
            free(amp);
            free(phi);
        }

        applyFFTW1D_backwardI2R(shot->seismogram, seismo_freq, nt, nrec);

        fftw_free(seismo_freq);
    }

    
    fftw_free(trace_freq);
    free(trace_time);
    free(amp_freq);
    free(phi_freq);

return;
}

// void PZEqualization(Dataset3D *dataset)
// {
 
//     Shot3D *shot  = &((dataset->shot)[0]);
//     int   nt = shot->nt;
//     float dt = shot->dt;
//     float ot = shot->ot;

//     float *trace_time = CPU_zaloc1F(nt);
    
//     int nshot  = dataset->nshot;
//     int ishot, it, nrec, irec, ifreq, nfreq;
//     float dfreq;
//     float velWater = 1490.0;
//     float offMin=0.0, offMax=3000.0;

//     // Forming trace representative of the dataset
//     for(ishot=0; ishot<nshot; ishot++)
//     {
//         Shot3D *shot  = &((dataset->shot)[ishot]);
//         int nrec = shot->nrec;        
//         buildSourceFromDirectArrival(trace_time, shot->seismogram, nt, dt, ot, nrec, \
//                                      shot->xs, shot->ys, shot->zs, shot->xr, shot->yr, shot->zr, \
//                                      offMin, offMax, 0.0, 0.0, nt*dt, 0.0, velWater);
//     }

//     // Transforming to frequency-domain
//     fftw_complex *trace_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt);
//     applyFFTW1D_forwardR2I(trace_freq, trace_time, nt, 1);
//     float *amp_freq = CPU_zaloc1F(nfreq);
//     float *phi_freq = CPU_zaloc1F(nfreq);
        
//     float RE, IM;
//     for(ifreq=0; ifreq<nfreq; ifreq++)
//     {
//         RE = trace_freq[ifreq][0];
//         IM = trace_freq[ifreq][1];
//         amp_freq[ifreq] = sqrtf(RE*RE + IM*IM);
//         phi_freq[ifreq] = atan2f(IM,RE);
//     }

//     // Deconvolution loop
//     for(ishot=0; ishot<nshot; ishot++)
//     {
//         Shot3D *shot  = &((dataset->shot)[ishot]);

//         nrec  = shot->nrec;
//         nfreq = nt;
//         dfreq = 1.0f/(nt*dt);

//         fftw_complex *seismo_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt * nrec);

//         applyFFTW1D_forwardR2I(seismo_freq, shot->seismogram, nt, nrec);
        
//         // Deconvolving seismogram
//         for(irec=0; irec<nrec; irec++)
//         {
//             float *amp = CPU_zaloc1F(nfreq);
//             float *phi = CPU_zaloc1F(nfreq);
//             for(ifreq=0; ifreq<nfreq; ifreq++)
//             {
//                 int i = ifreq + irec*nfreq;
//                 float RE = seismo_freq[i][0];
//                 float IM = seismo_freq[i][1];
                
//                 amp[ifreq] = sqrtf(RE*RE + IM*IM) / amp_freq[ifreq];
//                 phi[ifreq] = atan2f(IM,RE) - phi_freq[ifreq];

//                 seismo_freq[i][0] = amp[ifreq] * cosf(phi[ifreq]);
//                 seismo_freq[i][1] = amp[ifreq] * sinf(phi[ifreq]);
//             }
//             free(amp);
//             free(phi);
//         }

//         applyFFTW1D_backwardI2R(shot->seismogram, seismo_freq, nt, nrec);

//         fftw_free(seismo_freq);
//     }

    
//     fftw_free(trace_freq);
//     free(trace_time);
//     free(amp_freq);
//     free(phi_freq);

// return;
// }

long long int buildSourceFromDirectArrival(float *source, float *seismogram, int nt, float dt, float ot, int ntrc, \
                                 float *xs, float *ys, float *zs, float *xr, float *yr, float *zr, float offMin, float offMax, \
                                 float srcWin_ti, float srcWin_ti_taper, float srcWin_tf, float srcWin_tf_taper, float velWater)
{
    // Sinc interpolation parameters
    int r = 4;
    int r_ = r - 1;
    int diam = 2*r;

    int itrc, it;

    float t1 = srcWin_ti;
    float t2 = srcWin_ti + srcWin_ti_taper;
    float t3 = srcWin_tf - srcWin_tf_taper;
    float t4 = srcWin_tf;

    long long int n_contribution = 0;

    int step = 1;
    for(itrc=0; itrc<ntrc; itrc+=step)
    {
        float delta_x = xs[0] - xr[itrc];
        float delta_y = ys[0] - yr[itrc];
        float delta_z = zs[0] - zr[itrc];

        float offset = sqrtf(delta_x*delta_x + delta_y*delta_y);

        if( !( offMin<=offset  &&  offset<=offMax ) )    continue;

        float R = sqrtf(delta_x*delta_x + delta_y*delta_y + delta_z*delta_z);
        float T = R/velWater;

        n_contribution++;

        for(it=0; it<nt; it++)
        {
            float t = ot + it*dt;
            float t_shifted = t - T;
            if( t_shifted<t1  ||  t_shifted>t4 )    continue;
            float it0 = floorf((t_shifted-ot)/dt);
            float alpha = ( t_shifted - (it0*dt+ot) ) / dt;

            float sample = seismogram[itrc*nt + it];

            float *coef = interSinc1D_8P(alpha);
            int i;
            for(i=0; i<diam; i++)
            {
                int it_ = it0 - r_ + i;
                if(it_<0  ||  it_>=nt)    continue;
                source[it_] += coef[i] * sample;
            }
            free(coef);
        }
    }

    for(it=0; it<nt; it++)
    {       
        float t = ot + it*dt;
        float filter;

        if     ( t < t1  ||  t> t4 ) {
            filter = 0.0f;
        }   
        else if( t1<=t   &&  t< t2 ) {
            float arg = 0.5*M_PI * (t-t1)/(t2-t1);
            filter = sinf(arg) * sinf(arg);
        }    
        else if( t2<=t   &&  t<=t3 ) {
            filter = 1.0f;
        }   
        else if( t3< t   &&  t<=t4 ) {
            float arg = 0.5*M_PI * (t4-t)/(t4-t3);
            filter = sinf(arg) * sinf(arg);
        }   

        // source[it] *= filter;
    }
        

return n_contribution;
}
float* timeShiftDA(float *trace_inp, int nt, float dt, float ot, float xs, float ys, float zs, \
                   float xr, float yr, float zr, float offMin, float offMax, \
                   float srcWin_ti, float srcWin_ti_taper, float srcWin_tf, float srcWin_tf_taper, float velWater, int applyCosineFactor)
{
    // Sinc interpolation parameters
    int r = 4;
    int r_ = r - 1;
    int diam = 2*r;

    int itrc, it;

    float t1 = srcWin_ti;
    float t2 = srcWin_ti + srcWin_ti_taper;
    float t3 = srcWin_tf - srcWin_tf_taper;
    float t4 = srcWin_tf;

    long long int n_contribution = 0;

    float delta_x = xs - xr;
    float delta_y = ys - yr;
    float delta_z = zs - zr;

    float offset = sqrtf(delta_x*delta_x + delta_y*delta_y);

    if( !( offMin<=offset  &&  offset<=offMax ) )    return NULL;

    float R  = sqrtf(delta_x*delta_x + delta_y*delta_y + delta_z*delta_z);
    float T = R/velWater;

    float *trace_out = CPU_zaloc1F(nt);

    for(it=0; it<nt; it++)
    {
        float t = ot + it*dt;
        float t_shifted = t - T;
        if( t_shifted<t1  ||  t_shifted>t4 )    continue;
        float it0 = floorf((t_shifted-ot)/dt);
        float alpha = ( t_shifted - (it0*dt+ot) ) / dt;

        float sample = trace_inp[it];

        float *coef = interSinc1D_8P(alpha);
        int i;
        for(i=0; i<diam; i++)
        {
            int it_ = it0 - r_ + i;
            if(it_<0  ||  it_>=nt)    continue;
            float geometricalFactor;
            if(applyCosineFactor)    geometricalFactor = (1.0e-5*R*R) * (R/fabsf(delta_z));
            else                     geometricalFactor = (1.0e-5*R*R);
            trace_out[it_] += geometricalFactor * coef[i] * sample;
        }
        free(coef);
    }
    
    /*
    for(it=0; it<nt; it++)
    {       
        float t = ot + it*dt;
        float filter;

        if     ( t < t1  ||  t> t4 ) {
            filter = 0.0f;
        }   
        else if( t1<=t   &&  t< t2 ) {
            float arg = 0.5*M_PI * (t-t1)/(t2-t1);
            filter = sinf(arg) * sinf(arg);
        }    
        else if( t2<=t   &&  t<=t3 ) {
            filter = 1.0f;
        }   
        else if( t3< t   &&  t<=t4 ) {
            float arg = 0.5*M_PI * (t4-t)/(t4-t3);
            filter = sinf(arg) * sinf(arg);
        }   
        trace_out[it] *= filter;
    }
    //*/
        

return trace_out;
}
//*/

void ZeroUpOrDown(Shot3D *shot, int option)
{
    int   nt = shot->nt;
    float dt = shot->dt;
    float ot = shot->ot;
    int   irec, it;
    int   nrec  = shot->nrec;

    for(irec=0; irec<nrec; irec++)
    {
        if(option==1  &&  shot->zr[irec]>0.0)
        {
            for(it=0; it<nt; it++)
                (shot->seismogram)[it + irec*nt] = 0.0;
        }
        if(option==2  &&  shot->zr[irec]<0.0)
        {
            for(it=0; it<nt; it++)
                (shot->seismogram)[it + irec*nt] = 0.0;
        }
    }
        
return;
}
void UpOrDown_Coefficients(Shot3D *shot, float coeff_up, float coeff_down)
{
    int   nt = shot->nt;
    float dt = shot->dt;
    float ot = shot->ot;
    int   irec, it;
    int   nrec  = shot->nrec;

    for(irec=0; irec<nrec; irec++)
    {
        if(shot->zr[irec]>=0.0)
        {
            for(it=0; it<nt; it++)
                (shot->seismogram)[it + irec*nt] *= coeff_up;
        }
        if(shot->zr[irec]<0.0)
        {
            if(coeff_down==0.0) {
                shot->zr [irec] = 0.0;
                shot->izr[irec] = 0;
            }
            for(it=0; it<nt; it++)
                (shot->seismogram)[it + irec*nt]  *= coeff_down;
        }
    }
        
return;
}

void pickReceiver(Shot3D *shot, int irec_pick)
{
    int   nt = shot->nt;
    float dt = shot->dt;
    float ot = shot->ot;
    int   irec, it;
    int   nrec  = shot->nrec;

    for(irec=0; irec<nrec; irec++)
    {
        if(irec != irec_pick)
        {
            for(it=0; it<nt; it++)
                (shot->seismogram)[it + irec*nt] = 0.0;
        }
    }
        
return;
}

void finalizeDataset3DfromSEGY(Dataset3D *dataset, int shot_step)
{
    int ishot;
    for(ishot=0; ishot<dataset->nshot; ishot++)
    {

        Shot3D *shot = &(dataset->shot[ishot]);

        free(shot->xs);
        free(shot->ys);
        free(shot->zs);
        free(shot->xr);
        free(shot->yr);
        free(shot->zr);
    }

return;
}
int initShot3DfromFile(Dataset3D *dataset, int ishot)
{    
    int nshot = dataset->nshot;
    if(ishot>=nshot) {
        printf("\n  [initShot3DfromFile_metaData] The shot requested is not available (ishot=%d, nshot=%d). \n", ishot, dataset->nshot);
        return -1;
    }
        
    int nsou, nrec, nt;
    float ot, dt, maxFreq;
    long long int nStationsPreceding;
    readOdisseiaMetaFile_get_shotInfo(dataset->fMet1, ishot, &nsou, &nrec, &nt, &dt, &ot, &maxFreq, &nStationsPreceding);

    Shot3D *shot = &(dataset->shot[ishot]);

    shot->ishot        = ishot;
    shot->nshot        = nshot;
    shot->nsou         = nsou;
    shot->nrec         = nrec;
    shot->nt           = nt;
    shot->dt           = dt;
    shot->ot           = ot;
    shot->maxFrequency = maxFreq;

    shot->xs  = CPU_zaloc1F(nsou);
    shot->ys  = CPU_zaloc1F(nsou);
    shot->zs  = CPU_zaloc1F(nsou);
    shot->xr  = CPU_zaloc1F(nrec);
    shot->yr  = CPU_zaloc1F(nrec);
    shot->zr  = CPU_zaloc1F(nrec);

    shot->ixs  = CPU_zaloc1I(nsou);
    shot->iys  = CPU_zaloc1I(nsou);
    shot->izs  = CPU_zaloc1I(nsou);
    shot->ixr  = CPU_zaloc1I(nrec);
    shot->iyr  = CPU_zaloc1I(nrec);
    shot->izr  = CPU_zaloc1I(nrec);

    shot->alpha_xs  = CPU_zaloc1F(nsou);
    shot->alpha_ys  = CPU_zaloc1F(nsou);
    shot->alpha_zs  = CPU_zaloc1F(nsou);
    shot->coeff_sou = CPU_zaloc1F(24*nsou);
    shot->alpha_xr  = CPU_zaloc1F(nrec);
    shot->alpha_yr  = CPU_zaloc1F(nrec);
    shot->alpha_zr  = CPU_zaloc1F(nrec);
    shot->coeff_rec = CPU_zaloc1F(24*nrec);

    shot->nStationsPreceding = nStationsPreceding;

    readOdisseiaMetaFile_get_shotCoords(dataset->fMet1, dataset->fMet2, ishot, \
                                        nsou, nrec, nStationsPreceding, \
                                        shot->xs, shot->ys, shot->zs, \
                                        shot->xr, shot->yr, shot->zr);

return 0;
}
void initSeismogramAndSourceMemory(Dataset3D *dataset, int ishot)
{
    Shot3D *shot = &(dataset->shot[ishot]);
    
    int nsou = shot->nsou;
    int nrec = shot->nrec;
    int nt   = shot->nt;

    shot->source     = CPU_zaloc1F(nsou*nt);
    shot->seismogram = CPU_zaloc1F(nrec*nt);

return;
}
void initSeismogramMemory(Dataset3D *dataset, int ishot)
{
    Shot3D *shot = &(dataset->shot[ishot]);
    
    int nrec = shot->nrec;
    int nt   = shot->nt;
    shot->seismogram = CPU_zaloc1F(nrec*nt);

return;
}
void initSourceMemory(Dataset3D *dataset, int ishot)
{
    Shot3D *shot = &(dataset->shot[ishot]);
    
    int nsou = shot->nsou;
    int nt   = shot->nt;

    shot->source     = CPU_zaloc1F(nsou*nt);
return;
}
void loadSeismogram2Dataset(Dataset3D *dataset, int ishot, int normalizeInputDataAmplitudes)
{
    Shot3D *shot = &(dataset->shot[ishot]);
    
    int nsou = shot->nsou;
    int nrec = shot->nrec;
    int nt   = shot->nt;

    shot->source     = CPU_zaloc1F(nsou*nt);
    shot->seismogram = CPU_zaloc1F(nrec*nt);

    readSeismogramFromFile(dataset->fData, shot->source, shot->seismogram, ishot, nsou, nrec, nt, shot->nStationsPreceding);

    if(normalizeInputDataAmplitudes) 
    {
        double L2norm;
        long long int nelems;

        nelems = shot->nt * shot->nrec;
        L2norm = getArrayL2Norm(shot->seismogram, nelems);
        multiplyArrayByScalar( (1.0/L2norm), shot->seismogram, nelems);

        nelems = shot->nt * shot->nsou;
        L2norm = getArrayL2Norm(shot->source, nelems);
        multiplyArrayByScalar( (1.0/L2norm), shot->source, nelems);
    }

return;
}
int finalizeShot3DfromFile(Dataset3D *dataset, int ishot) 
{    
    if(ishot>=dataset->nshot)    return -1;
    Shot3D *shot = &(dataset->shot[ishot]);
    finalizeShot3DMetaInfo(shot);
return 0;
}
int finalizeSeismogramAndSource(Dataset3D *dataset, int ishot)
{    
    if(ishot>=dataset->nshot)    return -1;
    Shot3D *shot = &(dataset->shot[ishot]);
    finalizeShotSeismogramAndSource(shot);
return 0;
}
int finalizeShot3DMetaInfo(Shot3D *shot)
{    
    free(shot->xs);
    free(shot->ys);
    free(shot->zs);
    free(shot->xr);
    free(shot->yr);
    free(shot->zr);

    free(shot->ixs);
    free(shot->iys);
    free(shot->izs);
    free(shot->ixr);
    free(shot->iyr);
    free(shot->izr);

    free(shot->alpha_xs);
    free(shot->alpha_ys);
    free(shot->alpha_zs);
    free(shot->coeff_sou);
    free(shot->alpha_xr);
    free(shot->alpha_yr);
    free(shot->alpha_zr);
    free(shot->coeff_rec);

return 0;
}
int finalizeShotSeismogramAndSource(Shot3D *shot)
{    
    free(shot->seismogram);
    free(shot->source);

return 0;
}

double findAzimuth(float xs, float ys, float xr, float yr) {
    double deltaX = xs - xr;
    double deltaY = ys - yr;
    return atan2(deltaY, deltaX);
}

int getShotWithRestrictedAzimuths(Shot3D *shot_out, Shot3D *shot_inp, double maxAzimDegrees)
{            
    int    nsou = shot_inp->nsou;
    int    nrec = shot_inp->nrec;
    int    nt   = shot_inp->nt;
    float  ot   = shot_inp->ot;
    float  dt   = shot_inp->dt;

    // Find how many receivers are within the permitted azimuth range
    int irec, nrec_restricted=0;
    double xs  = shot_inp->xs[0];
    double ys  = shot_inp->ys[0];
    double zs  = shot_inp->zs[0];
    for(irec=0; irec<nrec; irec++) {
        double xr  = shot_inp->xr[irec];
        double yr  = shot_inp->yr[irec];
        double azm = (180.0/M_PI) * findAzimuth(xs, ys, xr, yr);
        if(fabs(azm) <= maxAzimDegrees)    nrec_restricted++;
    }
        
    // Initialize  arrays of the shot with restricted azimuths
    shot_out->source     = CPU_zaloc1F(nsou*nt);
    shot_out->seismogram = CPU_zaloc1F(nrec_restricted*nt);

    shot_out->ishot        = shot_inp->ishot;
    shot_out->nshot        = shot_inp->nshot;
    shot_out->nsou         = nsou;
    shot_out->nrec         = nrec_restricted;
    shot_out->nt           = nt;
    shot_out->dt           = dt;
    shot_out->ot           = ot;
    shot_out->maxFrequency = shot_inp->maxFrequency;

    shot_out->xs  = CPU_zaloc1F(nsou);
    shot_out->ys  = CPU_zaloc1F(nsou);
    shot_out->zs  = CPU_zaloc1F(nsou);
    shot_out->xr  = CPU_zaloc1F(nrec_restricted);
    shot_out->yr  = CPU_zaloc1F(nrec_restricted);
    shot_out->zr  = CPU_zaloc1F(nrec_restricted);

    shot_out->ixs  = CPU_zaloc1I(nsou);
    shot_out->iys  = CPU_zaloc1I(nsou);
    shot_out->izs  = CPU_zaloc1I(nsou);
    shot_out->ixr  = CPU_zaloc1I(nrec_restricted);
    shot_out->iyr  = CPU_zaloc1I(nrec_restricted);
    shot_out->izr  = CPU_zaloc1I(nrec_restricted);

    shot_out->alpha_xs  = CPU_zaloc1F(nsou);
    shot_out->alpha_ys  = CPU_zaloc1F(nsou);
    shot_out->alpha_zs  = CPU_zaloc1F(nsou);
    shot_out->coeff_sou = CPU_zaloc1F(24*nsou);
    shot_out->alpha_xr  = CPU_zaloc1F(nrec_restricted);
    shot_out->alpha_yr  = CPU_zaloc1F(nrec_restricted);
    shot_out->alpha_zr  = CPU_zaloc1F(nrec_restricted);
    shot_out->coeff_rec = CPU_zaloc1F(24*nrec_restricted);

    shot_out->nStationsPreceding = shot_inp->nStationsPreceding;


    // Copying values from old array
    xs  = shot_inp->xs[0];
    ys  = shot_inp->ys[0];
    zs  = shot_inp->zs[0];
    shot_out->xs[0] = shot_inp->xs[0];
    shot_out->ys[0] = shot_inp->ys[0];
    shot_out->zs[0] = shot_inp->zs[0];
    shot_out->ixs[0] = shot_inp->ixs[0];
    shot_out->iys[0] = shot_inp->iys[0];
    shot_out->izs[0] = shot_inp->izs[0];
    shot_out->alpha_xs[0]  = shot_inp->alpha_xs[0];
    shot_out->alpha_ys[0]  = shot_inp->alpha_ys[0];
    shot_out->alpha_zs[0]  = shot_inp->alpha_zs[0];
    shot_out->coeff_sou[0] = shot_inp->coeff_sou[0];
    memcpy(shot_out->source, shot_inp->source, nt*sizeof(float));
    
    int irec_restricted;
    int it;
    for(irec=0; irec<nrec; irec++) 
    {
        double xr  = shot_inp->xr[irec];
        double yr  = shot_inp->yr[irec];
        double azm = (180.0/M_PI) * findAzimuth(xs, ys, xr, yr);
        if(fabs(azm) <= maxAzimDegrees) 
        {
            shot_out->xr[irec_restricted]  = shot_inp->xr[irec];
            shot_out->yr[irec_restricted]  = shot_inp->yr[irec];
            shot_out->zr[irec_restricted]  = shot_inp->zr[irec];
            shot_out->ixr[irec_restricted] = shot_inp->ixr[irec];
            shot_out->iyr[irec_restricted] = shot_inp->iyr[irec];
            shot_out->izr[irec_restricted] = shot_inp->izr[irec];

            shot_out->alpha_xr[irec_restricted]  = shot_inp->alpha_xr[irec];
            shot_out->alpha_yr[irec_restricted]  = shot_inp->alpha_yr[irec];
            shot_out->alpha_zr[irec_restricted]  = shot_inp->alpha_zr[irec];
            shot_out->coeff_rec[irec_restricted] = shot_inp->coeff_rec[irec];

            long long int shift_inp = nt * irec;
            long long int shift_out = nt * irec_restricted;
            float *src = shot_inp->seismogram;
            float *dst = shot_out->seismogram;
            memcpy(dst+shift_out, src+shift_inp, nt*sizeof(float));

            irec_restricted++;
        }
    }

return 0;
}

void initDataset3D(Dataset3D *dataset, MSH *msh, Parms *parms)
{    
    initGeomParms(&(dataset->acqGeom), parms);

    dataset->nshot = (dataset->acqGeom).nshot;
    dataset->shot  = CPU_zaloc1Shot3D((dataset->acqGeom).nshot);

    initShotFilteringParms(dataset, parms);
    
    initShot3D(dataset, msh);


return;
}

void initShotFilteringParms(Dataset3D *dataset, Parms *parms)
{
    int ishot;
    for(ishot=0; ishot<dataset->nshot; ishot++)
    {
        // Parameters for direct arrival attenuation
        (dataset->shot)[ishot].DA_removeDA        = parms->DA_removeDA;
        (dataset->shot)[ishot].DA_waterSurfaceVel = parms->DA_waterSurfaceVel;
        (dataset->shot)[ishot].DA_offMin          = parms->DA_offMin;
        (dataset->shot)[ishot].DA_offMax          = parms->DA_offMax;

        // Parameters for top mute on shot gather
        (dataset->shot)[ishot].TopMute_apply = parms->TopMute_apply;
        (dataset->shot)[ishot].TopMute_t0    = parms->TopMute_t0;
        (dataset->shot)[ishot].TopMute_t1    = parms->TopMute_t1;
        (dataset->shot)[ishot].TopMute_off0  = parms->TopMute_off0;
        (dataset->shot)[ishot].TopMute_off1  = parms->TopMute_off1;
        (dataset->shot)[ishot].TopMute_ramp  = parms->TopMute_ramp;

        // Parameters for bottom mute on shot gather
        (dataset->shot)[ishot].BottomMute_apply = parms->BottomMute_apply;
        (dataset->shot)[ishot].BottomMute_t0    = parms->BottomMute_t0;
        (dataset->shot)[ishot].BottomMute_t1    = parms->BottomMute_t1;
        (dataset->shot)[ishot].BottomMute_off0  = parms->BottomMute_off0;
        (dataset->shot)[ishot].BottomMute_off1  = parms->BottomMute_off1;
        (dataset->shot)[ishot].BottomMute_ramp  = parms->BottomMute_ramp;
    }
return;
}

void initGeomParms(AcquisitionGeom *acqGeom, Parms *parms)
{
    acqGeom->shotOX = parms->shotOX;
    acqGeom->shotOY = parms->shotOY;
    acqGeom->shotOZ = parms->shotOZ;
    acqGeom->shotDX = parms->shotDX;
    acqGeom->shotDY = parms->shotDY;
    acqGeom->shotDZ = parms->shotDZ;
    acqGeom->shotNX = parms->shotNX;
    acqGeom->shotNY = parms->shotNY;
    acqGeom->shotNZ = parms->shotNZ;

    acqGeom->recOX  = parms->recOX;
    acqGeom->recOY  = parms->recOY;
    acqGeom->recOZ  = parms->recOZ;
    acqGeom->recDX  = parms->recDX;
    acqGeom->recDY  = parms->recDY;
    acqGeom->recDZ  = parms->recDZ;
    acqGeom->recNX  = parms->recNX;
    acqGeom->recNY  = parms->recNY;
    acqGeom->recNZ  = parms->recNZ;

    acqGeom->f1  = parms->f1;
    acqGeom->f2  = parms->f2;
    acqGeom->f3  = parms->f3;
    acqGeom->f4  = parms->f4;

    acqGeom->nshot = acqGeom->shotNX * acqGeom->shotNY * acqGeom->shotNZ;
    acqGeom->nrec  = acqGeom->recNX  * acqGeom->recNY  * acqGeom->recNZ;

return;
}

void finalizeDataset3D(Dataset3D *dataset)
{
    Shot3D *shot;
    int ishot, nshot;
    nshot = dataset->nshot;
/*
    for(ishot=0; ishot<nshot; ishot++)
    {
        shot = &((dataset->shot)[ishot]);

        free(shot->source);
        free(shot->seismogram);

        free(shot->ixs);
        free(shot->iys);
        free(shot->izs);
        free(shot->xs);
        free(shot->ys);
        free(shot->zs);
        free(shot->alpha_xs);
        free(shot->alpha_ys);
        free(shot->alpha_zs);
        // free(shot->coeff_xs);
        // free(shot->coeff_ys);
        // free(shot->coeff_zs);
        free(shot->coeff_sou);
        
        free(shot->ixr);
        free(shot->iyr);
        free(shot->izr);
        free(shot->xr);
        free(shot->yr);
        free(shot->zr);
        free(shot->alpha_xr);
        free(shot->alpha_yr);
        free(shot->alpha_zr);
        // free(shot->coeff_xr);
        // free(shot->coeff_yr);
        // free(shot->coeff_zr);
        free(shot->coeff_rec);
    }
*/

    free((dataset->shot));

return;
}

void initShot3D(Dataset3D *dataset, MSH *msh)
{
    Shot3D *shot;
    AcquisitionGeom *acqGeom = &(dataset->acqGeom);

    int ishot, nshot = dataset->nshot;
    long long int nStationsPreceding = 0;
    
    for(ishot=0; ishot<nshot; ishot++)
    {
        shot = &(dataset->shot[ishot]);
        shot->ishot        = ishot;
        shot->nshot        = nshot;
        shot->nsou         = 1;
        shot->nt           = msh->nt;
        shot->dt           = msh->dt;
        shot->ot           = msh->ot;
        shot->maxFrequency = acqGeom->f4;
        shot->nrec         = acqGeom->nrec;
        shot->seismogram   = CPU_zaloc1F(shot->nrec*msh->nt);
        shot->source       = CPU_zaloc1F(shot->nsou*msh->nt);
        
        shot->nStationsPreceding = nStationsPreceding;
        nStationsPreceding += (shot->nsou + shot->nrec);
    }

    initCoordinatesAndWavelet(dataset, msh, acqGeom);
}

Shot3D* initCopyShot3D(Dataset3D *dataset_inp, int ishot)
{
    Shot3D *shot_out = CPU_zaloc1Shot3D(1);
    Shot3D *shot_inp = &(dataset_inp->shot[ishot]);

    int nshot = dataset_inp->nshot;
    int nsou  = shot_inp->nsou;
    int nrec  = shot_inp->nrec;
    int   nt  = shot_inp->nt;
    float dt  = shot_inp->dt;
    float ot  = shot_inp->ot;
    
    shot_out->ishot        = ishot;
    shot_out->nshot        = nshot;
    shot_out->nsou         = 1;
    shot_out->nt           = nt;
    shot_out->dt           = dt;
    shot_out->ot           = ot;
    shot_out->maxFrequency = shot_inp->maxFrequency;
    shot_out->nrec         = nrec;
    shot_out->seismogram   = CPU_zaloc1F(nrec*nt);
    shot_out->source       = CPU_zaloc1F(nsou*nt);
        
    // Initializing arrays
    shot_out->ixs       = CPU_zaloc1I(nsou);
    shot_out->iys       = CPU_zaloc1I(nsou);
    shot_out->izs       = CPU_zaloc1I(nsou);
    shot_out->xs        = CPU_zaloc1F(nsou);
    shot_out->ys        = CPU_zaloc1F(nsou);
    shot_out->zs        = CPU_zaloc1F(nsou);
    shot_out->alpha_xs  = CPU_zaloc1F(nsou);
    shot_out->alpha_ys  = CPU_zaloc1F(nsou);
    shot_out->alpha_zs  = CPU_zaloc1F(nsou);
    shot_out->coeff_sou = CPU_zaloc1F(24*nsou);

    shot_out->ixr       = CPU_zaloc1I(nrec);
    shot_out->iyr       = CPU_zaloc1I(nrec);
    shot_out->izr       = CPU_zaloc1I(nrec);
    shot_out->xr        = CPU_zaloc1F(nrec);
    shot_out->yr        = CPU_zaloc1F(nrec);
    shot_out->zr        = CPU_zaloc1F(nrec);
    shot_out->alpha_xr  = CPU_zaloc1F(nrec);
    shot_out->alpha_yr  = CPU_zaloc1F(nrec);
    shot_out->alpha_zr  = CPU_zaloc1F(nrec);
    shot_out->coeff_rec = CPU_zaloc1F(24*nrec);

    int irec, isou, it, ic;
    for(irec=0; irec<nrec; irec++)
    {
        for(it=0; it<nt; it++)
            shot_out->seismogram[it + irec*nt] = shot_inp->seismogram[it + irec*nt];

        shot_out->xr[irec] = shot_inp->xr[irec];
        shot_out->yr[irec] = shot_inp->yr[irec];
        shot_out->zr[irec] = shot_inp->zr[irec];

        shot_out->ixr[irec] = shot_inp->ixr[irec];
        shot_out->iyr[irec] = shot_inp->iyr[irec];
        shot_out->izr[irec] = shot_inp->izr[irec];

        shot_out->alpha_xr[irec] = shot_inp->alpha_xr[irec];
        shot_out->alpha_yr[irec] = shot_inp->alpha_yr[irec];
        shot_out->alpha_zr[irec] = shot_inp->alpha_zr[irec];

        for(ic=0; ic<24; ic++)
            shot_out->coeff_rec[irec*24 + ic] = shot_inp->coeff_rec[irec*24 + ic];
    }
    for(isou=0; isou<nsou; isou++) 
    {
        for(it=0; it<nt; it++) 
            shot_out->source[it + isou*nt] = shot_inp->source[it + isou*nt];

        shot_out->xs[isou] = shot_inp->xs[isou];
        shot_out->ys[isou] = shot_inp->ys[isou];
        shot_out->zs[isou] = shot_inp->zs[isou];

        shot_out->ixs[isou] = shot_inp->ixs[isou];
        shot_out->iys[isou] = shot_inp->iys[isou];
        shot_out->izs[isou] = shot_inp->izs[isou];

        shot_out->alpha_xs[isou] = shot_inp->alpha_xs[isou];
        shot_out->alpha_ys[isou] = shot_inp->alpha_ys[isou];
        shot_out->alpha_zs[isou] = shot_inp->alpha_zs[isou];

        for(ic=0; ic<24; ic++)
            shot_out->coeff_sou[isou*24 + ic] = shot_inp->coeff_sou[isou*24 + ic];
    }

return shot_out;
}

void initCoordinatesAndWavelet(Dataset3D *dataset, MSH *msh, AcquisitionGeom *acqGeom) {
        
    int   nx  = msh->nx;
    int   ny  = msh->ny;
    int   nz  = msh->nz;
    int   nt  = msh->nt;
    int   nxb = msh->nxb;
    int   nyb = msh->nyb;
    int   nzb = msh->nzb;
    float dx  = msh->dx;
    float dy  = msh->dy;
    float dz  = msh->dz;
    float dt  = msh->dt;
    float ox  = msh->ox;
    float oy  = msh->oy;
    float oz  = msh->oz;
    float ot  = msh->ot;

    int nshot = acqGeom->nshot;
    int nrec  = acqGeom->nrec;

    float t0 = 0.0f;
    float fm = 5.0f;
    
    int Six, Siy, Siz;
    int Rix, Riy, Riz;

    for(Siz=0; Siz<acqGeom->shotNZ; Siz++)
        for(Siy=0; Siy<acqGeom->shotNY; Siy++)
            for(Six=0; Six<acqGeom->shotNX; Six++)
            {
                int ishot = Six + Siy * acqGeom->shotNX + Siz * acqGeom->shotNX*acqGeom->shotNY;
                
                Shot3D *shot = &(dataset->shot[ishot]);

                // Initializing arrays
                shot->ixs       = CPU_zaloc1I(1);
                shot->iys       = CPU_zaloc1I(1);
                shot->izs       = CPU_zaloc1I(1);
                shot->xs        = CPU_zaloc1F(1);
                shot->ys        = CPU_zaloc1F(1);
                shot->zs        = CPU_zaloc1F(1);
                shot->alpha_xs  = CPU_zaloc1F(1);
                shot->alpha_ys  = CPU_zaloc1F(1);
                shot->alpha_zs  = CPU_zaloc1F(1);
                shot->coeff_sou = CPU_zaloc1F(24);

                shot->ixr       = CPU_zaloc1I(nrec);
                shot->iyr       = CPU_zaloc1I(nrec);
                shot->izr       = CPU_zaloc1I(nrec);
                shot->xr        = CPU_zaloc1F(nrec);
                shot->yr        = CPU_zaloc1F(nrec);
                shot->zr        = CPU_zaloc1F(nrec);
                shot->alpha_xr  = CPU_zaloc1F(nrec);
                shot->alpha_yr  = CPU_zaloc1F(nrec);
                shot->alpha_zr  = CPU_zaloc1F(nrec);
                shot->coeff_rec = CPU_zaloc1F(24*nrec);

                // Initializing source coordinates
                shot->xs[0] = Six * acqGeom->shotDX + acqGeom->shotOX;
                shot->ys[0] = Siy * acqGeom->shotDY + acqGeom->shotOY;
                shot->zs[0] = Siz * acqGeom->shotDZ + acqGeom->shotOZ;

                // Initializing receiver coordinates
                for(Riz=0; Riz<acqGeom->recNZ; Riz++)
                    for(Riy=0; Riy<acqGeom->recNY; Riy++)
                        for(Rix=0; Rix<acqGeom->recNX; Rix++)
                        {
                            int irec = Rix + Riy * acqGeom->recNX + Riz * acqGeom->recNX*acqGeom->recNY;
                            shot->xr[irec] = Rix * acqGeom->recDX + acqGeom->recOX + shot->xs[0];
                            shot->yr[irec] = Riy * acqGeom->recDY + acqGeom->recOY + shot->ys[0];
                            shot->zr[irec] = Riz * acqGeom->recDZ + acqGeom->recOZ + shot->zs[0];
                        }

                //shot->source = getRicker(nt, dt, fm, t0);
                shot->source = getBandpassPulse3D(nt, dt, acqGeom->f1, acqGeom->f2, acqGeom->f3, acqGeom->f4, t0);
            }
}



void integrateRecTime3D(Shot3D *shot)
{
    int    ir, nrec, it, nt;
    float  dt, t, t1, taper, taperPerc=0.1;
    nrec = shot->nrec;
    nt = shot->nt;
    dt = shot->dt;
    t1 = (taperPerc*nt) * dt;
    for(ir=0; ir<nrec; ir++)
    {
        for(it=0; it<shot->nt; it++) 
        {
            t = it*dt;
            if(t<t1)    taper = t/t1;
            else        taper = 1.0;
            shot->seismogram[(it+1)+ir*shot->nt] += taper*dt*shot->seismogram[it+ir*shot->nt];
        }

    }

}

void AGC3D(Dataset3D *dataset, float time_window) {

    Shot3D  *shot;
    float  dt;
    int    ishot, nshot, ir, nrec, it, nt;
    int    its, nts, nsamp;

    shot = &((dataset->shot)[0]);

    nshot = shot->nshot;
    nrec  = shot->nrec;
    nt    = shot->nt;
    dt    = shot->dt;

    nts = ((int) ceilf(time_window/dt))/2;

    double *amps  = CPU_zaloc1D(nt);
    float *trace;

    for(ishot=0; ishot<nshot; ishot++)
    {
        shot = &((dataset->shot)[ishot]);
        
        for(ir=0; ir<nrec; ir++)
        {
            trace  = (shot->seismogram)+(ir*nt);
            
            // compute amplitudes            
            memset(amps, 0, nt*sizeof(double));
            
            for(it=0; it<nt; it++)
            {   
                nsamp = 0;
                int its0 = max(0,it-nts);
                int its1 = min(nt-1,it+nts);
                for(its=its0; its<=its1; its++)
                {
                    nsamp++;
                    amps[it] += trace[its]*trace[its];
                }
                amps[it] /= nsamp;
            }
            

            // checks if there is any amplitude too low
            double avg = 0.0;
            for(it=0; it<nt; it++)
            {
                avg += amps[it]/nt;
            }
            if(avg==0.0)
            {
                for(it=0; it<nt; it++)
                    amps[it] = 1.0e-6;
            }
            else
            {
                for(it=0; it<nt; it++) {
                    if(amps[it]<1.0e-6*avg)    amps[it] = 1.0e-6*avg;
                }
            }
            
            // Normalizes data
            for(it=0; it<nt; it++) {
                trace[it] /= amps[it];
            }
        }
        
    }

return;
}

void applyRecTaper3D(Shot3D *shot, float perc, float offMax, int sensorSide)
{
    int   it, ir, nrec;
    float xs, ys, zs, xr, yr, zr;
    float offset, aoffset, taper;
    float off0, off1, maxOff;

    if(shot->nsou>1)    return;

    xs   = shot->xs[0];
    ys   = shot->ys[0];
    zs   = shot->zs[0];
    nrec = shot->nrec;

    maxOff = 0.0;
    for(ir=0; ir<nrec; ir++)
    {
        xr =  shot->xr[ir];
        yr =  shot->yr[ir];
        aoffset = sqrtf( (xs-xr)*(xs-xr) + (ys-yr)*(ys-yr) );
        if(maxOff<aoffset)    maxOff = aoffset;
    }
    if(offMax>maxOff)    offMax = maxOff;

    off1 = offMax;
    off0 = (1.0-perc) * offMax;

    for(ir=0; ir<nrec; ir++)
    {
        xr =  shot->xr[ir];
        yr =  shot->yr[ir];
        float offset_x = (xs-xr);
        aoffset = sqrtf( (xs-xr)*(xs-xr) + (ys-yr)*(ys-yr) );
        
        if(off0>aoffset)
            taper = 1.0f;
        else if(off0<=aoffset && aoffset<off1)
            taper = 1.0 - (aoffset-off0)/(off1-off0);
        else
            taper = 0.0f;
        
        if( (fsignf(offset_x)!=sensorSide)  &&  sensorSide!=0 )
            taper = 0.0;
        
        for(it=0; it<shot->nt; it++)
            shot->seismogram[it+ir*shot->nt] *= taper;
    }
}


void finalizeShot3D(Shot3D *shot) {
    
    free(shot->source);
    free(shot->seismogram);
    
    free(shot->ixr);
    free(shot->iyr);
    free(shot->izr);
    free(shot->xr);
    free(shot->yr);
    free(shot->zr);
    free(shot->ixs);
    free(shot->iys);
    free(shot->izs);
    free(shot->xs);
    free(shot->ys);
    free(shot->zs);

    free(shot->alpha_xs);
    free(shot->alpha_ys);
    free(shot->alpha_zs);
    free(shot->coeff_sou);
    free(shot->alpha_xr);
    free(shot->alpha_yr);
    free(shot->alpha_zr);
    free(shot->coeff_rec);

return;
}

double addObjFunc3D(Shot3D *shot) {
    int it, ir;
    int nt = shot->nt;
    int nrec = shot->nrec;
    double sum = 0.0;
    for(ir=0; ir<nrec; ir++) {
        for(it=0; it<nt; it++) {
            sum += (shot->seismogram)[it+ir*nt] * (shot->seismogram)[it+ir*nt];
        }
    }
return sum;
}
void removeDA3D(Shot3D *shot, Shot3D *shot_da) {
    int it, ir;
    int nt = shot->nt;
    int nrec = shot->nrec;
    for(ir=0; ir<nrec; ir++) {
        for(it=0; it<nt; it++) {
            (shot->seismogram)[it+ir*nt] -= (shot_da->seismogram)[it+ir*nt];
        }
    }
return;
}
void copyData3D(Shot3D *shot_inp, Shot3D *shot_out) {
    int it, ir;
    int nt = shot_inp->nt;
    int nrec = shot_inp->nrec;
    for(ir=0; ir<nrec; ir++) {
        for(it=0; it<nt; it++) {
            (shot_out->seismogram)[it+ir*nt] = (shot_inp->seismogram)[it+ir*nt];
        }
    }
return;
}

void applyHalfDerivative2Source3D(Shot3D *shot) {
    applyHalfDerivative(shot->source, shot->nt, shot->dt, +1.0);
return;
}
void applyHalfDerivative2Data3D(Shot3D *shot) {
    int ir;
    for(ir=0; ir<shot->nrec; ir++) {
        applyHalfDerivative( (shot->seismogram)+ir*shot->nt, shot->nt, shot->dt, -1.0);
    }
return;
}


float *getBandpassPulse3D(int nt, float dt, float f1, float f2, float f3, float f4, float t0) 
{
    int it;
    float *pulse = CPU_zaloc1F(nt);
    t0 += 1.0f/(0.5*(f2+f1));
    int it0 = (int) (0.5f + t0/dt);
    pulse[it0] = 1.0f;
    bandpass(pulse, nt, dt, f1, f2, f3, f4);
    for(it=0; it<nt; it++) 
    {
        float t = it*dt;
        if(t<t0/3)
        {
            pulse[it]       *= t/(t0/3);
            pulse[2*it0-it] *= t/(t0/3);
        }
        else if(it>2*it0)    pulse[it] = 0.0f;
    }
return pulse;
}

void applyBandpassToTrace3D(Shot3D *shot, float f1, float f2, float f3, float f4) {
    int ir;
    for(ir=0; ir<shot->nrec; ir++)
        bandpass((shot->seismogram)+(ir*shot->nt), shot->nt, shot->dt, f1, f2, f3, f4);
return;
}
void applyBandpassToSourceAndSeismogram(Shot3D *shot, float f1, float f2, float f3, float f4) {
    int ir, is;
    for(is=0; is<shot->nsou; is++)    bandpass((shot->source    )+(is*shot->nt), shot->nt, shot->dt, f1, f2, f3, f4);
    for(ir=0; ir<shot->nrec; ir++)    bandpass((shot->seismogram)+(ir*shot->nt), shot->nt, shot->dt, f1, f2, f3, f4);
return;
}
void multiplySeismogramByTime3D(Shot3D *shot, float power) {
    int ir;
    for(ir=0; ir<shot->nrec; ir++)
        multiplyTraceByTime3D((shot->seismogram)+(ir*shot->nt), shot->nt, shot->dt, power);
return;
}
void multiplyTraceByTime3D(float *trace, int nt, float dt, float power) {
    int it;
    for(it=0; it<nt; it++)
        trace[it] *= pow(it*dt, power);
}

void demultiplexSeismogram3D(float *seismogram, int nt, int nrec) {
    float *seisTemp = CPU_zaloc1F(nt*nrec);
    int ir, it, i;
    for(ir=0; ir<nrec; ir++)
        for(it=0; it<nt; it++) {
            seisTemp[it+ir*nt] = seismogram[ir+it*nrec];
    }
    for(i=0; i<nrec*nt; i++)    seismogram[i] = seisTemp[i];

    free(seisTemp);
}


void makeEncodedDataset3D(int idev, Dataset3D *data_inp, Dataset3D *data_tmp, Dataset3D *data_enc, Dataset3D *data_res, int codelength, int codestep, int realizations)
{  
    Shot3D *shot_inp, *shot_tmp, *shot_enc, *shot_res;
    int   ishot, irec, it, nshot, nsou, nrec, nt, instance, seed;
    float dt;

    data_enc->nshot = realizations;
    data_res->nshot = realizations;
    data_enc->shot  = CPU_zaloc1Shot3D(realizations);
    data_res->shot  = CPU_zaloc1Shot3D(realizations);

    nt = shotEncoder3D(idev, data_inp, NULL, codelength, codestep, 0);
    dt = (data_inp->shot)[0].dt;

    // Total number of shots in the input data
    nshot = data_inp->nshot;
    nsou = nshot;
    
    // Number of receivers in the encoded data (sum of all receivers from all shots in the input data)
    nrec = 0;
    for(ishot=0; ishot<nshot; ishot++)
    {
        shot_tmp = &((data_tmp->shot)[ishot]);
        nrec += shot_tmp->nrec;
    }

    for(instance=0; instance<realizations; instance++)
    {
        printf("\n Encoding: Instance=%d \n", instance);

        seed = 10*(instance+1);
        shotEncoder3D(idev, data_inp, data_tmp, codelength, codestep, seed);
        
        shot_enc = &((data_enc->shot)[instance]);
        shot_res = &((data_res->shot)[instance]);

        shot_enc->nsou = nshot;
        shot_enc->nrec = nrec;
        shot_enc->nt   = nt;
        shot_enc->dt   = dt;
        shot_res->nsou = nshot;
        shot_res->nrec = nrec;
        shot_res->nt   = nt;
        shot_res->dt   = dt;

        shot_enc->ixs = CPU_zaloc1I(nshot);
        shot_enc->iys = CPU_zaloc1I(nshot);
        shot_enc->izs = CPU_zaloc1I(nshot);
        shot_enc->xs  = CPU_zaloc1F(nshot);
        shot_enc->ys  = CPU_zaloc1F(nshot);
        shot_enc->zs  = CPU_zaloc1F(nshot);

        shot_res->ixs = CPU_zaloc1I(nshot);
        shot_res->iys = CPU_zaloc1I(nshot);
        shot_res->izs = CPU_zaloc1I(nshot);
        shot_res->xs  = CPU_zaloc1F(nshot);
        shot_res->ys  = CPU_zaloc1F(nshot);
        shot_res->zs  = CPU_zaloc1F(nshot);

        shot_enc->ixr = CPU_zaloc1I(nrec);
        shot_enc->iyr = CPU_zaloc1I(nrec);
        shot_enc->izr = CPU_zaloc1I(nrec);
        shot_enc->xr  = CPU_zaloc1F(nrec);
        shot_enc->yr  = CPU_zaloc1F(nrec);
        shot_enc->zr  = CPU_zaloc1F(nrec);

        shot_res->ixr = CPU_zaloc1I(nrec);
        shot_res->iyr = CPU_zaloc1I(nrec);
        shot_res->izr = CPU_zaloc1I(nrec);
        shot_res->xr  = CPU_zaloc1F(nrec);
        shot_res->yr  = CPU_zaloc1F(nrec);
        shot_res->zr  = CPU_zaloc1F(nrec);

        shot_enc->source     = CPU_zaloc1F(nshot*nt);
        shot_res->source     = CPU_zaloc1F(nshot*nt);
        shot_enc->seismogram = CPU_zaloc1F(nrec*nt);
        shot_res->seismogram = CPU_zaloc1F(nrec*nt);

        // Closing and opening source and seismogram arrays with new time length
        int nr_ = 0; // total number of receivers in the precedent shots
        for(ishot=0; ishot<nshot; ishot++)
        {
            shot_tmp = &((data_tmp->shot)[ishot]);

            int shift = ishot*nt;
            size_t bytes = nt*sizeof(float);
            memcpy((shot_enc->source)+shift, (shot_tmp->source), bytes);
            memcpy((shot_res->source)+shift, (shot_tmp->source), bytes);

            shift = nr_ * nt;
            bytes = shot_tmp->nrec * nt * sizeof(float);
            memcpy((shot_enc->seismogram)+shift, (shot_tmp->seismogram), bytes);
            memcpy((shot_res->seismogram)+shift, (shot_tmp->seismogram), bytes);

            shift = nr_;
            bytes = shot_tmp->nrec * sizeof(int);
            memcpy((shot_enc->ixr)+shift, (shot_tmp->ixr), bytes);
            memcpy((shot_enc->iyr)+shift, (shot_tmp->iyr), bytes);
            memcpy((shot_enc->izr)+shift, (shot_tmp->izr), bytes);
            memcpy((shot_res->ixr)+shift, (shot_tmp->ixr), bytes);
            memcpy((shot_res->iyr)+shift, (shot_tmp->iyr), bytes);
            memcpy((shot_res->izr)+shift, (shot_tmp->izr), bytes);

            shift = nr_;
            bytes = shot_tmp->nrec * sizeof(float);
            memcpy((shot_enc->xr)+shift, (shot_tmp->xr), bytes);
            memcpy((shot_enc->yr)+shift, (shot_tmp->yr), bytes);
            memcpy((shot_enc->zr)+shift, (shot_tmp->zr), bytes);
            memcpy((shot_res->xr)+shift, (shot_tmp->xr), bytes);
            memcpy((shot_res->yr)+shift, (shot_tmp->yr), bytes);
            memcpy((shot_res->zr)+shift, (shot_tmp->zr), bytes);

            shot_enc->ixs[ishot] = shot_tmp->ixs[0];
            shot_enc->iys[ishot] = shot_tmp->iys[0];
            shot_enc->izs[ishot] = shot_tmp->izs[0];
            shot_enc->xs[ishot]  = shot_tmp->xs[0];
            shot_enc->ys[ishot]  = shot_tmp->ys[0];
            shot_enc->zs[ishot]  = shot_tmp->zs[0];

            shot_res->ixs[ishot] = shot_tmp->ixs[0];
            shot_res->iys[ishot] = shot_tmp->iys[0];
            shot_res->izs[ishot] = shot_tmp->izs[0];
            shot_res->xs[ishot]  = shot_tmp->xs[0];
            shot_res->ys[ishot]  = shot_tmp->ys[0];
            shot_res->zs[ishot]  = shot_tmp->zs[0];

            nr_ += shot_tmp->nrec;
        }
    }

return;
}


void windowData3D(int idev, Dataset3D *dataset, float tmin, float tmax, float dt)
{
    int ishot, it, it0, isou, irec;
    int nshot = dataset->nshot;
    for(ishot=0; ishot<nshot; ishot++)
    {
        // Resamp source
        Shot3D *shot = &(dataset->shot[ishot]);
        int   nsou = shot->nsou;
        int   nrec = shot->nrec;
        int   nt0 = shot->nt;
        float dt0 = shot->dt;

        float *source = CPU_zaloc1F(nt0*nsou);
        float *seismo = CPU_zaloc1F(nt0*nrec);
        memcpy(source, shot->source    , nt0*nsou*sizeof(float));
        memcpy(seismo, shot->seismogram, nt0*nrec*sizeof(float));

        int nt  = floorf(nt0*dt0/dt);

        shot->nt = nt;
        shot->dt = dt;
        free(shot->source);
        free(shot->seismogram);
        shot->source     = CPU_zaloc1F(nt*nsou);
        shot->seismogram = CPU_zaloc1F(nt*nrec);

        if(dt==dt0)
        {
            for(isou=0; isou<nsou; isou++)
            {
                for(it=0; it<nt; it++)
                    shot->source[isou*nt+it] = source[isou*nt+it];
            }
            for(irec=0; irec<nrec; irec++)
            {
                for(it=0; it<nt; it++)
                    shot->seismogram[irec*nt+it] = seismo[irec*nt+it];
            }
        }

        if(idev<0)
        {
            for(isou=0; isou<nsou; isou++)    resamp(source+isou*nt0, (shot->source    )+isou*nt, nt0, dt0, nt, dt);
            for(irec=0; irec<nrec; irec++)    resamp(seismo+irec*nt0, (shot->seismogram)+irec*nt, nt0, dt0, nt, dt);
        }
        else
        {
            // GPU_WRAP_SeismogramResamp(idev, source, shot->source    , nsou, nt0, dt0, nt, dt);
            // GPU_WRAP_SeismogramResamp(idev, seismo, shot->seismogram, nrec, nt0, dt0, nt, dt);
        }
        free(source);
        free(seismo);
    }

return;
}


void copyDataset3D(Dataset3D *dst, Dataset3D *src)
{  
    Shot3D *shot_inp, *shot_tmp, *shot_enc, *shot_res;
    int ishot, irec, it, nshot, nsou, nrec, nt, instance, seed;

    dst->nshot = src->nshot;
    nshot      = dst->nshot;

    dst->shot  = CPU_zaloc1Shot3D(nshot);

    int ns_ = 0;
    int nr_ = 0;
    for(ishot=0; ishot<nshot; ishot++)
    {
        
        Shot3D *shot_dst = &((dst->shot)[ishot]);
        Shot3D *shot_src = &((src->shot)[ishot]);

        shot_dst->ishot = shot_src->ishot;
        shot_dst->nshot = shot_src->nshot;
        shot_dst->nsou  = shot_src->nsou;
        shot_dst->nrec  = shot_src->nrec;
        shot_dst->nt    = shot_src->nt;
        shot_dst->ixMin = shot_src->ixMin;
        shot_dst->ixMax = shot_src->ixMax;
        shot_dst->dt    = shot_src->dt;
        shot_dst->ot    = shot_src->ot;

        nsou = shot_dst->nsou;
        nrec = shot_dst->nrec;
        nt   = shot_dst->nt;

        shot_dst->ixs = CPU_zaloc1I(nsou);
        shot_dst->iys = CPU_zaloc1I(nsou);
        shot_dst->izs = CPU_zaloc1I(nsou);
        shot_dst->xs  = CPU_zaloc1F(nsou);
        shot_dst->ys  = CPU_zaloc1F(nsou);
        shot_dst->zs  = CPU_zaloc1F(nsou);
        shot_dst->ixr = CPU_zaloc1I(nrec);
        shot_dst->iyr = CPU_zaloc1I(nrec);
        shot_dst->izr = CPU_zaloc1I(nrec);
        shot_dst->xr  = CPU_zaloc1F(nrec);
        shot_dst->yr  = CPU_zaloc1F(nrec);
        shot_dst->zr  = CPU_zaloc1F(nrec);

        shot_dst->source     = CPU_zaloc1F(nsou*nt);
        shot_dst->seismogram = CPU_zaloc1F(nrec*nt);

        size_t bytes = nsou * nt * sizeof(float);
        memcpy((shot_dst->source), (shot_src->source), bytes);

        bytes = nrec * nt * sizeof(float);
        memcpy((shot_dst->seismogram), (shot_src->seismogram), bytes);

        bytes = nsou * sizeof(int);
        memcpy((shot_dst->ixs), (shot_src->ixs), bytes);
        memcpy((shot_dst->iys), (shot_src->iys), bytes);
        memcpy((shot_dst->izs), (shot_src->izs), bytes);
        bytes = nsou * sizeof(float);
        memcpy((shot_dst->xs ), (shot_src->xs ), bytes);
        memcpy((shot_dst->ys ), (shot_src->ys ), bytes);
        memcpy((shot_dst->zs ), (shot_src->zs ), bytes);

        bytes = nrec * sizeof(int);
        memcpy((shot_dst->ixr), (shot_src->ixr), bytes);
        memcpy((shot_dst->iyr), (shot_src->iyr), bytes);
        memcpy((shot_dst->izr), (shot_src->izr), bytes);
        bytes = nrec * sizeof(float);
        memcpy((shot_dst->xr ), (shot_src->xr ), bytes);
        memcpy((shot_dst->yr ), (shot_src->yr ), bytes);
        memcpy((shot_dst->zr ), (shot_src->zr ), bytes);        
    }

return;
}

int shotEncoder3D(int idev, Dataset3D *data_in, Dataset3D *data_ou, int codelength, int codestep, int seed)
{  
    float *out;
    Shot3D *shot;
    int ishot, irec, it, nshot, nrec, nt;
    nshot = data_in->nshot;
    nt    = (data_in->shot)[0].nt;

    int ncdl = (codelength-1) * codestep;
    int  nt_ = nt + ncdl - 1;
    float *codes = codeGen(ncdl, codestep, nshot, seed);

    if(data_ou==NULL)
    {
        printf("\n\n outputing nt_=%d   codelength=%d  codestep=%d   ncdl=%d    ntOrig=%d  \n\n", nt_, codelength, codestep, ncdl, nt);
        return nt_;
    }

    outputSmart2d("codes", codes, ncdl, (data_in->shot)[0].dt, 0, nshot, 1, 0);

    // Closing and opening source and seismogram arrays with new time length
    for(ishot=0; ishot<nshot; ishot++)
    {
        nrec = (data_ou->shot)[ishot].nrec;
        shot = &((data_ou->shot)[ishot]);
        shot->nt = nt_;
        free(shot->source);
        free(shot->seismogram);
        shot->source = CPU_zaloc1F(nt_);
        shot->seismogram = CPU_zaloc1F(nt_*nrec);
    }

    // Producing the encoded source and seismogram
    for(ishot=0; ishot<nshot; ishot++)
    {
        float *codeForThisShot = &(codes[ishot*ncdl]);
        shot = &((data_ou->shot)[ishot]);
        nrec = (data_ou->shot)[ishot].nrec;

        float *trc = (data_in->shot)[ishot].source;
        float *out = convolution(trc, codeForThisShot, nt, ncdl);
        // float *out = GPU_WRAP_convolution(idev, trc, codeForThisShot, nt, ncdl);
        for(it=0; it<nt_; it++)        ((data_ou->shot)[ishot].source)[it] = out[it];
        free(out);

        if(idev<0)
        {
            for(irec=0; irec<nrec; irec++)
            {
                float *trc = ((data_in->shot)[ishot].seismogram)+irec*nt;
                float *out = convolution(trc, codeForThisShot, nt, ncdl);
                // float *out = GPU_WRAP_convolution(idev, trc, codeForThisShot, nt, ncdl);
                for(it=0; it<nt_; it++)    ((data_ou->shot)[ishot].seismogram)[it+irec*nt_] = out[it];
                free(out);
            }
        }
        else
        {
            trc = (data_in->shot)[ishot].seismogram;
            // out = GPU_WRAP_SeismogramCoding(idev, trc, codeForThisShot, nrec, nt, ncdl);
            memcpy((data_ou->shot)[ishot].seismogram, out, nrec*nt_*sizeof(float));
            free(out);
        }
    }
    free(codes);
return nt_;
}

void correlateWavelets3D(Dataset3D *data)
{
    float *out;
    Shot3D *shot;
    Shot3D *shot_temp;
    int ishot, ishot_, it, nshot, nrec, nt;
    float dt;
    nshot = data->nshot;
    nt    = (data->shot)[0].nt;
    dt    = (data->shot)[0].dt;

    float *stack = CPU_zaloc1F(2*nt-1);

    // Closing and opening source and seismogram arrays with new time length
    FILE *fSou_Hdr, *fSou_Bin;
    initFilesSmart3d("xcorr_wavelets", &fSou_Hdr, &fSou_Bin, 2*nt-1, dt, 0,   nshot, 1, 0,   nshot, 1, 0);
    printf("\n correlating wavelets 001 \n");
    for(ishot=0; ishot<nshot; ishot++)
    {
        printf("\n *** ishot=%d \n  ishot_ = ", ishot);
        shot = &((data->shot)[ishot]);
        for(ishot_=0; ishot_<nshot; ishot_++)
        {
            printf(" %d ", ishot_);
            shot_temp = &((data->shot)[ishot_]);
            float *corr = correlation(shot->source, shot_temp->source, nt, nt);
            for(it=0;it<2*nt-1; it++)    stack[it] += corr[it];
            writeSamples(corr, fSou_Bin, 2*nt-1);
            free(corr);
        }
    }
    fclose(fSou_Bin);
    fclose(fSou_Hdr);
    outputSmart1d("stack", stack, 2*nt-1, dt, 0);
    
return;
}

void writeSeisData2File3D(Dataset3D *dataset, int jshot, int iproc, char *fileNameSeismos, char *fileNameSources)
{
    int   ishot, nshot=dataset->nshot;
    int   nsou, nrec, nt;
    float dt;

    nrec = (dataset->shot)[0].nrec;
    nsou = (dataset->shot)[0].nsou;
    nt   = (dataset->shot)[0].nt;
    dt   = (dataset->shot)[0].dt;

    // Save seismogram with no DA to file
    FILE *fSeis_Hdr, *fSeis_Bin;
    FILE *fSou_Hdr, *fSou_Bin;
    if(iproc==0)
    {
        initFilesSmart3d(fileNameSeismos, &fSeis_Hdr, &fSeis_Bin, nt, dt, 0, nrec, 1, 0, dataset->nshot/jshot, 1, 0);

        if(fileNameSources!=NULL) initFilesSmart3d(fileNameSources, &fSou_Hdr, &fSou_Bin, nt, dt, 0, nsou, 1, 0, dataset->nshot/jshot, 1, 0);

        for(ishot=0; ishot<nshot; ishot+=jshot)
        {
            Shot3D *shot = &((dataset->shot)[ishot]);
            writeSamples(shot->seismogram, fSeis_Bin, nt*shot->nrec);
            if(fileNameSources!=NULL)
                writeSamples(shot->source, fSou_Bin, nsou*nt);
        }
        fclose(fSeis_Bin);
        fclose(fSeis_Hdr);
        if(fileNameSources!=NULL)
        {
            fclose(fSou_Bin);
            fclose(fSou_Hdr);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

return;
}

int dumpShotCoordinates3D(MSH *msh, Shot3D *shot)
{
    int is, ir;
    int ixs, iys, izs;
    int ixr, iyr, izr;

    printf("\n\n !!! Dumping Shot Coordinates !!! \n\n");

    FILE *fout = NULL;
    fout = fopen("shotCoords_debug", "w");
    if(fout==NULL) printf("\n\n CANNOT OPEN shotCoords_debug \n\n");
    fprintf(fout, "\n Source coordinates \n");
    
    for(is=0; is<shot->nsou; is++)
    {
        ixs = (int) ( (shot->xs[is] - msh->ox) / msh->dx );
        iys = (int) ( (shot->ys[is] - msh->oy) / msh->dy );
        izs = (int) ( (shot->zs[is] - msh->oz) / msh->dz );

        fprintf(fout, "xs=%f  ys=%f  zs=%f        ixs=%d  iys=%d  izs=%d \n", shot->xs[is], shot->ys[is], shot->zs[is], ixs, iys, izs);
    }

    fprintf(fout, "\n Receiver coordinates \n");
    for(ir=0; ir<shot->nrec; ir++)
    {
        ixr = (int) ( (shot->xr[ir] - msh->ox) / msh->dx );
        iyr = (int) ( (shot->yr[ir] - msh->oy) / msh->dy );
        izr = (int) ( (shot->zr[ir] - msh->oz) / msh->dz );

        fprintf(fout, "xr=%f  yr=%f  zr=%f        ixr=%d  iyr=%d  izr=%d \n", shot->xr[ir], shot->yr[ir], shot->zr[ir], ixr, iyr, izr);
    }

    fclose(fout);

return 0;
}

int checkShotBoundaries3D(MSH *msh, Shot3D *shot)
{
    int is, ir;
    int ixs, iys, izs;
    int ixr, iyr, izr;
    for(is=0; is<shot->nsou; is++)
    {
        ixs = (int) ( (shot->xs[is] - msh->ox) / msh->dx );
        iys = (int) ( (shot->ys[is] - msh->oy) / msh->dy );
        izs = (int) ( (shot->zs[is] - msh->oz) / msh->dz );

        if( !( 0.0<=ixs  &&  ixs<msh->nx )  )    { printf("\n is=%d xs[is]=%f ixs[is]=%d nx=%d dx=%f ox=%f\n", is, shot->xs[is], shot->ixs[is], msh->nx, msh->dx, msh->ox); return 1; }
        if( !( 0.0<=iys  &&  iys<msh->ny )  )    { printf("\n is=%d ys[is]=%f iys[is]=%d ny=%d dy=%f oy=%f\n", is, shot->ys[is], shot->iys[is], msh->ny, msh->dy, msh->oy); return 1; }
        if( !( 0.0<=izs  &&  izs<msh->nz )  )    { printf("\n is=%d zs[is]=%f izs[is]=%d nz=%d dz=%f oz=%f\n", is, shot->zs[is], shot->izs[is], msh->nz, msh->dz, msh->oz); return 1; }
    }
    for(ir=0; ir<shot->nrec; ir++)
    {
        ixr = (int) ( (shot->xr[ir] - msh->ox) / msh->dx );
        iyr = (int) ( (shot->yr[ir] - msh->oy) / msh->dy );
        izr = (int) ( (shot->zr[ir] - msh->oz) / msh->dz );

        if( !( 0.0<=ixr  &&  ixr<msh->nx )  )    { printf("\n ir=%d xr[ir]=%f ixr[ir]=%d nx=%d dx=%f ox=%f\n", ir, shot->xr[ir], shot->ixr[ir], msh->nx, msh->dx, msh->ox); return 1; }
        if( !( 0.0<=iyr  &&  iyr<msh->ny )  )    { printf("\n ir=%d yr[ir]=%f iyr[ir]=%d ny=%d dy=%f oy=%f\n", ir, shot->yr[ir], shot->iyr[ir], msh->ny, msh->dy, msh->oy); return 1; }
        if( !( 0.0<=izr  &&  izr<msh->nz )  )    { printf("\n ir=%d zr[ir]=%f izr[ir]=%d nz=%d dz=%f oz=%f\n", ir, shot->zr[ir], shot->izr[ir], msh->nz, msh->dz, msh->oz); return 1; }
    }

return 0;
}

int gridShot3D(MSH *msh, Shot3D *shot)
{
    int is, ir;
    for(is=0; is<shot->nsou; is++)
    {
        shot->ixs[is] = (int) ( (shot->xs[is] - msh->ox) / msh->dx );
        shot->iys[is] = (int) ( (shot->ys[is] - msh->oy) / msh->dy );
        shot->izs[is] = (int) ( (shot->zs[is] - msh->oz) / msh->dz );

        shot->alpha_xs[is] = (shot->xs[is] - (shot->ixs[is] * msh->dx + msh->ox)) / msh->dx;
        shot->alpha_ys[is] = (shot->ys[is] - (shot->iys[is] * msh->dy + msh->oy)) / msh->dy;
        shot->alpha_zs[is] = (shot->zs[is] - (shot->izs[is] * msh->dz + msh->oz)) / msh->dz;

        interSinc3D_8P(shot->coeff_sou+24*is+8, shot->coeff_sou+24*is+16, shot->coeff_sou+24*is+0, shot->alpha_xs[is], shot->alpha_ys[is], shot->alpha_zs[is]);

        int out_of_grid_limits = 0;
        if( !( 0<=shot->ixs[is]  &&  shot->ixs[is]<msh->nx )  )    { printf("\n ishot=%d  is=%d xs[is]=%f ixs[is]=%d nx=%d dx=%f ox=%f\n", shot->ishot, is, shot->xs[is], shot->ixs[is], msh->nx, msh->dx, msh->ox); out_of_grid_limits = 1; }
        if( !( 0<=shot->iys[is]  &&  shot->iys[is]<msh->ny )  )    { printf("\n ishot=%d  is=%d ys[is]=%f iys[is]=%d ny=%d dy=%f oy=%f\n", shot->ishot, is, shot->ys[is], shot->iys[is], msh->ny, msh->dy, msh->oy); out_of_grid_limits = 1; }
        if( !( 0<=shot->izs[is]  &&  shot->izs[is]<msh->nz )  )    { printf("\n ishot=%d  is=%d zs[is]=%f izs[is]=%d nz=%d dz=%f oz=%f\n", shot->ishot, is, shot->zs[is], shot->izs[is], msh->nz, msh->dz, msh->oz); out_of_grid_limits = 1; }
        if(out_of_grid_limits)    return 1; 
    }
    for(ir=0; ir<shot->nrec; ir++)
    {
        shot->ixr[ir] = (int) ( (shot->xr[ir] - msh->ox) / msh->dx );
        shot->iyr[ir] = (int) ( (shot->yr[ir] - msh->oy) / msh->dy );
        shot->izr[ir] = (int) ( (shot->zr[ir] - msh->oz) / msh->dz );

        shot->alpha_xr[ir] = (shot->xr[ir] - (shot->ixr[ir] * msh->dx + msh->ox)) / msh->dx;
        shot->alpha_yr[ir] = (shot->yr[ir] - (shot->iyr[ir] * msh->dy + msh->oy)) / msh->dy;
        shot->alpha_zr[ir] = (shot->zr[ir] - (shot->izr[ir] * msh->dz + msh->oz)) / msh->dz;

        interSinc3D_8P(shot->coeff_rec+24*ir+8, shot->coeff_rec+24*ir+16, shot->coeff_rec+24*ir+0, shot->alpha_xr[ir], shot->alpha_yr[ir], shot->alpha_zr[ir]);

        int out_of_grid_limits = 0;
        if( !( 0<=shot->ixr[ir]  &&  shot->ixr[ir]<msh->nx )  )    { printf("\n ishot=%d  ir=%d xr[ir]=%f ixr[ir]=%d nx=%d dx=%f ox=%f\n", shot->ishot, ir, shot->xr[ir], shot->ixr[ir], msh->nx, msh->dx, msh->ox); out_of_grid_limits = 1; }
        if( !( 0<=shot->iyr[ir]  &&  shot->iyr[ir]<msh->ny )  )    { printf("\n ishot=%d  ir=%d yr[ir]=%f iyr[ir]=%d ny=%d dy=%f oy=%f\n", shot->ishot, ir, shot->yr[ir], shot->iyr[ir], msh->ny, msh->dy, msh->oy); out_of_grid_limits = 1; }
        if( !( 0<=shot->izr[ir]  &&  shot->izr[ir]<msh->nz )  )    { printf("\n ishot=%d  ir=%d zr[ir]=%f izr[ir]=%d nz=%d dz=%f oz=%f\n", shot->ishot, ir, shot->zr[ir], shot->izr[ir], msh->nz, msh->dz, msh->oz); out_of_grid_limits = 1; }
        if(out_of_grid_limits)    return 1; 
    }

    resampleSourceAndRec(shot, msh->dt, msh->nt);

    // Save source in file
    // FILE *fSou_Hdr, *fSou_Bin;
    // initFilesSmart1d("source_resamp", &fSou_Hdr, &fSou_Bin, msh->nt, msh->dt, 0.0);
    // writeSamples(shot->source, fSou_Bin, msh->nt);
    // fclose(fSou_Bin);
    // fclose(fSou_Hdr);

return 0;
}

void computeDataBoundingBox(Dataset3D *dataset, float p1[], float p2[], float p3[], float p4[], \
                            float skirtLeft, float skirtRight, float skirtUp, float skirtDown, \
                            int applyCoordsShift)
{
    int ishot, irec;
    double LARGEDOUBLE = 1.0e+12;
    double  xmin=+LARGEDOUBLE, xmax=-LARGEDOUBLE, ymin=+LARGEDOUBLE, ymax=-LARGEDOUBLE;
    
    int nshot = dataset->nshot; printf("\n\n In [computeDataBoundingBox] nshot=%d \n\n", nshot);

    // Computes max and min x and y.
    for(ishot=0; ishot<nshot; ishot++)
    {
        Shot3D *shot = &((dataset->shot)[ishot]);

        if( xmin > shot->xs[0] )    xmin = shot->xs[0];
        if( xmax < shot->xs[0] )    xmax = shot->xs[0];
        if( ymin > shot->ys[0] )    ymin = shot->ys[0];
        if( ymax < shot->ys[0] )    ymax = shot->ys[0];

        // if(ishot%20==0)  printf("\n  shot->xs[0]=%lf   xmin=%lf    shot->ys[0]=%lf   ymin=%lf ", shot->xs[0], xmin, shot->ys[0], ymin);

        int nrec = shot->nrec;
        int irec;
        for(irec=0; irec<nrec; irec++)
        {
            if( xmin > shot->xr[irec] )    xmin = shot->xr[irec];
            if( xmax < shot->xr[irec] )    xmax = shot->xr[irec];
            if( ymin > shot->yr[irec] )    ymin = shot->yr[irec];
            if( ymax < shot->yr[irec] )    ymax = shot->yr[irec];
        }
    }

    xmin -= skirtLeft;
    xmax += skirtRight;
    ymin -= skirtDown;
    ymax += skirtUp;

    p1[0] = xmin;
    p1[1] = ymin;
    p2[0] = xmin;
    p2[1] = ymax;
    p3[0] = xmax;
    p3[1] = ymin;
    p4[0] = xmax;
    p4[1] = ymax;

    printf("\n Found dataset bounding box:  \n");
    printf("\n       p1x = %f    p1y = %f  ", p1[0], p1[1]);
    printf("\n       p2x = %f    p2y = %f  ", p2[0], p2[1]);
    printf("\n       p3x = %f    p3y = %f  ", p3[0], p3[1]);
    printf("\n       p4x = %f    p4y = %f  ", p4[0], p4[1]);

    printf("\n Relative dataset bounding box:  \n");
    printf("\n       p1x = %f    p1y = %f  ", p1[0]-xmin, p1[1]-ymin);
    printf("\n       p2x = %f    p2y = %f  ", p2[0]-xmin, p2[1]-ymin);
    printf("\n       p3x = %f    p3y = %f  ", p3[0]-xmin, p3[1]-ymin);
    printf("\n       p4x = %f    p4y = %f  ", p4[0]-xmin, p4[1]-ymin);

    if(applyCoordsShift)
    {
        for(ishot=0; ishot<nshot; ishot++)
        {
            Shot3D *shot = &((dataset->shot)[ishot]);

            shot->xs[0] -= xmin;
            shot->ys[0] -= ymin;
            // shot->ys[0] -= ymax;

            int nrec = shot->nrec;
            int irec;
            for(irec=0; irec<nrec; irec++)
            {
                shot->xr[irec] -= xmin;
                shot->yr[irec] -= ymin;
                // shot->yr[irec] -= ymax;
            }
        }
    }
    

return;
}

//*
void findDataGeometry_getShotLineDir(Dataset3D *dataset, float Vec_shotLine[]) 
{
    Vec_shotLine[0] = 0.0f;
    Vec_shotLine[1] = 0.0f;

    int    ishot, irec;
    double delta_x, delta_y;
    double delta_x_avg, delta_y_avg;
    double sign_x, sign_y;
    int    n_pos_x=0, n_neg_x=0;
    int    n_pos_y=0, n_neg_y=0;
    int    isPos_x, isPos_y;
    
    int nshot = dataset->nshot;

    for(ishot=0; ishot<nshot-1; ishot++)
    {
        Shot3D *shot_a = &((dataset->shot)[ishot  ]);
        Shot3D *shot_b = &((dataset->shot)[ishot+1]);

        if( (shot_b->xs)[0] - (shot_a->xs)[0] > 0.0 )    n_pos_x++;
        else                                             n_neg_x++;
        
        if( (shot_b->ys)[0] - (shot_a->ys)[0] > 0.0 )    n_pos_y++;
        else                                             n_neg_y++;

        delta_x += fabsf( (shot_b->xs)[0] - (shot_a->xs)[0] );
        delta_y += fabsf( (shot_b->ys)[0] - (shot_a->ys)[0] );
    }

    if      (n_pos_x>n_neg_x)    sign_x = +1;
    else if (n_pos_x<n_neg_x)    sign_x = -1;
    else                         sign_x =  0;

    if      (n_pos_y>n_neg_y)    sign_y = +1;
    else if (n_pos_y<n_neg_y)    sign_y = -1;
    else                         sign_y =  0;

    delta_x_avg = delta_x / (nshot-1);
    delta_y_avg = delta_y / (nshot-1);

    int count = 0;
    for(ishot=0; ishot<nshot-1; ishot++)
    {
        Shot3D *shot_a = &((dataset->shot)[ishot  ]);
        Shot3D *shot_b = &((dataset->shot)[ishot+1]);

        delta_x = (shot_b->xs)[0] - (shot_a->xs)[0];
        delta_y = (shot_b->ys)[0] - (shot_a->ys)[0];

        // int cond_dir_x = ( delta_x/fabsf(delta_x) == sign_x );
        // int cond_dir_y = ( delta_y/fabsf(delta_y) == sign_y );
        // int cond_dis_x = ( fabsf(delta_x) < 3*delta_x_avg );
        // int cond_dis_y = ( fabsf(delta_y) < 3*delta_y_avg );

        int cond_dir_x, cond_dir_y, cond_dis_x, cond_dis_y; 

        if(delta_x != 0.0)    cond_dir_x = ( delta_x/fabsf(delta_x) == sign_x );
        else                  cond_dir_x = 1;
        if(delta_y != 0.0)    cond_dir_y = ( delta_y/fabsf(delta_y) == sign_y );
        else                  cond_dir_y = 1;

        cond_dis_x = ( fabsf(delta_x) < 3*delta_x_avg );
        cond_dis_y = ( fabsf(delta_y) < 3*delta_y_avg );

        int cond = cond_dir_x && cond_dis_x  &&  cond_dir_y && cond_dis_y;

        if(cond) 
        {
            Vec_shotLine[0] += delta_x;
            Vec_shotLine[1] += delta_y;
            count++;
        }
    }
    Vec_shotLine[0] /= count;
    Vec_shotLine[1] /= count;

}
void findDataGeometry_getRecvLineDir(Dataset3D *dataset, float Vec_recvLine[])
{
    Vec_recvLine[0] = 0.0;
    Vec_recvLine[1] = 0.0;

    int    ishot, irec;
    double delta_x, delta_y;
    double delta_x_avg, delta_y_avg;
    double sign_x, sign_y;
    int    n_pos_x=0, n_neg_x=0;
    int    n_pos_y=0, n_neg_y=0;
    int    isPos_x, isPos_y;
    int    nrecTotal;

    int nshot = dataset->nshot;
    nrecTotal = 0;
    for(ishot=0; ishot<nshot; ishot++)
    {
        Shot3D *shot = &((dataset->shot)[ishot]);
        int nrec = shot->nrec;
        nrecTotal += nrec;
        for(irec=0; irec<nrec-1; irec++) 
        {
            if( (shot->xr)[irec+1] - (shot->xr)[irec] > 0.0 )    n_pos_x++;
            else                                                 n_neg_x++;
        
            if( (shot->yr)[irec+1] - (shot->yr)[irec] > 0.0 )    n_pos_y++;
            else                                                 n_neg_y++;

            delta_x += fabsf( (shot->xr)[irec+1] - (shot->xr)[irec] );
            delta_y += fabsf( (shot->yr)[irec+1] - (shot->yr)[irec] );
        }
    }

    if      (n_pos_x>n_neg_x)    sign_x = +1;
    else if (n_pos_x<n_neg_x)    sign_x = -1;
    else                         sign_x =  0;

    if      (n_pos_y>n_neg_y)    sign_y = +1;
    else if (n_pos_y<n_neg_y)    sign_y = -1;
    else                         sign_y =  0;

    delta_x_avg = delta_x / (nrecTotal-1);
    delta_y_avg = delta_y / (nrecTotal-1);

    int count = 0;
    for(ishot=0; ishot<nshot; ishot++)
    {
        Shot3D *shot = &((dataset->shot)[ishot]);
        int nrec = shot->nrec;
        for(irec=0; irec<nrec-1; irec++) 
        {
            delta_x = (shot->xr)[irec+1] - (shot->xr)[irec];
            delta_y = (shot->yr)[irec+1] - (shot->yr)[irec];

            int cond_dir_x, cond_dir_y, cond_dis_x, cond_dis_y; 

            if(delta_x != 0.0)    cond_dir_x = ( delta_x/fabsf(delta_x) == sign_x );
            else                  cond_dir_x = 1;
            if(delta_y != 0.0)    cond_dir_y = ( delta_y/fabsf(delta_y) == sign_y );
            else                  cond_dir_y = 1;

            cond_dis_x = ( fabsf(delta_x) < 3*delta_x_avg );
            cond_dis_y = ( fabsf(delta_y) < 3*delta_y_avg );

            int cond = cond_dir_x && cond_dis_x  &&  cond_dir_y && cond_dis_y;

            if(cond) 
            {
                Vec_recvLine[0] += delta_x;
                Vec_recvLine[1] += delta_y;
                count++;
            }
        }
    }
    Vec_recvLine[0] /= count;
    Vec_recvLine[1] /= count;

}
void changeStationHorizontalCoordinates(Dataset3D *dataset, double P0[], double angle)
{
    int ishot, isou, irec;
    int nshot = dataset->nshot;    
    for(ishot=0; ishot<nshot; ishot++)
    {
        Shot3D *shot = &((dataset->shot)[ishot]);
        
        int nrec = shot->nrec;
        int nsou = shot->nsou;

        for(isou=0; isou<nsou; isou++)
        {
            double xs = shot->xs[isou];
            double ys = shot->ys[isou];
            transformCoords_2D(&xs, &ys, P0, angle);
            (shot->xs)[isou] = xs;
            (shot->ys)[isou] = ys;
        }
        for(irec=0; irec<nrec; irec++) 
        {
            double xr = (shot->xr)[irec];
            double yr = (shot->yr)[irec];
            transformCoords_2D(&xr, &yr, P0, angle);
            (shot->xr)[irec] = xr;
            (shot->yr)[irec] = yr;
        }
    }

}


void changeStationHorizontalCoordinates2NewFrame(Dataset3D *dataset, double p1[], double p2[], double p3[], double p4[], double Pref_x, double Pref_y)
{
    double vec_x[2], vec_y[2], origin[2];
    
    // >>>>>>>>> Computing Origin and Iline/Xline vectors <<<<<<<<<<
    // origin
    origin[0] = p1[0]; 
    origin[1] = p1[1];
    
    // iline vector
    vec_x[0] = p2[0] - p1[0]; 
    vec_x[1] = p2[1] - p1[1];
    double Lx = sqrt( vec_x[0]*vec_x[0] + vec_x[1]*vec_x[1] );
    vec_x[0] /= Lx;
    vec_x[1] /= Lx;
    
    // xline vector
    vec_y[0] = p3[0] - p1[0];
    vec_y[1] = p3[1] - p1[1];
    double Ly = sqrt( vec_y[0]*vec_y[0] + vec_y[1]*vec_y[1] );
    vec_y[0] /= Ly;
    vec_y[1] /= Ly;
    
    // >>>>>>>>> Transforming coordinates <<<<<<<<<<
    int ishot, isou, irec;
    int nshot = dataset->nshot;
    for(ishot=0; ishot<nshot; ishot++)
    {
        Shot3D *shot = &((dataset->shot)[ishot]);
        
        int nrec = shot->nrec;
        int nsou = shot->nsou;

        for(isou=0; isou<nsou; isou++)
        {
            // transform shot coordinates
            double xs = shot->xs[isou];
            double ys = shot->ys[isou];
            transformCoordsReferenceSystem_2D(&xs, &ys, vec_x, vec_y, origin);
            (shot->xs)[isou] = xs;
            (shot->ys)[isou] = ys;
        }
        for(irec=0; irec<nrec; irec++)
        {
            double xr = (shot->xr)[irec];
            double yr = (shot->yr)[irec];
            transformCoordsReferenceSystem_2D(&xr, &yr, vec_x, vec_y, origin);
            (shot->xr)[irec] = xr;
            (shot->yr)[irec] = yr;
        }
    }

    printf("\n Final absolute grid:  \n");
    printf("\n       p1x = %f    p1y = %f  ", p1[0]+Pref_x, p1[1]+Pref_y);
    printf("\n       p2x = %f    p2y = %f  ", p2[0]+Pref_x, p2[1]+Pref_y);
    printf("\n       p3x = %f    p3y = %f  ", p3[0]+Pref_x, p3[1]+Pref_y);
    printf("\n       p4x = %f    p4y = %f  ", p4[0]+Pref_x, p4[1]+Pref_y);
    
    printf("\n Final relative grid:  \n");
    printf("\n       p1x = %f    p1y = %f  ", 0.0, 0.0);
    printf("\n       p2x = %f    p2y = %f  ",  Lx, 0.0);
    printf("\n       p3x = %f    p3y = %f  ", 0.0,  Ly);
    printf("\n       p4x = %f    p4y = %f  ",  Lx,  Ly);
}



void resampleSourceAndRec(Shot3D *shot, float dt, int nt)
{
    int it, isou, irec;
    float t;
    int nsou = shot->nsou;
    int nrec = shot->nrec;

    float *source     = CPU_zaloc1F(nt*nsou);
    float *seismogram = CPU_zaloc1F(nt*nrec);

    int r    = 4;
    int r_   = r - 1;
    int diam = 2*r;
    int ic;

    // Resampling source
    for(isou=0; isou<nsou; isou++)
    {
        int shiftShot = isou*shot->nt;
        for(it=0; it<nt; it++)
        {
            t = it*dt;
            int it_ = floorf(t/shot->dt);
            if (it_>=shot->nt-1)    it_ = shot->nt - 2;
            float alpha = (t - it_*shot->dt) / shot->dt;
            float *coef = interSinc1D_8P(alpha);
            for(ic=0; ic<diam; ic++)
            {
                int it_sinc = it_ + ic - r_;
                if(0<=it_sinc  &&  it_sinc<shot->nt)
                    source[it+isou*nt] += coef[ic] * shot->source[it_sinc + shiftShot];
            }
            free(coef);
            // float p1 = 1.0f - alpha; 
            // float p2 = alpha;
            // source[it+isou*nt] = p1 * shot->source[it_ + shiftShot] + p2 * shot->source[it_+1 + shiftShot];
        }
    }

    // Resampling receiver
    for(irec=0; irec<nrec; irec++)
    {
        int shiftRec = irec*shot->nt;
        for(it=0; it<nt; it++)
        {
            t = it*dt;
            int it_ = floorf(t/shot->dt);
            if (it_>=shot->nt-1)    it_ = shot->nt - 2;
            float alpha = (t - it_*shot->dt) / shot->dt;
            float *coef = interSinc1D_8P(alpha);
            for(ic=0; ic<diam; ic++)
            {
                int it_sinc = it_ + ic - r_;
                if(0<=it_sinc  &&  it_sinc<shot->nt)
                    seismogram[it+irec*nt] += coef[ic] * shot->seismogram[it_sinc + shiftRec];
            }
            free(coef);

            // float p1 = 1.0f - alpha;
            // float p2 = alpha;
            // seismogram[it+irec*nt] = p1 * shot->seismogram[it_ + shiftRec] + p2 * shot->seismogram[it_+1 + shiftRec];
        }
    }

    free(shot->source);
    free(shot->seismogram);
    shot->source     = source;
    shot->seismogram = seismogram;
    shot->dt         = dt;
    shot->nt         = nt;

return;
}

void applyTopMute_3D(Shot3D *shot)
{
    int is, ir, it;
    double xs, ys, zs;
    double xr, yr, zr;
    double R;
    double delta_x, delta_y, delta_z, delta_t;

    double t0 = shot->TopMute_t0;
    double t1 = shot->TopMute_t1;
    double off0 = shot->TopMute_off0;
    double off1 = shot->TopMute_off1;
    double ramp = shot->TopMute_ramp;
    double dtdR = (t1 - t0) / (off1 - off0);
    
    int    nt = shot->nt;
    double dt = shot->dt;
    double ot = shot->ot;
    double t, ta, tb, filter;

    for(is=0; is<shot->nsou; is++)
    {
        xs = shot->xs[is];
        ys = shot->ys[is];
        zs = shot->zs[is];
        for(ir=0; ir<shot->nrec; ir++)
        {
            xr = shot->xr[ir];
            yr = shot->yr[ir];
            zr = shot->zr[ir];
            delta_x = xs-xr;
            delta_y = ys-yr;
            delta_z = zs-zr;

            R  = sqrt(delta_x*delta_x + delta_y*delta_y + delta_z*delta_z);
            
            if     (R<off0              ) ta = t0;
            else if(R>off1              ) ta = t1;
            else if(off0<=R  &&  R<=off1) ta = t0 + (R-off0) * dtdR;
            tb = ta - ramp;

            // printf("\n R=%f  off0=%f  off1=%f  tb=%f  ta=%f\n", R, off0, off1, tb, ta);

            for(it=0; it<shot->nt; it++)
            {
                t = ot + it*dt;
                if     (t<tb            )    filter = 0.0;
                else if(t>ta            )    filter = 1.0;
                else if(tb<=t  &&  t<=ta)    filter = (t-tb)/(ta-tb);

                shot->seismogram[it+ir*nt] *= filter;
            }
        }
    }    
return;
}

void attenuateDA_3D(Shot3D *shot)
{

    float *seisDA = CPU_zaloc1F(shot->nt*shot->nrec);
    modelDA_3D(shot,seisDA,+1);
    // modelDA(shot,seisDA,-1);
    if( (removeDA_3D(shot,seisDA)) == 1 )  {
        printf("\n  Could not remove DA (direct arrival) from shot [%d], because offset range was not available \n", shot->ishot);
    }    
    free(seisDA);
return;
}

void modelDA_3D(Shot3D *shot, float *seisDA, int ghostSign)
{
    int is, ir, it;
    double xs, ys, zs;
    double xr, yr, zr;
    double R, Rh;
    double geometrical_spreading;
    double delta_x, delta_y, delta_z, delta_t;

    for(is=0; is<shot->nsou; is++)
    {
        xs = shot->xs[is];
        ys = shot->ys[is];
        zs = shot->zs[is] * ghostSign;
        for(ir=0; ir<shot->nrec; ir++)
        {
            xr = shot->xr[ir];
            yr = shot->yr[ir];
            zr = shot->zr[ir];
            delta_x = xs-xr;
            delta_y = ys-yr;
            delta_z = zs-zr;

            R  = sqrt(delta_x*delta_x + delta_y*delta_y + delta_z*delta_z);
            Rh = sqrt(delta_x*delta_x + delta_y*delta_y);
            geometrical_spreading = 1.0 / R;
            delta_t = R / ((double) shot->DA_waterSurfaceVel);

            for(it=0; it<shot->nt; it++)
            {
                int it_rec = it + floorf(delta_t/shot->dt);
                float alpha = delta_t/shot->dt - floorf(delta_t/shot->dt);
                float *coef = interSinc1D_8P(alpha);
                int r    = 4;
                int r_   = r - 1;
                int diam = 2*r;
                int ic;
                for(ic=0; ic<diam; ic++)
                {
                    int it0 = it_rec - r_ + ic;
                    if(0<=it0  &&  it0<shot->nt)
                        seisDA[it0+ir*shot->nt] += coef[ic] * shot->source[it+is*shot->nt] * geometrical_spreading * ghostSign;
                }
                free(coef);
            }
        }
    }

return;
}

int removeDA_3D(Shot3D *shot, float *seisDA)
{
    int is, ir, it;
    double xs, ys, zs;
    double xr, yr, zr;
    double R;
    double delta_x, delta_y, delta_z, delta_t;
    double ratio = 0;
    int ntraces = 0;

    // Computing optimal equalization factor (ratio)
    for(is=0; is<shot->nsou; is++)
    {
        // Find time of peak
        int it_max = 0;
        double maxAmp = 0.0;
        for(it=0; it<shot->nt; it++)
        {
            if( maxAmp < fabsf(shot->source[it+is*shot->nt]) ) 
            {
                maxAmp = fabsf(shot->source[it+is*shot->nt]);
                it_max = it;
            }
        }

        // Compute the ratio
        xs = shot->xs[is];
        ys = shot->ys[is];
        zs = shot->zs[is];
        for(ir=0; ir<shot->nrec; ir++)
        {
            xr = shot->xr[ir];
            yr = shot->yr[ir];
            zr = shot->zr[ir];
            delta_x = xs-xr;
            delta_y = ys-yr;
            delta_z = zs-zr;

            R  = sqrt(delta_x*delta_x + delta_y*delta_y + delta_z*delta_z);
            delta_t = R / ((double) shot->DA_waterSurfaceVel);
            int it_rec   = it_max + floorf(0.5f + delta_t/shot->dt);

            if(shot->DA_offMin<=R  &&  R<=shot->DA_offMax  &&  0<=it_rec  &&  it_rec<shot->nt)
            {
                ratio += shot->seismogram[it_rec+ir*shot->nt]/seisDA[it_rec+ir*shot->nt];
                ntraces++;
            }            
        }
    }
    if(ntraces>0)    ratio /= ntraces;
    else             return 1;

    printf("\n REMOVING DA ratio=%f  ntraces=%d\n", ratio, ntraces);

    // Removing direct arrival
    for(ir=0; ir<shot->nrec; ir++)
    {
        for(it=0; it<shot->nt; it++) {
            shot->seismogram[it+ir*shot->nt] -= ratio*seisDA[it+ir*shot->nt];
        }
    }

return 0;
}


void removeNoisyTraces(Shot3D *shot, float t0, float t1) {
    int it, ir;
    float  ot   = shot->ot;
    float  dt   = shot->dt;
    int    nt   = shot->nt;
    int    nrec = shot->nrec;

    for(ir=0; ir<nrec; ir++) 
    {
        double rx = shot->xr[ir] - shot->xs[0];
        double ry = shot->yr[ir] - shot->ys[0];
        double rz = shot->zr[ir] - shot->zs[0];
        double R = sqrt(rx*rx + ry*ry + rz*rz);
        double T = R / 1500.0;
        double taperT = t1-t0;

        // Compute trace statistics
        double max_amp = -1.0e6;
        double norm_trace_begin = 0;
        int n = 0;
        for(it=0; it<nt; it++) 
        {
            float t = ot + it*dt;
            if(t0<=t  &&  t<=(t1+T))  {
                norm_trace_begin += fabsf((shot->seismogram)[it+ir*nt]);
                n++;
            }
            else if( max_amp<fabsf((shot->seismogram)[it+ir*nt]) )    max_amp = fabsf((shot->seismogram)[it+ir*nt]);
        }
        if(n>0)    norm_trace_begin /= n;
        for(it=0; it<nt; it++)
        {
            float t = ot + it*dt;
            if(norm_trace_begin>=0.2*max_amp  ||  t<T) {
                (shot->seismogram)[it+ir*nt] = 0.0;
            }            
            if(T<=t  &&  t<T+taperT) {
                float taper = (t-T)/taperT; 
                (shot->seismogram)[it+ir*nt] *= taper;
            }
        }
    }

return;
}