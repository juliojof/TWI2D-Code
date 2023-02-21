void initDataset(Dataset *dataset, int nshot, float dx_rec, \
                 int nt, int nx, int nz, float dt, float dx, float dz, \
                 float offMax, float f1, float f2, float f3, float f4) {
    Shot *shot;
    int ishot;
    dataset->nshot = nshot;
    dataset->shot = CPU_zaloc1Shot(nshot);

    // printf("\n  Checking Shot:  nshot=%d  \n", nshot);

    for(ishot=0; ishot<nshot; ishot++)
    {
        shot = &((dataset->shot)[ishot]);
        initShot(shot, ishot, nshot, dx_rec, nt, nx, nz, dt, dx, dz, offMax, f1, f2, f3, f4);
    }

return;
}

//*
void initDatasetFromSEGY(FILE *fp, char *filename, Dataset *dataset) {

    int ntraces, nt;
    float dt;
    int *shotIndexes, *tracesPerShot;

    int nshot = indexShotGathersFromSEGY(fp, filename, &shotIndexes, &tracesPerShot, &nt, &dt, &ntraces);

    dataset->nshot = nshot;
    dataset->shot = CPU_zaloc1Shot(nshot);

    Shot *shot;
    int ishot;
    /*
    for(ishot=0; ishot<nshot; ishot++)
    {
        shot = &((dataset->shot)[ishot]);
        initShot(shot, ishot, nshot, dx_rec, nt, nx, nz, dt, dx, dz, offMax, f1, f2, f3, f4);
    }
    */

return;
}
//*/

// (Odisseia Shot) BEGIN
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
void initTimeLapseDataset3DfromFile(Dataset3D *dataset_base, Dataset3D *dataset_moni, Parms *parms)
{   
    int nshot;
    
    // Base
    printf("\n Input seismic file base: <%s> \n", parms->inputSeismicFileName_base);
    openOdisseiaFiles(&(dataset_base->fMet0), &(dataset_base->fMet1), &(dataset_base->fMet2), \
                      &(dataset_base->fData), parms->inputSeismicFileName_base, \
                      (dataset_base->OdisseiaFileType), "r");
    readOdisseiaMetaFile_get_nshot(dataset_base->fMet1, &nshot);
    dataset_base->nshot = nshot;
    dataset_base->shot  = CPU_zaloc1Shot3D(nshot);

    initShotFilteringParms(dataset_base, parms);

    // Monitor
    printf("\n Input seismic file moni: <%s> \n", parms->inputSeismicFileName_moni);
    openOdisseiaFiles(&(dataset_moni->fMet0), &(dataset_moni->fMet1), &(dataset_moni->fMet2), \
                      &(dataset_moni->fData), parms->inputSeismicFileName_moni, \
                      (dataset_moni->OdisseiaFileType), "r");
    readOdisseiaMetaFile_get_nshot(dataset_moni->fMet1, &nshot);
    dataset_moni->nshot = nshot;
    dataset_moni->shot  = CPU_zaloc1Shot3D(nshot);

    initShotFilteringParms(dataset_moni, parms);

return;
}
void initDataset3DfromFile(Dataset3D *dataset, Parms *parms)
{   
    // int OdisseiaFileType;
    printf("\n Input seismic file: <%s> \n", parms->inputSeismicFileName); fflush(stdout);
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

    acqGeom->recCoordRefFromShot  = parms->recCoordRefFromShot;

    // printf("\n  Debug:  acqGeom->shotNZ=%d  acqGeom->shotDZ=%f  \n ", acqGeom->shotNZ, acqGeom->shotDZ);
    // printf("\n  Debug:  acqGeom->recNZ=%d   acqGeom->recDZ=%f   \n ", acqGeom->recNZ , acqGeom->recDZ);

    acqGeom->f1  = parms->f1;
    acqGeom->f2  = parms->f2;
    acqGeom->f3  = parms->f3;
    acqGeom->f4  = parms->f4;

    acqGeom->nshot = acqGeom->shotNX * acqGeom->shotNY * acqGeom->shotNZ;
    acqGeom->nrec  = acqGeom->recNX  * acqGeom->recNY  * acqGeom->recNZ;

return;
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

    int recCoordRefFromShot = acqGeom->recCoordRefFromShot;


    // printf("\n  Debug:  recCoordRefFromShot=%d  \n ", recCoordRefFromShot);
    // printf("\n  Debug:  acqGeom->shotNX=%d  acqGeom->shotDX=%f  \n ", acqGeom->shotNX, acqGeom->shotDX);
    // printf("\n  Debug:  acqGeom->recNX=%d   acqGeom->recDX=%f   \n ", acqGeom->recNX , acqGeom->recDX);
    // printf("\n  Debug:  acqGeom->shotNZ=%d  acqGeom->shotDZ=%f  \n ", acqGeom->shotNZ, acqGeom->shotDZ);
    // printf("\n  Debug:  acqGeom->recNZ=%d   acqGeom->recDZ=%f   \n ", acqGeom->recNZ , acqGeom->recDZ);

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

                // printf("\n xs=%f  ys=%f  zs=%f  \n", shot->xs[0], shot->zs[0], shot->zs[0]);

                // Initializing receiver coordinates
                for(Riz=0; Riz<acqGeom->recNZ; Riz++)
                    for(Riy=0; Riy<acqGeom->recNY; Riy++)
                        for(Rix=0; Rix<acqGeom->recNX; Rix++)
                        {
                            int irec = Rix + Riy * acqGeom->recNX + Riz * acqGeom->recNX*acqGeom->recNY;
                            shot->xr[irec] = Rix * acqGeom->recDX + acqGeom->recOX + recCoordRefFromShot*shot->xs[0];
                            shot->yr[irec] = Riy * acqGeom->recDY + acqGeom->recOY + recCoordRefFromShot*shot->ys[0];
                            shot->zr[irec] = Riz * acqGeom->recDZ + acqGeom->recOZ + recCoordRefFromShot*shot->zs[0];
                            // printf("\n xr=%f  yr=%f  zr=%f   \n", shot->xr[irec], shot->zr[irec], shot->zr[irec]);
                        }

                //shot->source = getRicker(nt, dt, fm, t0);
                shot->source = getBandpassPulse3D(nt, dt, acqGeom->f1, acqGeom->f2, acqGeom->f3, acqGeom->f4, t0);
            }
}
void finalizeShot3D(Shot3D *shot) {
    finalizeShotSeismogramAndSource(shot);
    finalizeShot3DMetaInfo(shot);
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

void finalizeDataset3D(Dataset3D *dataset)
{
    Shot3D *shot;
    int ishot, nshot;
    nshot = dataset->nshot;
    free((dataset->shot));
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
// (Odisseia Shot)  END

void finalizeDataset(Dataset *dataset) {
    Shot *shot;
    int ishot, nshot;
    nshot = dataset->nshot;
    for(ishot=0; ishot<nshot; ishot++)
    {
        shot = &((dataset->shot)[ishot]);

        free(shot->source);
        free(shot->seismogram);

        free(shot->ixs);
        free(shot->izs);
        free(shot->xs);
        free(shot->zs);
        
        free(shot->ixr);
        free(shot->izr);
        free(shot->xr);
        free(shot->zr);
    }
    free((dataset->shot));
return;
}
/*
void initShotFromFile(Shot **shot, size_t itrace, size_t tracesInShot, int nt, float dt) 
{
    shot->nrec = tracesInSeismogram;
    shot->nshot = nshot;
    shot->nsou  = 1;
    shot->nt = nt;
    shot->dt = dt;
    shot->ot = 0.0;

    shot->seismogram = CPU_zaloc1F(shot->nrec*nt);


return;
}
*/
void initShot(Shot *shot, int ishot, int nshot, float dx_rec, \
              int nt, int nx, int nz, float dt, float dx, float dz, \
              float offMax, float f1, float f2, float f3, float f4) {

    shot->ishot = ishot;
    shot->nshot = nshot;
    shot->nsou  = 1;
    shot->nt    = nt;
    shot->dt    = dt;
    shot->ot    = 0.0;

    initSou(shot, ishot, nshot, nt, nx, nz, dt, dx, dz, offMax, f1, f2, f3, f4);
    initRec(shot, ishot, nshot, dx_rec, nx, nz, dx, dz);
        
    shot->seismogram = CPU_zaloc1F(shot->nrec*nt);

    // Save source in file
    FILE *fSou_Hdr, *fSou_Bin;
    initFilesSmart1d("source", &fSou_Hdr, &fSou_Bin, nt, dt, 0.0);
    writeSamples(shot->source, fSou_Bin, nt);
    fclose(fSou_Bin);
    fclose(fSou_Hdr);

    // applyHalfDerivative2Source(shot);
}

void initSou(Shot *shot, int ishot, int nshot, int nt, int nx, int nz, float dt, float dx, float dz, float offMax, float f1, float f2, float f3, float f4) {

    float t0 = 0.0f;
    float fm = 5.0f;

    float Lacq    = (nx-2*nxb-1)*dx;
    float dx_shot = floorf(0.5+Lacq/(nshot-1));
    float Lshots  = dx_shot*(nshot-1);
    float shift   = 0*0.5f * (Lacq - Lshots);
    float xs      = ishot * dx_shot + shift;
    int   ixs     = floorf(0.5f + xs/dx) + nxb;
    int   izs     = nzb + 5;

    if     (ixs>=nx) ixs = nx-1;
    else if(ixs<0  ) ixs = 0;
    if     (izs>=nz) izs = nz-1;
    else if(izs<0  ) izs = nzb;
    shot->ixs = CPU_zaloc1I(1);
    shot->izs = CPU_zaloc1I(1);
    shot->xs  = CPU_zaloc1F(1);
    shot->zs  = CPU_zaloc1F(1);
    shot->ixs[0] = ixs;
    shot->izs[0] = izs;
    shot->xs[0]  = (ixs-nxb) * dx;
    shot->zs[0]  = (izs-nzb) * dz;

    // printf("\n  Checking Shot:  ishot=%d  dx_shot=%f  dx=%f  xs=%f  ixs=%d\n", ishot, dx_shot, dx, xs, ixs);

    shot->ixMin = max( ((ixs-nxb)*dx - offMax)/dx - 1 + nxb,      nxb);
    shot->ixMax = min( ((ixs-nxb)*dx + offMax)/dx + 1 + nxb, nx-nxb-1);

    //shot->source = getRicker(nt, dt, fm, t0);
    shot->source = getBandpassPulse(nt, dt, f1, f2, f3, f4, t0);
    //applyShotTaper(shot, nx);
    applyHalfDerivative2Source(shot);
}


void initRec(Shot *shot, int ishot, int nshot, float dx_rec, int nx, int nz, float dx, float dz) {

    int ir, nrec;
    float Lacq   = (nx-2*nxb-1)*dx;
    nrec         = floorf(Lacq/dx_rec) + 1;
    shot->nrec   = nrec;
    float shift  = 0*0.5f*(Lacq-dx_rec*(nrec-1));
    shot->ixr    = CPU_zaloc1I(nrec);
    shot->izr    = CPU_zaloc1I(nrec);
    shot->xr     = CPU_zaloc1F(nrec);
    shot->zr     = CPU_zaloc1F(nrec);

    for(ir=0; ir<nrec; ir++) 
    {
        float xr = ir*dx_rec + shift;
        int  ixr = floorf(0.5f + xr/dx) + nxb;
        if(ixr>=nx) ixr=nx-1;
        shot->ixr[ir] = ixr;
        shot->izr[ir] = nzb;
        shot->xr[ir]  = (shot->ixr[ir]-nxb)*dx;
        shot->zr[ir]  = (shot->izr[ir]-nzb)*dz;
    }
}

void integrateRecTime(Shot *shot) {
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

void normalizeDatasetAmplitudes(Dataset *dataset) {

    Shot  *shot;
    float  dt;
    int    ishot, nshot, ir, nrec, it, nt;
    shot = &((dataset->shot)[0]);
    nshot = shot->nshot;
    nrec  = shot->nrec;
    nt    = shot->nt;
    dt    = shot->dt;
    // Find maximum amplitude
    double maxAmp = 0.0;
    for(ishot=0; ishot<nshot; ishot++) {
        shot = &((dataset->shot)[ishot]);
        for(ir=0; ir<nrec; ir++)
            for(it=0; it<nt; it++)
                if(maxAmp < fabs(shot->seismogram[ir*nt+it]))    maxAmp = fabs(shot->seismogram[ir*nt+it]);
    }
    // Normalize data
    for(ishot=0; ishot<nshot; ishot++) {
        shot = &((dataset->shot)[ishot]);
        for(ir=0; ir<nrec; ir++)
            for(it=0; it<nt; it++)
                shot->seismogram[ir*nt+it] /= maxAmp;
    }

return;
}

void AGC(Dataset *dataset, float time_window) {

    Shot  *shot;
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
            
            
            // computes amplitudes
            
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


void finalizeShot(Shot *shot) {
    free(shot->source);
    free(shot->ixr);
    free(shot->izr);
    free(shot->xr);
    free(shot->zr);
    free(shot->ixs);
    free(shot->izs);
    free(shot->xs);
    free(shot->zs);
    free(shot->seismogram);
return;
}

double addObjFunc(Shot *shot) {
    int it, ir;
    int nt = shot->nt;
    int nrec = shot->nrec;
    double sum = 0.0;
    for(ir=0; ir<nrec; ir++) {
        for(it=0; it<nt; it++) {
            sum += (double) ( (shot->seismogram)[it+ir*nt] * (shot->seismogram)[it+ir*nt] );
        }
    }
return sum;
}
void addShots(Shot *shot_F, float coef_A, Shot *shot_A, float coef_B, Shot *shot_B) {
    int it, ir;
    int nt = shot_A->nt;
    int nrec = shot_A->nrec;
    for(ir=0; ir<nrec; ir++)
        for(it=0; it<nt; it++) {
            (shot_F->seismogram)[it+ir*nt] = coef_A * (shot_A->seismogram)[it+ir*nt] + coef_B * (shot_B->seismogram)[it+ir*nt];
        }

return;
}
void removeDA(Shot *shot, Shot *shot_da) {
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
void removeDA_3D(Shot3D *shot, Shot3D *shot_da) {
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
void copyData(Shot *shot_inp, Shot *shot_out) {
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

void applyHalfDerivative2Source(Shot *shot) {
    applyHalfDerivative(shot->source, shot->nt, shot->dt, +1.0);
return;
}
void applyHalfDerivative2Source3D(Shot3D *shot) {
    applyHalfDerivative(shot->source, shot->nt, shot->dt, +1.0);
return;
}
void applyHalfDerivative2Data(Shot *shot) {
    int ir;
    for(ir=0; ir<shot->nrec; ir++) {
        applyHalfDerivative( (shot->seismogram)+ir*shot->nt, shot->nt, shot->dt, -1.0);
    }
return;
}
void applyHalfDerivative2Shot(Shot3D *shot) {
    int ir;
    for(ir=0; ir<shot->nrec; ir++) {
        applyHalfDerivative( (shot->seismogram)+ir*shot->nt, shot->nt, shot->dt, -1.0);
    }
return;
}
void applyHalfDerivative2Dataset(Dataset3D *dataset) {
    int ishot;
    for(ishot=0; ishot<dataset->nshot; ishot++)
        applyHalfDerivative2Shot(&(dataset->shot[ishot]));
}

void applyShotTaper(Shot *shot, int nx) {
    if(shot->nsou>1)    return;
    int it;
    float ixs = shot->ixs[0] - nxb;
    float nx0 = nx-2*nxb;
    float perc = 0.25;
    if(ixs<perc*nx0) {
        float taper = 1.0f - 0.85*(perc*nx0 - ixs)/(perc*nx0);
        for(it=0; it<shot->nt; it++)
            shot->source[it] *= taper;
    }
    if(ixs>(1-perc)*nx0) {
        float taper = 1.0f - 0.85*(ixs - (1-perc)*nx0)/(perc*nx0);
        for(it=0; it<shot->nt; it++)
            shot->source[it] *= taper;
    }
}

void applyRecTaper(Shot3D *shot, float perc, float offMax, int sensorSide)
{
    int   it, ir, nrec;
    float xs, zs, xr, zr;
    float offset, aoffset, taper;
    float off0, off1, maxOff;

    if(shot->nsou>1)    return;

    xs   = shot->xs[0];
    zs   = shot->zs[0];
    nrec = shot->nrec;

    maxOff = 0.0;
    for(ir=0; ir<nrec; ir++)
    {
        xr      =  shot->xr[ir];
        aoffset = fabsf(xs-xr);
        if(maxOff<aoffset)    maxOff = aoffset;
    }
    if(offMax>maxOff)    offMax = maxOff;

    off1 = offMax;
    off0 = (1.0-perc) * offMax;

    for(ir=0; ir<nrec; ir++)
    {
        xr =  shot->xr[ir];
        offset  = xs-xr;
        aoffset = fabsf(offset);
        
        if(off0>aoffset)
            taper = 1.0f;
        else if(off0<=aoffset && aoffset<off1)
            taper = 1.0 - (aoffset-off0)/(off1-off0);
        else
            taper = 0.0f;
        
        if( (fsignf(offset)!=sensorSide)  &&  sensorSide!=0 )
            taper = 0.0;
        
        for(it=0; it<shot->nt; it++)
            shot->seismogram[it+ir*shot->nt] *= taper;
    }
}


void offsetMaxMute(Shot3D *shot, float perc, float offMax, int sensorSide)
{
    int   it, ir, nrec;
    float xs, zs, xr, zr;
    float offset, aoffset, taper;
    float off0, off1;

    if(shot->nsou>1)    return;

    xs   = shot->xs[0];
    zs   = shot->zs[0];
    nrec = shot->nrec;

    off1 = offMax;
    off0 = (1.0-perc) * offMax;

    for(ir=0; ir<nrec; ir++)
    {
        xr =  shot->xr[ir];
        offset  = xs-xr;
        aoffset = fabsf(offset);
        
        if(off0>aoffset)
            taper = 1.0f;

        else if(off0<=aoffset && aoffset<off1)
            taper = 1.0 - (aoffset-off0)/(off1-off0);
            
        else
            taper = 0.0f;
        
        if( (fsignf(offset)!=sensorSide)  &&  sensorSide!=0 )
            taper = 0.0;
        
        for(it=0; it<shot->nt; it++)
            shot->seismogram[it+ir*shot->nt] *= taper;
    }
}

void offsetMinMute(Shot3D *shot, float perc, float offMin, int sensorSide)
{
    int   it, ir, nrec;
    float xs, zs, xr, zr;
    float offset, aoffset, taper;
    float off0, off1;

    if(shot->nsou>1)    return;

    xs   = shot->xs[0];
    zs   = shot->zs[0];
    nrec = shot->nrec;

    off1 = offMin;
    off0 = (1.0+perc) * offMin;

    for(ir=0; ir<nrec; ir++)
    {
        xr =  shot->xr[ir];
        offset  = xs-xr;
        aoffset = fabsf(offset);
        
        if(off0<aoffset)
            taper = 1.0f;

        else if(off0>=aoffset && aoffset>off1)
            taper = 1.0 - (aoffset-off1)/(off0-off1);
            
        else
            taper = 0.0f;
        
        if( (fsignf(offset)!=sensorSide)  &&  sensorSide!=0 )
            taper = 0.0;
        
        for(it=0; it<shot->nt; it++)
            shot->seismogram[it+ir*shot->nt] *= taper;
    }
}

void applyWaterBottomMute(Shot *shot, float waterBottomDepth, float taper)
{
    int   it, ir, nrec;
    float xs, zs, xr, zr;
    float aoffset;

    if(shot->nsou>1)    return;

    xs   = shot->xs[0];
    zs   = shot->zs[0];
    nrec = shot->nrec;

    float t0 = waterBottomDepth / 1500.0;
    for(ir=0; ir<nrec; ir++)
    {
        xr      =  shot->xr[ir];
        aoffset = fabsf(xs-xr);

        float h = aoffset/2.0;        
        float t2 = t0*t0 + h*h / (1500.0*1500.0);
        float waterBottomTraveltime = 2.0 * sqrtf(t2);

        for(it=0; it<shot->nt; it++)
        {
            float t = (it * shot->dt) - shot->ot;
            float mute;
            if(t < (waterBottomTraveltime-taper) )    mute = 0.0;
            else if( t < waterBottomTraveltime )      mute = ( t - (waterBottomTraveltime-taper) ) / taper;
            else mute = 1.0;
            shot->seismogram[it+ir*shot->nt] *= mute;
        }

    }

return;
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

return 0;
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

float *getBandpassPulse(int nt, float dt, float f1, float f2, float f3, float f4, float t0) {
    int it;
    float *pulse = CPU_zaloc1F(nt);
    t0 += 1.0f/(0.5*(f2+f1));
    int it0 = (int) (0.5f + t0/dt); 
    pulse[it0] = 1.0f;
    bandpass(pulse, nt, dt, f1, f2, f3, f4);
    for(it=0; it<nt; it++) {
        float t = it*dt;
        if(t<t0/3) {
            pulse[it]       *= t/(t0/3);
            pulse[2*it0-it] *= t/(t0/3);
        }
        else if(it>2*it0)    pulse[it] = 0.0f;
    }
return pulse;
}

void applyBandpassToTrace(Shot *shot, float f1, float f2, float f3, float f4) {
    int ir;
    for(ir=0; ir<shot->nrec; ir++)
        bandpass((shot->seismogram)+(ir*shot->nt), shot->nt, shot->dt, f1, f2, f3, f4);
return;
}
void multiplySeismogramByTime(Shot *shot, float power) {
    int ir;
    for(ir=0; ir<shot->nrec; ir++)
        multiplyTraceByTime((shot->seismogram)+(ir*shot->nt), shot->nt, shot->dt, power);
return;
}
void multiplyTraceByTime(float *trace, int nt, float dt, float power) {
    int it;
    for(it=0; it<nt; it++)
        trace[it] *= pow(it*dt, power);
}

void demultiplexSeismogram(float *seismogram, int nt, int nrec) {
    float *seisTemp = CPU_zaloc1F(nt*nrec);
    int ir, it, i;
    for(ir=0; ir<nrec; ir++)
        for(it=0; it<nt; it++) {
            seisTemp[it+ir*nt] = seismogram[ir+it*nrec];
    }
    for(i=0; i<nrec*nt; i++)    seismogram[i] = seisTemp[i];

    free(seisTemp);
}


void makeEncodedDataset(int idev, Dataset *data_inp, Dataset *data_tmp, Dataset *data_enc, Dataset *data_res, int codelength, int codestep, int realizations)
{  
    Shot *shot_inp, *shot_tmp, *shot_enc, *shot_res;
    int   ishot, irec, it, nshot, nsou, nrec, nt, instance, seed;
    float dt;

    data_enc->nshot = realizations;
    data_res->nshot = realizations;
    data_enc->shot  = CPU_zaloc1Shot(realizations);
    data_res->shot  = CPU_zaloc1Shot(realizations);

    nt = shotEncoder(idev, data_inp, NULL, codelength, codestep, 0);
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
        shotEncoder(idev, data_inp, data_tmp, codelength, codestep, seed);
        
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
        shot_enc->izs = CPU_zaloc1I(nshot);
        shot_enc->xs  = CPU_zaloc1F(nshot);
        shot_enc->zs  = CPU_zaloc1F(nshot);

        shot_res->ixs = CPU_zaloc1I(nshot);
        shot_res->izs = CPU_zaloc1I(nshot);
        shot_res->xs  = CPU_zaloc1F(nshot);
        shot_res->zs  = CPU_zaloc1F(nshot);

        shot_enc->ixr = CPU_zaloc1I(nrec);
        shot_enc->izr = CPU_zaloc1I(nrec);
        shot_enc->xr  = CPU_zaloc1F(nrec);
        shot_enc->zr  = CPU_zaloc1F(nrec);

        shot_res->ixr = CPU_zaloc1I(nrec);
        shot_res->izr = CPU_zaloc1I(nrec);
        shot_res->xr  = CPU_zaloc1F(nrec);
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
            memcpy((shot_enc->izr)+shift, (shot_tmp->izr), bytes);
            memcpy((shot_res->ixr)+shift, (shot_tmp->ixr), bytes);
            memcpy((shot_res->izr)+shift, (shot_tmp->izr), bytes);

            shift = nr_;
            bytes = shot_tmp->nrec * sizeof(float);
            memcpy((shot_enc->xr)+shift, (shot_tmp->xr), bytes);
            memcpy((shot_enc->zr)+shift, (shot_tmp->zr), bytes);
            memcpy((shot_res->xr)+shift, (shot_tmp->xr), bytes);
            memcpy((shot_res->zr)+shift, (shot_tmp->zr), bytes);

            shot_enc->ixs[ishot] = shot_tmp->ixs[0];
            shot_enc->izs[ishot] = shot_tmp->izs[0];
            shot_enc->xs[ishot]  = shot_tmp->xs[0];
            shot_enc->zs[ishot]  = shot_tmp->zs[0];

            shot_res->ixs[ishot] = shot_tmp->ixs[0];
            shot_res->izs[ishot] = shot_tmp->izs[0];
            shot_res->xs[ishot]  = shot_tmp->xs[0];
            shot_res->zs[ishot]  = shot_tmp->zs[0];

            nr_ += shot_tmp->nrec;
        }
    }

return;
}


void windowData(int idev, Dataset *dataset, float tmin, float tmax, float dt)
{
    int ishot, it, it0, isou, irec;
    int nshot = dataset->nshot;
    for(ishot=0; ishot<nshot; ishot++)
    {
        // Resamp source
        Shot *shot = &(dataset->shot[ishot]);
        int   nsou = shot->nsou;
        int   nrec = shot->nrec;
        int   nt0 = shot->nt;
        float dt0 = shot->dt;

        float *source = CPU_zaloc1F(nt0*nsou);
        float *seismo = CPU_zaloc1F(nt0*nrec);
        memcpy(source, shot->source    , nt0*nsou*sizeof(float));
        memcpy(seismo, shot->seismogram, nt0*nrec*sizeof(float));

        // if(ishot==0)
        // for(irec=0; irec<nrec; irec++)
        // {
        //     for(it=0; it<nt0; it+=100)
        //         printf("\n irec=%d  it=%d  seismo[irec*nt0+it]=%f  shot->seismogram[irec*nt0+it]=%f", irec, it, seismo[irec*nt0+it], shot->seismogram[irec*nt0+it]);
        // }

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
            GPU_WRAP_SeismogramResamp(idev, source, shot->source    , nsou, nt0, dt0, nt, dt);
            GPU_WRAP_SeismogramResamp(idev, seismo, shot->seismogram, nrec, nt0, dt0, nt, dt);
        }
        free(source);
        free(seismo);
    }

return;
}


void copyDataset(Dataset *dst, Dataset *src)
{  
    Shot *shot_inp, *shot_tmp, *shot_enc, *shot_res;
    int ishot, irec, it, nshot, nsou, nrec, nt, instance, seed;

    dst->nshot = src->nshot;
    nshot      = dst->nshot;

    dst->shot  = CPU_zaloc1Shot(nshot);

    // printf("\n >>>>>>> nshot=%d \n", nshot);

    int ns_ = 0;
    int nr_ = 0;
    for(ishot=0; ishot<nshot; ishot++)
    {
        
        Shot *shot_dst = &((dst->shot)[ishot]);
        Shot *shot_src = &((src->shot)[ishot]);

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
        shot_dst->izs = CPU_zaloc1I(nsou);
        shot_dst->xs  = CPU_zaloc1F(nsou);
        shot_dst->zs  = CPU_zaloc1F(nsou);
        shot_dst->ixr = CPU_zaloc1I(nrec);
        shot_dst->izr = CPU_zaloc1I(nrec);
        shot_dst->xr  = CPU_zaloc1F(nrec);
        shot_dst->zr  = CPU_zaloc1F(nrec);

        shot_dst->source     = CPU_zaloc1F(nsou*nt);
        shot_dst->seismogram = CPU_zaloc1F(nrec*nt);


        size_t bytes = nsou * nt * sizeof(float);
        memcpy((shot_dst->source), (shot_src->source), bytes);

        bytes = nrec * nt * sizeof(float);
        memcpy((shot_dst->seismogram), (shot_src->seismogram), bytes);

        bytes = nsou * sizeof(int);
        memcpy((shot_dst->ixs), (shot_src->ixs), bytes);
        memcpy((shot_dst->izs), (shot_src->izs), bytes);
        bytes = nsou * sizeof(float);
        memcpy((shot_dst->xs ), (shot_src->xs ), bytes);
        memcpy((shot_dst->zs ), (shot_src->zs ), bytes);

        bytes = nrec * sizeof(int);
        memcpy((shot_dst->ixr), (shot_src->ixr), bytes);
        memcpy((shot_dst->izr), (shot_src->izr), bytes);
        bytes = nrec * sizeof(float);
        memcpy((shot_dst->xr ), (shot_src->xr ), bytes);
        memcpy((shot_dst->zr ), (shot_src->zr ), bytes);        
    }

return;
}

int shotEncoder(int idev, Dataset *data_in, Dataset *data_ou, int codelength, int codestep, int seed)
{  
    float *out;
    Shot *shot;
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
            out = GPU_WRAP_SeismogramCoding(idev, trc, codeForThisShot, nrec, nt, ncdl);
            memcpy((data_ou->shot)[ishot].seismogram, out, nrec*nt_*sizeof(float));
            free(out);
        }
    }
    free(codes);
    
return nt_;
}

void correlateWavelets(Dataset *data)
{  
    float *out;
    Shot *shot;
    Shot *shot_temp;
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



void writeSeisData2File(Dataset *dataset, int jshot, int iproc, char *fileNameSeismos, char *fileNameSources)
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
            Shot *shot = &((dataset->shot)[ishot]);
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