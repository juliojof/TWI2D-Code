int accumulateModel(float *model_out, MSH *msh_out, float *model_inp, MSH *msh_inp);
int accumulateExtendedModel(float *model_out, MSH *msh_out, float *model_inp, MSH *msh_inp, long long int lx0);
int getShotExtendedModel(float *model_ref, MSH *msh_ref, float *model_sht, MSH *msh_sht, long long int lx0);
int getShotGrid(MSH *msh_out, float *model_inp, MSH *msh_inp, Shot3D *shot, Parms *parms);
int getShotModel(float **model_out, MSH *msh_out, float *model_inp, MSH *msh_inp, Shot3D *shot, Parms *parms);
int getShotFrame(Frame *shotFrame, MSH *mshShot, MSH *mshRef, Shot3D *shot, float AppertureExt_x);
void findMinMaxWithinFrame(float *min, float *max, Frame *frame, MSH *msh, float *array);
void interpShotGrid_GPUCompliant(MSH *ou, MSH *in, Parms *parms, int strideX, int strideY, int strideZ, float spaceSampling, float vmax);
float getSpaceSampling(float maxFrequency, float vmin);
float getStableTimeSampling(float spaceSampling, float vmax);
void stackLargeArray_MPI(float *array, int inode, int inode_manager, long long int N, MPI_Datatype datatype, MPI_Comm comm);
void stackArray_MPI(float *array, int inode, int inode_manager, long long int N, MPI_Datatype datatype, MPI_Comm comm);

void AcousticDataModeling(Dataset3D *dataset, Dataset3D *dataset_da, float *vp, float *vp_da, \
                          MSH *msh, Parms *parms, int iproc, int nproc, int readInputData) {
    
    int   ishot, nshot=dataset->nshot;
    time_t time_now, time_start;
    time_start = time(NULL);

    FILE *fSeis_inp_Hdr, *fSeis_inp_Bin;
    FILE *fSeis_da_Hdr, *fSeis_da_Bin;

    for(ishot=iproc; ishot<nshot; ishot+=nproc)
    {
        // Initializing seismogram memoery
        initSeismogramMemory(dataset, ishot);
        Shot3D *shot = &(dataset->shot)[ishot];
        applyHalfDerivative2Source3D(shot);

        MSH msh_shot;
        int processThisShot;
        if(getShotGrid(&msh_shot, vp, msh, shot, parms))    processThisShot = 0;
        else                                                processThisShot = 1;

        // Opening and saving files with input and EBM shots
        if(iproc==0  &&  ishot==0) initFilesSmart3d("IN_Seismogram", &fSeis_inp_Hdr, &fSeis_inp_Bin, shot->nt, shot->dt, shot->ot, shot->nrec, 1, 0, dataset->nshot/nproc, 1, 0);
        if(iproc==0  &&  ishot==0) initFilesSmart3d("DA_Seismogram", &fSeis_da_Hdr, &fSeis_da_Bin, shot->nt, shot->dt, shot->ot, shot->nrec, 1, 0, dataset->nshot/nproc, 1, 0);

        if(processThisShot)
        {
            float *vpShot;
            getShotModel(&vpShot, &msh_shot, vp, msh, shot, parms);
            
            if(iproc==0)    printf(" (%3d ; ", ishot);
            modelShot_GPU(procGPU, vpShot, &(dataset->shot)[ishot], msh_shot.jt, msh_shot.nt, \
                          msh_shot.nx, msh_shot.nz, msh_shot.nxb, msh_shot.nzb, \
                          msh_shot.dx, msh_shot.dz, msh_shot.dt, \
                          msh_shot.ox, msh_shot.oz, msh_shot.ot, 1, 0);
            if(ishot%nproc==0)    writeSamples(shot->seismogram   , fSeis_inp_Bin, shot->nt    * shot->nrec);
            time_now = time(NULL);
            if(iproc==0)    printf("%3ld)", time_now-time_start);
            free(vpShot);

            if(vp_da!=NULL  &&  dataset_da!=NULL)
            {
                initSeismogramMemory(dataset_da, ishot);
                Shot3D *shot_da = &(dataset_da->shot)[ishot];
                applyHalfDerivative2Source3D(shot_da);

                float *vpShot;
                getShotModel(&vpShot, &msh_shot, vp_da, msh, shot_da, parms);
                initSeismogramMemory(dataset_da, ishot);
                modelShot_GPU(procGPU, vpShot, &(dataset_da->shot)[ishot], msh_shot.jt, msh_shot.nt, \
                          msh_shot.nx, msh_shot.nz, msh_shot.nxb, msh_shot.nzb, \
                          msh_shot.dx, msh_shot.dz, msh_shot.dt, \
                          msh_shot.ox, msh_shot.oz, msh_shot.ot, 0, 0);
                if(ishot%nproc==0)    writeSamples(shot_da->seismogram, fSeis_da_Bin, shot_da->nt * shot_da->nrec);
                removeDA_3D(shot, shot_da);
                finalizeSeismogramAndSource(dataset_da, ishot);
                free(vpShot);
            }
            // applyRecTaper(&((dataset->shot)[ishot]), parms->offPerc, parms->offMax, parms->acquisitionGeom);
        }

        saveDataset3D_binaryData(dataset, FILE_POINTERS_INPUT, ishot);
        finalizeSeismogramAndSource(dataset, ishot);
    }

    if(iproc==0)
    {
        fclose(fSeis_inp_Hdr);
        fclose(fSeis_inp_Bin);
        fclose(fSeis_da_Hdr);
        fclose(fSeis_da_Bin);
    }

    MPI_Barrier(MPI_COMM_WORLD);

return;
}

void migration_TAU(Dataset3D *dataset, float *image, float *ilumin, float *Tcig, float *vp_bg,\
                   MSH *msh, Eimage_MSH *eimsh, Parms *parms, int ism, int recordMovie, int iproc, int nproc) {
    
    int ret = 0;

    long long int   nt = (dataset->shot)[0].nt;
    float           dt = (dataset->shot)[0].dt;
    float           ot = ((dataset->shot)[0]).ot;
    
    long long int ishot, nshot;
    nshot = dataset->nshot;

    time_t time_now, time_start;
    time_start = time(NULL);

    if(iproc==0)
        printf("\n\n Starting migration nshot=%lld \n (Shot,time): \n", nshot);    fflush(stdout);

    for(ishot=iproc; ishot<nshot; ishot+=nproc)
    {
        initShot3DfromFile(dataset, ishot);
        loadSeismogram2Dataset(dataset, ishot, parms->normalizeInputDataAmplitudes);
        Shot3D *shot = &(dataset->shot)[ishot];
        applyHalfDerivative2Source3D(shot);
        applyHalfDerivative2Shot(shot);
        
        MSH msh_shot;
        float *vpShot;
        int processThisShot;
        if(getShotGrid(&msh_shot, vp_bg, msh, shot, parms))    processThisShot = 0;
        else                                                   processThisShot = 1;

        if(processThisShot)
        {
            getShotModel(&vpShot, &msh_shot, vp_bg, msh, shot, parms);

            long long int size_xz     = msh_shot.nx * msh_shot.nz;
            long long int size_xz_tau = size_xz * ((long long int) eimsh->ntau);

            float *image_shot  = CPU_zaloc1F(size_xz);
            float *ilumin_shot = CPU_zaloc1F(size_xz);
            float *Tcig_shot  = CPU_zaloc1F(size_xz_tau);
            emigrateShot_TAU_Mem_GPU(procGPU, image_shot, ilumin_shot, Tcig_shot, vpShot, \
                                     &((dataset->shot)[ishot]), &msh_shot, eimsh, (recordMovie && ism==ishot), ism);
            
            accumulateModel( image, msh,  image_shot, &msh_shot);
            accumulateModel(ilumin, msh, ilumin_shot, &msh_shot);
            accumulateExtendedModel(Tcig, msh, Tcig_shot, &msh_shot, eimsh->tau0);

            free(image_shot);
            free(ilumin_shot);
            free(Tcig_shot);
            free(vpShot);
        }
        ret = finalizeShot3DfromFile(dataset, ishot);
        if(ret!=0)    printf("\n CANNOT FINALIZE SHOT %lld \n", ishot);
        ret = finalizeSeismogramAndSource(dataset, ishot);
        if(ret!=0)    printf("\n CANNOT FINALIZE SHOT %lld \n", ishot);

        time_now = time(NULL);
        if(iproc==0)    printf(" (%3lld ; %3ld) ", ishot, time_now-time_start);
    } // loop over shots
    MPI_Barrier(MPI_COMM_WORLD);

    if(iproc==0) {
        printf("\n\n Migrated shots \n\n");
        printf("\n Elapsed time: %3ld \n", time(NULL)-time_start);
    }

    long long int ntau  = eimsh->ntau;
    long long int nx    = msh->nx;
    long long int nz    = msh->nz;
    // Reducing results into manager and redistributing back to every process
    stackLargeArray_MPI(Tcig  , iproc, 0, nx*nz*ntau, MPI_FLOAT, MPI_COMM_WORLD);
    stackLargeArray_MPI(image , iproc, 0, nx*nz     , MPI_FLOAT, MPI_COMM_WORLD);
    stackLargeArray_MPI(ilumin, iproc, 0, nx*nz     , MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

return;
}

void migration_HOCIG_VOCIG(Dataset3D *dataset, float *image, float *Hocig, float *Vocig, float *vp_bg,\
                           MSH *msh, Eimage_MSH *eimsh, Parms *parms, int ism, int recordMovie, int iproc, int nproc) {
    
    int ret = 0;

    long long int   nt = (dataset->shot)[0].nt;
    float           dt = (dataset->shot)[0].dt;
    float           ot = ((dataset->shot)[0]).ot;

    float *ilumin = CPU_zaloc1F(msh->nx * msh->nz);
    
    long long int ishot, nshot;
    nshot = dataset->nshot;

    time_t time_now, time_start;
    time_start = time(NULL);

    if(iproc==0)
        printf("\n\n Starting migration nshot=%lld \n (Shot,time): \n", nshot);    fflush(stdout);

    for(ishot=iproc; ishot<nshot; ishot+=nproc)
    {
        initShot3DfromFile(dataset, ishot);
        loadSeismogram2Dataset(dataset, ishot, parms->normalizeInputDataAmplitudes);
        Shot3D *shot = &(dataset->shot)[ishot];
        applyHalfDerivative2Source3D(shot);
        applyHalfDerivative2Shot(shot);
        
        MSH msh_shot;
        float *vpShot;
        int processThisShot;
        if(getShotGrid(&msh_shot, vp_bg, msh, shot, parms))    processThisShot = 0;
        else                                                   processThisShot = 1;

        if(processThisShot)
        {
            getShotModel(&vpShot, &msh_shot, vp_bg, msh, shot, parms);

            long long int size_xz    = msh_shot.nx * msh_shot.nz;
            long long int size_xz_lx = size_xz * eimsh->nlx;
            long long int size_xz_lz = size_xz * eimsh->nlz;

            float *image_shot  = CPU_zaloc1F(size_xz);
            float *ilumin_shot = CPU_zaloc1F(size_xz);
            float *Hocig_shot  = CPU_zaloc1F(size_xz_lx);
            float *Vocig_shot  = CPU_zaloc1F(size_xz_lz);
            emigrateShot_HOCIG_VOCIG_Mem_GPU(procGPU, image_shot, ilumin_shot, Hocig_shot, Vocig_shot, vpShot, \
                                             &((dataset->shot)[ishot]), &msh_shot, eimsh, (recordMovie && ism==ishot), ism);
            
            accumulateModel( image, msh,  image_shot, &msh_shot);
            accumulateModel(ilumin, msh, ilumin_shot, &msh_shot);
            accumulateExtendedModel(Hocig, msh, Hocig_shot, &msh_shot, eimsh->lx0);
            accumulateExtendedModel(Vocig, msh, Vocig_shot, &msh_shot, eimsh->lz0);

            free(image_shot);
            free(ilumin_shot);
            free(Hocig_shot);
            free(Vocig_shot);
            free(vpShot);
        }
        ret = finalizeShot3DfromFile(dataset, ishot);
        if(ret!=0)    printf("\n CANNOT FINALIZE SHOT %lld \n", ishot);
        ret = finalizeSeismogramAndSource(dataset, ishot);
        if(ret!=0)    printf("\n CANNOT FINALIZE SHOT %lld \n", ishot);

        time_now = time(NULL);
        if(iproc==0)    printf(" (%3lld ; %3ld) ", ishot, time_now-time_start);
    } // loop over shots
    MPI_Barrier(MPI_COMM_WORLD);

    if(iproc==0) {
        printf("\n\n Migrated shots \n\n");
        printf("\n Elapsed time: %3ld \n", time(NULL)-time_start);
    }

    long long int nlx   = eimsh->nlx;
    long long int nlz   = eimsh->nlz;
    long long int nx    = msh->nx;
    long long int nz    = msh->nz;
    // Reducing results into manager and redistributing back to every process
    stackLargeArray_MPI(Hocig , iproc, 0, nx*nz*nlx, MPI_FLOAT, MPI_COMM_WORLD);
    stackLargeArray_MPI(Vocig , iproc, 0, nx*nz*nlz, MPI_FLOAT, MPI_COMM_WORLD);
    stackLargeArray_MPI(image , iproc, 0, nx*nz    , MPI_FLOAT, MPI_COMM_WORLD);
    stackLargeArray_MPI(ilumin, iproc, 0, nx*nz    , MPI_FLOAT, MPI_COMM_WORLD);
    long long int rx=100, rz=5;
    modelSmooth(ilumin, ilumin, nx, nz, rx, rz);
    // normalizeArray(ilumin,     nx*nz);
    normalizeArray(image,     nx*nz);
    // normalizeArray( Hocig, nlx*nx*nz);
    // normalizeArray( Vocig, nlz*nx*nz);
    float epsilon = 0.005;
    // applyIlum2EImage(Hocig, Hocig, image, image, ilumin, nx, nz, dx, lx0, epsilon);
    // applyIlum2EImage(Vocig, Vocig, image, image, ilumin, nx, nz, dx, lz0, epsilon);
    free(ilumin);

    // normalizeArray( image,     nx*nz);
    // normalizeArray( Hocig, nlx*nx*nz);
    // normalizeArray( Vocig, nlz*nx*nz);
    MPI_Barrier(MPI_COMM_WORLD);

return;
}

//*
void migration_TLCIG(Dataset *dataset, float *image, float *TLcig, float *vp_bg,\
                     long long int nx, long long int nz, \
                     float dx, float dz, float ox, float oz, long long int lt0, \
                     long long int jt, int ism, int recordMovie, int iproc, int nproc) {
    
    float *ilumin = CPU_zaloc1F(nx*nz);
    
    long long int   nt = (dataset->shot)[0].nt;
    float dt = (dataset->shot)[0].dt;
    float ot = ((dataset->shot)[0]).ot;

    long long int nlt = 2*lt0 + 1;
    
    long long int ishot, nshot;
    nshot = dataset->nshot;

    time_t time_now, time_start;
    time_start = time(NULL);

    if(iproc==0)
        printf("\n\n Starting migration nshot=%lld \n (Shot,time): \n", nshot);

    
    int lx0=0, lz0=0;
    GPU_OPER_initEnvironment(procGPU, nx, nz, nt, dx, dz, dt, ((dataset->shot)[0]).nsou, ((dataset->shot)[0]).nrec, lx0, lz0, lt0, nxb, nzb, MODE_CIGS_FULL);
    GPU_OPER_alocArrays_acoustic();

    for(ishot=iproc; ishot<nshot; ishot+=nproc)
    {
        if(ishot==iproc)    recordMovie = 1;
        else                recordMovie = 0;
        recordMovie = 0;
        emigrateShot_TLCIG_Mem_GPU(procGPU, image, ilumin, TLcig, vp_bg, &((dataset->shot)[ishot]), \
                                   jt, nt, nx, nz, dx, dz, dt, ox, oz, ot, lt0, recordMovie, ishot);
        time_now = time(NULL);
        if(iproc==0)    printf(" (%3lld ; %3ld) ", ishot, time_now-time_start);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    GPU_OPER_copyTLCIG_GPU2CPU(TLcig);
    GPU_OPER_copyIlumin_GPU2CPU(ilumin);
    GPU_OPER_copyImage_GPU2CPU(image);
    GPU_OPER_freeArrays_acoustic();

    // Reducing results into manager and redistributing back to every process
    float *image_tmp;
    float *TLcig_tmp;
    float *ilumin_tmp;
    if(iproc==0) {
        image_tmp  = CPU_zaloc1F(nx*nz);
        TLcig_tmp  = CPU_zaloc1F(nx*nz*nlt);
        ilumin_tmp = CPU_zaloc1F(nx*nz);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce( image,  image_tmp, nx*nz    , MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce( TLcig,  TLcig_tmp, nx*nz*nlt, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ilumin, ilumin_tmp, nx*nz    , MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(iproc==0)
    {
        memcpy( image,  image_tmp, nx*nz    *sizeof(float));
        memcpy( TLcig,  TLcig_tmp, nx*nz*nlt*sizeof(float));
        memcpy(ilumin, ilumin_tmp, nx*nz    *sizeof(float));
        free(image_tmp);
        free(TLcig_tmp);
        free(ilumin_tmp);
    }
    else {
        memset( image, 0, nx*nz    *sizeof(float));
        memset( TLcig, 0, nx*nz*nlt*sizeof(float));
        memset(ilumin, 0, nx*nz    *sizeof(float));
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast( image, nx*nz    , MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast( TLcig, nx*nz*nlt, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(ilumin, nx*nz    , MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    long long int rx=1000, rz=5;
    modelSmooth(ilumin, ilumin, nx, nz, rx, rz);
    
    // normalizeArray(ilumin,     nx*nz);
    normalizeArray( image,     nx*nz);
    // normalizeArray( TLcig, nlt*nx*nz);
    
    float epsilon = 0.005;
    // applyIlum2EImage(TLcig, TLcig, image, image, ilumin, nx, nz, dx, lt0, epsilon);
    free(ilumin); 

    // normalizeArray( image,     nx*nz);
    // normalizeArray( TLcig, nlt*nx*nz);

    MPI_Barrier(MPI_COMM_WORLD);    

return;
}
//*/

double extendedBornModeling(Dataset *dataset, float *eimage, float *Hocig, float *Vocig, float *vp_bg, \
                            long long int nx, long long int nz, long long int lx0, long long int lz0, long long int lt0, \
                            float dx, float dz, float ox, float oz, \
                            long long int jt, float perc, float offMax, int acquisitionGeom, \
                            char *data_FileName, int iproc, int nproc) {

    long long int   nt = ((dataset->shot)[0]).nt;
    float dt = ((dataset->shot)[0]).dt;
    float ot = ((dataset->shot)[0]).ot;
    
    float power = 1.2;

    Shot *shot = dataset->shot;

    double objFunc = 0;
    long long int ishot, nshot;
    nshot = ((dataset->shot)[0]).nshot;

    time_t time_now, time_start;
    time_start = time(NULL);

    FILE *fSeis_Hdr, *fSeis_Bin;

    if(data_FileName!=NULL)
        initFilesSmart3d(data_FileName, &fSeis_Hdr, &fSeis_Bin, nt, dt, ot, shot[0].nrec, 1, 0, dataset->nshot/10, 1, 0);
    
    int OPTIM = 1;
    int mode_cigs;
    if(Hocig!=NULL  &&  Vocig!=NULL)    mode_cigs = MODE_CIGS_ORTH;
    else                                mode_cigs = MODE_CIGS_FULL;
    if(OPTIM)
    {
        GPU_OPER_initEnvironment(procGPU, nx, nz, nt, dx, dz, dt, ((dataset->shot)[0]).nsou, ((dataset->shot)[0]).nrec, lx0, lz0, 0, nxb, nzb, mode_cigs);
        GPU_OPER_alocArrays_acoustic();
        GPU_OPER_alocArrays_acoustic_ExtBornSct();
        if(mode_cigs == MODE_CIGS_FULL) {
            GPU_OPER_copyEImage_CPU2GPU(eimage);
        }
        else if(mode_cigs == MODE_CIGS_ORTH) {
            GPU_OPER_copyOCIGS_CPU2GPU(Hocig, Vocig);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(iproc==0)  printf("\n\n Starting extended Born modeling of %lld shots: \n  (Shot,time): \n", nshot);  fflush(stdout);
    for(ishot=iproc; ishot<nshot; ishot+=nproc)
    {
        // printf("\n  nx=%lld  nxb=%d  iproc=%d  ishot=%lld  ixs=%d  \n", nx, nxb, iproc, ishot, (dataset->shot)[ishot].ixs[0]);  fflush(stdout);
        // if(iproc!=0)    continue;
        modelShotBorn_lx_lz_lt_GPU(procGPU, eimage, Hocig, Vocig, vp_bg, &((dataset->shot)[ishot]), \
                                   jt, nt, nx, nz, dx, dz, dt, ox, oz, ot, lx0, lz0, lt0, 1, 0, OPTIM);

        // applyRecTaper(&((dataset->shot)[ishot]), perc, offMax, acquisitionGeom);
        
        time_now = time(NULL);
        if(iproc==0)    printf(" (%3lld ; %3ld) ", ishot, time_now-time_start);  fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(OPTIM)
    {
        GPU_OPER_freeArrays_acoustic_ExtBornSct();
        GPU_OPER_freeArrays_acoustic();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(iproc==0) {
        printf("\n\n Extended Born modeled shots concluded \n\n");
        printf("\n Elapsed time: %3ld \n", time(NULL)-time_start);  fflush(stdout);
    }
    for(ishot=iproc; ishot<nshot; ishot+=nproc)
    {        
        if(ishot%10==0  &&  data_FileName!=NULL)
            writeSamples((dataset->shot)[ishot].seismogram, fSeis_Bin, nt * (dataset->shot)[0].nrec);
    }
    // Sending seismograms to manager process
    for(ishot=0; ishot<nshot; ishot++)
    {
        long long int nt   = (dataset->shot)[ishot].nt;
        long long int nrec = (dataset->shot)[ishot].nrec;
        if((ishot%nproc)!=0)
        {
            if(iproc==ishot%nproc)
            {
                MPI_Send((dataset->shot)[ishot].seismogram, nt*nrec, MPI_FLOAT,      0, 0, MPI_COMM_WORLD);
            }
            if(iproc==0)
            {
                int sender = ishot%nproc;
                MPI_Recv((dataset->shot)[ishot].seismogram, nt*nrec, MPI_FLOAT, sender, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast((dataset->shot)[ishot].seismogram, nt*nrec, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if(iproc==0)
    {
        printf("\n\n Finishing transfering data around \n\n");
        printf("\n Elapsed time: %3ld \n", time(NULL)-time_start);  fflush(stdout);
    }

    objFunc = 0.0;
    for(ishot=0; ishot<nshot; ishot++)
    {
        objFunc += addObjFunc(&((dataset->shot)[ishot]));
    }

    MPI_Barrier(MPI_COMM_WORLD);

return objFunc;
}

void ExtendedBornModeling_TAU(Dataset3D *dataset, float *vp_bg, float *Tcig, MSH *msh, Eimage_MSH *eimsh, \
                              Parms *parms, int ism, int recordMovie, int iproc, int nproc) {
    
    int ret = 0;

    float         dx   = msh->dx;
    float         dz   = msh->dz;
    float         dt   = msh->dt;
    float         ox   = msh->ox;
    float         oz   = msh->oz;
    float         ot   = msh->ot;
    long long int nx   = msh->nx;
    long long int nz   = msh->nz;
    long long int nt   = msh->nt;
    long long int ntau = eimsh->ntau;
    long long int tau0 = eimsh->tau0;

    
    long long int ishot, nshot;
    nshot = dataset->nshot;

    FILE *fSeis_inp_Hdr, *fSeis_inp_Bin;
    FILE *fSeis_ebm_Hdr, *fSeis_ebm_Bin;

    time_t time_now, time_start;
    time_start = time(NULL);
    for(ishot=iproc; ishot<nshot; ishot+=nproc)
    {
        initShot3DfromFile(dataset, ishot);
        loadSeismogram2Dataset(dataset, ishot, parms->normalizeInputDataAmplitudes);
        Shot3D *shot = &(dataset->shot)[ishot];
        applyHalfDerivative2Source3D(shot);

        // Saving file with dimensions (nz, nx, ntau, nlx, nlz) of the dataset, 
        // to make it easier to read from a Jupyter notebook.
        if(iproc==0  &&  ishot==0)
        {
            char dimsFileName[2048];
            sprintf(dimsFileName, "%s/EBM_dims", parms->executionDirectory);
            printf("DEBUG IN EBM: dimsFileName = <%s>", dimsFileName);
            FILE *fp = fopen(dimsFileName,"w");
            float EBM_nshot = dataset->nshot/nproc;
            float EBM_nrec  = shot->nrec;
            float EBM_nt    = shot->nt;
            fwrite(&EBM_nt   , sizeof(float), 1, fp);
            fwrite(&EBM_nrec , sizeof(float), 1, fp);
            fwrite(&EBM_nshot, sizeof(float), 1, fp);        
            fclose(fp);
        }
    
        MSH msh_shot;
        float *vpShot;
        int processThisShot;
        if(getShotGrid(&msh_shot, vp_bg, msh, shot, parms))    processThisShot = 0;
        else                                                   processThisShot = 1;

        // Opening and saving files with input and EBM shots
        if(iproc==0  &&  ishot==0) initFilesSmart3d("INP_Seismogram", &fSeis_inp_Hdr, &fSeis_inp_Bin, shot->nt, shot->dt, shot->ot, shot->nrec, 1, 0, dataset->nshot/nproc, 1, 0);
        if(iproc==0  &&  ishot==0) initFilesSmart3d("EBM_Seismogram", &fSeis_ebm_Hdr, &fSeis_ebm_Bin, shot->nt, shot->dt, shot->ot, shot->nrec, 1, 0, dataset->nshot/nproc, 1, 0);

        if(processThisShot)
        {
            long long int nx  = msh_shot.nx;
            long long int nz  = msh_shot.nz;
            float         dx  = msh_shot.dx;
            float         dz  = msh_shot.dz;
            float         ox  = msh_shot.ox;
            float         oz  = msh_shot.oz;

            float     *Tcig_shot = CPU_zaloc1F(nx*nz*eimsh->ntau);
            getShotModel(&vpShot, &msh_shot, vp_bg, msh, shot, parms);
            getShotExtendedModel(Tcig, msh, Tcig_shot, &msh_shot, eimsh->tau0);

            Shot3D *shotEBM = initCopyShot3D(dataset, ishot);
            ExtendedBornModeling_TAU_Mem_GPU(procGPU, Tcig_shot, vpShot, shotEBM, &msh_shot, eimsh, ism, (recordMovie && ishot==ism));

            if(ishot%nproc==0)    writeSamples(shot->seismogram   , fSeis_inp_Bin, shot->nt    * shot->nrec); 
            if(ishot%nproc==0)    writeSamples(shotEBM->seismogram, fSeis_ebm_Bin, shotEBM->nt * shotEBM->nrec);

            finalizeShot3D(shotEBM);

            time_now = time(NULL);
            if(iproc==0)    printf(" (%3lld ; %3ld) ", ishot, time_now-time_start);  fflush(stdout);

            free(Tcig_shot);
        }

        ret = finalizeShot3DfromFile(dataset, ishot);
        if(ret!=0)    printf("\n CANNOT FINALIZE SHOT %lld \n", ishot);
        ret = finalizeSeismogramAndSource(dataset, ishot);
        if(ret!=0)    printf("\n CANNOT FINALIZE SHOT %lld \n", ishot);
    }
    if(iproc==0)
    {
        fclose(fSeis_inp_Hdr);
        fclose(fSeis_inp_Bin);
        fclose(fSeis_ebm_Hdr);
        fclose(fSeis_ebm_Bin);
    }

return;
}

void BornGradient_DataDomain_TAU(Dataset3D *dataset, float *gradient, float *ilumin, float *vp_bg, \
                                 float *Tcig, GRDProc *grdProc, MSH *msh, Eimage_MSH *eimsh, \
                                 Parms *parms, int ism, int recordMovie, int iter, int iproc, int nproc) {
    
    int ret = 0;

    float         dx   = msh->dx;
    float         dz   = msh->dz;
    float         dt   = msh->dt;
    float         ox   = msh->ox;
    float         oz   = msh->oz;
    float         ot   = msh->ot;
    long long int nx   = msh->nx;
    long long int nz   = msh->nz;
    long long int nt   = msh->nt;
    long long int ntau = eimsh->ntau;
    long long int tau0 = eimsh->tau0;
    
    float *gradient1      = CPU_zaloc1F(nx*nz);
    float *gradient2      = CPU_zaloc1F(nx*nz);
    float *ilumin1        = CPU_zaloc1F(nx*nz);
    float *ilumin2        = CPU_zaloc1F(nx*nz);
        
    char grad1ShotName[1024];
    char grad2ShotName[1024];

    long long int ishot, nshot;
    nshot = dataset->nshot;

    int mode_cigs = MODE_CIGS_TIME;

    MPI_Barrier(MPI_COMM_WORLD);

    time_t time_now, time_start;
    time_start = time(NULL);
    for(ishot=iproc; ishot<nshot; ishot+=nproc)
    {
        initShot3DfromFile(dataset, ishot);
        loadSeismogram2Dataset(dataset, ishot, parms->normalizeInputDataAmplitudes);
        Shot3D *shot = &(dataset->shot)[ishot];
        applyHalfDerivative2Source3D(shot);
        applyHalfDerivative2Shot(shot);

        MSH msh_shot;
        float *vpShot;
        int processThisShot;
        if(getShotGrid(&msh_shot, vp_bg, msh, shot, parms))    processThisShot = 0;
        else                                                   processThisShot = 1;

        if(processThisShot)
        {
            long long int nx  = msh_shot.nx;
            long long int nz  = msh_shot.nz;    
            float         dx  = msh_shot.dx;
            float         dz  = msh_shot.dz;
            float         ox  = msh_shot.ox;
            float         oz  = msh_shot.oz;

            float *gradient1_shot = CPU_zaloc1F(nx*nz);
            float *gradient2_shot = CPU_zaloc1F(nx*nz);
            float   *ilumin1_shot = CPU_zaloc1F(nx*nz);
            float   *ilumin2_shot = CPU_zaloc1F(nx*nz);
            float     *Tcig_shot  = CPU_zaloc1F(nx*nz*eimsh->ntau);

            getShotModel(&vpShot, &msh_shot, vp_bg, msh, shot, parms);

            getShotExtendedModel(Tcig, msh, Tcig_shot, &msh_shot, eimsh->tau0);

            //migrateExtBorn_lx_lz_lt_Mem_GPU
            migrateExtBorn_TAU_Mem_GPU(procGPU, gradient1_shot, gradient2_shot, ilumin1_shot, ilumin2_shot, \
                                       Tcig_shot, vpShot, shot, &msh_shot, eimsh, ism, (recordMovie && ishot==ism) );
            
            accumulateModel(gradient1, msh, gradient1_shot, &msh_shot);
            accumulateModel(gradient2, msh, gradient2_shot, &msh_shot);
            accumulateModel(  ilumin1, msh,   ilumin1_shot, &msh_shot);
            accumulateModel(  ilumin2, msh,   ilumin2_shot, &msh_shot);

            time_now = time(NULL);
            if(iproc==0)    printf(" (%3lld ; %3ld) ", ishot, time_now-time_start);  fflush(stdout);

            char gradient1ShotName[1024];    sprintf(gradient1ShotName, "gradient1_shot%lld", ishot);
            char gradient2ShotName[1024];    sprintf(gradient2ShotName, "gradient2_shot%lld", ishot);
            // outputSmart2d(gradient1ShotName, gradient1_shot, nz, dz, oz, nx, dx, ox);
            // outputSmart2d(gradient2ShotName, gradient2_shot, nz, dz, oz, nx, dx, ox);

            free(gradient1_shot);
            free(gradient2_shot);
            free(  ilumin1_shot);
            free(  ilumin2_shot);
            free(Tcig_shot);
        }
        ret = finalizeShot3DfromFile(dataset, ishot);
        if(ret!=0)    printf("\n CANNOT FINALIZE SHOT %lld \n", ishot);
        ret = finalizeSeismogramAndSource(dataset, ishot);
        if(ret!=0)    printf("\n CANNOT FINALIZE SHOT %lld \n", ishot);
    }

    stackLargeArray_MPI(     Tcig, iproc, 0, nx*nz*ntau, MPI_FLOAT, MPI_COMM_WORLD);
    stackLargeArray_MPI(gradient1, iproc, 0, nx*nz     , MPI_FLOAT, MPI_COMM_WORLD);
    stackLargeArray_MPI(gradient2, iproc, 0, nx*nz     , MPI_FLOAT, MPI_COMM_WORLD);

    normalizeArray(gradient1, nx*nz);
    normalizeArray(gradient2, nx*nz);
    // addArrays(gradient, 1.0, gradient1, 1.0, gradient2, nx*nz);
    
    if(grdProc->applySmoothing)
    {
        int rx = grdProc->halfLengthX/dx;
        int rz = grdProc->halfLengthZ/dz;
        modelSmooth(gradient, gradient, nx, nz, rx, rz);
    }
    if(grdProc->applyXTaper)
        modelTaperingX(gradient, gradient, nx, nz, dx, grdProc->x1, grdProc->x2, grdProc->x3, grdProc->x4);

    if(grdProc->applyZTaper) {
        modelTaperingZ(gradient, gradient, nx, nz, dz, grdProc->z1, grdProc->z2);
        modelTaperingZ(gradient, gradient, nx, nz, dz, grdProc->z4, grdProc->z3);
    }

    takeSquareRootArrays(ilumin1, ilumin1, nx*nz);
    takeSquareRootArrays(ilumin2, ilumin2, nx*nz);
    // int rx = 200/dx;
    // int rz = 50/dz;
    int rx = 5000/dx;
    int rz = 100/dz;
    modelSmooth(ilumin1, ilumin1, nx, nz, rx, rz);
    modelSmooth(ilumin2, ilumin2, nx, nz, rx, rz);
    addArrays(ilumin, 1.0, ilumin1, 1.0, ilumin2, nx*nz);
    normalizeArray(ilumin , nx*nz);
    // normalizeArray(ilumin1, nx*nz);
    // normalizeArray(ilumin2, nx*nz);
    float epsilon1 = 0.1;
    float epsilon2 = 0.1;
    // applyIlum2EImage(NULL, NULL, gradient, gradient, ilumin, nx, nz, dx, tau0, epsilon1);
    // normalizeArray(gradient, nx*nz);

    applyIlum2EImage(NULL, NULL, gradient1, gradient1, ilumin1, nx, nz, dx, tau0, epsilon1);
    normalizeArray(gradient1, nx*nz);
    applyIlum2EImage(NULL, NULL, gradient2, gradient2, ilumin2, nx, nz, dx, tau0, epsilon2);
    normalizeArray(gradient2, nx*nz);
    addArrays(gradient, 1.0, gradient1, 1.0, gradient2, nx*nz);
    normalizeArray(gradient, nx*nz);


    char gradient1Name[1024];
    char gradient2Name[1024];
    char ilumin1Name[1024];
    char ilumin2Name[1024];
    char iluminName[1024];
    sprintf(gradient1Name, "gradient1_iter%d", iter);
    sprintf(gradient2Name, "gradient2_iter%d", iter);
    sprintf(ilumin1Name,     "ilumin1_iter%d", iter);
    sprintf(ilumin2Name,     "ilumin2_iter%d", iter);
    sprintf(iluminName,       "ilumin_iter%d", iter);
    outputSmart2d(gradient1Name, gradient1, nz, dz, oz, nx, dx, ox);
    outputSmart2d(gradient2Name, gradient2, nz, dz, oz, nx, dx, ox);
    outputSmart2d(ilumin1Name,     ilumin1, nz, dz, oz, nx, dx, ox);
    outputSmart2d(ilumin2Name,     ilumin2, nz, dz, oz, nx, dx, ox);
    outputSmart2d(iluminName ,      ilumin, nz, dz, oz, nx, dx, ox);

    free(gradient1); 
    free(gradient2);
    free(ilumin1);
    free(ilumin2);

return;
}

void BornGradient_DataDomain(Dataset3D *dataset, float *gradient, float *vp_bg, \
                             float *eimage, float *Heref, float *Veref, GRDProc *grdProc, MSH *msh, Eimage_MSH *eimsh, \
                             Parms *parms, int ism, int recordMovie, int iter, int iproc, int nproc) {
    
    int ret = 0;

    float         dx  = msh->dx;
    float         dz  = msh->dz;
    float         dt  = msh->dt;
    float         ox  = msh->ox;
    float         oz  = msh->oz;
    float         ot  = msh->ot;
    long long int nx  = msh->nx;
    long long int nz  = msh->nz;
    long long int nt  = msh->nt;
    long long int nlx = eimsh->nlx;
    long long int nlz = eimsh->nlz;
    long long int lx0 = eimsh->lx0;
    long long int lz0 = eimsh->lz0;
    
    float *gradient1      = CPU_zaloc1F(nx*nz);
    float *gradient2      = CPU_zaloc1F(nx*nz);
    float *ilumin         = CPU_zaloc1F(nx*nz);
    float *ilumin1        = CPU_zaloc1F(nx*nz);
    float *ilumin2        = CPU_zaloc1F(nx*nz);
    
    // long long int   nt = ((dataset->shot)[0]).nt;
    // float dt = ((dataset->shot)[0]).dt;
    // float ot = ((dataset->shot)[0]).ot;
    
    char grad1ShotName[1024];
    char grad2ShotName[1024];


    long long int ishot, nshot;
    nshot = dataset->nshot;

    int mode_cigs;
    if(Heref!=NULL  &&  Veref!=NULL)    mode_cigs = MODE_CIGS_ORTH;
    else                                mode_cigs = MODE_CIGS_FULL;

    time_t time_now, time_start;
    time_start = time(NULL);
    for(ishot=iproc; ishot<nshot; ishot+=nproc)
    {
        initShot3DfromFile(dataset, ishot);
        loadSeismogram2Dataset(dataset, ishot, parms->normalizeInputDataAmplitudes);
        Shot3D *shot = &(dataset->shot)[ishot];
        applyHalfDerivative2Source3D(shot);
        applyHalfDerivative2Shot(shot);

        MSH msh_shot;
        float *vpShot;
        int processThisShot;
        if(getShotGrid(&msh_shot, vp_bg, msh, shot, parms))    processThisShot = 0;
        else                                                   processThisShot = 1;

        if(processThisShot)
        {
            long long int nx  = msh_shot.nx;
            long long int nz  = msh_shot.nz;
            float         dx  = msh_shot.dx;
            float         dz  = msh_shot.dz;
            float         ox  = msh_shot.ox;
            float         oz  = msh_shot.oz;

            float *gradient1_shot = CPU_zaloc1F(nx*nz);
            float *gradient2_shot = CPU_zaloc1F(nx*nz);
            float   *ilumin1_shot = CPU_zaloc1F(nx*nz);
            float   *ilumin2_shot = CPU_zaloc1F(nx*nz);
            float    *eimage_shot  = CPU_zaloc1F(nx*nz*eimsh->nlx);
            float     *Heref_shot = CPU_zaloc1F(nx*nz*eimsh->nlx);
            float     *Veref_shot = CPU_zaloc1F(nx*nz*eimsh->nlz);

            getShotModel(&vpShot, &msh_shot, vp_bg, msh, shot, parms);

            getShotExtendedModel( Heref, msh,  Heref_shot, &msh_shot, eimsh->lx0);
            getShotExtendedModel( Veref, msh,  Veref_shot, &msh_shot, eimsh->lz0);
            if(eimage!=NULL)
                getShotExtendedModel(eimage, msh, eimage_shot, &msh_shot, eimsh->lx0);

            migrateExtBorn_lx_lz_lt_Mem_GPU(procGPU, gradient1_shot, gradient2_shot, ilumin1_shot, ilumin2_shot, \
                                            eimage_shot, Heref_shot, Veref_shot, vpShot, shot, &msh_shot, eimsh, ism, (recordMovie && ishot==ism) );
            
            accumulateModel(gradient1, msh, gradient1_shot, &msh_shot);
            accumulateModel(gradient2, msh, gradient2_shot, &msh_shot);
            accumulateModel(  ilumin1, msh,   ilumin1_shot, &msh_shot);
            accumulateModel(  ilumin2, msh,   ilumin2_shot, &msh_shot);

            time_now = time(NULL);
            if(iproc==0)    printf(" (%3lld ; %3ld) ", ishot, time_now-time_start);  fflush(stdout);

            char gradient1ShotName[1024];    sprintf(gradient1ShotName, "gradient1_shot%lld", ishot);
            char gradient2ShotName[1024];    sprintf(gradient2ShotName, "gradient2_shot%lld", ishot);
            outputSmart2d(gradient1ShotName, gradient1_shot, nz, dz, oz, nx, dx, ox);
            outputSmart2d(gradient2ShotName, gradient2_shot, nz, dz, oz, nx, dx, ox);

            free(gradient1_shot);
            free(gradient2_shot);
            free(  ilumin1_shot);
            free(  ilumin2_shot);

            free(eimage_shot);
            free(Heref_shot);
            free(Veref_shot);
        }
        ret = finalizeShot3DfromFile(dataset, ishot);
        if(ret!=0)    printf("\n CANNOT FINALIZE SHOT %lld \n", ishot);
        ret = finalizeSeismogramAndSource(dataset, ishot);
        if(ret!=0)    printf("\n CANNOT FINALIZE SHOT %lld \n", ishot);
    }

    stackLargeArray_MPI(Heref    , iproc, 0, nx*nz*nlx, MPI_FLOAT, MPI_COMM_WORLD);
    stackLargeArray_MPI(Veref    , iproc, 0, nx*nz*nlz, MPI_FLOAT, MPI_COMM_WORLD);
    stackLargeArray_MPI(gradient1, iproc, 0, nx*nz    , MPI_FLOAT, MPI_COMM_WORLD);
    stackLargeArray_MPI(gradient2, iproc, 0, nx*nz    , MPI_FLOAT, MPI_COMM_WORLD);

    addArrays(gradient, 1.0, gradient1, 1.0, gradient2, nx*nz);
    normalizeArray(gradient1, nx*nz);
    normalizeArray(gradient2, nx*nz);
    
    if(grdProc->applySmoothing)
    {
        int rx = grdProc->halfLengthX/dx;
        int rz = grdProc->halfLengthZ/dz;
        modelSmooth(gradient, gradient, nx, nz, rx, rz);
    }
    if(grdProc->applyXTaper)
        modelTaperingX(gradient, gradient, nx, nz, dx, grdProc->x1, grdProc->x2, grdProc->x3, grdProc->x4);

    if(grdProc->applyZTaper) {
        modelTaperingZ(gradient, gradient, nx, nz, dz, grdProc->z1, grdProc->z2);
        modelTaperingZ(gradient, gradient, nx, nz, dz, grdProc->z4, grdProc->z3);
    }

    int rx = 5000/dx;
    int rz = 100/dz;
    modelSmooth(ilumin1, ilumin1, nx, nz, rx, rz);
    modelSmooth(ilumin2, ilumin2, nx, nz, rx, rz);
    ilumin = CPU_zaloc1F(nx*nz);
    addArrays(ilumin, 1.0, ilumin1, 1.0, ilumin2, nx*nz);
    normalizeArray(ilumin, nx*nz);
    normalizeArray(ilumin1, nx*nz);
    normalizeArray(ilumin2, nx*nz);
    float epsilon = 0.002;
    applyIlum2EImage(NULL, NULL, gradient, gradient, ilumin, nx, nz, dx, lx0, epsilon);
    normalizeArray(gradient, nx*nz);

    char gradient1Name[1024];
    char gradient2Name[1024];
    char ilumin1Name[1024];
    char ilumin2Name[1024];
    char iluminName[1024];
    sprintf(gradient1Name, "gradient1_iter%d", iter);
    sprintf(gradient2Name, "gradient2_iter%d", iter);
    sprintf(ilumin1Name,     "ilumin1_iter%d", iter);
    sprintf(ilumin2Name,     "ilumin2_iter%d", iter);
    sprintf(iluminName,       "ilumin_iter%d", iter);
    outputSmart2d(gradient1Name, gradient1, nz, dz, oz, nx, dx, ox);
    outputSmart2d(gradient2Name, gradient2, nz, dz, oz, nx, dx, ox);
    outputSmart2d(ilumin1Name,   ilumin1, nz, dz, oz, nx, dx, ox);
    outputSmart2d(ilumin2Name,   ilumin2, nz, dz, oz, nx, dx, ox);
    outputSmart2d(iluminName ,   ilumin, nz, dz, oz, nx, dx, ox);

    free(gradient1); 
    free(gradient2);
    free(ilumin1);
    free(ilumin2);
    free(ilumin);

return;
}
/*
void BornGradient_DataDomain_backup(Dataset3D *dataset, float *gradient, float *vp_bg, \
                             float *eimage, float *Heref, float *Veref, GRDProc *grdProc, \
                             long long int nx, long long int nz, \
                             long long int lx0, long long int lz0, long long int lt0, \
                             float dx, float dz, float ox, float oz, \
                             long long int jt, long long int jt_EIC, int iter, int imageDomain, \
                             int iproc, int nproc) {
    
    long long int   nt = ((dataset->shot)[0]).nt;
    float dt = ((dataset->shot)[0]).dt;
    float ot = ((dataset->shot)[0]).ot;

    float *gradient1_shot = CPU_zaloc1F(nx*nz);
    float *gradient2_shot = CPU_zaloc1F(nx*nz);
    float *gradient1 = CPU_zaloc1F(nx*nz);
    float *gradient2 = CPU_zaloc1F(nx*nz);
    float *ilumin    = CPU_zaloc1F(nx*nz);
    float *ilumin1   = CPU_zaloc1F(nx*nz);
    float *ilumin2   = CPU_zaloc1F(nx*nz);
    
    char grad1ShotName[1024];
    char grad2ShotName[1024];

    long long int nlx = 2*lx0 + 1;
    long long int nlz = 2*lz0 + 1;

    if(eimage!=NULL) {
        normalizeArray(eimage, nx*nz*nlx*nlz);
    }
    else if (Heref!=NULL  &&  Veref!=NULL) {
        // normalizeArray(Heref, nx*nz*nlx);
        // normalizeArray(Veref, nx*nz*nlz);
        // printf("\n\n Normalizing Heref and Veref \n\n");
        // normalizeArraysJointly(Heref, nx*nz*nlx, Veref, nx*nz*nlz);
    }

    long long int ishot, nshot;
    nshot = dataset->nshot;

    int OPTIM = 1;
    int mode_cigs;
    if(Heref!=NULL  &&  Veref!=NULL)    mode_cigs = MODE_CIGS_ORTH;
    else                                mode_cigs = MODE_CIGS_FULL;
    if(OPTIM)
    {
        GPU_OPER_initEnvironment(procGPU, nx, nz, nt, dx, dz, dt, ((dataset->shot)[0]).nsou, ((dataset->shot)[0]).nrec, lx0, lz0, 0, nxb, nzb, mode_cigs);
        GPU_OPER_alocArrays_acoustic();
        GPU_OPER_alocArrays_acoustic_ExtBornSct();
        if(mode_cigs == MODE_CIGS_FULL) {
            GPU_OPER_copyEImage_CPU2GPU(eimage);
        }
        else if(mode_cigs == MODE_CIGS_ORTH) {
            GPU_OPER_copyOCIGS_CPU2GPU(Heref, Veref);
        }
    }
    time_t time_now, time_start;
    time_start = time(NULL);
    for(ishot=iproc; ishot<nshot; ishot+=nproc)
    {
        long long int recMv=0, shtMv=nshot/2;
        if( ishot<0 )
        {
            recMv = 1;
            shtMv = ishot;
        }
        else    recMv = 1;

        // if(ishot==20)
        migrateExtBorn_lx_lz_lt_Mem_GPU_OPTIM(procGPU, gradient1, gradient2, gradient1_shot, gradient2_shot, ilumin1, ilumin2, eimage, Heref, Veref, \
                                              vp_bg, &((dataset->shot)[ishot]), jt, jt_EIC, nt, nx, nz, dx, dz, dt, ox, oz, ot, lx0, lz0, lt0, \
                                              shtMv, recMv, imageDomain);
        time_now = time(NULL);
        if(iproc==0)    printf(" (%3lld ; %3ld) ", ishot, time_now-time_start);  fflush(stdout);

        char gradient1ShotName[1024];    sprintf(gradient1ShotName, "gradient1_shot%lld", ishot);
        char gradient2ShotName[1024];    sprintf(gradient2ShotName, "gradient2_shot%lld", ishot);
        outputSmart2d(gradient1ShotName, gradient1, nz, dz, oz, nx, dx, ox);
        outputSmart2d(gradient2ShotName, gradient2, nz, dz, oz, nx, dx, ox);        
    }
    if(OPTIM)
    {
        GPU_OPER_copyImage_BRNSCT_GPU2CPU(gradient1, gradient2, 1);
        GPU_OPER_copyIlumin_BRNSCT_GPU2CPU(ilumin1, ilumin2, 1);
        GPU_OPER_freeArrays_acoustic_ExtBornSct();
        GPU_OPER_freeArrays_acoustic();
    }
    if(iproc==0)
    {
        outputSmart3d("Heref", Heref, nlx, dx, -lx0*dx,  nz, dz, oz, nx, dx, ox);
        outputSmart3d("Veref", Veref, nlz, dz, -lz0*dz,  nz, dz, oz, nx, dx, ox);
        outputSmart2d("grad1_iproc0", gradient1, nz, dz, oz, nx, dx, ox);
        outputSmart2d("grad2_iproc0", gradient2, nz, dz, oz, nx, dx, ox);
    }
    float *gradient1_tmp;
    float *gradient2_tmp;
    float *ilumin1_tmp;
    float *ilumin2_tmp;
    if(iproc==0)
    {
        gradient1_tmp  = CPU_zaloc1F(nx*nz);
        gradient2_tmp  = CPU_zaloc1F(nx*nz);
          ilumin1_tmp  = CPU_zaloc1F(nx*nz);
          ilumin2_tmp  = CPU_zaloc1F(nx*nz);
    }
    MPI_Reduce(gradient1, gradient1_tmp, nx*nz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(gradient2, gradient2_tmp, nx*nz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(  ilumin1,   ilumin1_tmp, nx*nz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(  ilumin2,   ilumin2_tmp, nx*nz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(iproc==0)
    {
        memcpy(gradient1, gradient1_tmp, nx*nz*sizeof(float));
        memcpy(gradient2, gradient2_tmp, nx*nz*sizeof(float));
        memcpy(  ilumin1, ilumin1_tmp, nx*nz*sizeof(float));
        memcpy(  ilumin2, ilumin2_tmp, nx*nz*sizeof(float));
        free(gradient1_tmp);
        free(gradient2_tmp);
        free(  ilumin1_tmp);
        free(  ilumin2_tmp);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(gradient1, nx*nz, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(gradient2, nx*nz, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(  ilumin1, nx*nz, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(  ilumin2, nx*nz, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);


    // this is used when the model parameter is velocity itself, instead of square slowness
    // untransformVel(vp_bg, nx, nz, dx, dt);
    // divideArraysPow(gradient1, 1.0, gradient1, 1.0, vp_bg, nx*nz);
    // divideArraysPow(gradient2, 1.0, gradient2, 1.0, vp_bg, nx*nz);
    // transformVel(vp_bg, nx, nz, dx, dt);

    // multiplyArrayByScalar(-1.0, gradient1, nx*nz);
    // multiplyArrayByScalar(-1.0, gradient2, nx*nz);

    addArrays(gradient, 1.0, gradient1, 1.0, gradient2, nx*nz);
    normalizeArray(gradient1, nx*nz);
    normalizeArray(gradient2, nx*nz);
    
    if(grdProc->applySmoothing)
    {
        int rx = grdProc->halfLengthX/dx;
        int rz = grdProc->halfLengthX/dz;
        modelSmooth(gradient, gradient, nx, nz, rx, rz);
    }
    if(grdProc->applyXTaper)
        modelTaperingX(gradient, gradient, nx, nz, dx, grdProc->x1, grdProc->x2, grdProc->x3, grdProc->x4);

    if(grdProc->applyZTaper) {
        modelTaperingZ(gradient, gradient, nx, nz, dz, grdProc->z1, grdProc->z2);
        modelTaperingZ(gradient, gradient, nx, nz, dz, grdProc->z4, grdProc->z3);
    }
        

    int rx = 5000/dx;
    int rz = 100/dz;
    modelSmooth(ilumin1, ilumin1, nx, nz, rx, rz);
    modelSmooth(ilumin2, ilumin2, nx, nz, rx, rz);
    ilumin = CPU_zaloc1F(nx*nz);
    addArrays(ilumin, 1.0, ilumin1, 1.0, ilumin2, nx*nz);
    normalizeArray(ilumin, nx*nz);
    normalizeArray(ilumin1, nx*nz);
    normalizeArray(ilumin2, nx*nz);
    float epsilon = 0.002;
    applyIlum2EImage(NULL, NULL, gradient, gradient, ilumin, nx, nz, dx, lx0, epsilon);
    normalizeArray(gradient, nx*nz);

    char gradient1Name[1024];
    char gradient2Name[1024];
    char ilumin1Name[1024];
    char ilumin2Name[1024];
    char iluminName[1024];
    sprintf(gradient1Name, "gradient1_iter%d", iter);
    sprintf(gradient2Name, "gradient2_iter%d", iter);
    sprintf(ilumin1Name,     "ilumin1_iter%d", iter);
    sprintf(ilumin2Name,     "ilumin2_iter%d", iter);
    sprintf(iluminName,       "ilumin_iter%d", iter);
    outputSmart2d(gradient1Name, gradient1, nz, dz, oz, nx, dx, ox);
    outputSmart2d(gradient2Name, gradient2, nz, dz, oz, nx, dx, ox);
    outputSmart2d(ilumin1Name,   ilumin1, nz, dz, oz, nx, dx, ox);
    outputSmart2d(ilumin2Name,   ilumin2, nz, dz, oz, nx, dx, ox);
    outputSmart2d(iluminName ,   ilumin, nz, dz, oz, nx, dx, ox);

    free(gradient1); 
    free(gradient2);
    free(gradient1_shot);
    free(gradient2_shot);
    free(ilumin1);
    free(ilumin2);
    free(ilumin);

return;
}
*/



// Accumulating model
int accumulateModel(float *model_out, MSH *msh_out, float *model_inp, MSH *msh_inp)
{
    float ox_inp = msh_inp->ox + msh_inp->nxb * msh_inp->dx;
    float ox_out = msh_out->ox + msh_out->nxb * msh_out->dx;
    float dx = msh_inp->dx;
    int nx_  = msh_inp->nx - 2 * msh_inp->nxb;
    int nz_  = msh_inp->nz - 2 * msh_inp->nzb;
    int nxb_inp = msh_inp->nxb;
    int nzb_inp = msh_inp->nzb;
    int nxb_out = msh_out->nxb;
    int nzb_out = msh_out->nzb;

    long long int ix0 = (int) ((ox_inp - ox_out) / dx);
    
    // printf("\n nxb_inp=%d  nzb_inp=%d  nxb_out=%d  nzb_out=%d  nx_=%d  nz_=%d  nx_out=%d  nz_out=%d  ix0=%lld\n", \
           nxb_inp, nzb_inp, nxb_out, nzb_out, nx_, nz_, msh_out->nx, msh_out->nz, ix0);
    
    long long int ix, iz;
    for(ix=0; ix<nx_; ix++)
        for(iz=0; iz<nz_; iz++) {
            long long int i_inp = (iz+nzb_inp) + (ix+nxb_inp)     * msh_inp->nz;
            long long int i_out = (iz+nzb_out) + (ix+nxb_out+ix0) * msh_out->nz;
            int cond_out = i_out < msh_out->nx*msh_out->nz  &&  i_out >= 0;
            int cond_inp = i_inp < msh_inp->nx*msh_inp->nz  &&  i_inp >= 0;
            if( cond_out &&  cond_inp )
                model_out[i_out] += model_inp[i_inp];
        }
return 0;
}
int accumulateExtendedModel(float *model_out, MSH *msh_out, float *model_inp, MSH *msh_inp, long long int lx0)
{
    long long int nlx = 2*lx0 + 1;

    float ox_inp = msh_inp->ox + msh_inp->nxb * msh_inp->dx;
    float ox_out = msh_out->ox + msh_out->nxb * msh_out->dx;
    float dx = msh_inp->dx;
    int nx_  = msh_inp->nx - 2 * msh_inp->nxb;
    int nz_  = msh_inp->nz - 2 * msh_inp->nzb;
    int nxb_inp = msh_inp->nxb;
    int nzb_inp = msh_inp->nzb;
    int nxb_out = msh_out->nxb;
    int nzb_out = msh_out->nzb;

    long long int ix0 = (int) ((ox_inp - ox_out) / dx);

    // printf("\n nxb_inp=%d  nzb_inp=%d  nxb_out=%d  nzb_out=%d  nx_=%d  nz_=%d  nx_out=%d  nz_out=%d  ix0=%lld\n", \
           nxb_inp, nzb_inp, nxb_out, nzb_out, nx_, nz_, msh_out->nx, msh_out->nz, ix0);
    
    long long int ix, iz, lx;
    for(ix=0; ix<nx_; ix++)
        for(iz=0; iz<nz_; iz++)
            for(lx=0; lx<nlx; lx++)
            {
                long long int i_inp = lx + (iz+nzb_inp) * nlx + (ix+nxb_inp)     * msh_inp->nz * nlx;
                long long int i_out = lx + (iz+nzb_out) * nlx + (ix+nxb_out+ix0) * msh_out->nz * nlx;
                
                int cond_out = i_out < msh_out->nx*msh_out->nz*nlx  &&  i_out >= 0;
                int cond_inp = i_inp < msh_inp->nx*msh_inp->nz*nlx  &&  i_inp >= 0;
                if( cond_out &&  cond_inp )
                    model_out[i_out] += model_inp[i_inp];
        }
return 0;
}
int getShotExtendedModel(float *model_ref, MSH *msh_ref, float *model_sht, MSH *msh_sht, long long int lx0)
{
    long long int nlx = 2*lx0 + 1;

    float ox_sht = msh_sht->ox + msh_sht->nxb * msh_sht->dx;
    float ox_ref = msh_ref->ox + msh_ref->nxb * msh_ref->dx;
    float dx     = msh_sht->dx;
    int nx_      = msh_sht->nx - 2 * msh_sht->nxb;
    int nz_      = msh_sht->nz - 2 * msh_sht->nzb;
    int nxb_sht  = msh_sht->nxb;
    int nzb_sht  = msh_sht->nzb;
    int nxb_ref  = msh_ref->nxb;
    int nzb_ref  = msh_ref->nzb;

    long long int ix0 = (int) ((ox_sht - ox_ref) / dx);
    
    long long int ix, iz, lx;
    for(ix=0; ix<nx_; ix++)
        for(iz=0; iz<nz_; iz++)
            for(lx=0; lx<nlx; lx++)
            {
                long long int i_sht = lx + (iz+nzb_sht)*nlx + (ix+nxb_sht)     * msh_sht->nz * nlx;
                long long int i_ref = lx + (iz+nzb_ref)*nlx + (ix+nxb_ref+ix0) * msh_ref->nz * nlx;
                
                int cond_ref = i_ref < msh_ref->nx*msh_ref->nz*nlx  &&  i_ref >= 0;
                int cond_sht = i_sht < msh_sht->nx*msh_sht->nz*nlx  &&  i_sht >= 0;
                if( cond_ref &&  cond_sht )
                    model_sht[i_sht] = model_ref[i_ref];
            }
return 0;
}


// Finding shot grid and model
int getShotGrid(MSH *msh_out, float *model_inp, MSH *msh_inp, Shot3D *shot, Parms *parms)
{
    int ret;

    Frame shotFrame;
    ret = getShotFrame(&shotFrame, msh_out, msh_inp, shot, parms->AppertureExtension_x);
    if(ret)    return 1;

    float vmax, vmin;
    findMinMaxWithinFrame(&vmin, &vmax, &shotFrame, msh_inp, model_inp);
    float spaceSampling = msh_inp->dx;
    // float spaceSampling = getSpaceSampling(shot->maxFrequency, vmin);

    int strideX=32, strideY=1, strideZ=1;
    interpShotGrid_GPUCompliant(msh_out, msh_out, parms, strideX, strideY, strideZ, spaceSampling, vmax);

return 0;
}

int getShotModel(float **model_out, MSH *msh_out, float *model_inp, MSH *msh_inp, Shot3D *shot, Parms *parms)
{
    float scalar = 1.0f;
    *model_out = extendResampInputModel(model_inp, msh_inp->nx, msh_inp->nz, \
                                        msh_inp->dx, msh_inp->dz, \
                                        msh_inp->ox, msh_inp->oz, \
                                        msh_out->nx, msh_out->nz, \
                                        msh_out->dx, msh_out->dz, \
                                        msh_out->ox+msh_out->nxb*msh_out->dx, \
                                        msh_out->oz+msh_out->nzb*msh_out->dz, \
                                        msh_out->nxb, msh_out->nzb, scalar);
    int ret = gridShot3D(msh_out, shot);
    if(ret)
    {
        printf("\n Function gridShot3D returned %d, which means \"unable to grid this shot\". Function [getShotModel] will return. \n", ret);
        return 1;
    }

return 0;
}

int getShotFrame(Frame *shotFrame, MSH *mshShot, MSH *mshRef, Shot3D *shot, float AppertureExt_x)
{
    int ir, is;

    if(gridShot3D(mshRef, shot))    return 1;

    shotFrame->ixMin = shot->ixr[0];
    shotFrame->ixMax = shot->ixr[0];
    shotFrame->iyMin = shot->iyr[0];
    shotFrame->iyMax = shot->iyr[0];

    for(ir=1; ir<shot->nrec; ir++)
    {
        if(shotFrame->ixMin > shot->ixr[ir])    shotFrame->ixMin = shot->ixr[ir];
        if(shotFrame->ixMax < shot->ixr[ir])    shotFrame->ixMax = shot->ixr[ir];

        if(shotFrame->iyMin > shot->iyr[ir])    shotFrame->iyMin = shot->iyr[ir];
        if(shotFrame->iyMax < shot->iyr[ir])    shotFrame->iyMax = shot->iyr[ir];
    }

    for(is=0; is<shot->nsou; is++)
    {
        if(shotFrame->ixMin > shot->ixs[is])    shotFrame->ixMin = shot->ixs[is];
        if(shotFrame->ixMax < shot->ixs[is])    shotFrame->ixMax = shot->ixs[is];

        if(shotFrame->iyMin > shot->iys[is])    shotFrame->iyMin = shot->iys[is];
        if(shotFrame->iyMax < shot->iys[is])    shotFrame->iyMax = shot->iys[is];
    }
    
    shotFrame->ixMin = max(shotFrame->ixMin - (int) floorf(AppertureExt_x/mshRef->dx),            0);
    shotFrame->ixMax = min(shotFrame->ixMax + (int) floorf(AppertureExt_x/mshRef->dx), mshRef->nx-1);

    shotFrame->ixMin = max(shotFrame->ixMin-5,            0);
    shotFrame->ixMax = min(shotFrame->ixMax+5, mshRef->nx-1);

    mshShot->dx = mshRef->dx;
    mshShot->dy = mshRef->dy;
    mshShot->dz = mshRef->dz;
    mshShot->dt = mshRef->dt;

    // P1 coords
    shotFrame->x1 = shotFrame->ixMin * mshRef->dx + mshRef->ox;
    shotFrame->y1 = shotFrame->iyMin * mshRef->dy + mshRef->oy;

    // P2 coords
    shotFrame->x2 = shotFrame->ixMax * mshRef->dx + mshRef->ox;
    shotFrame->y2 = shotFrame->iyMin * mshRef->dy + mshRef->oy;

    // P3 coords
    shotFrame->x3 = shotFrame->ixMin * mshRef->dx + mshRef->ox;
    shotFrame->y3 = shotFrame->iyMax * mshRef->dy + mshRef->oy;

    // P4 coords
    shotFrame->x4 = shotFrame->ixMax * mshRef->dx + mshRef->ox;
    shotFrame->y4 = shotFrame->iyMax * mshRef->dy + mshRef->oy;

    mshShot->nx = shotFrame->ixMax - shotFrame->ixMin + 1;
    mshShot->ny = shotFrame->iyMax - shotFrame->iyMin + 1;
    mshShot->nz = mshRef->nz;
    mshShot->nt = mshRef->nt;

    mshShot->nxb = 0;
    mshShot->nyb = 0;
    mshShot->nzb = 0;

    mshShot->ox = shotFrame->x1;
    mshShot->oy = shotFrame->y1;
    mshShot->oz = mshRef->oz;
    mshShot->ot = mshRef->ot;

    mshShot->LX = (mshShot->nx - 1) * mshShot->dx;
    mshShot->LY = (mshShot->ny - 1) * mshShot->dy;
    mshShot->LZ = mshRef->LZ;
    mshShot->LT = mshRef->LT;

    mshShot->jt     = mshRef->jt;
    mshShot->jt_lt  = mshRef->jt_lt;
    mshShot->jt_EIC = mshRef->jt_EIC;

    /*
    printf("\n AppertureExt_x=%f  \n", AppertureExt_x);
    printf("\n mshRef: nx=%d ny=%d nz=%d dx=%f dy=%f dz=%f ox=%f oy=%f oz=%f LX=%f LY=%f LZ=%f  jt=%d  jt_lt=%d  jt_EIC=%d \n", \
            mshRef->nx, mshRef->ny, mshRef->nz, mshRef->dx, mshRef->dy, mshRef->dz, mshRef->ox, mshRef->oy, mshRef->oz, \
            mshRef->LX, mshRef->LY, mshRef->LZ, mshRef->jt, mshRef->jt_lt, mshRef->jt_EIC);
    printf("\n mshShot: nx=%d ny=%d nz=%d dx=%f dy=%f dz=%f ox=%f oy=%f oz=%f LX=%f LY=%f LZ=%f  jt=%d  jt_lt=%d  jt_EIC=%d \n", \
            mshShot->nx, mshShot->ny, mshShot->nz, mshShot->dx, mshShot->dy, mshShot->dz, \
            mshShot->ox, mshShot->oy, mshShot->oz, mshShot->LX, mshShot->LY, mshShot->LZ, mshShot->jt, mshShot->jt_lt, mshShot->jt_EIC);
    printf("\n shotFrame: ixMin=%d  ixMax=%d   iyMin=%d  iyMax=%d\n", \
            shotFrame->ixMin, shotFrame->ixMax, shotFrame->iyMin, shotFrame->iyMax);
    //*/

return 0;
}

void findMinMaxWithinFrame(float *min, float *max, Frame *frame, MSH *msh, float *array)
{
    int ix, iy, iz;

    int iz0 = 0;
    int ix0 = frame->ixMin;
    int iy0 = frame->iyMin;
    int i0 = iz0 + ix0*msh->nz + iy0*msh->nz*msh->nx;

    // printf("\n  ixMin=%d  ixMax=%d    iyMin=%d  iyMax=%d\n", frame->ixMin, frame->ixMax, frame->iyMin, frame->iyMax); fflush(stdout);
    // printf("\n  nx=%d ny=%d nz=%d  \n", msh->nx, msh->ny, msh->nz); fflush(stdout);

    *min = array[i0];
    *max = array[i0];

    for(iy=frame->iyMin; iy<=frame->iyMax; iy++)
        for(ix=frame->ixMin; ix<=frame->ixMax; ix++)
            for(iz=0; iz<msh->nz; iz++)
            {
                int i = iz + ix*msh->nz + iy*msh->nz*msh->nx;
                if(*min>array[i])    *min = array[i];
                if(*max<array[i])    *max = array[i];
            }
return;
}
void interpShotGrid_GPUCompliant(MSH *ou, MSH *in, Parms *parms, int strideX, int strideY, int strideZ, float spaceSampling, float vmax)
{
    ou->dx = spaceSampling;
    ou->dy = spaceSampling;
    ou->dz = spaceSampling;
    // ou->dt = getStableTimeSampling(spaceSampling, vmax);

    ou->ox = in->ox;
    ou->oy = in->oy;
    ou->oz = in->oz;
    ou->ot = in->ot;

    ou->nxb = parms->nxb;
    ou->nyb = parms->nyb;
    ou->nzb = parms->nzb;

    ou->LX = in->LX;
    ou->LY = in->LY;
    ou->LZ = in->LZ;
    ou->LT = in->LT;

    ou->nx = (int) (1 + ceilf(ou->LX/ou->dx) + 2*ou->nxb );
    ou->ny = (int) (1 + ceilf(ou->LY/ou->dy) + 2*ou->nyb );
    ou->nz = (int) (1 + ceilf(ou->LZ/ou->dz) + 2*ou->nzb );
    ou->nt = (int) (1 + ceilf(ou->LT/ou->dt) );

    if( (ou->nx)%strideX != 0 )    ou->nx = ( 1 + (ou->nx)/strideX ) * strideX;
    if( (ou->ny)%strideY != 0 )    ou->ny = ( 1 + (ou->ny)/strideY ) * strideY;
    if( (ou->nz)%strideZ != 0 )    ou->nz = ( 1 + (ou->nz)/strideZ ) * strideZ;

    ou->LX = (ou->nx - 2*ou->nxb - 1) * ou->dx;
    ou->LY = (ou->ny - 2*ou->nyb - 1) * ou->dy;
    ou->LZ = (ou->nz - 2*ou->nzb - 1) * ou->dz;

    ou->ox -= ou->nxb * ou->dx;
    ou->oy -= ou->nyb * ou->dy;
    ou->oz -= ou->nzb * ou->dz;

    ou->jt     = in->jt;
    ou->jt_lt  = in->jt_lt;
    ou->jt_EIC = in->jt_EIC;

    // printf("\n msh GPU compliant: nx=%d ny=%d nz=%d dx=%f dy=%f dz=%f ox=%f oy=%f oz=%f LX=%f LY=%f LZ=%f jt=%d jt_lt=%d jt_EIC=%d\n", \
             ou->nx, ou->ny, ou->nz, ou->dx, ou->dy, ou->dz, ou->ox, ou->oy, ou->oz, ou->LX, ou->LY, ou->LZ, ou->jt, ou->jt_lt, ou->jt_EIC);

return;
}

float getSpaceSampling(float maxFrequency, float vmin)
{
    float samplesPerWavelength = 3.2;
    float minWavelength = vmin / maxFrequency;
    float spaceSampling = minWavelength / samplesPerWavelength;

return spaceSampling;
}

float getStableTimeSampling(float spaceSampling, float vmax)
{
    float coef = 3.0f * sqrtf(3.0f);
    float timeSampling = spaceSampling / (coef * vmax);

return timeSampling;
}


void stackArray_MPI(float *array, int inode, int inode_manager, long long int N, MPI_Datatype datatype, MPI_Comm comm)
{
    float *array_tmp;
    
    if(inode==inode_manager)    array_tmp  = CPU_zaloc1F(N);
    MPI_Barrier(comm);

    MPI_Reduce(array, array_tmp, N, datatype, MPI_SUM, 0, comm);
    MPI_Barrier(comm);

    if(inode==inode_manager) {
        memcpy( array,  array_tmp, N*sizeof(float));
        free(array_tmp);
    }
    else {
        memset( array, 0, N*sizeof(float));
    }
    MPI_Barrier(comm);

    MPI_Bcast(array, N, datatype, inode_manager, comm);
    MPI_Barrier(comm);

return;
}
void stackLargeArray_MPI(float *array, int inode, int inode_manager, long long int N, MPI_Datatype datatype, MPI_Comm comm)
{
    long long int maxInt = 1073741824; // equal to 2^30. The maximum in is actually 2ˆ31, but I'm leaving a little safety margin just in case.

    if(N<=maxInt)    stackArray_MPI(array, inode, inode_manager, N, datatype, comm);
    else
    {
        long long int N_;
        long long int remainder;
        long long int nturns, iturn;

        N_ = min(N, maxInt);
        remainder = N%maxInt;
        nturns = 1 + N/maxInt;

        for(iturn=0; iturn<nturns; iturn++)
        {
            long long int chunk;
            if(iturn<(nturns-1))    chunk = N_;
            else                    chunk = remainder;

            size_t shiftMem = N_ * iturn;
            stackArray_MPI(array+shiftMem, inode, inode_manager, chunk, datatype, comm);
        }
    }

return;
}