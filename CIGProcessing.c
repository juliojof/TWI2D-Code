void get_HVocig_ExtExpRef(float *Hocig_background, float *Vocig_background, \
                          float *Hocig_residual,   float *Vocig_residual, \
                          float *Hocig_background_phaseShift, float *Hocig_background_original, \
                          float *Hocig_residual_phaseShift, float *Hocig_residual_original, \
                          float *Vocig_background_phaseShift, float *Vocig_background_original, \
                          float *Vocig_residual_phaseShift, float *Vocig_residual_original, \
                          long long int lx0, long long int lz0, long long int nx, long long int nz, long long int nt, \
                          float dx, float dz, float dt, float ox, float oz, float ot, int iproc, int nproc, float *vp_bg)
{
    
    long long int nlx = 2*lx0 + 1;
    long long int nlz = 2*lz0 + 1;

    time_t time_now, time_start;
    time_start = time(NULL);

    float *vel = CPU_zaloc1F(nx*nz);
    transp(vel, vp_bg, nz, nx);
    untransformVel(vel, nx, nz, dx, dt);    
    int iproc_bck_org, iproc_bck_phs, iproc_res_org, iproc_res_phs;
    if(nproc>=4)
    {
        iproc_bck_org = 0;
        iproc_bck_phs = 1;
        iproc_res_org = 2;
        iproc_res_phs = 3;
    }
    else if(nproc>=2)
    {
        iproc_bck_org = 0;
        iproc_bck_phs = 1;
        iproc_res_org = 0;
        iproc_res_phs = 1;
    }
    else
    {
        iproc_bck_org = 0;
        iproc_bck_phs = 0;
        iproc_res_org = 0;
        iproc_res_phs = 0;
    }

    if(iproc==iproc_bck_org) {
        printf("\n iproc %d performing Hocig_background_original  procGPU=%d", iproc, procGPU);
        time_now = time(NULL);
        printf("\n\n before ..... in phase-shift elapsed time %ld seconds iproc=%d \n\n", time_now-time_start, iproc);
        GPU_WRAP_HocigPhaseShift(procGPU, Hocig_background_original  , Hocig_background, vel, nx, nz, nt/50, nxb, nzb, dx, dz, dt, lx0, 1);
        time_now = time(NULL);
        printf("\n\n after ..... in phase-shift elapsed time %ld seconds iproc=%d \n\n", time_now-time_start, iproc);
    }   
    if(iproc==iproc_bck_phs) {
        printf("\n iproc %d performing Hocig_background_phaseShift  procGPU=%d", iproc, procGPU);
        time_now = time(NULL);
        printf("\n\n before ..... in phase-shift elapsed time %ld seconds iproc=%d \n\n", time_now-time_start, iproc);
        GPU_WRAP_HocigPhaseShift(procGPU, Hocig_background_phaseShift, Hocig_background, vel, nx, nz, nt/50, nxb, nzb, dx, dz, dt, lx0, 2);
        time_now = time(NULL);
        printf("\n\n after ..... in phase-shift elapsed time %ld seconds iproc=%d \n\n", time_now-time_start, iproc);
    }   
    if(iproc==iproc_res_org) {
        printf("\n iproc %d performing Hocig_residual_original  procGPU=%d", iproc, procGPU);
        time_now = time(NULL);
        printf("\n\n before ..... in phase-shift elapsed time %ld seconds iproc=%d \n\n", time_now-time_start, iproc);
        GPU_WRAP_HocigPhaseShift(procGPU, Hocig_residual_original    , Hocig_residual  , vel, nx, nz, nt/50, nxb, nzb, dx, dz, dt, lx0, 1);
        time_now = time(NULL);
        printf("\n\n after ..... in phase-shift elapsed time %ld seconds iproc=%d \n\n", time_now-time_start, iproc);
    }   
    if(iproc==iproc_res_phs) {
        printf("\n iproc %d performing Hocig_residual_phaseShift  procGPU=%d", iproc, procGPU);
        time_now = time(NULL);
        printf("\n\n before ..... in phase-shift elapsed time %ld seconds iproc=%d \n\n", time_now-time_start, iproc);
        GPU_WRAP_HocigPhaseShift(procGPU, Hocig_residual_phaseShift  , Hocig_residual  , vel, nx, nz, nt/50, nxb, nzb, dx, dz, dt, lx0, 2);
        time_now = time(NULL);
        printf("\n\n after ..... in phase-shift elapsed time %ld seconds iproc=%d \n\n", time_now-time_start, iproc);
    }   
    MPI_Barrier(MPI_COMM_WORLD);

    free(vel);

    float *Hocig_background_phaseShift_tmp, *Hocig_residual_original_tmp, *Hocig_residual_phaseShift_tmp;
    if(nproc>=2)
    {
        Hocig_background_phaseShift_tmp  = CPU_zaloc1F(nx*nz*nlx);
        Hocig_residual_original_tmp      = CPU_zaloc1F(nx*nz*nlx);
        Hocig_residual_phaseShift_tmp    = CPU_zaloc1F(nx*nz*nlx);

        MPI_Reduce(Hocig_background_phaseShift, Hocig_background_phaseShift_tmp, nx*nz*nlx, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(Hocig_residual_original    , Hocig_residual_original_tmp    , nx*nz*nlx, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(Hocig_residual_phaseShift  , Hocig_residual_phaseShift_tmp  , nx*nz*nlx, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        if(iproc==0)
        {
            memcpy(Hocig_background_phaseShift, Hocig_background_phaseShift_tmp, nx*nz*nlx*sizeof(float));
            memcpy(Hocig_residual_original    , Hocig_residual_original_tmp    , nx*nz*nlx*sizeof(float));
            memcpy(Hocig_residual_phaseShift  , Hocig_residual_phaseShift_tmp  , nx*nz*nlx*sizeof(float));
        }
        MPI_Barrier(MPI_COMM_WORLD);
        free(Hocig_background_phaseShift_tmp);
        free(Hocig_residual_original_tmp    );
        free(Hocig_residual_phaseShift_tmp  );

        MPI_Bcast(Hocig_background_original  , nx*nz*nlx, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Bcast(Hocig_background_phaseShift, nx*nz*nlx, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Bcast(Hocig_residual_original    , nx*nz*nlx, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Bcast(Hocig_residual_phaseShift  , nx*nz*nlx, MPI_FLOAT, 0, MPI_COMM_WORLD);
        
        MPI_Barrier(MPI_COMM_WORLD);
    }

    time_now = time(NULL);
    printf("\n\n >>>>>> Performed phase-shift in %ld seconds \n\n", time_now-time_start);

    if(iproc==0)
    {
        outputSmart3d("Hocig_background"           , Hocig_background           , nlx, dx, -lx0*dx,  nz, dz, oz,  nx, dx, ox);
        outputSmart3d("Hocig_background_original"  , Hocig_background_original  , nlx, dx, -lx0*dx,  nz, dz, oz,  nx, dx, ox);
        outputSmart3d("Hocig_background_phaseShift", Hocig_background_phaseShift, nlx, dx, -lx0*dx,  nz, dz, oz,  nx, dx, ox);
        outputSmart3d("Hocig_residual"             , Hocig_residual             , nlx, dx, -lx0*dx,  nz, dz, oz,  nx, dx, ox);
        outputSmart3d("Hocig_residual_original"    , Hocig_residual_original    , nlx, dx, -lx0*dx,  nz, dz, oz,  nx, dx, ox);
        outputSmart3d("Hocig_residual_phaseShift"  , Hocig_residual_phaseShift  , nlx, dx, -lx0*dx,  nz, dz, oz,  nx, dx, ox);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);


return;
}


void get_HVocig_ExtExpRef_fromDocig(float *Docig_background, float *Docig_residual, \
                                    float *Docig_background_phaseShift, float *Docig_background_original, \
                                    float *Docig_residual_phaseShift, float *Docig_residual_original, \
                                    long long int lx0, long long int ndip, float dipMin, float ddip, \
                                    long long int nx, long long int nz, long long int nt, \
                                    float dx, float dz, float dt, float ox, float oz, float ot, int iproc, int nproc, float *vp_bg)
{
    dipMin *= M_PI/180.0;
    ddip   *= M_PI/180.0;

    long long int nlx = 2*lx0 + 1;

    time_t time_now, time_start;
    time_start = time(NULL);

    float *vel = CPU_zaloc1F(nx*nz);
    transp(vel, vp_bg, nz, nx);
    untransformVel(vel, nx, nz, dx, dt);    
    int iproc_bck_org, iproc_bck_phs, iproc_res_org, iproc_res_phs;
    if(nproc>=4)
    {
        iproc_bck_org = 0;
        iproc_bck_phs = 1;
        iproc_res_org = 2;
        iproc_res_phs = 3;
    }
    else if(nproc>=2)
    {
        iproc_bck_org = 0;
        iproc_bck_phs = 1;
        iproc_res_org = 0;
        iproc_res_phs = 1;
    }
    else
    {
        iproc_bck_org = 0;
        iproc_bck_phs = 0;
        iproc_res_org = 0;
        iproc_res_phs = 0;
    }

    long long int idip;
    float dip;
    // long long int nt_ = max(nt/50,10);
    long long int nt_ = 20;
    for(idip=0; idip<ndip; idip++)
    {
        dip = dipMin + idip*ddip;
        long long int shiftDip = idip * nlx*nx*nz;
        if(iproc==iproc_bck_org)
        {
            GPU_WRAP_DocigPhaseShift(procGPU, Docig_background_original+shiftDip  , Docig_background+shiftDip, \
                                     vel, dip, nx, nz, nt_, nxb, nzb, dx, dz, dt, lx0, 1);
        }   
        if(iproc==iproc_bck_phs)
        {
            GPU_WRAP_DocigPhaseShift(procGPU, Docig_background_phaseShift+shiftDip, Docig_background+shiftDip, \
                                     vel, dip, nx, nz, nt_, nxb, nzb, dx, dz, dt, lx0, 2);
        }   
        if(iproc==iproc_res_org)
        {
            GPU_WRAP_DocigPhaseShift(procGPU, Docig_residual_original+shiftDip    , Docig_residual+shiftDip  , \
                                     vel, dip, nx, nz, nt_, nxb, nzb, dx, dz, dt, lx0, 1);
        }   
        if(iproc==iproc_res_phs)
        {
            GPU_WRAP_DocigPhaseShift(procGPU, Docig_residual_phaseShift+shiftDip  , Docig_residual+shiftDip  , \
                                     vel, dip, nx, nz, nt_, nxb, nzb, dx, dz, dt, lx0, 2);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    free(vel);

    float *Docig_background_phaseShift_tmp, *Docig_residual_original_tmp, *Docig_residual_phaseShift_tmp;
    if(nproc>=2)
    {
        Docig_background_phaseShift_tmp  = CPU_zaloc1F(nx*nz*nlx*ndip);
        Docig_residual_original_tmp      = CPU_zaloc1F(nx*nz*nlx*ndip);
        Docig_residual_phaseShift_tmp    = CPU_zaloc1F(nx*nz*nlx*ndip);

        MPI_Reduce(Docig_background_phaseShift, Docig_background_phaseShift_tmp, nx*nz*nlx*ndip, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(Docig_residual_original    , Docig_residual_original_tmp    , nx*nz*nlx*ndip, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(Docig_residual_phaseShift  , Docig_residual_phaseShift_tmp  , nx*nz*nlx*ndip, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        if(iproc==0)
        {
            memcpy(Docig_background_phaseShift, Docig_background_phaseShift_tmp, nx*nz*nlx*ndip*sizeof(float));
            memcpy(Docig_residual_original    , Docig_residual_original_tmp    , nx*nz*nlx*ndip*sizeof(float));
            memcpy(Docig_residual_phaseShift  , Docig_residual_phaseShift_tmp  , nx*nz*nlx*ndip*sizeof(float));
        }
        MPI_Barrier(MPI_COMM_WORLD);
        free(Docig_background_phaseShift_tmp);
        free(Docig_residual_original_tmp    );
        free(Docig_residual_phaseShift_tmp  );

        MPI_Bcast(Docig_background_original  , nx*nz*nlx*ndip, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Bcast(Docig_background_phaseShift, nx*nz*nlx*ndip, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Bcast(Docig_residual_original    , nx*nz*nlx*ndip, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Bcast(Docig_residual_phaseShift  , nx*nz*nlx*ndip, MPI_FLOAT, 0, MPI_COMM_WORLD);
        
        MPI_Barrier(MPI_COMM_WORLD);
    }

    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);  if(iproc==0)  printf("\n\n >>>>>> Performed phase-shift in %ld seconds \n\n", time_now-time_start);

    /*
    if(iproc==0)
    {
        outputSmart3d("Docig_background"           , Docig_background           , nlx, dx, -lx0*dx,  nz, dz, oz,  nx, dx, ox);
        outputSmart3d("Docig_background_original"  , Docig_background_original  , nlx, dx, -lx0*dx,  nz, dz, oz,  nx, dx, ox);
        outputSmart3d("Docig_background_phaseShift", Docig_background_phaseShift, nlx, dx, -lx0*dx,  nz, dz, oz,  nx, dx, ox);
        outputSmart3d("Docig_residual"             , Docig_residual             , nlx, dx, -lx0*dx,  nz, dz, oz,  nx, dx, ox);
        outputSmart3d("Docig_residual_original"    , Docig_residual_original    , nlx, dx, -lx0*dx,  nz, dz, oz,  nx, dx, ox);
        outputSmart3d("Docig_residual_phaseShift"  , Docig_residual_phaseShift  , nlx, dx, -lx0*dx,  nz, dz, oz,  nx, dx, ox);
    }
    //*/
    
    MPI_Barrier(MPI_COMM_WORLD);


return;
}