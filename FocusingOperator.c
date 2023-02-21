void FocusingOperator_lx_lz_lt(int adj, float *eimage_diff, float *eimage_contracted, float *eimage_contracted_ref, float *eimage, \
                               float zMin, float zMax, float maxDepth, long long int nx, long long int nz, float dx, float dz, float ox, float oz, \
                               int lx0, int lz0, int lt0, int ltheta0, float dtheta, float otheta, \
                               float dipMin, float dipMax, float ddip, float lambdaMax, float dlambda, int applyADFocusing, \
                               int iter, float mu, float epsFocus, float lambdaNull, int outputs, int iproc, int nproc) {
    
    time_t time_now, time_start, time_begin, time_elapsed;
    
    int nlx    = 2*lx0+1;
    int nlz    = 2*lz0+1;
    int nlt    = 2*lt0+1;
    int ntheta = 2*ltheta0+1;

    float dtanth = dtheta * M_PI/180.0;
    int ltanth0 = (int) ceilf(tanf(ltheta0*dtheta*M_PI/180.0f)/dtanth);
    int ntanth = 2*ltanth0+1;

    time_start = time(NULL);
    if(iproc==0)    printf("\n\n Applying Demoveout .....  lz0=%d    mu=%f    eps=%f   lambdaMax=%f", lz0, mu, epsFocus, lambdaMax);
    float *eimage_contracted_tmp     = CPU_zaloc1F(nx*nz*nlx*nlz*nlt);
    float *eimage_contracted_tmp_ref = CPU_zaloc1F(nx*nz*nlx*nlz*nlt);
    
    int FreeContractedRef = 0;
    if(eimage_contracted_ref==NULL)
    {
        eimage_contracted_ref = CPU_zaloc1F(nx*nz*nlx*nlz*nlt);
        FreeContractedRef = 1;
    }

    float xMin = 0;
    float xMax = (nx-2*nxb)*dx;

    // Apply focusing operator in the angle-domain
    if(applyADFocusing)
    {
    /*
    float *adcig                  = CPU_zaloc1F(nx*nz*ntheta);
    float *adcig_tanTheta         = CPU_zaloc1F(nx*nz*ntanth);
    float *adcig_ADcontracted     = CPU_zaloc1F(nx*nz*ntheta);
    float *adcig_ADcontracted_ref = CPU_zaloc1F(nx*nz*ntheta);
    float *adcig_TDcontracted     = CPU_zaloc1F(nx*nz*ntanth);
    float *adcig_TDcontracted_ref = CPU_zaloc1F(nx*nz*ntanth);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Transforming ODCIG to ADCIG to apply focusing operator in the angle-domain ....");
    GPU_WRAP_lx2theta(0, procGPU, 0, eimage, adcig, nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);
    
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Transforming from theta to tangent: ntanth=%d  dtheta=%f....", ntanth, dtheta);
    GPU_WRAP_adcig_theta2tanTheta(0, procGPU, 0, adcig, adcig_tanTheta, nx, dx, nz, dz, ltheta0, dtheta, ltanth0, dtanth, nxb, nzb, 0, (nx-2*nxb)*dx, 0, (nz-2*nzb)*dz);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf(" done \n");


    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Applying focusing operator in the angle-domain ....");
    // GPU_WRAP_Demoveout_AD(procGPU, 0, adj, adcig, adcig_ADcontracted    , nx, dx, nz, dz, lx0, ltheta0, dtheta, nxb, nzb, xMin, xMax, zMin, zMax, 1.00-mu, lx0*dx);
    // GPU_WRAP_Demoveout_AD(procGPU, 0, adj, adcig, adcig_ADcontracted_ref, nx, dx, nz, dz, lx0, ltheta0, dtheta, nxb, nzb, xMin, xMax, zMin, zMax, 1.00   , lx0*dx);
    GPU_WRAP_Demoveout_AD(procGPU, 0, adj, adcig_tanTheta, adcig_TDcontracted    , nx, dx, nz, dz, lx0, ltanth0, dtanth, nxb, nzb, xMin, xMax, zMin, zMax, 1.00-mu, lx0*dx);
    GPU_WRAP_Demoveout_AD(procGPU, 0, adj, adcig_tanTheta, adcig_TDcontracted_ref, nx, dx, nz, dz, lx0, ltanth0, dtanth, nxb, nzb, xMin, xMax, zMin, zMax, 1.00   , lx0*dx);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf(" done \n");

    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Transforming from tangent to theta....");
    GPU_WRAP_adcig_theta2tanTheta(1, procGPU, 0, adcig_ADcontracted_ref, adcig_TDcontracted_ref, nx, dx, nz, dz, ltheta0, dtheta, ltanth0, dtanth, nxb, nzb, 0, (nx-2*nxb)*dx, 0, (nz-2*nzb)*dz);
    GPU_WRAP_adcig_theta2tanTheta(1, procGPU, 0, adcig_ADcontracted    , adcig_TDcontracted    , nx, dx, nz, dz, ltheta0, dtheta, ltanth0, dtanth, nxb, nzb, 0, (nx-2*nxb)*dx, 0, (nz-2*nzb)*dz);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf(" done \n");

    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Recovering ODCIG from contracted ADCIG....");
    float *eimage_ADcontracted     = cg_theta2lx_gpu(procGPU, adcig_ADcontracted    , nx, nz, dx, dz, lx0, ntheta, dtheta, otheta, 5);
    float *eimage_ADcontracted_ref = cg_theta2lx_gpu(procGPU, adcig_ADcontracted_ref, nx, nz, dx, dz, lx0, ntheta, dtheta, otheta, 5);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf(" done \n");
    if(outputs)
    {
        char adcig_input_FileName[1024];
        char adcig_ADcontracted_FileName[1024];
        char adcig_ADcontracted_ref_FileName[1024];
        char adcig_TD_FileName[1024];
        char adcig_TDcontracted_FileName[1024];
        char adcig_TDcontracted_ref_FileName[1024];
        char eimage_ADcontracted_FileName[1024];
        char eimage_ADcontracted_ref_FileName[1024];
        char adcig_FileName[1024];
        sprintf(adcig_input_FileName                      , "adcig_input_Iter%d"            , iter);
        sprintf(eimage_ADcontracted_FileName              , "eimage_ADcontracted_Iter%d"    , iter);
        sprintf(eimage_ADcontracted_ref_FileName          , "eimage_ADcontracted_ref_Iter%d", iter);
        sprintf(adcig_ADcontracted_FileName               , "adcigs_ADcontracted_Iter%d"    , iter);
        sprintf(adcig_ADcontracted_ref_FileName           , "adcigs_ADcontracted_ref_Iter%d", iter);
        sprintf(adcig_TD_FileName                         , "adcigs_TD_Iter%d"              , iter);
        sprintf(adcig_TDcontracted_FileName               , "adcigs_TDcontracted_Iter%d"    , iter);
        sprintf(adcig_TDcontracted_ref_FileName           , "adcigs_TDcontracted_ref_Iter%d", iter);
        outputSmart3d(adcig_input_FileName                , adcig                           , ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
        outputSmart3d(adcig_ADcontracted_ref_FileName     , adcig_ADcontracted_ref          , ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
        outputSmart3d(adcig_ADcontracted_FileName         , adcig_ADcontracted              , ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
        outputSmart3d(adcig_TD_FileName                   , adcig_tanTheta                  , ntanth, dtanth, -ltanth0*dtanth, nz, dz, oz, nx, dx, ox);
        outputSmart3d(adcig_TDcontracted_ref_FileName     , adcig_TDcontracted_ref          , ntanth, dtanth, -ltanth0*dtanth, nz, dz, oz, nx, dx, ox);
        outputSmart3d(adcig_TDcontracted_FileName         , adcig_TDcontracted              , ntanth, dtanth, -ltanth0*dtanth, nz, dz, oz, nx, dx, ox);
        outputSmart3d(eimage_ADcontracted_ref_FileName    , eimage_ADcontracted_ref         ,    nlx,     dx,         -lx0*dx, nz, dz, oz, nx, dx, ox);
        outputSmart3d(eimage_ADcontracted_FileName        , eimage_ADcontracted             ,    nlx,     dx,         -lx0*dx, nz, dz, oz, nx, dx, ox);
    }
    free(adcig                  );
    free(adcig_ADcontracted     );
    free(adcig_ADcontracted_ref );
    free(eimage_ADcontracted    );
    free(eimage_ADcontracted_ref);
    //*/
    }

    if(lx0>0  &&  lz0==0  &&  lt0==0)
    {
        GPU_WRAP_Demoveout(procGPU, 0, adj, eimage, eimage_contracted_tmp    , nx, dx, nz, dz, lx0, nxb, nzb, xMin, xMax, zMin, zMax, 1.00-mu, 0*lambdaNull);
        GPU_WRAP_Demoveout(procGPU, 0, adj, eimage, eimage_contracted_tmp_ref, nx, dx, nz, dz, lx0, nxb, nzb, xMin, xMax, zMin, zMax, 1.00   , 0*lambdaNull);
        // GPU_WRAP_Demoveout(procGPU, 0, adj, eimage, eimage_contracted_tmp    , nx, dx, nz, dz, lx0, nxb, nzb, xMin, xMax, zMin, zMax, 1.00-mu, lx0*dx);
        // GPU_WRAP_Demoveout(procGPU, 0, adj, eimage, eimage_contracted_tmp_ref, nx, dx, nz, dz, lx0, nxb, nzb, xMin, xMax, zMin, zMax, 1.00   , lx0*dx);
    }
    if(lx0>0  &&  lz0>0  &&  lt0==0)
    {
        GPU_WRAP_Demoveout_lx_lz(procGPU, 0, adj, eimage, eimage_contracted_tmp    , nx, dx, nz, dz, lx0, lz0, \
                                 dipMin, dipMax, ddip, lambdaMax, dlambda, nxb, nzb, xMin, xMax, zMin, zMax, 1.00-mu, lambdaMax);

        GPU_WRAP_Demoveout_lx_lz(procGPU, 0, adj, eimage, eimage_contracted_tmp_ref, nx, dx, nz, dz, lx0, lz0, \
                                 dipMin, dipMax, ddip, lambdaMax, dlambda, nxb, nzb, xMin, xMax, zMin, zMax, 1.00   , lambdaMax);
    }
    if(lx0>0  &&  lz0>0  &&  lt0>0)
    {
        // GPU_WRAP_Demoveout_lx_lz_lt(procGPU, 0, adj, eimage, eimage_contracted_tmp    , nx, dx, nz, dz, lx0, lz0, \
        //                             dipMin, dipMax, ddip, lambdaMax, dlambda, nxb, nzb, xMin, xMax, zMin, zMax, 1.00-mu, lambdaMax);

        // GPU_WRAP_Demoveout_lx_lz_lt(procGPU, 0, adj, eimage, eimage_contracted_tmp_ref, nx, dx, nz, dz, lx0, lz0, lt0, \
        //                             dipMin, dipMax, ddip, lambdaMax, dlambda, nxb, nzb, xMin, xMax, zMin, zMax, 1.00   , lambdaMax);
    }
    time_now = time(NULL);
    if(iproc==0)
        printf("  Concluded in %ld seconds\n", time_now-time_start);

    time_start = time(NULL);
    if(iproc==0)
        printf("\n Applying contraction .....");
    // float stretchingFactor = (1.0/(1.0-mu)) - epsFocus;
    // if(adj) {
    //     stretchingFactor = (1.00-mu) + epsFocus;
    // }
    float stretchingFactor = 1.0 + epsFocus;
    if(adj) {
          stretchingFactor = 1.0 - epsFocus;
    }
    


    // lx2Contraction_forward(0, eimage_contracted_tmp    , eimage_contracted    , nx, dx,  nz, lx0, stretchingFactor);
    // lx2Contraction_forward(0, eimage_contracted_tmp_ref, eimage_contracted_ref, nx, dx,  nz, lx0,  1.000          );
    if(lz0==0)
    {
        GPU_WRAP_Contraction_lx(procGPU, 0, eimage_contracted_tmp    , eimage_contracted    , nx, dx, nz, dz, lx0, \
                                nxb, nzb, xMin, xMax, zMin, zMax, stretchingFactor);
        GPU_WRAP_Contraction_lx(procGPU, 0, eimage_contracted_tmp_ref, eimage_contracted_ref, nx, dx, nz, dz, lx0, \
                                nxb, nzb, xMin, xMax, zMin, zMax, 1.00);
    }
    else
    {
        // memcpy(eimage_contracted_ref, eimage_contracted_tmp_ref, nx*nz*nlx*nlz*sizeof(float));
        // memcpy(eimage_contracted    , eimage_contracted_tmp    , nx*nz*nlx*nlz*sizeof(float));

        GPU_WRAP_Contraction_lx_lz(procGPU, 0, eimage_contracted_tmp    , eimage_contracted    , nx, dx, nz, dz, lx0, lz0, \
                                   dipMin, dipMax, ddip, lambdaMax, dlambda, \
                                   nxb, nzb, xMin, xMax, zMin, zMax, stretchingFactor);
        GPU_WRAP_Contraction_lx_lz(procGPU, 0, eimage_contracted_tmp_ref, eimage_contracted_ref, nx, dx, nz, dz, lx0, lz0, \
                                   dipMin, dipMax, ddip, lambdaMax, dlambda, \
                                   nxb, nzb, xMin, xMax, zMin, zMax, 1.0);
    }
    
    time_now = time(NULL);
    if(iproc==0)    printf("  Concluded in %ld seconds \n\n", time_now-time_start);

    float *adcig_contracted     = CPU_zaloc1F(nx*nz*ntheta);
    float *adcig_contracted_ref = CPU_zaloc1F(nx*nz*ntheta);
    float *adcig_diff           = CPU_zaloc1F(nx*nz*ntheta);
    // lx2theta_forward(0, eimage_contracted    , adcig_contracted    , nx, dx, nz,  dz, lx0, ntheta, dtheta, otheta);
    // lx2theta_forward(0, eimage_contracted_ref, adcig_contracted_ref, nx, dx, nz,  dz, lx0, ntheta, dtheta, otheta);
    GPU_WRAP_lx2theta(0, procGPU, 0, eimage_contracted    , adcig_contracted    , nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);
    GPU_WRAP_lx2theta(0, procGPU, 0, eimage_contracted_ref, adcig_contracted_ref, nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);

    addArrays(eimage_diff, 1.0, eimage_contracted, -1.0, eimage_contracted_ref, nx*nz*nlx*nlz);
    addArrays( adcig_diff, 1.0,  adcig_contracted, -1.0,  adcig_contracted_ref, nx*nz*ntheta);
    
   
    // Save contracted image to file
    if(outputs)
    {
        char eimage_contracted_FileName[1024];
        char eimage_contracted_ref_FileName[1024];
        char eimage_contracted_tmp_FileName[1024];
        char eimage_contracted_tmp_ref_FileName[1024];
        char adcigs_contracted_FileName[1024];
        char adcigs_contracted_ref_FileName[1024];
        char eimage_diff_FileName[1024];
        char adcigs_diff_FileName[1024];

        sprintf(eimage_contracted_FileName         , "eimage_contracted_Iter%d", iter);
        sprintf(eimage_contracted_tmp_FileName     , "eimage_contracted_tmp_Iter%d", iter);
        sprintf(adcigs_contracted_FileName         , "adcigs_contracted_Iter%d", iter);
        sprintf(eimage_contracted_ref_FileName     , "eimage_contracted_ref_Iter%d", iter);
        sprintf(eimage_contracted_tmp_ref_FileName , "eimage_contracted_tmp_ref_Iter%d", iter);
        sprintf(adcigs_contracted_ref_FileName     , "adcigs_contracted_ref_Iter%d", iter);
        sprintf(eimage_diff_FileName               , "eimage_diff_Iter%d", iter);
        sprintf(adcigs_diff_FileName               , "adcigs_diff_Iter%d", iter);
        
        outputSmart3d(adcigs_contracted_ref_FileName    , adcig_contracted_ref      , ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
        outputSmart3d(adcigs_contracted_FileName        , adcig_contracted          , ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
        outputSmart3d(adcigs_diff_FileName              , adcig_diff                , ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
        outputSmart3d(eimage_contracted_ref_FileName    , eimage_contracted_ref     ,    nlx,     dx,         -lx0*dx, nz, dz, oz, nx, dx, ox);
        outputSmart3d(eimage_contracted_FileName        , eimage_contracted         ,    nlx,     dx,         -lx0*dx, nz, dz, oz, nx, dx, ox);
        outputSmart3d(eimage_contracted_tmp_ref_FileName, eimage_contracted_tmp_ref ,    nlx,     dx,         -lx0*dx, nz, dz, oz, nx, dx, ox);
        outputSmart3d(eimage_contracted_tmp_FileName    , eimage_contracted_tmp     ,    nlx,     dx,         -lx0*dx, nz, dz, oz, nx, dx, ox);
        outputSmart3d(eimage_diff_FileName              , eimage_diff               ,    nlx,     dx,         -lx0*dx, nz, dz, oz, nx, dx, ox);
    }

    free(eimage_contracted_tmp);
    free(eimage_contracted_tmp_ref);
    free(adcig_contracted);
    free(adcig_contracted_ref);
    free(adcig_diff);
    if(FreeContractedRef)
        free(eimage_contracted_ref);

return;  
}

/*
void FocusingOperator_TimeLag(int adj, float *eimage_diff, float *eimage_contracted, float *eimage_contracted_ref, float *eimage, \
                              float zMin, float zMax, long long int nx, long long int nz, float dx, float dz, float ox, float oz, \
                              int lt0, int iter, float mu, float lt0_null, int outputs, int iproc, int nproc) {
    
    time_t time_now, time_start, time_begin, time_elapsed;
    
    int nlt    = 2*lt0+1;

    time_start = time(NULL);
    if(iproc==0)    printf("\n\n Applying Demoveout .....  lz0=%d    mu=%f    eps=%f   lambdaMax=%f", lz0, mu, epsFocus, lambdaMax);
    
    int FreeContractedRef = 0;
    if(eimage_contracted_ref==NULL) {
        eimage_contracted_ref = CPU_zaloc1F(nx*nz*nlt);
        FreeContractedRef = 1;
    }

    float xMin = 0;
    float xMax = (nx-2*nxb)*dx;


    GPU_WRAP_Contraction_lt(procGPU, 0, eimage_contracted    , eimage, nx, dx, nz, dz, lt0, \
                            nxb, nzb, xMin, xMax, zMin, zMax, 1.0-mu);
    GPU_WRAP_Contraction_lt(procGPU, 0, eimage_contracted_ref, eimage, nx, dx, nz, dz, lt0, \
                            nxb, nzb, xMin, xMax, zMin, zMax, 1.00);


    addArrays(eimage_diff, 1.0, eimage_contracted, -1.0, eimage_contracted_ref, nx*nz*nlx*nlz);
    
   
    // Save contracted image to file
    if(outputs)
    {
        char eimage_contracted_FileName[1024];
        char eimage_contracted_ref_FileName[1024];
        char eimage_diff_FileName[1024];

        sprintf(eimage_contracted_FileName         , "eimage_contracted_Iter%d", iter);
        sprintf(eimage_contracted_ref_FileName     , "eimage_contracted_ref_Iter%d", iter);
        sprintf(eimage_diff_FileName               , "eimage_diff_Iter%d", iter);
        
        outputSmart3d(eimage_contracted_ref_FileName    , eimage_contracted_ref     ,    nlx,     dx,         -lx0*dx, nz, dz, oz, nx, dx, ox);
        outputSmart3d(eimage_contracted_FileName        , eimage_contracted         ,    nlx,     dx,         -lx0*dx, nz, dz, oz, nx, dx, ox);
        outputSmart3d(eimage_diff_FileName              , eimage_diff               ,    nlx,     dx,         -lx0*dx, nz, dz, oz, nx, dx, ox);
    }

    if(FreeContractedRef)
        free(eimage_contracted_ref);

return;  
}
//*/


void FocusingOperator_AD(int adj, float *eimage_diff, float *eimage_contracted, float *eimage_contracted_ref, float *eimage, \
                         float zMin, float zMax, float maxDepth, \
                         long long int nx, long long int nz, float dx, float dz, float ox, float oz, \
                         int lx0, int ltheta0, float dtheta, float otheta, \
                         int iter, float mu, float epsFocus, float lambdaNull, int outputs, int iproc, int nproc) {
    
    time_t time_now, time_start, time_begin, time_elapsed;
    
    int nlx    = 2*lx0+1;
    int ntheta = 2*ltheta0+1;

    float dtanth = dtheta * M_PI/180.0;
    int ltanth0 = (int) ceilf(tanf(ltheta0*dtheta*M_PI/180.0f)/dtanth);
    int ntanth = 2*ltanth0+1;

    time_start = time(NULL);
    if(iproc==0)    printf("\n\n Applying Demoveout .....  mu=%f    eps=%f   lambdaMax=%f", mu, epsFocus, lx0*dx);
    float *eimage_contracted_tmp     = CPU_zaloc1F(nx*nz*nlx);
    float *eimage_contracted_tmp_ref = CPU_zaloc1F(nx*nz*nlx);
    

    float xMin = 0;
    float xMax = (nx-2*nxb)*dx;

    // Apply focusing operator in the angle-domain

    int niterCG = 10;

    float *adcig                  = CPU_zaloc1F(nx*nz*ntheta);
    float *adcig_tanTheta         = CPU_zaloc1F(nx*nz*ntanth);
    float *adcig_recover          = CPU_zaloc1F(nx*nz*ntheta);
    float *adcig_adjoint          = CPU_zaloc1F(nx*nz*ntheta);
    float *adcig_ADcontracted     = CPU_zaloc1F(nx*nz*ntheta);
    float *adcig_ADcontracted_ref = CPU_zaloc1F(nx*nz*ntheta);
    float *adcig_TDcontracted     = CPU_zaloc1F(nx*nz*ntanth);
    float *adcig_TDcontracted_ref = CPU_zaloc1F(nx*nz*ntanth);

    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Transforming ODCIG to ADCIG to apply focusing operator in the angle-domain ....");
    GPU_WRAP_lx2theta(0, procGPU, 0, eimage, adcig, nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);
    
    /*
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Transforming from theta to tangent: ntanth=%d  dtheta=%f....", ntanth, dtheta);
    GPU_WRAP_adcig_theta2tanTheta(0, procGPU, 0, adcig, adcig_tanTheta, nx, dx, nz, dz, ltheta0, dtheta, ltanth0, dtanth, nxb, nzb, 0, (nx-2*nxb)*dx, 0, (nz-2*nzb)*dz);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf(" done \n");

    cg_tanTheta2theta(adcig_tanTheta, adcig_recover, nx, nz, dx, dz, 0, (nx-2*nxb)*dx, 0, (nz-2*nzb)*dz, ltheta0, dtheta, ltanth0, dtanth, niterCG, iproc);
    cg_tanTheta2theta(adcig_tanTheta, adcig_adjoint, nx, nz, dx, dz, 0, (nx-2*nxb)*dx, 0, (nz-2*nzb)*dz, ltheta0, dtheta, ltanth0, dtanth, 1, iproc);

    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Applying focusing operator in the angle-domain ....");
    GPU_WRAP_Demoveout_AD(procGPU, 0, adj, adcig_tanTheta, adcig_TDcontracted    , nx, dx, nz, dz, lx0, ltanth0, dtanth, nxb, nzb, xMin, xMax, zMin, zMax, 1.00-mu, lx0*dx);
    GPU_WRAP_Demoveout_AD(procGPU, 0, adj, adcig_tanTheta, adcig_TDcontracted_ref, nx, dx, nz, dz, lx0, ltanth0, dtanth, nxb, nzb, xMin, xMax, zMin, zMax, 1.00   , lx0*dx);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf(" done \n");

    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Transforming from tangentTheta to theta....");
    // GPU_WRAP_adcig_theta2tanTheta(1, procGPU, 0, adcig_ADcontracted_ref, adcig_TDcontracted_ref, nx, dx, nz, dz, ltheta0, dtheta, ltanth0, dtanth, nxb, nzb, 0, (nx-2*nxb)*dx, 0, (nz-2*nzb)*dz);
    // GPU_WRAP_adcig_theta2tanTheta(1, procGPU, 0, adcig_ADcontracted    , adcig_TDcontracted    , nx, dx, nz, dz, ltheta0, dtheta, ltanth0, dtanth, nxb, nzb, 0, (nx-2*nxb)*dx, 0, (nz-2*nzb)*dz);
    
    cg_tanTheta2theta(adcig_TDcontracted_ref, adcig_ADcontracted_ref, nx, nz, dx, dz, 0, (nx-2*nxb)*dx, 0, (nz-2*nzb)*dz, ltheta0, dtheta, ltanth0, dtanth, niterCG, iproc);
    cg_tanTheta2theta(adcig_TDcontracted    , adcig_ADcontracted    , nx, nz, dx, dz, 0, (nx-2*nxb)*dx, 0, (nz-2*nzb)*dz, ltheta0, dtheta, ltanth0, dtanth, niterCG, iproc);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf(" done \n");

    //*/

    printf("\n\n");
    time_start = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Applying focusing operator in the angle-domain ....");
    // GPU_WRAP_Demoveout_AD(procGPU, 0, adj, adcig, adcig_ADcontracted    , nx, dx, nz, dz, lx0, ltheta0, dtheta, nxb, nzb, xMin, xMax, zMin, zMax, 1.00-mu, 0.0*lx0*dx);
    // GPU_WRAP_Demoveout_AD(procGPU, 0, adj, adcig, adcig_ADcontracted_ref, nx, dx, nz, dz, lx0, ltheta0, dtheta, nxb, nzb, xMin, xMax, zMin, zMax, 1.00   , 0.0*lx0*dx);

    GPU_WRAP_Contraction_theta(procGPU, 0, adcig, adcig_ADcontracted    , nx, dx, nz, dz, ltheta0, dtheta, nxb, nzb, xMin, xMax, zMin, zMax,  mu, 60.0);
    GPU_WRAP_Contraction_theta(procGPU, 0, adcig, adcig_ADcontracted_ref, nx, dx, nz, dz, ltheta0, dtheta, nxb, nzb, xMin, xMax, zMin, zMax, 0.0, 60.0);

    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf(" done in %ld seconds \n\n", time_now-time_start);
    // MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf(" done \n");

    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Recovering ODCIG from contracted ADCIG.... ");
    float *eimage_ADcontracted     = cg_theta2lx_gpu(procGPU, adcig_ADcontracted    , nx, nz, dx, dz, lx0, ntheta, dtheta, otheta, 5);
    float *eimage_ADcontracted_ref = cg_theta2lx_gpu(procGPU, adcig_ADcontracted_ref, nx, nz, dx, dz, lx0, ntheta, dtheta, otheta, 5);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf(" done \n");

    if(outputs)
    {
        char adcig_input_FileName[1024];
        char adcig_recover_FileName[1024];
        char adcig_adjoint_FileName[1024];
        char adcig_ADcontracted_FileName[1024];
        char adcig_ADcontracted_ref_FileName[1024];
        char adcig_TD_FileName[1024];
        char adcig_TDcontracted_FileName[1024];
        char adcig_TDcontracted_ref_FileName[1024];
        char eimage_ADcontracted_FileName[1024];
        char eimage_ADcontracted_ref_FileName[1024];
        char adcig_FileName[1024];
        sprintf(adcig_input_FileName                      , "adcigs_input_Iter%d"           , iter);
        sprintf(adcig_recover_FileName                    , "adcigs_recover_Iter%d"         , iter);
        sprintf(adcig_adjoint_FileName                    , "adcigs_adjoint_Iter%d"         , iter);
        sprintf(eimage_ADcontracted_FileName              , "eimage_ADcontracted_Iter%d"    , iter);
        sprintf(eimage_ADcontracted_ref_FileName          , "eimage_ADcontracted_ref_Iter%d", iter);
        sprintf(adcig_ADcontracted_FileName               , "adcigs_ADcontracted_Iter%d"    , iter);
        sprintf(adcig_ADcontracted_ref_FileName           , "adcigs_ADcontracted_ref_Iter%d", iter);
        sprintf(adcig_TD_FileName                         , "adcigs_TD_Iter%d"              , iter);
        sprintf(adcig_TDcontracted_FileName               , "adcigs_TDcontracted_Iter%d"    , iter);
        sprintf(adcig_TDcontracted_ref_FileName           , "adcigs_TDcontracted_ref_Iter%d", iter);
        outputSmart3d(adcig_input_FileName                , adcig                           , ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
        outputSmart3d(adcig_recover_FileName              , adcig_recover                   , ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
        outputSmart3d(adcig_adjoint_FileName              , adcig_adjoint                   , ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
        outputSmart3d(adcig_ADcontracted_ref_FileName     , adcig_ADcontracted_ref          , ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
        outputSmart3d(adcig_ADcontracted_FileName         , adcig_ADcontracted              , ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
        outputSmart3d(adcig_TD_FileName                   , adcig_tanTheta                  , ntanth, dtanth, -ltanth0*dtanth, nz, dz, oz, nx, dx, ox);
        outputSmart3d(adcig_TDcontracted_ref_FileName     , adcig_TDcontracted_ref          , ntanth, dtanth, -ltanth0*dtanth, nz, dz, oz, nx, dx, ox);
        outputSmart3d(adcig_TDcontracted_FileName         , adcig_TDcontracted              , ntanth, dtanth, -ltanth0*dtanth, nz, dz, oz, nx, dx, ox);
        outputSmart3d(eimage_ADcontracted_ref_FileName    , eimage_ADcontracted_ref         ,    nlx,     dx,         -lx0*dx, nz, dz, oz, nx, dx, ox);
        outputSmart3d(eimage_ADcontracted_FileName        , eimage_ADcontracted             ,    nlx,     dx,         -lx0*dx, nz, dz, oz, nx, dx, ox);
    }
    free(adcig                  );
    free(adcig_recover          );
    free(adcig_adjoint          );
    free(adcig_ADcontracted     );
    free(adcig_ADcontracted_ref );

    time_now = time(NULL);
    if(iproc==0)
        printf("  Concluded in %ld seconds\n", time_now-time_start);

    time_start = time(NULL);
    if(iproc==0)
        printf("\n Applying contraction .....");
    float stretchingFactor = 1.0 + epsFocus;
    if(adj) {
          stretchingFactor = 1.0 - epsFocus;
    }

    GPU_WRAP_Contraction_lx(procGPU, 0, eimage_ADcontracted    , eimage_contracted    , nx, dx, nz, dz, lx0, nxb, nzb, xMin, xMax, zMin, zMax, stretchingFactor);
    GPU_WRAP_Contraction_lx(procGPU, 0, eimage_ADcontracted_ref, eimage_contracted_ref, nx, dx, nz, dz, lx0, nxb, nzb, xMin, xMax, zMin, zMax, 1.00);
    
    time_now = time(NULL);
    if(iproc==0)    printf("  Concluded in %ld seconds \n\n", time_now-time_start);

    float *adcig_contracted     = CPU_zaloc1F(nx*nz*ntheta);
    float *adcig_contracted_ref = CPU_zaloc1F(nx*nz*ntheta);
    float *adcig_diff           = CPU_zaloc1F(nx*nz*ntheta);
    GPU_WRAP_lx2theta(0, procGPU, 0, eimage_contracted    , adcig_contracted    , nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);
    GPU_WRAP_lx2theta(0, procGPU, 0, eimage_contracted_ref, adcig_contracted_ref, nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);

    addArrays(eimage_diff, 1.0, eimage_contracted, -1.0, eimage_contracted_ref, nx*nz*nlx);
    addArrays( adcig_diff, 1.0,  adcig_contracted, -1.0,  adcig_contracted_ref, nx*nz*ntheta);
   
    // Save contracted image to file
    if(outputs)
    {
        char eimage_contracted_FileName[1024];
        char eimage_contracted_ref_FileName[1024];
        char eimage_contracted_tmp_FileName[1024];
        char eimage_contracted_tmp_ref_FileName[1024];
        char adcigs_contracted_FileName[1024];
        char adcigs_contracted_ref_FileName[1024];
        char eimage_diff_FileName[1024];
        char adcigs_diff_FileName[1024];

        sprintf(eimage_contracted_FileName         , "eimage_contracted_Iter%d", iter);
        sprintf(eimage_contracted_tmp_FileName     , "eimage_contracted_tmp_Iter%d", iter);
        sprintf(adcigs_contracted_FileName         , "adcigs_contracted_Iter%d", iter);
        sprintf(eimage_contracted_ref_FileName     , "eimage_contracted_ref_Iter%d", iter);
        sprintf(eimage_contracted_tmp_ref_FileName , "eimage_contracted_tmp_ref_Iter%d", iter);
        sprintf(adcigs_contracted_ref_FileName     , "adcigs_contracted_ref_Iter%d", iter);
        sprintf(eimage_diff_FileName               , "eimage_diff_Iter%d", iter);
        sprintf(adcigs_diff_FileName               , "adcigs_diff_Iter%d", iter);
        
        outputSmart3d(adcigs_contracted_ref_FileName    , adcig_contracted_ref      , ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
        outputSmart3d(adcigs_contracted_FileName        , adcig_contracted          , ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
        outputSmart3d(adcigs_diff_FileName              , adcig_diff                , ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
        outputSmart3d(eimage_contracted_ref_FileName    , eimage_contracted_ref     ,    nlx,     dx,         -lx0*dx, nz, dz, oz, nx, dx, ox);
        outputSmart3d(eimage_contracted_FileName        , eimage_contracted         ,    nlx,     dx,         -lx0*dx, nz, dz, oz, nx, dx, ox);
        outputSmart3d(eimage_contracted_tmp_ref_FileName, eimage_contracted_tmp_ref ,    nlx,     dx,         -lx0*dx, nz, dz, oz, nx, dx, ox);
        outputSmart3d(eimage_contracted_tmp_FileName    , eimage_contracted_tmp     ,    nlx,     dx,         -lx0*dx, nz, dz, oz, nx, dx, ox);
        outputSmart3d(eimage_diff_FileName              , eimage_diff               ,    nlx,     dx,         -lx0*dx, nz, dz, oz, nx, dx, ox);
    }

    free(eimage_ADcontracted    );
    free(eimage_ADcontracted_ref);
    free(adcig_contracted);
    free(adcig_contracted_ref);
    free(adcig_diff);
    if(eimage_contracted_ref==NULL)
        free(eimage_contracted_ref);

return;  
}




void FocusingOperator_Gocig(int adj, float *eimage_diff, float *eimage_contracted, float *eimage_contracted_ref, float *eimage, \
                             float zMin, float zMax, float maxDepth, \
                             long long int nx, long long int nz, float dx, float dz, float ox, float oz, long long int nxb, long long int nzb, \
                             long long int lx0, long long int lz0, long long int ltheta0, float dtheta, float otheta, \
                             float lambdaMax, float dipMin, float dipMax, float ddip, \
                             int iter, float muFocus, float epsFocus, float lambdaNull, \
                             int outputs, int iproc, int nproc) {
        
    long long int nlx    = 2*lx0+1;
    long long int nlz    = 2*lz0+1;
    
    float *eimage_contracted_tmp     = CPU_zaloc1F(nx*nz*nlx*nlz);
    float *eimage_contracted_tmp_ref = CPU_zaloc1F(nx*nz*nlx*nlz);

    float    xMin = 0;
    float    xMax = nx*dx;
    float   shift = 1.00 - muFocus;
    if(adj) shift = 1.00 + muFocus;

    GPU_WRAP_Demoveout_Gocig(procGPU, 0, adj, eimage, eimage_contracted_tmp    , \
                             nx, dx, nz, dz, nxb, nzb, \
                             lx0, lz0, lambdaMax, dipMin, dipMax, ddip, \
                             xMin, xMax, zMin, zMax, shift, lambdaNull);

    GPU_WRAP_Demoveout_Gocig(procGPU, 0, adj, eimage, eimage_contracted_tmp_ref, \
                             nx, dx, nz, dz, nxb, nzb, \
                             lx0, lz0, lambdaMax, dipMin, dipMax, ddip, \
                             xMin, xMax, zMin, zMax, 1.0, lambdaNull);

    float    stretchingFactor = 1.0 + epsFocus;
    if(adj)  stretchingFactor = 1.0 - epsFocus;
    GPU_WRAP_Contraction_lx_lz(procGPU, 0, eimage_contracted_tmp    , eimage_contracted    , \
                               nx, dx, nz, dz, lx0, lz0, \
                               dipMin, dipMax, ddip, lambdaMax, dx, \
                               nxb, nzb, xMin, xMax, zMin, zMax, stretchingFactor);

    GPU_WRAP_Contraction_lx_lz(procGPU, 0, eimage_contracted_tmp_ref, eimage_contracted_ref, \
                               nx, dx, nz, dz, lx0, lz0, \
                               dipMin, dipMax, ddip, lambdaMax, dx, \
                               nxb, nzb, xMin, xMax, zMin, zMax, 1.0);

    addArrays(eimage_diff, 1.0, eimage_contracted, -1.0, eimage_contracted_ref, nx*nz*nlx*nlz);    
   
    free(eimage_contracted_tmp);    
    free(eimage_contracted_tmp_ref);

return;  
}

void FocusingOperator_Tcig(int adj, float *Tcig_foc, float *Tcig_ref, float *Tcig, float zMin, float zMax, \
                           long long int nx, long long int nz, float dx, float dz, float ox, float oz, \
                           long long int nxb, long long int nzb, long long int tau0, \
                           float muFocus, float tauNull, int outputs, int iproc, int nproc) {
        
    long long int ntau = 2*tau0+1; 

    float   xMin = 0;
    float   xMax = nx*dx;
    float   shift = 1.00 - muFocus;
    if(adj) shift = 1.00 + muFocus;

    GPU_WRAP_Contraction_lx(procGPU, 0, Tcig, Tcig_foc, nx, dx, nz, dz, tau0, nxb, nzb, xMin, xMax, zMin, zMax, shift);
    GPU_WRAP_Contraction_lx(procGPU, 0, Tcig, Tcig_ref, nx, dx, nz, dz, tau0, nxb, nzb, xMin, xMax, zMin, zMax,   1.0);

    // GPU_WRAP_Demoveout(procGPU, 0, 0, Tcig, Tcig_foc, nx, dx, nz, dz, tau0, nxb, nzb, xMin, xMax, zMin, zMax, shift, tauNull);
    // GPU_WRAP_Demoveout(procGPU, 0, 0, Tcig, Tcig_ref, nx, dx, nz, dz, tau0, nxb, nzb, xMin, xMax, zMin, zMax,   1.0, tauNull);
    
    // GPU_WRAP_Demoveout_Tcig(procGPU, 0, adj, Tcig, Tcig_foc, nx, dx, nz, dz, tau0, nxb, nzb, xMin, xMax, zMin, zMax, shift, tauNull);
    // GPU_WRAP_Demoveout_Tcig(procGPU, 0, adj, Tcig, Tcig_ref, nx, dx, nz, dz, tau0, nxb, nzb, xMin, xMax, zMin, zMax,   1.0, tauNull);

    MPI_Barrier(MPI_COMM_WORLD);

return;  
}

void FocusingOperator_Docig(int adj, float *eimage_diff, float *eimage_contracted, float *eimage_contracted_ref, float *eimage, \
                            float zMin, float zMax, long long int nx, long long int nz, float dx, float dz, float ox, float oz, \
                            long long int nxb, long long int nzb, long long int lx0, long long int ndip, float dipMin, float ddip, \
                            float muFocus, float epsFocus, float lambdaNull, int outputs, int iproc, int nproc) {
        
    long long int idip;
    long long int nlx = 2*lx0+1;

    float   xMin = 0;
    float   xMax = nx*dx;
    float   shift = 1.00 - muFocus;
    if(adj) shift = 1.00 + muFocus;

    float    stretchingFactor = 1.0 + epsFocus;
    if(adj)  stretchingFactor = 1.0 - epsFocus;

    /*
    for(idip=0; idip<ndip; idip++)
    {
        float *eimage_contracted_tmp     = CPU_zaloc1F(nx*nz*nlx);
        float *eimage_contracted_ref_tmp = CPU_zaloc1F(nx*nz*nlx);

        float dip = dipMin + idip*ddip;
        size_t shiftMem = idip*nx*nz*nlx;

        GPU_WRAP_Demoveout_Docig(procGPU, 0, adj, eimage+shiftMem, eimage_contracted_tmp    , nx, dx, nz, dz, nxb, nzb, lx0, dip, xMin, xMax, zMin, zMax, shift, lambdaNull);
        GPU_WRAP_Demoveout_Docig(procGPU, 0, adj, eimage+shiftMem, eimage_contracted_ref_tmp, nx, dx, nz, dz, nxb, nzb, lx0, dip, xMin, xMax, zMin, zMax,   1.0, lambdaNull);
        
        GPU_WRAP_Contraction_lx(procGPU, 0, eimage_contracted_tmp    , eimage_contracted+shiftMem    , nx, dx, nz, dz, lx0, nxb, nzb, xMin, xMax, zMin, zMax, stretchingFactor);
        GPU_WRAP_Contraction_lx(procGPU, 0, eimage_contracted_ref_tmp, eimage_contracted_ref+shiftMem, nx, dx, nz, dz, lx0, nxb, nzb, xMin, xMax, zMin, zMax, 1.00);   

        free(eimage_contracted_tmp);
        free(eimage_contracted_ref_tmp);
    }
    //*/
    //*
    for(idip=iproc; idip<ndip; idip+=nproc)
    {
        float *eimage_contracted_tmp     = CPU_zaloc1F(nx*nz*nlx);
        float *eimage_contracted_ref_tmp = CPU_zaloc1F(nx*nz*nlx);

        float dip = dipMin + idip*ddip;
        size_t shiftMem = idip*nx*nz*nlx;

        GPU_WRAP_Demoveout_Docig(procGPU, 0, adj, eimage+shiftMem, eimage_contracted_tmp    , nx, dx, nz, dz, nxb, nzb, lx0, dip, xMin, xMax, zMin, zMax, shift, lambdaNull);
        GPU_WRAP_Demoveout_Docig(procGPU, 0, adj, eimage+shiftMem, eimage_contracted_ref_tmp, nx, dx, nz, dz, nxb, nzb, lx0, dip, xMin, xMax, zMin, zMax,   1.0, lambdaNull);
        
        GPU_WRAP_Contraction_lx(procGPU, 0, eimage_contracted_tmp    , eimage_contracted+shiftMem    , nx, dx, nz, dz, lx0, nxb, nzb, xMin, xMax, zMin, zMax, stretchingFactor);
        GPU_WRAP_Contraction_lx(procGPU, 0, eimage_contracted_ref_tmp, eimage_contracted_ref+shiftMem, nx, dx, nz, dz, lx0, nxb, nzb, xMin, xMax, zMin, zMax, 1.00);   

        free(eimage_contracted_tmp);
        free(eimage_contracted_ref_tmp);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // Sending seismograms to manager process
    for(idip=0; idip<ndip; idip++)
    {
        size_t shiftMem = idip*nx*nz*nlx;
        if((idip%nproc)!=0)
        {
            if(iproc==idip%nproc)
                MPI_Send(eimage_contracted_ref+shiftMem, nx*nz*nlx, MPI_FLOAT,      0, 0, MPI_COMM_WORLD);
            if(iproc==0) {
                int sender = idip%nproc;
                MPI_Recv(eimage_contracted_ref+shiftMem, nx*nz*nlx, MPI_FLOAT, sender, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if(iproc==idip%nproc)
                MPI_Send(eimage_contracted    +shiftMem, nx*nz*nlx, MPI_FLOAT,      0, 0, MPI_COMM_WORLD);
            if(iproc==0) {
                int sender = idip%nproc;
                MPI_Recv(eimage_contracted    +shiftMem, nx*nz*nlx, MPI_FLOAT, sender, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(eimage_contracted_ref+shiftMem, nx*nz*nlx, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Bcast(eimage_contracted    +shiftMem, nx*nz*nlx, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    //*/

    addArrays(eimage_diff, 1.0, eimage_contracted, -1.0, eimage_contracted_ref, nx*nz*nlx*ndip);    

return;  
}

void FocusingOperator_TL_Docig(int adj, float *eimage_contracted, float *eimage, \
                               float zMin, float zMax, long long int nx, long long int nz, float dx, float dz, float ox, float oz, \
                               long long int nxb, long long int nzb, long long int lx0, long long int ndip, float dipMin, float ddip, \
                               float muFocus, float epsFocus, float lambdaNull, int outputs, int iproc, int nproc) {
        
    long long int idip;
    long long int nlx = 2*lx0+1;

    float   xMin = 0;
    float   xMax = nx*dx;
    float   shift = 1.00 - muFocus;
    if(adj) shift = 1.00 + muFocus;

    float    stretchingFactor = 1.0 + epsFocus;
    if(adj)  stretchingFactor = 1.0 - epsFocus;

    for(idip=iproc; idip<ndip; idip+=nproc)
    {
        float *eimage_contracted_tmp     = CPU_zaloc1F(nx*nz*nlx);

        float dip = dipMin + idip*ddip;
        size_t shiftMem = idip*nx*nz*nlx;

        GPU_WRAP_Demoveout_Docig(procGPU, 0, adj, eimage+shiftMem, eimage_contracted_tmp    , nx, dx, nz, dz, nxb, nzb, lx0, dip, xMin, xMax, zMin, zMax, shift, lambdaNull);

        GPU_WRAP_Contraction_lx(procGPU, 0, eimage_contracted_tmp    , eimage_contracted+shiftMem    , nx, dx, nz, dz, lx0, nxb, nzb, xMin, xMax, zMin, zMax, stretchingFactor);

        free(eimage_contracted_tmp);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // Sending seismograms to manager process
    for(idip=0; idip<ndip; idip++)
    {
        size_t shiftMem = idip*nx*nz*nlx;
        if((idip%nproc)!=0)
        {
            if(iproc==idip%nproc)
                MPI_Send(eimage_contracted    +shiftMem, nx*nz*nlx, MPI_FLOAT,      0, 0, MPI_COMM_WORLD);
            if(iproc==0) {
                int sender = idip%nproc;
                MPI_Recv(eimage_contracted    +shiftMem, nx*nz*nlx, MPI_FLOAT, sender, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(eimage_contracted    +shiftMem, nx*nz*nlx, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
    }
return;  
}

void FocusingOperator_lambda(int adj, float *eimage_diff, float *eimage_contracted, float *eimage_contracted_ref, float *eimage, \
                             float zMin, float zMax, float maxDepth, long long int nx, long long int nz, float dx, float dz, float ox, float oz, \
                             int nxb, int nzb, int lx0, int ltheta0, float dtheta, float otheta, \
                             int iter, float muFocus, float epsFocus, float lambdaNull, \
                             int outputs, int iproc, int nproc) {
    
    time_t time_now, time_start, time_begin, time_elapsed;
    
    int nlx    = 2*lx0+1;
    int ntheta = 2*ltheta0+1;

    time_start = time(NULL);
    // if(iproc==0)    printf("\n\n Applying Focusing Operator with muFocus=%f, eps=%f, lambdaMax=%f ...... ", muFocus, epsFocus, lx0*dx);
    // printf("\n\n iproc=%d  >>>>  Applying Focusing Operator with muFocus=%f, eps=%f, lambdaMax=%f  eimage_contracted_ref=%p  eimage_contracted=%p ...... ", \
    //        iproc, muFocus, epsFocus, lx0*dx, eimage_contracted_ref, eimage_contracted);
    
    float *eimage_contracted_tmp     = CPU_zaloc1F(nx*nz*nlx);
    float *eimage_contracted_tmp_ref = CPU_zaloc1F(nx*nz*nlx);
    if(eimage_contracted_ref==NULL){
        printf("\n\n ERROR: iproc %d had no memory allocated for array eimage_contracted_ref. Program will abort. \n\n", iproc);
        exit(-1);
    }
    if(eimage_contracted    ==NULL) {
        printf("\n\n ERROR: iproc %d had no memory allocated for array eimage_contracted. Program will abort. \n\n", iproc);
        exit(-1);
    }

    float    xMin = 0;
    float    xMax = nx*dx;
    float   shift = 1.00 - muFocus;
    if(adj) shift = 1.00 + muFocus;
    GPU_WRAP_Demoveout(procGPU, 0, adj, eimage, eimage_contracted_tmp    , nx, dx, nz, dz, lx0, nxb, nzb, xMin, xMax, zMin, zMax, shift, 0*lambdaNull);
    GPU_WRAP_Demoveout(procGPU, 0, adj, eimage, eimage_contracted_tmp_ref, nx, dx, nz, dz, lx0, nxb, nzb, xMin, xMax, zMin, zMax, 1.00 , 0*lambdaNull);
    // GPU_WRAP_Demoveout(procGPU, 0, adj, eimage, eimage_contracted_tmp    , nx, dx, nz, dz, lx0, nxb, nzb, xMin, xMax, zMin, zMax, shift, lx0*dx);
    // GPU_WRAP_Demoveout(procGPU, 0, adj, eimage, eimage_contracted_tmp_ref, nx, dx, nz, dz, lx0, nxb, nzb, xMin, xMax, zMin, zMax, 1.00 , lx0*dx);
    
    time_now = time(NULL);
    // if(iproc==0)    printf("  concluded in %ld seconds\n", time_now-time_start);
    // printf("  concluded in %ld seconds\n", time_now-time_start);

    time_start = time(NULL);
    // if(iproc==0)    printf("\n Applying contraction .....");
    float    stretchingFactor = 1.0 + epsFocus;
    if(adj)  stretchingFactor = 1.0 - epsFocus;
    GPU_WRAP_Contraction_lx(procGPU, 0, eimage_contracted_tmp    , eimage_contracted    , nx, dx, nz, dz, lx0, nxb, nzb, xMin, xMax, zMin, zMax, stretchingFactor);
    GPU_WRAP_Contraction_lx(procGPU, 0, eimage_contracted_tmp_ref, eimage_contracted_ref, nx, dx, nz, dz, lx0, nxb, nzb, xMin, xMax, zMin, zMax, 1.00);   
    time_now = time(NULL);
    // if(iproc==0)    printf("  concluded in %ld seconds \n\n", time_now-time_start);

    addArrays(eimage_diff, 1.0, eimage_contracted, -1.0, eimage_contracted_ref, nx*nz*nlx);    
   
    // Save contracted image to file
    if(outputs)
    {
        float *adcig_contracted     = CPU_zaloc1F(nx*nz*ntheta);
        float *adcig_contracted_ref = CPU_zaloc1F(nx*nz*ntheta);
        float *adcig_diff           = CPU_zaloc1F(nx*nz*ntheta);
        GPU_WRAP_lx2theta(0, procGPU, 0, eimage_contracted    , adcig_contracted    , nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);
        GPU_WRAP_lx2theta(0, procGPU, 0, eimage_contracted_ref, adcig_contracted_ref, nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);
        addArrays( adcig_diff, 1.0,  adcig_contracted, -1.0,  adcig_contracted_ref, nx*nz*ntheta);

        char eimage_contracted_FileName[1024];
        char eimage_contracted_ref_FileName[1024];
        char eimage_contracted_tmp_FileName[1024];
        char eimage_contracted_tmp_ref_FileName[1024];
        char adcigs_contracted_FileName[1024];
        char adcigs_contracted_ref_FileName[1024];
        char eimage_diff_FileName[1024];
        char adcigs_diff_FileName[1024];

        sprintf(eimage_contracted_FileName         , "eimage_contracted_Iter%d", iter);
        sprintf(eimage_contracted_tmp_FileName     , "eimage_contracted_tmp_Iter%d", iter);
        sprintf(adcigs_contracted_FileName         , "adcigs_contracted_Iter%d", iter);
        sprintf(eimage_contracted_ref_FileName     , "eimage_contracted_ref_Iter%d", iter);
        sprintf(eimage_contracted_tmp_ref_FileName , "eimage_contracted_tmp_ref_Iter%d", iter);
        sprintf(adcigs_contracted_ref_FileName     , "adcigs_contracted_ref_Iter%d", iter);
        sprintf(eimage_diff_FileName               , "eimage_diff_Iter%d", iter);
        sprintf(adcigs_diff_FileName               , "adcigs_diff_Iter%d", iter);
        
        outputSmart3d(adcigs_contracted_ref_FileName    , adcig_contracted_ref      , ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
        outputSmart3d(adcigs_contracted_FileName        , adcig_contracted          , ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
        outputSmart3d(adcigs_diff_FileName              , adcig_diff                , ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
        outputSmart3d(eimage_contracted_ref_FileName    , eimage_contracted_ref     ,    nlx,     dx,         -lx0*dx, nz, dz, oz, nx, dx, ox);
        outputSmart3d(eimage_contracted_FileName        , eimage_contracted         ,    nlx,     dx,         -lx0*dx, nz, dz, oz, nx, dx, ox);
        outputSmart3d(eimage_contracted_tmp_ref_FileName, eimage_contracted_tmp_ref ,    nlx,     dx,         -lx0*dx, nz, dz, oz, nx, dx, ox);
        outputSmart3d(eimage_contracted_tmp_FileName    , eimage_contracted_tmp     ,    nlx,     dx,         -lx0*dx, nz, dz, oz, nx, dx, ox);
        outputSmart3d(eimage_diff_FileName              , eimage_diff               ,    nlx,     dx,         -lx0*dx, nz, dz, oz, nx, dx, ox);

        free(adcig_contracted);
        free(adcig_contracted_ref);
        free(adcig_diff);
    }
    free(eimage_contracted_tmp);    
    free(eimage_contracted_tmp_ref);

return;  
}




void FocusingOperator_simple_lx(int adj, float *eimage_contracted, float *eimage_input, float zMin, float zMax, \
                                float maxDepth, long long int nx, long long int nz, float dx, float dz, float ox, float oz, \
                                int lx0, int ltheta0, float dtheta, float otheta, int iter, float mu, float epsFocus, int iproc, int nproc) {
    
    int nlx = 2*lx0+1;
    int ntheta = 2*ltheta0 + 1;

    if(iproc==0)
        printf("\n\n Applying Demoveout  .....");
    float *eimage_contracted_tmp     = CPU_zaloc1F(nx*nz*nlx);

    // float xMin = x1;
    // float xMax = x4;
    float xMin = 0;
    float xMax = (nx-2*nxb)*dx;
    GPU_WRAP_Demoveout(procGPU, 0, adj, eimage_input, eimage_contracted_tmp, 
                       nx, dx, nz, dz, lx0, nxb, nzb, xMin, xMax, zMin, zMax, 1.00-mu, lx0*dx);
    // float stretchingFactor = (1.000/(1.00-mu)) + epsFocus;
    // if(adj) {
    //     stretchingFactor = (1.00-mu) - epsFocus;
    // }
    float stretchingFactor = 1.0 + epsFocus;
    if(adj) {
          stretchingFactor = 1.0 - epsFocus;
    }
    // lx2Contraction_forward(0, eimage_contracted_tmp, eimage_contracted, nx, dx,  nz, lx0, stretchingFactor);
    GPU_WRAP_Contraction_lx(procGPU, 0, eimage_contracted_tmp, eimage_contracted, nx, dx, nz, dz, lx0, \
                            nxb, nzb, xMin, xMax, zMin, zMax, stretchingFactor);

    if(iproc==0)
        printf("  Concluded\n\n");

return;
}








// void FocusingOperator_lx_lz(int adj, float *eimage_diff, float *eimage_contracted, float *eimage_contracted_ref, float *eimage, \
//                             float zMin, float zMax, \
//                             float maxDepth, long long int nx, long long int nz, float dx, float dz, float ox, float oz, \
//                             int lx0, int ltheta0, float dtheta, float otheta, int iter, float mu, float epsFocus, \
//                             int outputs, int iproc, int nproc) {
    
//     time_t time_now, time_start, time_begin, time_elapsed;
    
//     int nlx = 2*lx0+1;
//     int ntheta = 2*ltheta0 + 1;

    
//     float *eimage_contracted_tmp     = CPU_zaloc1F(nx*nz*nlx);
//     float *eimage_contracted_tmp_ref = CPU_zaloc1F(nx*nz*nlx);
    
//     if(eimage_contracted_ref==NULL)
//         eimage_contracted_ref = CPU_zaloc1F(nx*nz*nlx);

//     time_start = time(NULL);
//     if(iproc==0)    printf("\n\n Applying Demoveout .....");

//     float xMin = 0;
//     float xMax = (nx-2*nxb)*dx;
//     GPU_WRAP_Demoveout(procGPU, 0, adj, eimage, eimage_contracted_tmp    , nx, dx, nz, dz, lx0, nxb, nzb, xMin, xMax, zMin, zMax, 1.00-mu, lx0*dx);
//     GPU_WRAP_Demoveout(procGPU, 0, adj, eimage, eimage_contracted_tmp_ref, nx, dx, nz, dz, lx0, nxb, nzb, xMin, xMax, zMin, zMax, 1.00   , lx0*dx);
    
//     time_now = time(NULL);
//     if(iproc==0)    printf("  Concluded in %ld seconds\n", time_now-time_start);

//     time_start = time(NULL);
//     if(iproc==0)    printf("\n Applying contraction .....");

//     float stretchingFactor;
//     if(adj==0)    stretchingFactor = (1.000/(1.00-mu)) + epsFocus;
//     else          stretchingFactor = (1.00-mu) - epsFocus;
    
//     // lx2Contraction_forward(0, eimage_contracted_tmp    , eimage_contracted    , nx, dx,  nz, lx0, stretchingFactor);
//     // lx2Contraction_forward(0, eimage_contracted_tmp_ref, eimage_contracted_ref, nx, dx,  nz, lx0,  1.000          );
//     GPU_WRAP_Contraction_lx(procGPU, 0, eimage_contracted_tmp    , eimage_contracted    , nx, dx, nz, dz, lx0, \
//                             nxb, nzb, xMin, xMax, zMin, zMax, stretchingFactor);
//     GPU_WRAP_Contraction_lx(procGPU, 0, eimage_contracted_tmp_ref, eimage_contracted_ref, nx, dx, nz, dz, lx0, \
//                             nxb, nzb, xMin, xMax, zMin, zMax, 1.00);
    
//     time_now = time(NULL);
//     if(iproc==0)
//         printf("  Concluded in %ld seconds \n\n", time_now-time_start);

//     float *adcig_contracted     = CPU_zaloc1F(nx*nz*ntheta);
//     float *adcig_contracted_ref = CPU_zaloc1F(nx*nz*ntheta);
//     float *adcig_diff           = CPU_zaloc1F(nx*nz*ntheta);
//     // lx2theta_forward(0, eimage_contracted    , adcig_contracted    , nx, dx, nz,  dz, lx0, ntheta, dtheta, otheta);
//     // lx2theta_forward(0, eimage_contracted_ref, adcig_contracted_ref, nx, dx, nz,  dz, lx0, ntheta, dtheta, otheta);
//     GPU_WRAP_lx2theta(0, procGPU, 0, eimage_contracted    , adcig_contracted    , nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);
//     GPU_WRAP_lx2theta(0, procGPU, 0, eimage_contracted_ref, adcig_contracted_ref, nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);

//     addArrays(eimage_diff, 1.0, eimage_contracted, -1.0, eimage_contracted_ref, nx*nz*nlx);
//     addArrays( adcig_diff, 1.0,  adcig_contracted, -1.0,  adcig_contracted_ref, nx*nz*ntheta);
    
    
//     // Save contracted image to file
//     if(outputs)
//     {
//         char eimage_contracted_FileName[1024];
//         char eimage_contracted_ref_FileName[1024];
//         char eimage_contracted_tmp_FileName[1024];
//         char eimage_contracted_tmp_ref_FileName[1024];
//         char adcigs_contracted_FileName[1024];
//         char adcigs_contracted_ref_FileName[1024];
//         char eimage_diff_FileName[1024];
//         char adcigs_diff_FileName[1024];

//         sprintf(eimage_contracted_FileName         , "eimage_contracted_Iter%d", iter);
//         sprintf(eimage_contracted_tmp_FileName     , "eimage_contracted_tmp_Iter%d", iter);
//         sprintf(adcigs_contracted_FileName         , "adcigs_contracted_Iter%d", iter);
//         sprintf(eimage_contracted_ref_FileName     , "eimage_contracted_ref_Iter%d", iter);
//         sprintf(eimage_contracted_tmp_ref_FileName , "eimage_contracted_tmp_ref_Iter%d", iter);
//         sprintf(adcigs_contracted_ref_FileName     , "adcigs_contracted_ref_Iter%d", iter);
//         sprintf(eimage_diff_FileName               , "eimage_diff_Iter%d", iter);
//         sprintf(adcigs_diff_FileName               , "adcigs_diff_Iter%d", iter);
        
//         outputSmart3d(adcigs_contracted_ref_FileName, adcig_contracted_ref  , ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
//         outputSmart3d(adcigs_contracted_FileName    , adcig_contracted      , ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
//         outputSmart3d(adcigs_diff_FileName          , adcig_diff            , ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
//         outputSmart3d(eimage_contracted_ref_FileName, eimage_contracted_ref , nlx, dx, -lx0*dx, nz, dz, oz, nx, dx, ox);
//         outputSmart3d(eimage_contracted_FileName    , eimage_contracted     , nlx, dx, -lx0*dx, nz, dz, oz, nx, dx, ox);
//         outputSmart3d(eimage_contracted_tmp_ref_FileName, eimage_contracted_tmp_ref , nlx, dx, -lx0*dx, nz, dz, oz, nx, dx, ox);
//         outputSmart3d(eimage_contracted_tmp_FileName, eimage_contracted_tmp     , nlx, dx, -lx0*dx, nz, dz, oz, nx, dx, ox);
//         outputSmart3d(eimage_diff_FileName          , eimage_diff           , nlx, dx, -lx0*dx, nz, dz, oz, nx, dx, ox);
//     }

//     free(eimage_contracted_tmp);
//     free(eimage_contracted_tmp_ref);
//     free(adcig_contracted);
//     free(adcig_contracted_ref);
//     free(adcig_diff);
//     if(eimage_contracted_ref==NULL)
//         free(eimage_contracted_ref);

// return;  
// }
