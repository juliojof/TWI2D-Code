float computeResiduals(Dataset *dataset_input, Dataset *dataset_residuals, Dataset *dataset_reference, Dataset *dataset_misfit, \
                       float *eimage_contracted, float *eimage_contracted_ref, float *eimage_diff_diff, float *eimage_residuals_contracted, 
                       MSH *msh, float *vp_bg, EIProc *eiproc, Eimage_MSH *eimage_msh, float offMax, int acquisitionGeom, float perc, \
                       float zMin, float zMax, float maxDepth, float muFocus, float epsFocus, float lambdaNull, int iter, int ism, int outputs, \
                       int imageDomain, int iproc, int nproc) {

    int nxb         = msh->nxb;
    int nzb         = msh->nzb;
    int nx          = msh->nx;
    int nz          = msh->nz;
    int nt          = msh->nt;
    float dx        = msh->dx;
    float dz        = msh->dz;
    float dt        = msh->dt;
    float ox        = msh->ox;
    float oz        = msh->oz;
    float ot        = msh->ot;
    int jt          = msh->jt;
    int jt_EIC      = msh->jt_EIC;
    int jt_lt       = msh->jt_lt;
    int ltheta0     = eimage_msh->ltheta0;
    int ntheta      = eimage_msh->ntheta;
    float dtheta    = eimage_msh->dtheta;
    float otheta    = eimage_msh->otheta;
    int lx0         = eimage_msh->lx0;
    int lz0         = eimage_msh->lz0;
    int lv0         = eimage_msh->lv0;
    int lt0         = eimage_msh->lt0;
    int nlx         = eimage_msh->nlx;
    int nlz         = eimage_msh->nlz;
    int nlv         = eimage_msh->nlv;
    int nlt         = eimage_msh->nlt;
    int lambda0     = eimage_msh->nlambda/2;
    int   nxLD      = eimage_msh->nx_LambdaDip;
    int   nzLD      = eimage_msh->nz_LambdaDip;
    float oxLD      = eimage_msh->ox_LambdaDip;
    float ozLD      = eimage_msh->oz_LambdaDip;
    int   nlambda   = eimage_msh->nlambda;
    float lambdaMax = eimage_msh->lambdaMax;
    float dlambda   = eimage_msh->dlambda;
    float dipMin    = eimage_msh->dipMin;
    float dipMax    = eimage_msh->dipMax;
    float ddip      = eimage_msh->ddip;
    int   ndip      = eimage_msh->ndip;
    

    char refDataFileName[1024];
    char misDataFileName[1024];
    char resDataFileName[1024];
    sprintf( refDataFileName , "seismogram_referenceData_Iter%d", iter);
    sprintf( misDataFileName , "seismogram_misfitData_Iter%d", iter);
    sprintf( resDataFileName , "seismogram_residuals_Iter%d", iter);
    
    float objFunc = 0.0;

    zeroArray(eimage_contracted, nx*nz*nlx*nlz*nlt);
    if(imageDomain)
        zeroArray(eimage_diff_diff , nx*nz*nlx*nlz*nlt);
    else if(eimage_residuals_contracted!=NULL)
        zeroArray(eimage_residuals_contracted, nx*nz*nlx*nlz*nlt);

    if(iproc==0)
        printf("\n\n  nlx=%d  nlz=%d  nlt=%d  nlv=%d  jt_lt=%d  ntheta=%d  \n\n", nlx, nlz, nlt, nlv, jt_lt, ntheta);

    float  *image, *eimage, *eimage_2B_Focused, *eimage_taper_angle, *eimage_residuals_taper_angle;
    float *eimage_diff_contracted;
    float *eimage_diff;
    image   = CPU_zaloc1F(nx*nz);
    eimage  = CPU_zaloc1F(nx*nz*nlx*nlz*nlt);

    // >>>>>>>>>>> Migrate shots <<<<<<<<<<<
    if(lv0==0)
    {
        float *eimage_theta             = CPU_zaloc1F(nx*nz*ntheta);
        float *eimage_theta_filterECIG  = CPU_zaloc1F(nx*nz*ntheta);
        float *eimage_theta_taper_angle = CPU_zaloc1F(nx*nz*ntheta);

        MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n\n Migrating data  lv0=%d ...... ", lv0);
        migration(dataset_input, image, eimage, vp_bg, nx, nz, lx0, lz0, lt0, jt_lt, ltheta0, dx, dz, \
                  dtheta, ox, oz, otheta, jt, ism, 0, iproc, nproc);
        MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n\n Migrating data  lv0=%d ...... concluded", lv0);

    /*
    if(outputs)
    {
        outputModel2d("image0_debug", image , nz, dz, 0.0, nx, dx, 0.0, nzb, nxb);
        outputSmart3d("eimage_debug", eimage, nlx, dx, -lx0*dx,  nz, dz, oz,  nx, dx, ox);        
    }
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n\n STOPPING PROGRAM ABRUPTLY - DEBUG \n\n");
    if(1)    exit(-1);
    //*/

        GPU_WRAP_lx2theta(0, procGPU, 0, eimage, eimage_theta           , nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);

        if(outputs) {
            // outputSmart2d(image_FileName            ,  image                                                  , nz, dz, oz, nx, dx, ox);
            outputModel2d(image_FileName            ,  image                                                  , nz, dz, 0.0, nx, dx, 0.0, nzb, nxb);
            outputSmart3d(eimage_FileName           , eimage                 , nlx   , dx    , -lx0*dx        , nz, dz, oz, nx, dx, ox);
            outputSmart3d(adcigs_FileName           , eimage_theta           , ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
        }

        filterECIG(eimage, eiproc, eimage_msh, msh);
        GPU_WRAP_lx2theta(0, procGPU, 0, eimage, eimage_theta_filterECIG, nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);

        if(outputs) {
            outputSmart3d(eimage_filterECIG_FileName, eimage                 , nlx   , dx    , -lx0*dx        , nz, dz, oz, nx, dx, ox);
            outputSmart3d(adcigs_filterECIG_FileName, eimage_theta_filterECIG, ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
        }

        eimage_taper_angle = filterECIG_Angle(eimage, eiproc, nx, nz, nxb, nzb, \
                                              dx, dz, ox, oz, ltheta0, dtheta, lx0, iproc, procGPU);

        
        if(eimage_taper_angle!=NULL)  {

            GPU_WRAP_lx2theta(0, procGPU, 0, eimage_taper_angle, eimage_theta_taper_angle, nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);
            if(outputs) {
                outputSmart3d(eimage_taper_angle_FileName, eimage_taper_angle, nlx , dx    , -lx0*dx        , nz, dz, oz, nx, dx, ox);
                outputSmart3d(adcigs_taper_angle_FileName, eimage_theta_taper_angle, ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
            }
            memcpy(eimage, eimage_taper_angle, nx*nz*nlx*sizeof(float));
            free(eimage_taper_angle);
        }
        
        free(eimage_theta);
        free(eimage_theta_filterECIG);
        free(eimage_theta_taper_angle);

        //  >>>>> FOCUSING IMAGE <<<<<<
        eimage_diff = CPU_zaloc1F(nx*nz*nlx*nlz);
        MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Applying Focusing Operator ...... ");
        FocusingOperator_lx_lz_lt(0, eimage_diff, eimage_contracted, eimage_contracted_ref, eimage, zMin, zMax, maxDepth,  \
                                  nx, nz, dx, dz, ox, oz, lx0, lz0, lt0, ltheta0, dtheta, otheta, \
                                  dipMin, dipMax, ddip, lambdaMax, dlambda, eiproc->applyADFocusing, \
                                  iter, muFocus, epsFocus, lambdaNull, outputs, iproc, nproc);
        MPI_Barrier(MPI_COMM_WORLD);  if(iproc==0)    printf("\n done ");

        MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Filtering focused image ...... ");
        // filterECIG(eimage_diff          , eiproc, eimage_msh, msh);
        // filterECIG(eimage_contracted    , eiproc, eimage_msh, msh);
        // filterECIG(eimage_contracted_ref, eiproc, eimage_msh, msh);
        MPI_Barrier(MPI_COMM_WORLD);  if(iproc==0)    printf("\n done ");
    }
    else
    {
        MPI_Barrier(MPI_COMM_WORLD);  
        if(iproc==0)    printf("\n\n CODE SHOULD NOT PASS HERE -> PROGRAM WILL ABORT ... Goodbye! \n\n"); 
        MPI_Barrier(MPI_COMM_WORLD);
        exit(-1);
    }

    // applyLinearWeight_ODCIG(eimage_diff          , eimage_diff      , dx, nx, nz, lx0);
    // applyLinearWeight_ODCIG(eimage_contracted    , eimage_contracted, dx, nx, nz, lx0);
    // applyLinearWeight_ODCIG(eimage_contracted_ref, eimage_contracted, dx, nx, nz, lx0);
    MPI_Barrier(MPI_COMM_WORLD);

    if(dataset_reference != NULL)
    {
        if(iproc==0)    printf("\n Modeling reference data \n");
        extendedBornModeling(dataset_reference, eimage_contracted_ref, NULL, NULL, vp_bg, nx, nz, lx0, lz0, lt0, \
                             dx, dz, ox, oz, jt, perc, offMax, acquisitionGeom, NULL, iproc, nproc);
        int jshot = 10;
        writeSeisData2File(dataset_reference, jshot, iproc, refDataFileName, NULL);
        if(iproc==0)    printf("\n Modeling reference data  ----->>>>> concluded\n");
    }
    if(dataset_misfit    != NULL)
    {
        if(iproc==0)    printf("\n Modeling misfit data \n");
        extendedBornModeling(dataset_misfit, eimage_contracted, NULL, NULL, vp_bg, nx, nz, lx0, lz0, lt0, \
                            dx, dz, ox, oz, jt, perc, offMax, acquisitionGeom, NULL, iproc, nproc);
        int jshot = 10;
        writeSeisData2File(dataset_misfit   , jshot, iproc, misDataFileName, NULL);
        if(iproc==0)    printf("\n Modeling misfit data  ----->>>>> concluded\n");
    }

    //  >>>>> IMAGE RESIDUALS <<<<<<
    if(imageDomain)
    {
        eimage_diff_contracted         = CPU_zaloc1F(nx*nz*nlx*nlz);
        FocusingOperator_lx_lz_lt(1, eimage_diff_diff, eimage_diff_contracted, NULL, eimage_diff, zMin, zMax, maxDepth, 
                                  nx, nz, dx, dz, ox, oz, lx0, lz0, lt0, ltheta0, dtheta, otheta, \
                                  eimage_msh->dipMin, eimage_msh->dipMax, eimage_msh->ddip, eimage_msh->lambdaMax, eimage_msh->dlambda, \
                                  eiproc->applyADFocusing, iter, muFocus, epsFocus, lambdaNull, outputs, iproc, nproc);
        if(eiproc->applyBandass2EImage)
        {
            lxApplyBandpass(eimage_diff_diff       , nx, nz, dz, lx0, eiproc->lambda1, eiproc->lambda2, eiproc->lambda3, eiproc->lambda4);
            lxApplyBandpass(eimage_diff_contracted , nx, nz, dz, lx0, eiproc->lambda1, eiproc->lambda2, eiproc->lambda3, eiproc->lambda4);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        lxTapering_lx_varDepth(eimage_diff_diff, eimage_diff_diff, nx, nz, dx, lx0, \
                               eiproc->lx1_i, eiproc->lx2_i, eiproc->lx3_i, eiproc->lx4_i, \
                               eiproc->lx1_f, eiproc->lx2_f, eiproc->lx3_f, eiproc->lx4_f, eiproc->z_i, eiproc->z_f);
        objFunc = computeL2NormArray(eimage_diff, nx*nz*nlx*nlz);
    }
    MPI_Barrier(MPI_COMM_WORLD);


    //  >>>>> DATA RESIDUALS <<<<<<  
    if(imageDomain==0)
    {
        // lxNullTapering(eimage_diff, eimage_diff, nx, nz, dx, lx0, lambdaNull);
        if(iproc==0)    printf("\n\n Computing data residuals. \n (Shot,time): \n");
        objFunc = extendedBornModeling(dataset_residuals, eimage_diff, NULL, NULL, vp_bg, nx, nz, lx0, lz0, lt0, \
                                       dx, dz, ox, oz, jt, perc, offMax, acquisitionGeom, NULL, iproc, nproc);
        MPI_Barrier(MPI_COMM_WORLD);

        // // // // Second-term migrated residuals - begin // // // //
        if(eimage_residuals_contracted != NULL)
        {
            float  *image_residuals, *eimage_residuals, *eimage_residuals_theta;
            image_residuals             = CPU_zaloc1F(nx*nz);
            eimage_residuals            = CPU_zaloc1F(nx*nz*nlx*nlz);
            eimage_residuals_theta      = CPU_zaloc1F(nx*nz*ntheta);

            // Migrate shots
            if(iproc==0)    printf("\n\n Migrating data residuals - second term. \n");
            migration(dataset_residuals, image_residuals, eimage_residuals, vp_bg, nx, nz, lx0, lz0, lt0, jt_lt, \
                      ltheta0, dx, dz, dtheta, ox, oz, otheta, jt, ism, 0, iproc, nproc);
            if(iproc==0)    printf("\n Migrating data residuals - second term -->> concluded. \n");
            if(outputs && iter==0)
            {
                // outputSmart2d("image_residuals",   image_residuals, nz, dz, oz, nx, dx, ox);
                outputModel2d("image_residuals",   image_residuals, nz, dz, 0.0, nx, dx, 0.0, nzb, nxb);
                outputSmart3d("eimage_residuals", eimage_residuals, nlx, dx, -lx0*dx, nz, dz, oz, nx, dx, ox);
            }
            MPI_Barrier(MPI_COMM_WORLD);

            if(eiproc->applyXTaper)
                lxTapering_x(eimage_residuals, eimage_residuals, nx, nz, dx, lx0, eiproc->x1, eiproc->x2, eiproc->x3, eiproc->x4);
            if(eiproc->applyZTaper) {
                lxTapering(eimage_residuals, eimage_residuals, nx, nz, dz, lx0, lz0, eiproc->z0, eiproc->z1);
                lxTapering(eimage_residuals, eimage_residuals, nx, nz, dz, lx0, lz0, eiproc->z3, eiproc->z2);
            }
            if(eiproc->applyLXTaper)
                lxTapering_lx_varDepth(eimage_residuals, eimage_residuals, nx, nz, dx, lx0, \
                                       eiproc->lx1_i, eiproc->lx2_i, eiproc->lx3_i, eiproc->lx4_i, \
                                       eiproc->lx1_f, eiproc->lx2_f, eiproc->lx3_f, eiproc->lx4_f, eiproc->z_i, eiproc->z_f);

            if(eiproc->applyTaperAngle)
            {
                GPU_WRAP_lx2theta(0, procGPU, 0, eimage_residuals, eimage_residuals_theta, nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);

                // Applying cosine to ADCIG (Jacobian of the Offset to Angle transformation)
                // applyCosADCIG(eimage_residuals_theta, eimage_residuals_theta, nx, nz, ntheta, dtheta);

                // Tapering angles
                taperADCIG(eimage_residuals_theta, eimage_residuals_theta, nx, nz, ntheta, dtheta, eiproc->theta2, eiproc->theta1, 1.0);
                taperADCIG(eimage_residuals_theta, eimage_residuals_theta, nx, nz, ntheta, dtheta, eiproc->theta3, eiproc->theta4, 1.0);
                taperADCIG(eimage_residuals_theta, eimage_residuals_theta, nx, nz, ntheta, dtheta, eiproc->theta4, +dtheta*ltheta0, 0.0);
                taperADCIG(eimage_residuals_theta, eimage_residuals_theta, nx, nz, ntheta, dtheta, eiproc->theta1, -dtheta*ltheta0, 0.0);
                if(iproc==0)
                    printf("\n Recovering eimage_residuals_theta after angle tapering");

                MPI_Barrier(MPI_COMM_WORLD);
                eimage_residuals_taper_angle = cg_theta2lx_gpu(procGPU, eimage_residuals_theta, nx, nz, dx, dz, lx0, ntheta, dtheta, otheta, 5);
                MPI_Barrier(MPI_COMM_WORLD);
                eimage_2B_Focused = eimage_residuals_taper_angle;
            }
            else {
                eimage_2B_Focused = eimage_residuals;
            }


            FocusingOperator_simple_lx(1, eimage_residuals_contracted, eimage_2B_Focused, zMin, zMax, maxDepth, \
                                       nx, nz, dx, dz, ox, oz, lx0, ltheta0, dtheta, otheta, iter, muFocus, epsFocus, iproc, nproc);
            if(eiproc->applyBandass2EImage) {
                lxApplyBandpass(eimage_residuals_contracted , nx, nz, dz, lx0, eiproc->lambda1, eiproc->lambda2, eiproc->lambda3, eiproc->lambda4);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            // applyLinearWeight_ODCIG(eimage_residuals_contracted, eimage_residuals_contracted, dx, nx, nz, lx0);

            free(image_residuals);
            free(eimage_residuals);
            free(eimage_residuals_theta);
            if(eiproc->applyTaperAngle)    free(eimage_residuals_taper_angle);
            
            if(outputs && iter==0)
            {
                outputSmart3d("eimage_residuals_contracted", eimage_residuals, nlx, dx, -lx0*dx, nz, dz, oz, nx, dx, ox);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }


    }

    free(image);
    free(eimage);
    free(eimage_diff);
    if(imageDomain)
    {
        free(eimage_diff_contracted); 
    }
    

return objFunc;
}
float computeResiduals_AD(Dataset *dataset_input, Dataset *dataset_residuals, Dataset *dataset_reference, Dataset *dataset_misfit, \
                          float *Hocig_contracted, float *Hocig_contracted_ref, float *Vocig_contracted, float *Vocig_contracted_ref, \
                          float *Hocig_diff, float *Vocig_diff, \
                          MSH *msh, float *vp_bg, EIProc *eiproc, Eimage_MSH *eimage_msh, float offMax, int acquisitionGeom, float perc, \
                          float zMin, float zMax, float maxDepth, float muFocus, float epsFocus, float lambdaNull, int iter, int ism, int outputs, \
                          int imageDomain, int computingGradient, int iproc, int nproc) {

    long long int nxb         = msh->nxb;
    long long int nzb         = msh->nzb;
    long long int nx          = msh->nx;
    long long int nz          = msh->nz;
    long long int nt          = msh->nt;
    float dx                  = msh->dx;
    float dz                  = msh->dz;
    float dt                  = msh->dt;
    float ox                  = msh->ox;
    float oz                  = msh->oz;
    float ot                  = msh->ot;
    long long int jt          = msh->jt;
    long long int jt_EIC      = msh->jt_EIC;
    long long int jt_lt       = msh->jt_lt;
    long long int ltheta0     = eimage_msh->ltheta0;
    long long int ntheta      = eimage_msh->ntheta;
    float dtheta              = eimage_msh->dtheta;
    float otheta              = eimage_msh->otheta;
    long long int lt0         = eimage_msh->lt0;
    long long int lx0         = eimage_msh->lx0;
    long long int lv0         = eimage_msh->lv0;
    long long int nlx         = eimage_msh->nlx;
    long long int nlv         = eimage_msh->nlv;
  
    char refDataFileName[1024];
    char misDataFileName[1024];
    char resDataFileName[1024];
    sprintf( refDataFileName , "seismogram_referenceData_Iter%d", iter);
    sprintf( misDataFileName , "seismogram_misfitData_Iter%d", iter);
    sprintf( resDataFileName , "seismogram_residuals_Iter%d", iter);

    time_t time_start, time_now;

    if(iproc==0)    printf("\n\n  nlx=%lld  nlv=%lld  jt_lt=%lld  ntheta=%lld  \n\n", nlx, nlv, jt_lt, ntheta);

    float  *image = CPU_zaloc1F(nx*nz);

    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n\n Migrating data H&V OCIGS  lx0=%lld  lv0=%lld ...... ", lx0, lv0);
    time_start = time(NULL);
    float *Hocig = CPU_zaloc1F(nx*nz*nlx);
    float *Vocig = CPU_zaloc1F(nx*nz*nlv);
    migration_HOCIG_VOCIG(dataset_input, image, Hocig, Vocig, vp_bg, nx, nz, \
                          lx0, lv0, ltheta0, dx, dz, dtheta, ox, oz, \
                          otheta, jt, ism, 1, iproc, nproc);
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

    //  H&V CIG filterings  
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n\n Filtering H&V CIGS ...... ");
    time_start = time(NULL);
    filterECIG(Hocig, eiproc, eimage_msh, msh);
    eimage_msh->lx0 = lv0;
    filterECIG(Vocig, eiproc, eimage_msh, msh);
    eimage_msh->lx0 = lx0;
    float *Hocig_taper_angle = filterECIG_Angle(Hocig, eiproc, nx, nz, nxb, nzb, dx, dz, ox, oz, \
                                                ltheta0, dtheta, lx0, iproc, procGPU);
    float *Vocig_taper_angle = filterECIG_Angle(Vocig, eiproc, nx, nz, nxb, nzb, dx, dz, ox, oz, \
                                                ltheta0, dtheta, lv0, iproc, procGPU);
    if(Hocig_taper_angle!=NULL)
    {
        memcpy(Hocig, Hocig_taper_angle, nx*nz*nlx*sizeof(float));
        free(Hocig_taper_angle);
        memcpy(Vocig, Vocig_taper_angle, nx*nz*nlv*sizeof(float));
        free(Vocig_taper_angle);
    }
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

    // Saving adcig
    float *Hacig = CPU_zaloc1F(nx*nz*ntheta);
    GPU_WRAP_lx2theta(0, procGPU, 0, Hocig, Hacig, nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);

    // Applying focusing operator
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Applying focusing operator to Hocig ...... "); time_start = time(NULL);
    FocusingOperator_AD(0, Hocig_diff, Hocig_contracted, Hocig_contracted_ref, Hocig, zMin, zMax, maxDepth, \
                        nx, nz, dx, dz, ox, oz, lx0, ltheta0, dtheta, otheta, iter, muFocus, epsFocus, lambdaNull, outputs, iproc, nproc);
    MPI_Barrier(MPI_COMM_WORLD);    time_now = time(NULL); if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);


    // Saving Misfit and Reference datasets
    if(dataset_reference != NULL)
    {
        if(iproc==0)    printf("\n Modeling reference data \n");
        extendedBornModeling(dataset_reference, NULL, Hocig_contracted_ref, Vocig_contracted_ref, vp_bg, nx, nz, lx0, lv0, lt0, \
                             dx, dz, ox, oz, jt, perc, offMax, acquisitionGeom, NULL, iproc, nproc);

        int jshot = 10;
        writeSeisData2File(dataset_reference, jshot, iproc, refDataFileName, NULL);
        if(iproc==0)    printf("\n Modeling reference data  ----->>>>> concluded\n");
    }
    if(dataset_misfit    != NULL)
    {
        if(iproc==0)    printf("\n Modeling misfit data \n");
        extendedBornModeling(dataset_misfit, NULL, Hocig_contracted, Vocig_contracted, vp_bg, nx, nz, lx0, lv0, lt0, \
                             dx, dz, ox, oz, jt, perc, offMax, acquisitionGeom, NULL, iproc, nproc);
        int jshot = 10;
        writeSeisData2File(dataset_misfit   , jshot, iproc, misDataFileName, NULL);
        if(iproc==0)    printf("\n Modeling misfit data  ----->>>>> concluded\n");
    }

    //  Computing data residuals
    lxNullTapering(Hocig_diff, Hocig_diff, nx, nz, dx, lx0, lambdaNull);
    double objFuncH = getArrayL2Norm(Hocig_diff, nx*nz*nlx);
    double objFuncV = getArrayL2Norm(Vocig_diff, nx*nz*nlv);
    float objFunc = objFuncH + objFuncV;
    if(iproc==0)    printf("\n\n objFuncH=%g  objFuncV=%g  objFunc=%g \n\n", objFuncH, objFuncV, objFunc);
    multiplyArrayByScalar(0.0, Vocig_diff, nx*nz*nlv);
    if(computingGradient)
    {
        if(iproc==0)    printf("\n\n Computing data residuals. \n (Shot,time): \n");
        extendedBornModeling(dataset_residuals, NULL, Hocig_diff, Vocig_diff, vp_bg, nx, nz, lx0, lv0, lt0, \
                             dx, dz, ox, oz, jt, perc, offMax, acquisitionGeom, NULL, iproc, nproc);
    }
    // if(iproc==0)    printf("\n\n Computing data residuals. \n (Shot,time): \n");
    // objFunc = extendedBornModeling(dataset_residuals, NULL, Hocig_diff, Vocig_diff, vp_bg, nx, nz, lx0, lv0, lt0, \
    //                                dx, dz, ox, oz, jt, perc, offMax, acquisitionGeom, NULL, iproc, nproc);
    MPI_Barrier(MPI_COMM_WORLD);

    free(image);
    

return objFunc;
}
//*
float computeResiduals_TimeLag(Dataset *dataset_input, Dataset *dataset_residuals, Dataset *dataset_reference, Dataset *dataset_misfit, \
                               float *TLcig_contracted, float *TLcig_contracted_ref, float *TLcig_normalizer, \
                               MSH *msh, float *vp_bg, EIProc *eiproc, Eimage_MSH *eimage_msh, float offMax, int acquisitionGeom, float perc, \
                               float zMin, float zMax, float maxDepth, float muFocus, float epsFocus, float lambdaNull, int iter, int ism, int outputs, \
                               int imageDomain, int computingGradient, int iproc, int nproc) {

    float         dip;
    long long int idip;
    long long int nxb         = msh->nxb;
    long long int nzb         = msh->nzb;
    long long int nx          = msh->nx;
    long long int nz          = msh->nz;
    long long int nt          = msh->nt;
    float dx                  = msh->dx;
    float dz                  = msh->dz;
    float dt                  = msh->dt;
    float ox                  = msh->ox;
    float oz                  = msh->oz;
    float ot                  = msh->ot;
    long long int jt          = msh->jt;
    long long int jt_EIC      = msh->jt_EIC;
    long long int jt_lt       = msh->jt_lt;
    long long int ltheta0     = eimage_msh->ltheta0;
    long long int ntheta      = eimage_msh->ntheta;
    float dtheta              = eimage_msh->dtheta;
    float otheta              = eimage_msh->otheta;
    long long int lx0         = eimage_msh->lx0;
    long long int lz0         = eimage_msh->lz0;
    long long int lv0         = eimage_msh->lv0;
    long long int lt0         = eimage_msh->lt0;
    long long int lambda0     = eimage_msh->nlambda/2;
    long long int nlx         = eimage_msh->nlx;
    long long int nlz         = eimage_msh->nlz;
    long long int nlv         = eimage_msh->nlv;
    long long int nlt         = eimage_msh->nlt;
    long long int nlambda     = eimage_msh->nlambda;
    float lambdaMax           = eimage_msh->lambdaMax;
    float dlambda             = eimage_msh->dlambda;
    

    char refDataFileName[1024];
    char misDataFileName[1024];
    char resDataFileName[1024];
    sprintf( refDataFileName , "seismogram_referenceData_Iter%d", iter);
    sprintf( misDataFileName , "seismogram_misfitData_Iter%d", iter);
    sprintf( resDataFileName , "seismogram_residuals_Iter%d", iter);

    time_t time_start, time_now;

    if(iproc==0)
        printf("\n\n  nlx=%lld  nlz=%lld  nlt=%lld  nlv=%lld  jt_lt=%lld  ntheta=%lld  \n\n", nlx, nlz, nlt, nlv, jt_lt, ntheta);

    float  *image = CPU_zaloc1F(nx*nz);

    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n\n Migrating data H&V OCIGS  lv0=%lld ...... ", lv0);
    time_start = time(NULL);
    float *TLcig = CPU_zaloc1F(nx*nz*nlt);    
    migration_TLCIG(dataset_input, image, TLcig, vp_bg, nx, nz, dx, dz, ox, oz, lt0, jt, ism, 0, iproc, nproc);
    
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

    float objFunc = 0.0;
    
return objFunc;
}
//*/


///////////// Residuals with Popping CIGs - begin ////////////////////////
float computeResiduals_ORTH_PopCIG(Dataset *dataset_input, \
                                   float *Hocig_background_phaseShift, float *Hocig_background_original, \
                                   float *Hocig_residual_phaseShift, float *Hocig_residual_original, \
                                   float *Vocig_background_phaseShift, float *Vocig_background_original, \
                                   float *Vocig_residual_phaseShift, float *Vocig_residual_original, \
                                   MSH *msh, float *vp_bg, EIProc *eiproc, Eimage_MSH *eimage_msh, \
                                   float offMax, int acquisitionGeom, float perc, float zMin, float zMax, float maxDepth, \
                                   float muFocus, float epsFocus, float lambdaNull, int iter, int ism, int outputs, \
                                   int computingGradient, int focusingQC, int iproc, int nproc) {

    float         dip;
    long long int idip;
    long long int nxb         = msh->nxb;
    long long int nzb         = msh->nzb;
    long long int nx          = msh->nx;
    long long int nz          = msh->nz;
    long long int nt          = msh->nt;
    float dx                  = msh->dx;
    float dz                  = msh->dz;
    float dt                  = msh->dt;
    float ox                  = msh->ox;
    float oz                  = msh->oz;
    float ot                  = msh->ot;
    long long int jt          = msh->jt;
    long long int jt_EIC      = msh->jt_EIC;
    long long int jt_lt       = msh->jt_lt;
    long long int ltheta0     = eimage_msh->ltheta0;
    long long int ntheta      = eimage_msh->ntheta;
    float dtheta              = eimage_msh->dtheta;
    float otheta              = eimage_msh->otheta;
    long long int lx0         = eimage_msh->lx0;
    long long int lz0         = eimage_msh->lz0;
    long long int lv0         = eimage_msh->lv0;
    long long int lambda0     = eimage_msh->nlambda/2;
    long long int nlx         = eimage_msh->nlx;
    long long int nlz         = eimage_msh->nlz;
    long long int nlv         = eimage_msh->nlv;
    long long int nlambda     = eimage_msh->nlambda;
    float lambdaMax           = eimage_msh->lambdaMax;
    float dlambda             = eimage_msh->dlambda;
    float dipMin              = eimage_msh->dipMin;
    float dipMax              = eimage_msh->dipMax;
    float ddip                = eimage_msh->ddip;
    long long int ndip        = eimage_msh->ndip;

    int niterCG = 8;  // number of CG iterations to transform back from Docig to HVocig

    time_t time_start, time_now;

    float  *image = CPU_zaloc1F(nx*nz);

    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n\n Migrating data H&V OCIGS: nlx=%lld  nlz=%lld  nlv=%lld  jt_lt=%lld  ntheta=%lld  lv0=%lld ...... ", \
                                                            nlx, nlz, nlv, jt_lt, ntheta, lv0);
    time_start = time(NULL);
    float *Hocig = CPU_zaloc1F(nx*nz*nlx);
    float *Vocig = CPU_zaloc1F(nx*nz*nlv);
    migration_HOCIG_VOCIG(dataset_input, image, Hocig, Vocig, vp_bg, nx, nz, \
                          lx0, lv0, ltheta0, dx, dz, dtheta, ox, oz, \
                          otheta, jt, ism, 1, iproc, nproc);
    
    
    reciprocityAngle(Hocig, eiproc, acquisitionGeom, nx, nz, nxb, nzb, dx, dz, ox, oz, ltheta0, dtheta, lx0, iproc, procGPU);
    if(EIMSH_CIGS_VCIG) {
        reciprocityAngle(Vocig, eiproc, acquisitionGeom, nx, nz, nxb, nzb, dx, dz, ox, oz, ltheta0, dtheta, lv0, iproc, procGPU);
    }
    
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

    //  H&V CIG filterings  
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n\n Filtering H&V CIGS ...... ");
    time_start = time(NULL);
    filterECIG(Hocig, eiproc, eimage_msh, msh);
    if(EIMSH_CIGS_VCIG) {
        eimage_msh->lx0 = lv0;
        filterECIG(Vocig, eiproc, eimage_msh, msh);
        eimage_msh->lx0 = lx0;
    }

    float *Hocig_taper_angle = filterECIG_Angle(Hocig, eiproc, nx, nz, nxb, nzb, dx, dz, ox, oz, \
                                                ltheta0, dtheta, lx0, iproc, procGPU);
    float *Vocig_taper_angle = NULL;
    if(EIMSH_CIGS_VCIG) {
           Vocig_taper_angle = filterECIG_Angle(Vocig, eiproc, nx, nz, nxb, nzb, dx, dz, ox, oz, \
                                                ltheta0, dtheta, lv0, iproc, procGPU);
    }
    if(Hocig_taper_angle!=NULL) {
        memcpy(Hocig, Hocig_taper_angle, nx*nz*nlx*sizeof(float));
        free(Hocig_taper_angle);
    }
    if(Vocig_taper_angle!=NULL) {
        memcpy(Vocig, Vocig_taper_angle, nx*nz*nlv*sizeof(float));
        free(Vocig_taper_angle);
    }
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

    float *Hacig = CPU_zaloc1F(nx*nz*ntheta);
    GPU_WRAP_lx2theta(0, procGPU, 0, Hocig, Hacig, nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);


    // Transform to DOCIG
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    {printf("\n Transforming H&V CIGS to Docig ...... "); time_start = time(NULL);}
    float *Docig     = CPU_zaloc1F(nx*nz*nlambda*ndip);
    GPU_WRAP_VHocig2Docig(procGPU, Docig, Hocig, Vocig, nx, nz, dx, dz, ndip, ddip, dipMin, lx0, lv0, lambda0, iproc, EIMSH_CIGS_VCIG);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    {time_now = time(NULL); printf("done in %ld seconds \n", time_now-time_start);}


    float *Hocig_contracted      = CPU_zaloc1F(nx*nz*nlx);
    float *Hocig_contracted_ref  = CPU_zaloc1F(nx*nz*nlx);
    float *Vocig_contracted      = CPU_zaloc1F(nx*nz*nlv);
    float *Vocig_contracted_ref  = CPU_zaloc1F(nx*nz*nlv);
    float *Docig_contracted_ref  = CPU_zaloc1F(nx*nz*nlambda*ndip);
    float *Docig_contracted      = CPU_zaloc1F(nx*nz*nlambda*ndip);
    float *Docig_diff            = CPU_zaloc1F(nx*nz*nlambda*ndip);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Applying focusing operator to Docig ...... "); time_start = time(NULL);
    FocusingOperator_Docig(0, Docig_diff, Docig_contracted, Docig_contracted_ref, Docig, \
                           zMin, zMax, nx, nz, dx, dz, ox, oz, nxb, nzb, \
                           lambda0, ndip, dipMin, ddip, muFocus, epsFocus, lambdaNull, \
                           outputs, iproc, nproc);
    MPI_Barrier(MPI_COMM_WORLD);    time_now = time(NULL); if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

    float *Docig_contracted_phaseShift = CPU_zaloc1F(nx*nz*nlambda*ndip);
    float *Docig_contracted_original   = CPU_zaloc1F(nx*nz*nlambda*ndip);
    float *Docig_difference_phaseShift = CPU_zaloc1F(nx*nz*nlambda*ndip);
    float *Docig_difference_original   = CPU_zaloc1F(nx*nz*nlambda*ndip);

    if(computingGradient)  
    {
        get_HVocig_ExtExpRef_fromDocig(Docig_contracted, Docig_diff, Docig_contracted_phaseShift, Docig_contracted_original, \
                                       Docig_difference_phaseShift, Docig_difference_original, lambda0, ndip, dipMin, ddip, \
                                       nx, nz, nt, dx, dz, dt, ox, oz, ot, iproc, nproc, vp_bg);
        MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Transforming background and residual Docigs to HVocig ...... "); time_start = time(NULL);
        cg_VHocig2Docig(procGPU, Docig_contracted_original  , Hocig_background_original  , Vocig_background_original  , nx, nz, dx, dz, lx0, lv0, lambda0, ndip, ddip, dipMin, niterCG, iproc);
        cg_VHocig2Docig(procGPU, Docig_difference_original  , Hocig_residual_original    , Vocig_residual_original    , nx, nz, dx, dz, lx0, lv0, lambda0, ndip, ddip, dipMin, niterCG, iproc);
        cg_VHocig2Docig(procGPU, Docig_contracted_phaseShift, Hocig_background_phaseShift, Vocig_background_phaseShift, nx, nz, dx, dz, lx0, lv0, lambda0, ndip, ddip, dipMin, niterCG, iproc);
        cg_VHocig2Docig(procGPU, Docig_difference_phaseShift, Hocig_residual_phaseShift  , Vocig_residual_phaseShift  , nx, nz, dx, dz, lx0, lv0, lambda0, ndip, ddip, dipMin, niterCG, iproc);
        MPI_Barrier(MPI_COMM_WORLD);    time_now = time(NULL); if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);
    }
    if(focusingQC  ||  !computingGradient)
    {
        MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Transforming focused and reference Docigs to HVocig ...... "); time_start = time(NULL);
        cg_VHocig2Docig(procGPU, Docig_contracted    , Hocig_contracted    , Vocig_contracted    , nx, nz, dx, dz, lx0, lv0, lambda0, ndip, ddip, dipMin, niterCG, iproc);
        cg_VHocig2Docig(procGPU, Docig_contracted_ref, Hocig_contracted_ref, Vocig_contracted_ref, nx, nz, dx, dz, lx0, lv0, lambda0, ndip, ddip, dipMin, niterCG, iproc);
        MPI_Barrier(MPI_COMM_WORLD);    time_now = time(NULL); if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);
    }
    if(!computingGradient)
    {
        memcpy(Hocig_background_original, Hocig_contracted, nx*nz*nlx);
        memcpy(Vocig_background_original, Vocig_contracted, nx*nz*nlv);
        addArrays(Hocig_residual_original, 1.0, Hocig_contracted, -1.0, Hocig_contracted_ref, nx*nz*nlx);
        addArrays(Vocig_residual_original, 1.0, Vocig_contracted, -1.0, Vocig_contracted_ref, nx*nz*nlz);
    }
    MPI_Barrier(MPI_COMM_WORLD);    time_now = time(NULL); if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

    filterECIG(Hocig_background_original  , eiproc, eimage_msh, msh);
    filterECIG(Hocig_residual_original    , eiproc, eimage_msh, msh);
    if(computingGradient)  
    {
        filterECIG(Hocig_residual_phaseShift  , eiproc, eimage_msh, msh);
        filterECIG(Hocig_background_phaseShift, eiproc, eimage_msh, msh);
    }
    if(EIMSH_CIGS_VCIG)
    {
        eimage_msh->lx0 = lv0;
        filterECIG(Vocig_background_original  , eiproc, eimage_msh, msh);
        filterECIG(Vocig_residual_original    , eiproc, eimage_msh, msh);
        if(computingGradient)    
        {
            filterECIG(Vocig_background_phaseShift, eiproc, eimage_msh, msh);
            filterECIG(Vocig_residual_phaseShift  , eiproc, eimage_msh, msh);
        }
        eimage_msh->lx0 = lx0;
    }

    //  >>>>> DATA RESIDUALS <<<<<<
    lxNullTapering(Hocig_residual_original, Hocig_residual_original, nx, nz, dx, lx0, lambdaNull);
    lxNullTapering(Vocig_residual_original, Vocig_residual_original, nx, nz, dx, lv0, lambdaNull);
    double objFuncH = getArrayL2Norm(Hocig_residual_original, nx*nz*nlx);
    double objFuncV = getArrayL2Norm(Vocig_residual_original, nx*nz*nlv);
    float objFunc = objFuncH;
    if(EIMSH_CIGS_VCIG) objFunc += objFuncV;
    if(iproc==0)    printf("\n\n objFuncH=%g  objFuncV=%g  objFunc=%g \n\n", objFuncH, objFuncV, objFunc);


    float *Hacig_contracted     = CPU_zaloc1F(nx*nz*ntheta);
    float *Hacig_contracted_ref = CPU_zaloc1F(nx*nz*ntheta);
    if(outputs && focusingQC)
    {
        GPU_WRAP_lx2theta(0, procGPU, 0, Hocig_contracted_ref, Hacig_contracted_ref, nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);
        GPU_WRAP_lx2theta(0, procGPU, 0, Hocig_contracted    , Hacig_contracted    , nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);
    }

    // Saving GOCIG, H&V CIGs
    
    time_start = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);  if(iproc==0)  printf("\n Saving H&V&D OCIGs ...... ");
    if(outputs)
    {
        outputModel2d(image_FileName, image, nz, dz, 0.0, nx, dx, 0.0, nzb, nxb);
        outputSmart3d(Hacig_FileName, Hacig                       , ntheta , dtheta, -ltheta0*dtheta,  nz, dz, oz,  nx, dx, ox);
        outputSmart3d(Hocig_FileName, Hocig                       , nlx    ,     dx, -lx0*dx,  nz, dz, oz,  nx, dx, ox);
        outputSmart3d(Vocig_FileName, Vocig                       , nlv    ,     dz, -lv0*dz,  nz, dz, oz,  nx, dx, ox);
        outputSmart4d(Docig_FileName, Docig                       , nlambda,     dx, -lambda0*dx, nz, dz, oz, nx, dx, ox, ndip, ddip, dipMin);
        outputSmart3d(Hocig_diff_FileName, Hocig_residual_original, nlx    ,     dx, -lx0*dx,  nz, dz, oz, nx, dx, ox);
        // outputSmart3d(Hacig_diff_FileName, Hacig_diff, ntheta, dtheta, -ltheta0*dtheta,  nz, dz, oz, nx, dx, ox);
    }
    if(outputs && focusingQC)
    {
        outputSmart4d(Docig_contracted_FileName     , Docig_contracted    , nlambda, dx, -lambda0*dx, nz, dz, oz, nx, dx, ox, ndip, ddip, dipMin);
        outputSmart4d(Docig_contracted_ref_FileName , Docig_contracted_ref, nlambda, dx, -lambda0*dx, nz, dz, oz, nx, dx, ox, ndip, ddip, dipMin);
        outputSmart3d(Hocig_contracted_FileName     , Hocig_contracted    , nlx, dx, -lx0*dx,  nz, dz, oz, nx, dx, ox);
        outputSmart3d(Hocig_contracted_ref_FileName , Hocig_contracted_ref, nlx, dx, -lx0*dx,  nz, dz, oz, nx, dx, ox);
        outputSmart3d(Vocig_contracted_FileName     , Vocig_contracted    , nlv, dz, -lv0*dz,  nz, dz, oz, nx, dx, ox);
        outputSmart3d(Vocig_contracted_ref_FileName , Vocig_contracted_ref, nlv, dz, -lv0*dz,  nz, dz, oz, nx, dx, ox);
        outputSmart3d(Hacig_contracted_FileName     , Hacig_contracted    , ntheta, dtheta, -ltheta0*dtheta,  nz, dz, oz, nx, dx, ox);
        outputSmart3d(Hacig_contracted_ref_FileName , Hacig_contracted_ref, ntheta, dtheta, -ltheta0*dtheta,  nz, dz, oz, nx, dx, ox);
    }        
        time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);  if(iproc==0)  printf("done in %ld seconds \n", time_now-time_start);

    MPI_Barrier(MPI_COMM_WORLD);

    free(Hocig);
    free(Vocig);
    free(Hacig);
    free(Docig);
    free(Hocig_contracted);
    free(Hocig_contracted_ref);
    free(Vocig_contracted);
    free(Vocig_contracted_ref);
    free(Hacig_contracted);
    free(Hacig_contracted_ref);
    free(Docig_contracted_ref);
    free(Docig_contracted);
    free(Docig_diff);
    
    free(Docig_contracted_original);
    free(Docig_difference_original);
    free(Docig_contracted_phaseShift);
    free(Docig_difference_phaseShift);

    free(image);
    MPI_Barrier(MPI_COMM_WORLD);

    
return objFunc;
}


float computeResiduals_ORTH_PopDocig(Dataset *dataset_input, \
                                     float *Docig_background_phaseShift, float *Docig_background_original, \
                                     float *Docig_residual_phaseShift, float *Docig_residual_original, \
                                     MSH *msh, float *vp_bg, EIProc *eiproc, Eimage_MSH *eimage_msh, \
                                     float offMax, int acquisitionGeom, float perc, float zMin, float zMax, float maxDepth, \
                                     float muFocus, float epsFocus, float lambdaNull, int iter, int ism, int outputs, \
                                     int computingGradient, int focusingQC, int iproc, int nproc) {

    float         dip;
    long long int idip;
    long long int nxb         = msh->nxb;
    long long int nzb         = msh->nzb;
    long long int nx          = msh->nx;
    long long int nz          = msh->nz;
    long long int nt          = msh->nt;
    float dx                  = msh->dx;
    float dz                  = msh->dz;
    float dt                  = msh->dt;
    float ox                  = msh->ox;
    float oz                  = msh->oz;
    float ot                  = msh->ot;
    long long int jt          = msh->jt;
    long long int jt_EIC      = msh->jt_EIC;
    long long int jt_lt       = msh->jt_lt;
    long long int ltheta0     = eimage_msh->ltheta0;
    long long int ntheta      = eimage_msh->ntheta;
    float dtheta              = eimage_msh->dtheta;
    float otheta              = eimage_msh->otheta;
    long long int lx0         = eimage_msh->lx0;
    long long int lz0         = eimage_msh->lz0;
    long long int lv0         = eimage_msh->lv0;
    long long int lambda0     = eimage_msh->nlambda/2;
    long long int nlx         = eimage_msh->nlx;
    long long int nlz         = eimage_msh->nlz;
    long long int nlv         = eimage_msh->nlv;
    long long int nlambda     = eimage_msh->nlambda;
    float lambdaMax           = eimage_msh->lambdaMax;
    float dlambda             = eimage_msh->dlambda;
    float dipMin              = eimage_msh->dipMin;
    float dipMax              = eimage_msh->dipMax;
    float ddip                = eimage_msh->ddip;
    long long int ndip        = eimage_msh->ndip;

    int niterCG = 8;  // number of CG iterations to transform back from Docig to HVocig

    time_t time_start, time_now;

    float  *image = CPU_zaloc1F(nx*nz);

    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n\n Migrating data H&V OCIGS: nlx=%lld  nlz=%lld  nlv=%lld  jt_lt=%lld  ntheta=%lld  lv0=%lld  ndip=%lld ...... ", \
                                                            nlx, nlz, nlv, jt_lt, ntheta, lv0, ndip);
    time_start = time(NULL);
    float *Hocig = CPU_zaloc1F(nx*nz*nlx);
    float *Vocig = CPU_zaloc1F(nx*nz*nlv);
    migration_HOCIG_VOCIG(dataset_input, image, Hocig, Vocig, vp_bg, nx, nz, \
                          lx0, lv0, ltheta0, dx, dz, dtheta, ox, oz, \
                          otheta, jt, ism, 1, iproc, nproc);
    
    
    reciprocityAngle(Hocig, eiproc, acquisitionGeom, nx, nz, nxb, nzb, dx, dz, ox, oz, ltheta0, dtheta, lx0, iproc, procGPU);
    if(EIMSH_CIGS_VCIG) {
        reciprocityAngle(Vocig, eiproc, acquisitionGeom, nx, nz, nxb, nzb, dx, dz, ox, oz, ltheta0, dtheta, lv0, iproc, procGPU);
    }
    
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

    //  H&V OCIG filterings  
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n\n Filtering H&V CIGS ...... ");
    time_start = time(NULL);
    filterECIG(Hocig, eiproc, eimage_msh, msh);
    if(EIMSH_CIGS_VCIG) {
        eimage_msh->lx0 = lv0;
        filterECIG(Vocig, eiproc, eimage_msh, msh);
        eimage_msh->lx0 = lx0;
    }

    //  H&V ACIG filterings  
    float *Hocig_taper_angle = filterECIG_Angle(Hocig, eiproc, nx, nz, nxb, nzb, dx, dz, ox, oz, \
                                                ltheta0, dtheta, lx0, iproc, procGPU);
    float *Vocig_taper_angle = NULL;
    if(EIMSH_CIGS_VCIG) {
           Vocig_taper_angle = filterECIG_Angle(Vocig, eiproc, nx, nz, nxb, nzb, dx, dz, ox, oz, \
                                                ltheta0, dtheta, lv0, iproc, procGPU);
    }
    if(Hocig_taper_angle!=NULL) {
        memcpy(Hocig, Hocig_taper_angle, nx*nz*nlx*sizeof(float));
        free(Hocig_taper_angle);
    }
    if(Vocig_taper_angle!=NULL) {
        memcpy(Vocig, Vocig_taper_angle, nx*nz*nlv*sizeof(float));
        free(Vocig_taper_angle);
    }
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

    float *Hacig = CPU_zaloc1F(nx*nz*ntheta);
    GPU_WRAP_lx2theta(0, procGPU, 0, Hocig, Hacig, nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);


    // Transform to DOCIG
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    {printf("\n Transforming H&V CIGS to Docig ...... "); time_start = time(NULL);}
    float *Docig     = CPU_zaloc1F(nx*nz*nlambda*ndip);
    GPU_WRAP_VHocig2Docig(procGPU, Docig, Hocig, Vocig, nx, nz, dx, dz, ndip, ddip, dipMin, lx0, lv0, lambda0, iproc, EIMSH_CIGS_VCIG);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    {time_now = time(NULL); printf("done in %ld seconds \n", time_now-time_start);}


    // Apply focusing operator
    float *Docig_contracted_ref  = CPU_zaloc1F(nx*nz*nlambda*ndip);
    float *Docig_contracted      = CPU_zaloc1F(nx*nz*nlambda*ndip);
    float *Docig_diff            = CPU_zaloc1F(nx*nz*nlambda*ndip);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Applying focusing operator to Docig ...... "); time_start = time(NULL);
    FocusingOperator_Docig(0, Docig_diff, Docig_contracted, Docig_contracted_ref, Docig, \
                           zMin, zMax, nx, nz, dx, dz, ox, oz, nxb, nzb, \
                           lambda0, ndip, dipMin, ddip, muFocus, epsFocus, lambdaNull, \
                           outputs, iproc, nproc);
    MPI_Barrier(MPI_COMM_WORLD);    time_now = time(NULL); if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

    // Apply phase-shift
    if(computingGradient)  
    {
        MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Applying phase-shift to Docig ...... "); time_start = time(NULL);
        get_HVocig_ExtExpRef_fromDocig(Docig_contracted, Docig_diff, Docig_background_phaseShift, Docig_background_original, \
                                       Docig_residual_phaseShift, Docig_residual_original, lambda0, ndip, dipMin, ddip, \
                                       nx, nz, nt, dx, dz, dt, ox, oz, ot, iproc, nproc, vp_bg);
        MPI_Barrier(MPI_COMM_WORLD);    time_now = time(NULL); if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);
    }    

    //  Compute data residuals
    float objFunc = getArrayL2Norm(Docig_diff, nx*nz*nlambda*ndip);
    if(iproc==0)    printf("\n\n objFunc=%g \n\n", objFunc);

    MPI_Barrier(MPI_COMM_WORLD);

    // Outputs for QC
    float *Hocig_contracted     = CPU_zaloc1F(nx*nz*nlx);
    float *Hocig_contracted_ref = CPU_zaloc1F(nx*nz*nlx);
    float *Hocig_diff           = CPU_zaloc1F(nx*nz*nlx);
    float *Vocig_contracted     = CPU_zaloc1F(nx*nz*nlv);
    float *Vocig_contracted_ref = CPU_zaloc1F(nx*nz*nlv);
    float *Vocig_diff           = CPU_zaloc1F(nx*nz*nlv);
    float *Hacig_contracted     = CPU_zaloc1F(nx*nz*ntheta);
    float *Hacig_contracted_ref = CPU_zaloc1F(nx*nz*ntheta);
    float *Hacig_diff           = CPU_zaloc1F(nx*nz*ntheta);
    time_start = time(NULL);
    if(outputs && focusingQC)
    {
        // Computing gathers to output
        if(iproc==0)    printf("\n Transforming focused and reference Docigs to HVocig ...... "); time_start = time(NULL);
        cg_VHocig2Docig(procGPU, Docig_contracted    , Hocig_contracted    , Vocig_contracted    , nx, nz, dx, dz, lx0, lv0, lambda0, ndip, ddip, dipMin, niterCG, iproc);
        cg_VHocig2Docig(procGPU, Docig_contracted_ref, Hocig_contracted_ref, Vocig_contracted_ref, nx, nz, dx, dz, lx0, lv0, lambda0, ndip, ddip, dipMin, niterCG, iproc);
        GPU_WRAP_lx2theta(0, procGPU, 0, Hocig_contracted_ref, Hacig_contracted_ref, nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);
        GPU_WRAP_lx2theta(0, procGPU, 0, Hocig_contracted    , Hacig_contracted    , nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);
        addArrays(Hocig_diff, 1.0, Hocig_contracted, -1.0, Hocig_contracted_ref, nx*nz*nlx);
        addArrays(Vocig_diff, 1.0, Vocig_contracted, -1.0, Vocig_contracted_ref, nx*nz*nlz);
        addArrays(Hacig_diff, 1.0, Hacig_contracted, -1.0, Hacig_contracted_ref, nx*nz*ntheta);
        time_now = time(NULL); if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

        // Writing gathers to file
        if(iproc==0)  printf("\n Saving H&V&D OCIGs ...... ");
        outputModel2d(image_FileName                , image                  , nz     , dz    , 0.0            , nx, dx,0.0, nzb, nxb);
        outputSmart3d(Hacig_FileName                , Hacig                  , ntheta , dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
        outputSmart3d(Hocig_FileName                , Hocig                  , nlx    , dx    , -lx0*dx        , nz, dz, oz, nx, dx, ox);
        outputSmart3d(Vocig_FileName                , Vocig                  , nlv    , dz    , -lv0*dz        , nz, dz, oz, nx, dx, ox);
        outputSmart4d(Docig_FileName                , Docig                  , nlambda, dx    , -lambda0*dx    , nz, dz, oz, nx, dx, ox, ndip, ddip, dipMin);
        outputSmart3d(Hocig_diff_FileName           , Hocig_diff, nlx    , dx    , -lx0*dx        , nz, dz, oz, nx, dx, ox);
        outputSmart3d(Hacig_diff_FileName           , Hacig_diff             , ntheta , dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
        outputSmart4d(Docig_contracted_FileName     , Docig_contracted       , nlambda, dx    , -lambda0*dx    , nz, dz, oz, nx, dx, ox, ndip, ddip, dipMin);
        outputSmart4d(Docig_contracted_ref_FileName , Docig_contracted_ref   , nlambda, dx    , -lambda0*dx    , nz, dz, oz, nx, dx, ox, ndip, ddip, dipMin);
        outputSmart3d(Hocig_contracted_FileName     , Hocig_contracted       , nlx    , dx    , -lx0*dx        , nz, dz, oz, nx, dx, ox);
        outputSmart3d(Hocig_contracted_ref_FileName , Hocig_contracted_ref   , nlx    , dx    , -lx0*dx        , nz, dz, oz, nx, dx, ox);
        outputSmart3d(Vocig_contracted_FileName     , Vocig_contracted       , nlv    , dz    , -lv0*dz        , nz, dz, oz, nx, dx, ox);
        outputSmart3d(Vocig_contracted_ref_FileName , Vocig_contracted_ref   , nlv    , dz    , -lv0*dz        , nz, dz, oz, nx, dx, ox);
        outputSmart3d(Hacig_contracted_FileName     , Hacig_contracted       , ntheta , dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
        outputSmart3d(Hacig_contracted_ref_FileName , Hacig_contracted_ref   , ntheta , dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);

        outputSmart4d("Docig_background_original"  , Docig_background_original  , nlambda, dx    , -lambda0*dx    , nz, dz, oz, nx, dx, ox, ndip, ddip, dipMin);
        outputSmart4d("Docig_background_phaseShift", Docig_background_phaseShift, nlambda, dx    , -lambda0*dx    , nz, dz, oz, nx, dx, ox, ndip, ddip, dipMin);
        outputSmart4d("Docig_residual_original"    , Docig_residual_original    , nlambda, dx    , -lambda0*dx    , nz, dz, oz, nx, dx, ox, ndip, ddip, dipMin);
        outputSmart4d("Docig_residual_phaseShift"  , Docig_residual_phaseShift  , nlambda, dx    , -lambda0*dx    , nz, dz, oz, nx, dx, ox, ndip, ddip, dipMin);

        time_now = time(NULL);    if(iproc==0)  printf("done in %ld seconds \n", time_now-time_start);
    }        
    MPI_Barrier(MPI_COMM_WORLD);

    // Freeing memory
    free(Hocig);
    free(Vocig);
    free(Hacig);
    free(Docig);
    free(Hocig_contracted);
    free(Hocig_contracted_ref);
    free(Vocig_contracted);
    free(Vocig_contracted_ref);
    free(Hacig_contracted);
    free(Hacig_contracted_ref);
    free(Docig_contracted_ref);
    free(Docig_contracted);
    free(Docig_diff);
    free(Hocig_diff);
    free(Vocig_diff);
    free(Hacig_diff);

    free(image);
    MPI_Barrier(MPI_COMM_WORLD);

    
return objFunc;
}
///////////// Residuals with Popping CIGs - end /////////////////////////

float computeResiduals_ORTH(Dataset *dataset_input, Dataset *dataset_residuals, Dataset *dataset_reference, Dataset *dataset_misfit, \
                            float *Hocig_contracted, float *Hocig_contracted_ref, float *Hocig_normalizer, \
                            float *Vocig_contracted, float *Vocig_contracted_ref, float *Vocig_normalizer, \
                            float *Hocig_diff, float *Vocig_diff, \
                            MSH *msh, float *vp_bg, EIProc *eiproc, Eimage_MSH *eimage_msh, float offMax, int acquisitionGeom, float perc, \
                            float zMin, float zMax, float maxDepth, float muFocus, float epsFocus, float lambdaNull, int iter, int ism, int outputs, \
                            int imageDomain, int computingGradient, int iproc, int nproc) {

    float         dip;
    long long int idip;
    long long int nxb         = msh->nxb;
    long long int nzb         = msh->nzb;
    long long int nx          = msh->nx;
    long long int nz          = msh->nz;
    long long int nt          = msh->nt;
    float dx                  = msh->dx;
    float dz                  = msh->dz;
    float dt                  = msh->dt;
    float ox                  = msh->ox;
    float oz                  = msh->oz;
    float ot                  = msh->ot;
    long long int jt          = msh->jt;
    long long int jt_EIC      = msh->jt_EIC;
    long long int jt_lt       = msh->jt_lt;
    long long int ltheta0     = eimage_msh->ltheta0;
    long long int ntheta      = eimage_msh->ntheta;
    float dtheta              = eimage_msh->dtheta;
    float otheta              = eimage_msh->otheta;
    long long int lx0         = eimage_msh->lx0;
    long long int lz0         = eimage_msh->lz0;
    long long int lv0         = eimage_msh->lv0;
    long long int lt0         = eimage_msh->lt0;
    long long int lambda0     = eimage_msh->nlambda/2;
    long long int nlx         = eimage_msh->nlx;
    long long int nlz         = eimage_msh->nlz;
    long long int nlv         = eimage_msh->nlv;
    long long int nlt         = eimage_msh->nlt;
    long long int nlambda     = eimage_msh->nlambda;
    long long int nxLD        = eimage_msh->nx_LambdaDip;
    long long int nzLD        = eimage_msh->nz_LambdaDip;
    float oxLD                = eimage_msh->ox_LambdaDip;
    float ozLD                = eimage_msh->oz_LambdaDip;
    float lambdaMax           = eimage_msh->lambdaMax;
    float dlambda             = eimage_msh->dlambda;
    float dipMin              = eimage_msh->dipMin;
    float dipMax              = eimage_msh->dipMax;
    float ddip                = eimage_msh->ddip;
    long long int ndip        = eimage_msh->ndip;

    int niterCG = 8;  // number of CG iterations to transform back from Docig to HVocig
    

    char refDataFileName[1024];
    char misDataFileName[1024];
    char resDataFileName[1024];
    sprintf( refDataFileName , "seismogram_referenceData_Iter%d", iter);
    sprintf( misDataFileName , "seismogram_misfitData_Iter%d", iter);
    sprintf( resDataFileName , "seismogram_residuals_Iter%d", iter);

    time_t time_start, time_now;

    if(iproc==0)
        printf("\n\n  nlx=%lld  nlz=%lld  nlt=%lld  nlv=%lld  jt_lt=%lld  ntheta=%lld  \n\n", nlx, nlz, nlt, nlv, jt_lt, ntheta);

    float  *image = CPU_zaloc1F(nx*nz);

    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n\n Migrating data H&V OCIGS  lv0=%lld ...... ", lv0);
    time_start = time(NULL);
    float *Hocig = CPU_zaloc1F(nx*nz*nlx);
    float *Vocig = CPU_zaloc1F(nx*nz*nlv);
    migration_HOCIG_VOCIG(dataset_input, image, Hocig, Vocig, vp_bg, nx, nz, \
                          lx0, lv0, ltheta0, dx, dz, dtheta, ox, oz, \
                          otheta, jt, ism, 1, iproc, nproc);
    
    
    reciprocityAngle(Hocig, eiproc, acquisitionGeom, nx, nz, nxb, nzb, dx, dz, ox, oz, ltheta0, dtheta, lx0, iproc, procGPU);
    if(EIMSH_CIGS_VCIG) {
        reciprocityAngle(Vocig, eiproc, acquisitionGeom, nx, nz, nxb, nzb, dx, dz, ox, oz, ltheta0, dtheta, lv0, iproc, procGPU);
    }
    
    /*
    if(outputs)
    {
        outputModel2d("image_debug2", image, nz, dz, 0.0, nx, dx, 0.0, nzb, nxb);
        outputSmart3d("Hocig_debug2", Hocig, nlx, dx, -lx0*dx,  nz, dz, oz,  nx, dx, ox);        
    }
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n\n STOPPING PROGRAM ABRUPTLY - DEBUG \n\n");
    if(1)    exit(-1);
    */
    
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

    //  H&V CIG filterings  
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n\n Filtering H&V CIGS ...... ");
    time_start = time(NULL);
    filterECIG(Hocig, eiproc, eimage_msh, msh);
    if(EIMSH_CIGS_VCIG) {
        eimage_msh->lx0 = lv0;
        filterECIG(Vocig, eiproc, eimage_msh, msh);
        eimage_msh->lx0 = lx0;
    }

    float *Hocig_taper_angle = filterECIG_Angle(Hocig, eiproc, nx, nz, nxb, nzb, dx, dz, ox, oz, \
                                                ltheta0, dtheta, lx0, iproc, procGPU);
    float *Vocig_taper_angle = NULL;
    if(EIMSH_CIGS_VCIG) {
           Vocig_taper_angle = filterECIG_Angle(Vocig, eiproc, nx, nz, nxb, nzb, dx, dz, ox, oz, \
                                                ltheta0, dtheta, lv0, iproc, procGPU);
    }
    if(Hocig_taper_angle!=NULL) {
        memcpy(Hocig, Hocig_taper_angle, nx*nz*nlx*sizeof(float));
        free(Hocig_taper_angle);
    }
    if(Vocig_taper_angle!=NULL) {
        memcpy(Vocig, Vocig_taper_angle, nx*nz*nlv*sizeof(float));
        free(Vocig_taper_angle);
    }
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

    float *Hacig = CPU_zaloc1F(nx*nz*ntheta);
    GPU_WRAP_lx2theta(0, procGPU, 0, Hocig, Hacig, nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);

    ///////////// QC - BEGINs ///////////////
    /*
    if(outputs)
    {
        time_start = time(NULL);
        printf("\n Saving H&V&D OCIGs ...... ");
        // outputSmart2d(image_FileName, image, nz, dz, oz, nx, dx, ox);
        outputModel2d(image_FileName, image, nz, dz, 0.0, nx, dx, 0.0, nzb, nxb);
        outputSmart3d(Hacig_FileName,  Hacig        , ntheta, dtheta, -ltheta0*dtheta,  nz, dz, oz,  nx, dx, ox);
        outputSmart3d(Hocig_FileName,  Hocig        , nlx, dx, -lx0*dx,  nz, dz, oz,  nx, dx, ox);
        outputSmart3d(Vocig_FileName,  Vocig        , nlv, dz, -lv0*dz,  nz, dz, oz,  nx, dx, ox);        
        time_now = time(NULL);
        printf("done in %ld seconds \n", time_now-time_start);
    }
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Stopping program abruptly for QC!! ");
    if(1) exit(-1);
    //*/
    /////////////// QC - END //////////////////

    // Storing the ocig normalizer
    if(computingGradient) {
        memcpy(Hocig_normalizer, Hocig, nx*nz*nlx*sizeof(float));
        memcpy(Vocig_normalizer, Vocig, nx*nz*nlv*sizeof(float));
    }

    // normalizeHocigByEnvelope(Hocig, Hocig_normalizer, lx0, nx, nz, dz);
    // normalizeVocigByEnvelope(Vocig, Vocig_normalizer, lv0, nx, nz, dx);

    // Transform to DOCIG
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Transforming H&V CIGS to Docig ...... ");
    time_start = time(NULL);
    float *Docig     = CPU_zaloc1F(nx*nz*nlambda*ndip);
    // VHocig2Docig(Docig, Hocig, Vocig, nx, nz, ndip, dx, dz, ddip, dipMin, lx0, lv0, lambda0, iproc);
    GPU_WRAP_VHocig2Docig(procGPU, Docig, Hocig, Vocig, nx, nz, dx, dz, ndip, ddip, dipMin, lx0, lv0, lambda0, iproc, EIMSH_CIGS_VCIG);
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

    // Debug only
    /*
    float *Hocig_recover = CPU_zaloc1F(nx*nz*nlx);
    float *Vocig_recover = CPU_zaloc1F(nx*nz*nlv);
    // GPU_WRAP_Docig2VHocig(procGPU, Docig, Hocig_recover, Vocig_recover, nx, nz, dx, dz, ndip, ddip, dipMin, lx0, lv0, lambda0, iproc, EIMSH_CIGS_VCIG);
    // Docig2VHocig(Docig, Hocig_recover, Vocig_recover, nx, nz, ndip, dx, dz, ddip, dipMin, lx0, lv0, lambda0, iproc);
    cg_VHocig2Docig(procGPU, Docig, Hocig_recover, Vocig_recover, nx, nz, dx, dz, lx0, lv0, lambda0, ndip, ddip, dipMin, niterCG, iproc);
    float *Hacig_recover = CPU_zaloc1F(nx*nz*ntheta);
    GPU_WRAP_lx2theta(0, procGPU, 0, Hocig_recover, Hacig_recover, nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);
    //*/

    float *Docig_contracted_ref  = CPU_zaloc1F(nx*nz*nlambda*ndip);
    float *Docig_contracted      = CPU_zaloc1F(nx*nz*nlambda*ndip);
    float *Docig_diff            = CPU_zaloc1F(nx*nz*nlambda*ndip);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Applying focusing operator to Docig ...... "); time_start = time(NULL);
    FocusingOperator_Docig(0, Docig_diff, Docig_contracted, Docig_contracted_ref, Docig, \
                           zMin, zMax, nx, nz, dx, dz, ox, oz, nxb, nzb, \
                           lambda0, ndip, dipMin, ddip, muFocus, epsFocus, lambdaNull, \
                           outputs, iproc, nproc);

    
    MPI_Barrier(MPI_COMM_WORLD);    time_now = time(NULL); if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);
    
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Transforming focused and reference Docigs to HVocig ...... "); time_start = time(NULL);
    cg_VHocig2Docig(procGPU, Docig_contracted    , Hocig_contracted    , Vocig_contracted    , nx, nz, dx, dz, lx0, lv0, lambda0, ndip, ddip, dipMin, niterCG, iproc);
    cg_VHocig2Docig(procGPU, Docig_contracted_ref, Hocig_contracted_ref, Vocig_contracted_ref, nx, nz, dx, dz, lx0, lv0, lambda0, ndip, ddip, dipMin, niterCG, iproc);
    // Docig2VHocig(Docig_contracted    , Hocig_contracted    , Vocig_contracted    , nx, nz, ndip, dx, dz, ddip, dipMin, lx0, lv0, lambda0, iproc);
    // Docig2VHocig(Docig_contracted_ref, Hocig_contracted_ref, Vocig_contracted_ref, nx, nz, ndip, dx, dz, ddip, dipMin, lx0, lv0, lambda0, iproc);
    MPI_Barrier(MPI_COMM_WORLD);    time_now = time(NULL); if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n normalizing HVocig by envelope ...... "); time_start = time(NULL);
    // normalizeHocigByEnvelope(Hocig_contracted    , lx0, nx, nz, dz);
    // normalizeHocigByEnvelope(Hocig_contracted_ref, lx0, nx, nz, dz);
    // normalizeVocigByEnvelope(Vocig_contracted    , lv0, nx, nz, dx);
    // normalizeVocigByEnvelope(Vocig_contracted_ref, lv0, nx, nz, dx);
    // normalizeHocigByEnvelope(Hocig_contracted, Hocig_normalizer, lx0, nx, nz, dz);
    // normalizeHocigByEnvelope(Hocig_diff      , Hocig_normalizer, lx0, nx, nz, dz);
    // normalizeVocigByEnvelope(Vocig_contracted, Vocig_normalizer, lv0, nx, nz, dx);
    // normalizeVocigByEnvelope(Vocig_diff      , Vocig_normalizer, lv0, nx, nz, dx);
    MPI_Barrier(MPI_COMM_WORLD);    time_now = time(NULL); if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

    filterECIG(Hocig_contracted_ref, eiproc, eimage_msh, msh);
    filterECIG(Hocig_contracted    , eiproc, eimage_msh, msh);

    //  >>>>> DATA RESIDUALS <<<<<<  
    // multiplyArrayByScalar(0.0, Hocig_diff, nx*nz*nlx);
    // multiplyArrayByScalar(0.0, Vocig_diff, nx*nz*nlv);
    addArrays(Hocig_diff, 1.0, Hocig_contracted, -1.0, Hocig_contracted_ref, nx*nz*nlx);
    addArrays(Vocig_diff, 1.0, Vocig_contracted, -1.0, Vocig_contracted_ref, nx*nz*nlv);
    lxNullTapering(Hocig_diff, Hocig_diff, nx, nz, dx, lx0, lambdaNull);
    double objFuncH = getArrayL2Norm(Hocig_diff, nx*nz*nlx);
    double objFuncV = getArrayL2Norm(Vocig_diff, nx*nz*nlv);
    float objFunc = objFuncH + objFuncV;
    if(iproc==0)    printf("\n\n objFuncH=%g  objFuncV=%g  objFunc=%g \n\n", objFuncH, objFuncV, objFunc);
    // multiplyArrayByScalar(0.0, Vocig_diff, nx*nz*nlv);

    // filterECIG(Hocig_diff          , eiproc, eimage_msh, msh);
    
    if(EIMSH_CIGS_VCIG)
    {
        eimage_msh->lx0 = lv0;
        filterECIG(Vocig_contracted_ref, eiproc, eimage_msh, msh);
        filterECIG(Vocig_contracted    , eiproc, eimage_msh, msh);
        filterECIG(Vocig_diff          , eiproc, eimage_msh, msh);
        eimage_msh->lx0 = lx0;
    }

    float *Hacig_contracted     = CPU_zaloc1F(nx*nz*ntheta);
    float *Hacig_contracted_ref = CPU_zaloc1F(nx*nz*ntheta);
    GPU_WRAP_lx2theta(0, procGPU, 0, Hocig_contracted_ref, Hacig_contracted_ref, nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);
    GPU_WRAP_lx2theta(0, procGPU, 0, Hocig_contracted    , Hacig_contracted    , nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);

    // Saving GOCIG, H&V CIGs
    if(outputs)
    {
        time_start = time(NULL);
        printf("\n Saving H&V&D OCIGs ...... ");
        // outputSmart2d(image_FileName, image, nz, dz, oz, nx, dx, ox);
        outputModel2d(image_FileName, image, nz, dz, 0.0, nx, dx, 0.0, nzb, nxb);
        outputSmart3d(Hacig_FileName,  Hacig        , ntheta, dtheta, -ltheta0*dtheta,  nz, dz, oz,  nx, dx, ox);
        outputSmart3d(Hocig_FileName,  Hocig        , nlx, dx, -lx0*dx,  nz, dz, oz,  nx, dx, ox);
        outputSmart3d(Vocig_FileName,  Vocig        , nlv, dz, -lv0*dz,  nz, dz, oz,  nx, dx, ox);
        // outputSmart3d("Hocig_recover", Hocig_recover, nlx, dx, -lx0*dx,  nz, dz, oz,  nx, dx, ox);
        // outputSmart3d("Hacig_recover", Hacig_recover, ntheta, dtheta, -ltheta0*dtheta,  nz, dz, oz,  nx, dx, ox);
        // outputSmart3d("Vocig_recover", Vocig_recover, nlv, dz, -lv0*dz,  nz, dz, oz,  nx, dx, ox);
        outputSmart4d(Docig_FileName, Docig, nlambda, dx, -lambda0*dx, nz, dz, oz, nx, dx, ox, ndip, ddip, dipMin);


        outputSmart4d(Docig_contracted_FileName     , Docig_contracted    , nlambda, dx, -lambda0*dx, nz, dz, oz, nx, dx, ox, ndip, ddip, dipMin);
        outputSmart4d(Docig_contracted_ref_FileName , Docig_contracted_ref, nlambda, dx, -lambda0*dx, nz, dz, oz, nx, dx, ox, ndip, ddip, dipMin);
        outputSmart3d(Hocig_contracted_FileName     , Hocig_contracted    , nlx, dx, -lx0*dx,  nz, dz, oz, nx, dx, ox);
        outputSmart3d(Hocig_contracted_ref_FileName , Hocig_contracted_ref, nlx, dx, -lx0*dx,  nz, dz, oz, nx, dx, ox);
        outputSmart3d(Vocig_contracted_FileName     , Vocig_contracted    , nlv, dz, -lv0*dz,  nz, dz, oz, nx, dx, ox);
        outputSmart3d(Vocig_contracted_ref_FileName , Vocig_contracted_ref, nlv, dz, -lv0*dz,  nz, dz, oz, nx, dx, ox);
        outputSmart3d(Hacig_contracted_FileName     , Hacig_contracted    , ntheta, dtheta, -ltheta0*dtheta,  nz, dz, oz, nx, dx, ox);
        outputSmart3d(Hacig_contracted_ref_FileName , Hacig_contracted_ref, ntheta, dtheta, -ltheta0*dtheta,  nz, dz, oz, nx, dx, ox);

        outputSmart3d(Hocig_diff_FileName, Hocig_diff, nlx, dx, -lx0*dx,  nz, dz, oz, nx, dx, ox);
        // outputSmart3d(Hacig_diff_FileName, Hacig_diff, ntheta, dtheta, -ltheta0*dtheta,  nz, dz, oz, nx, dx, ox);
        
        time_now = time(NULL);
        printf("done in %ld seconds \n", time_now-time_start);
    }

    // NEW CODE SEGMENT ///////////// !!!!!!!!!!!!!!!!!!!!!
    // Saving Misfit and Reference datasets
    if(dataset_reference != NULL)
    {
        long long int lv0_tmp;
        if(EIMSH_CIGS_VCIG)  lv0_tmp = lv0;
        else                 lv0_tmp = 0;  
        if(iproc==0)    printf("\n Modeling reference data \n");
        extendedBornModeling(dataset_reference, NULL, Hocig_contracted_ref, Vocig_contracted_ref, vp_bg, nx, nz, lx0, lv0_tmp, lt0, \
                             dx, dz, ox, oz, jt, perc, offMax, acquisitionGeom, NULL, iproc, nproc);

        int jshot = 10;
        writeSeisData2File(dataset_reference, jshot, iproc, refDataFileName, NULL);
        if(iproc==0)    printf("\n Modeling reference data  ----->>>>> concluded\n");
    }
    if(dataset_misfit    != NULL)
    {
        long long int lv0_tmp;
        if(EIMSH_CIGS_VCIG)  lv0_tmp = lv0;
        else                 lv0_tmp = 0;  
        if(iproc==0)    printf("\n Modeling misfit data \n");
        extendedBornModeling(dataset_misfit, NULL, Hocig_contracted, Vocig_contracted, vp_bg, nx, nz, lx0, lv0_tmp, lt0, \
                             dx, dz, ox, oz, jt, perc, offMax, acquisitionGeom, NULL, iproc, nproc);
        int jshot = 10;
        writeSeisData2File(dataset_misfit   , jshot, iproc, misDataFileName, NULL);
        if(iproc==0)    printf("\n Modeling misfit data  ----->>>>> concluded\n");
    }
    if(computingGradient)
    {
        long long int lv0_tmp;
        if(EIMSH_CIGS_VCIG)  lv0_tmp = lv0;
        else                 lv0_tmp = 0;  
        if(iproc==0)    printf("\n\n Computing data residuals. \n (Shot,time): \n");
        extendedBornModeling(dataset_residuals, NULL, Hocig_diff, Vocig_diff, vp_bg, nx, nz, lx0, lv0_tmp, lt0, \
                             dx, dz, ox, oz, jt, perc, offMax, acquisitionGeom, NULL, iproc, nproc);
        int jshot = 10;
        writeSeisData2File(dataset_residuals  , jshot, iproc, resDataFileName, NULL);
    }
    // if(iproc==0)    printf("\n\n Computing data residuals. \n (Shot,time): \n");
    // objFunc = extendedBornModeling(dataset_residuals, NULL, Hocig_diff, Vocig_diff, vp_bg, nx, nz, lx0, lv0, lt0, \
    //                                dx, dz, ox, oz, jt, perc, offMax, acquisitionGeom, NULL, iproc, nproc);
    MPI_Barrier(MPI_COMM_WORLD);

    free(Hocig);
    free(Vocig);
    // free(Hocig_recover);    // for debug
    // free(Vocig_recover);    // for debug
    // free(Hacig_recover);    // for debug
    MPI_Barrier(MPI_COMM_WORLD);
    free(Hacig_contracted);
    free(Hacig_contracted_ref);
    free(Hacig);
    free(Docig);
    free(Docig_contracted_ref);
    free(Docig_contracted);
    free(Docig_diff);
    free(image);
    
return objFunc;
}

float computeResiduals_ORTH_Gocig(Dataset *dataset_input, Dataset *dataset_residuals, Dataset *dataset_reference, Dataset *dataset_misfit, \
                            float *Hocig_contracted, float *Hocig_contracted_ref, float *Vocig_contracted, float *Vocig_contracted_ref, \
                            MSH *msh, float *vp_bg, EIProc *eiproc, Eimage_MSH *eimage_msh, float offMax, int acquisitionGeom, float perc, \
                            float zMin, float zMax, float maxDepth, float muFocus, float epsFocus, float lambdaNull, int iter, int ism, int outputs, \
                            int imageDomain, int computingGradient, int iproc, int nproc) {

    long long int nxb         = msh->nxb;
    long long int nzb         = msh->nzb;
    long long int nx          = msh->nx;
    long long int nz          = msh->nz;
    long long int nt          = msh->nt;
    float dx                  = msh->dx;
    float dz                  = msh->dz;
    float dt                  = msh->dt;
    float ox                  = msh->ox;
    float oz                  = msh->oz;
    float ot                  = msh->ot;
    long long int jt          = msh->jt;
    long long int jt_EIC      = msh->jt_EIC;
    long long int jt_lt       = msh->jt_lt;
    long long int ltheta0     = eimage_msh->ltheta0;
    long long int ntheta      = eimage_msh->ntheta;
    float dtheta              = eimage_msh->dtheta;
    float otheta              = eimage_msh->otheta;
    long long int lx0         = eimage_msh->lx0;
    long long int lz0         = eimage_msh->lz0;
    long long int lv0         = eimage_msh->lv0;
    long long int lt0         = eimage_msh->lt0;
    long long int nlx         = eimage_msh->nlx;
    long long int nlz         = eimage_msh->nlz;
    long long int nlv         = eimage_msh->nlv;
    long long int nlt         = eimage_msh->nlt;
    long long int lambda0     = eimage_msh->nlambda/2;
    long long int   nxLD      = eimage_msh->nx_LambdaDip;
    long long int   nzLD      = eimage_msh->nz_LambdaDip;
    float oxLD                = eimage_msh->ox_LambdaDip;
    float ozLD                = eimage_msh->oz_LambdaDip;
    long long int   nlambda   = eimage_msh->nlambda;
    float lambdaMax           = eimage_msh->lambdaMax;
    float dlambda             = eimage_msh->dlambda;
    float dipMin              = eimage_msh->dipMin;
    float dipMax              = eimage_msh->dipMax;
    float ddip                = eimage_msh->ddip;
    long long int   ndip      = eimage_msh->ndip;
    

    char refDataFileName[1024];
    char misDataFileName[1024];
    char resDataFileName[1024];
    sprintf( refDataFileName , "seismogram_referenceData_Iter%d", iter);
    sprintf( misDataFileName , "seismogram_misfitData_Iter%d", iter);
    sprintf( resDataFileName , "seismogram_residuals_Iter%d", iter);

    time_t time_start, time_now;

    if(iproc==0)
        printf("\n\n  nlx=%lld  nlz=%lld  nlt=%lld  nlv=%lld  jt_lt=%lld  ntheta=%lld  \n\n", nlx, nlz, nlt, nlv, jt_lt, ntheta);

    float  *image = CPU_zaloc1F(nx*nz);

    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n\n Migrating data H&V OCIGS  lv0=%lld ...... ", lv0);
    time_start = time(NULL);
    float *Hocig = CPU_zaloc1F(nx*nz*nlx);
    float *Vocig = CPU_zaloc1F(nx*nz*nlv);
    migration_HOCIG_VOCIG(dataset_input, image, Hocig, Vocig, vp_bg, nx, nz, \
                          lx0, lv0, ltheta0, dx, dz, dtheta, ox, oz, \
                          otheta, jt, ism, 1, iproc, nproc);
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

    //  H&V CIG filterings  
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n\n Filtering H&V CIGS ...... ");
    time_start = time(NULL);
    filterECIG(Hocig, eiproc, eimage_msh, msh);
    eimage_msh->lx0 = lv0;
    filterECIG(Vocig, eiproc, eimage_msh, msh);
    eimage_msh->lx0 = lx0;
    float *Hocig_taper_angle = filterECIG_Angle(Hocig, eiproc, nx, nz, nxb, nzb, dx, dz, ox, oz, \
                                                ltheta0, dtheta, lx0, iproc, procGPU);
    float *Vocig_taper_angle = filterECIG_Angle(Vocig, eiproc, nx, nz, nxb, nzb, dx, dz, ox, oz, \
                                                ltheta0, dtheta, lv0, iproc, procGPU);
    if(Hocig_taper_angle!=NULL)
    {
        memcpy(Hocig, Hocig_taper_angle, nx*nz*nlx*sizeof(float));
        free(Hocig_taper_angle);
        memcpy(Vocig, Vocig_taper_angle, nx*nz*nlv*sizeof(float));
        free(Vocig_taper_angle);
    }
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

    // Transform to Gocig
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Transforming H&V CIGS to Gocig ...... ");
    time_start = time(NULL);
    float *Gocig;
    Gocig = VHocig2Gocig(Hocig, Vocig, nx, nz, dx, dz, lx0, lv0, lv0, iproc);
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);
    
    // GOCIG to H&V CIG
    float *Hocig_debug = CPU_zaloc1F(nx*nz*nlx);
    float *Vocig_debug = CPU_zaloc1F(nx*nz*nlv);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n \n Recovering DEBUG H&V OCIGs ...... ");
    time_start = time(NULL);
    Gocig2VHocig(Gocig, Hocig_debug, Vocig_debug, nx, nz, dx, dz, lx0, lv0, lv0, iproc);
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n \n", time_now-time_start);
    
    // GOCIG focusing
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n \n Applying DEBUG Focusing Operator in GOCIGs ...... ");
    time_start = time(NULL);
    float *Gocig_contracted_ref_dbg = CPU_zaloc1F(nx*nz*nlx*nlv);
    float *Gocig_contracted_dbg     = CPU_zaloc1F(nx*nz*nlx*nlv);
    float *Gocig_diff_dbg           = CPU_zaloc1F(nx*nz*nlx*nlv);
    FocusingOperator_Gocig(0, Gocig_diff_dbg, Gocig_contracted_dbg, Gocig_contracted_ref_dbg, Gocig, zMin, zMax, maxDepth, \
                           nx, nz, dx, dz, ox, oz, nxb, nzb, lx0, lv0, ltheta0, dtheta, otheta, lambdaMax, dipMin, dipMax, ddip, \
                           iter, muFocus, epsFocus, lambdaNull, 0, iproc, nproc);
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n \n", time_now-time_start);

    // Saving Focused Gocig
    if(outputs)
    {
        time_start = time(NULL);
        printf("\n Saving DEBUG H&V OCIGs ...... ");
        if(iter==0)
        {
            outputSmart3d("Hocig_debug", Hocig_debug, nlx, dx, -lx0*dx, nz, dz, oz, nx, dx, ox);
            outputSmart3d("Vocig_debug", Vocig_debug, nlv, dz, -lv0*dz, nz, dz, oz, nx, dx, ox);
            outputSmart4d("Gocig_contracted_debug", Gocig_contracted_dbg, nlx, dx, -lx0*dx, nlv, dz, -lv0*dz, nz, dz, oz, nx, dx, ox);
            outputSmart4d("Gocig_contracted_ref_debug", Gocig_contracted_ref_dbg, nlx, dx, -lx0*dx, nlv, dz, -lv0*dz, nz, dz, oz, nx, dx, ox);
        }
        time_now = time(NULL);
        printf("done in %ld seconds \n", time_now-time_start);
    }
    free(Hocig_debug);
    free(Vocig_debug);

    // Saving Gocig, H&V CIGs
    if(outputs)
    {
        time_start = time(NULL);
        printf("\n Saving H&V OCIGs ...... ");
        // outputSmart2d(image_FileName, image, nz, dz, oz, nx, dx, ox);
        outputModel2d(image_FileName, image, nz, dz, 0.0, nx, dx, 0.0, nzb, nxb);
        outputSmart3d(Hocig_FileName, Hocig, nlx, dx, -lx0*dx,  nz, dz, oz,  nx, dx, ox);
        outputSmart3d(Vocig_FileName, Vocig, nlv, dz, -lv0*dz,  nz, dz, oz,  nx, dx, ox);
        outputSmart4d(Gocig_FileName, Gocig, nlx, dx, -lx0*dx, nlv, dz, -lv0*dz, nz, dz, oz, nx, dx, ox);
        time_now = time(NULL);
        printf("done in %ld seconds \n", time_now-time_start);
    }
    MPI_Barrier(MPI_COMM_WORLD);


    // Focused GOCIG to H&V CIGs
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Transforming Gocig (focused and reference) back to H&V CIGs ...... ");
    time_start = time(NULL);
    Gocig2VHocig(Gocig_contracted_dbg    , Hocig_contracted    , Vocig_contracted    , nx, nz, dx, dz, lx0, lv0, lv0, iproc);
    Gocig2VHocig(Gocig_contracted_ref_dbg, Hocig_contracted_ref, Vocig_contracted_ref, nx, nz, dx, dz, lx0, lv0, lv0, iproc);
    freeMem(Gocig_contracted_ref_dbg, nx*nz*nlx*nlv*sizeof(float));
    freeMem(Gocig_contracted_dbg    , nx*nz*nlx*nlv*sizeof(float));
    freeMem(Gocig_diff_dbg          , nx*nz*nlx*nlv*sizeof(float));
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n done %ld seconds \n", time_now-time_start);

    filterECIG(Hocig_contracted    , eiproc, eimage_msh, msh);
    filterECIG(Hocig_contracted_ref, eiproc, eimage_msh, msh);


    // Saving Focused H&V CIGs
    MPI_Barrier(MPI_COMM_WORLD);  if(iproc==0)    printf("\n\n Saving contracted H&V OCIGs ...... ");
    time_start = time(NULL);
    if(outputs)
    {
        outputSmart3d(Hocig_contracted_FileName    , Hocig_contracted    , nlx, dx, -lx0*dx,  nz, dz,      oz, nx, dx, ox);
        outputSmart3d(Hocig_contracted_ref_FileName, Hocig_contracted_ref, nlx, dx, -lx0*dx,  nz, dz,      oz, nx, dx, ox);
        outputSmart3d(Vocig_contracted_FileName    , Vocig_contracted    , nlv, dz, -lv0*dz,  nz, dz,      oz, nx, dx, ox);
        outputSmart3d(Vocig_contracted_ref_FileName, Vocig_contracted_ref, nlv, dz, -lv0*dz,  nz, dz,      oz, nx, dx, ox);
    }
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);
    free(Hocig);
    free(Vocig);

    MPI_Barrier(MPI_COMM_WORLD);

    // Saving Misfit and Reference datasets
    if(dataset_reference != NULL)
    {
        if(iproc==0)    printf("\n Modeling reference data \n");
        extendedBornModeling(dataset_reference, NULL, Hocig_contracted_ref, Vocig_contracted_ref, vp_bg, nx, nz, lx0, lv0, lt0, \
                             dx, dz, ox, oz, jt, perc, offMax, acquisitionGeom, NULL, iproc, nproc);

        int jshot = 10;
        writeSeisData2File(dataset_reference, jshot, iproc, refDataFileName, NULL);
        if(iproc==0)    printf("\n Modeling reference data  ----->>>>> concluded\n");
    }
    if(dataset_misfit    != NULL)
    {
        if(iproc==0)    printf("\n Modeling misfit data \n");
        extendedBornModeling(dataset_misfit, NULL, Hocig_contracted, Vocig_contracted, vp_bg, nx, nz, lx0, lv0, lt0, \
                             dx, dz, ox, oz, jt, perc, offMax, acquisitionGeom, NULL, iproc, nproc);
        int jshot = 10;
        writeSeisData2File(dataset_misfit   , jshot, iproc, misDataFileName, NULL);
        if(iproc==0)    printf("\n Modeling misfit data  ----->>>>> concluded\n");
    }

    //  >>>>> DATA RESIDUALS <<<<<<  
    float *Hocig_diff = CPU_zaloc1F(nx*nz*nlx);
    float *Vocig_diff = CPU_zaloc1F(nx*nz*nlv);
    addArrays(Hocig_diff, 1.0, Hocig_contracted, -1.0, Hocig_contracted_ref, nx*nz*nlx);
    addArrays(Vocig_diff, 1.0, Vocig_contracted, -1.0, Vocig_contracted_ref, nx*nz*nlv);
    lxNullTapering(Hocig_diff, Hocig_diff, nx, nz, dx, lx0, lambdaNull);
    double objFuncH = getArrayL2Norm(Hocig_diff, nx*nz*nlx);
    double objFuncV = getArrayL2Norm(Vocig_diff, nx*nz*nlv);
    float objFunc = objFuncH + 0*objFuncV;
    if(iproc==0)    printf("\n\n objFuncH=%g  objFuncV=%g  objFunc=%g \n\n", objFuncH, objFuncV, objFunc);
    multiplyArrayByScalar(0.0, Vocig_diff, nx*nz*nlv);
    if(computingGradient)
    {
        if(iproc==0)    printf("\n\n Computing data residuals. \n (Shot,time): \n");
        extendedBornModeling(dataset_residuals, NULL, Hocig_diff, Vocig_diff, vp_bg, nx, nz, lx0, lv0, lt0, \
                             dx, dz, ox, oz, jt, perc, offMax, acquisitionGeom, NULL, iproc, nproc);
    }
    // if(iproc==0)    printf("\n\n Computing data residuals. \n (Shot,time): \n");
    // objFunc = extendedBornModeling(dataset_residuals, NULL, Hocig_diff, Vocig_diff, vp_bg, nx, nz, lx0, lv0, lt0, \
    //                                dx, dz, ox, oz, jt, perc, offMax, acquisitionGeom, NULL, iproc, nproc);
    MPI_Barrier(MPI_COMM_WORLD);

    free(image);
    free(Hocig_diff);
    free(Vocig_diff);
    

return objFunc;
}

float computeResiduals_ORTH_old(Dataset *dataset_input, Dataset *dataset_residuals, Dataset *dataset_reference, Dataset *dataset_misfit, \
                                float *Hocig_contracted, float *Hocig_contracted_ref, float *Vocig_contracted, float *Vocig_contracted_ref, \
                                MSH *msh, float *vp_bg, EIProc *eiproc, Eimage_MSH *eimage_msh, float offMax, int acquisitionGeom, float perc, \
                                float zMin, float zMax, float maxDepth, float muFocus, float epsFocus, float lambdaNull, int iter, int ism, int outputs, \
                                int imageDomain, int computingGradient, int iproc, int nproc) {

    long long int nxb         = msh->nxb;
    long long int nzb         = msh->nzb;
    long long int nx          = msh->nx;
    long long int nz          = msh->nz;
    long long int nt          = msh->nt;
    float dx                  = msh->dx;
    float dz                  = msh->dz;
    float dt                  = msh->dt;
    float ox                  = msh->ox;
    float oz                  = msh->oz;
    float ot                  = msh->ot;
    long long int jt          = msh->jt;
    long long int jt_EIC      = msh->jt_EIC;
    long long int jt_lt       = msh->jt_lt;
    long long int ltheta0     = eimage_msh->ltheta0;
    long long int ntheta      = eimage_msh->ntheta;
    float dtheta              = eimage_msh->dtheta;
    float otheta              = eimage_msh->otheta;
    long long int lx0         = eimage_msh->lx0;
    long long int lz0         = eimage_msh->lz0;
    long long int lv0         = eimage_msh->lv0;
    long long int lt0         = eimage_msh->lt0;
    long long int nlx         = eimage_msh->nlx;
    long long int nlz         = eimage_msh->nlz;
    long long int nlv         = eimage_msh->nlv;
    long long int nlt         = eimage_msh->nlt;
    long long int lambda0     = eimage_msh->nlambda/2;
    long long int   nxLD      = eimage_msh->nx_LambdaDip;
    long long int   nzLD      = eimage_msh->nz_LambdaDip;
    float oxLD                = eimage_msh->ox_LambdaDip;
    float ozLD                = eimage_msh->oz_LambdaDip;
    long long int   nlambda   = eimage_msh->nlambda;
    float lambdaMax           = eimage_msh->lambdaMax;
    float dlambda             = eimage_msh->dlambda;
    float dipMin              = eimage_msh->dipMin;
    float dipMax              = eimage_msh->dipMax;
    float ddip                = eimage_msh->ddip;
    long long int   ndip      = eimage_msh->ndip;

    // printf("\n\n\n\n ================================== \n\n ORTH CIGs PARAMETER CHECK: \n");
    // printf("nx=%lld nz=%lld nt=%lld dx=%f dz=%f dt=%f ox=%f oz=%f ot=%f lx0=%lld lz0=%lld lv0=%lld nlx=%lld nlz=%lld nlv=%lld", \
    //         nx, nz, nt, dx, dz, dt, ox, oz, ot, lx0, lz0, lv0, nlx, nlz, nlv);
    // printf("\n\n LambdaDip Parms: nxLD=%lld  nzLD=%lld  oxLD=%f  ozLD=%f \n\n", nxLD, nzLD, oxLD, ozLD);
    

    char refDataFileName[1024];
    char misDataFileName[1024];
    char resDataFileName[1024];
    sprintf( refDataFileName , "seismogram_referenceData_Iter%d", iter);
    sprintf( misDataFileName , "seismogram_misfitData_Iter%d", iter);
    sprintf( resDataFileName , "seismogram_residuals_Iter%d", iter);

    time_t time_start, time_now;

    if(iproc==0)
        printf("\n\n  nlx=%lld  nlz=%lld  nlt=%lld  nlv=%lld  jt_lt=%lld  ntheta=%lld  \n\n", nlx, nlz, nlt, nlv, jt_lt, ntheta);

    float  *image = CPU_zaloc1F(nx*nz);

    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n\n Migrating data H&V OCIGS  lv0=%lld ...... ", lv0);
    time_start = time(NULL);
    float *Hocig = CPU_zaloc1F(nx*nz*nlx);
    float *Vocig = CPU_zaloc1F(nx*nz*nlv);
    migration_HOCIG_VOCIG(dataset_input, image, Hocig, Vocig, vp_bg, nx, nz, \
                          lx0, lv0, ltheta0, dx, dz, dtheta, ox, oz, \
                          otheta, jt, ism, 1, iproc, nproc);
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

    //  H&V CIG filterings  
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n\n Filtering H&V CIGS ...... ");
    time_start = time(NULL);
    filterECIG(Hocig, eiproc, eimage_msh, msh);
    eimage_msh->lx0 = lv0;
    filterECIG(Vocig, eiproc, eimage_msh, msh);
    eimage_msh->lx0 = lx0;
    float *Hocig_taper_angle = filterECIG_Angle(Hocig, eiproc, nx, nz, nxb, nzb, dx, dz, ox, oz, \
                                                ltheta0, dtheta, lx0, iproc, procGPU);
    float *Vocig_taper_angle = filterECIG_Angle(Vocig, eiproc, nx, nz, nxb, nzb, dx, dz, ox, oz, \
                                                ltheta0, dtheta, lv0, iproc, procGPU);
    if(Hocig_taper_angle!=NULL)
    {
        memcpy(Hocig, Hocig_taper_angle, nx*nz*nlx*sizeof(float));
        free(Hocig_taper_angle);
        // freeMem(Hocig_taper_angle, nx*nz*nlx*sizeof(float));
        memcpy(Vocig, Vocig_taper_angle, nx*nz*nlv*sizeof(float));
        free(Vocig_taper_angle);
        // freeMem(Vocig_taper_angle, nx*nz*nlx*sizeof(float));
    }
    // multiplyArrayByScalar(0.0, Vocig, nx*nz*nlv);
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

    // Transform to Gocig
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Transforming H&V CIGS to Gocig ...... ");
    time_start = time(NULL);
    float *Gocig;
    Gocig = VHocig2Gocig(Hocig, Vocig, nx, nz, dx, dz, lx0, lv0, lv0, iproc);
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);
    
    // ========= DEBUG - begin ===========
    // ========= DEBUG - begin ===========
    float *Hocig_debug = CPU_zaloc1F(nx*nz*nlx);
    float *Vocig_debug = CPU_zaloc1F(nx*nz*nlv);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n \n Recovering DEBUG H&V OCIGs ...... ");
    time_start = time(NULL);
    Gocig2VHocig(Gocig, Hocig_debug, Vocig_debug, nx, nz, dx, dz, lx0, lv0, lv0, iproc);
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n \n", time_now-time_start);
    
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n \n Applying DEBUG Focusing Operator in GOCIGs ...... ");
    time_start = time(NULL);
    float *Gocig_contracted_ref_dbg = CPU_zaloc1F(nx*nz*nlx*nlv);
    float *Gocig_contracted_dbg     = CPU_zaloc1F(nx*nz*nlx*nlv);
    float *Gocig_diff_dbg           = CPU_zaloc1F(nx*nz*nlx*nlv);
    FocusingOperator_Gocig(0, Gocig_diff_dbg, Gocig_contracted_dbg, Gocig_contracted_ref_dbg, Gocig, zMin, zMax, maxDepth, \
                           nx, nz, dx, dz, ox, oz, nxb, nzb, lx0, lv0, ltheta0, dtheta, otheta, lambdaMax, dipMin, dipMax, ddip, \
                           iter, muFocus, epsFocus, lambdaNull, 0, iproc, nproc);
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n \n", time_now-time_start);

    if(outputs)
    {
        time_start = time(NULL);
        printf("\n Saving DEBUG H&V OCIGs ...... ");
        if(iter==0)
        {
            outputSmart3d("Hocig_debug", Hocig_debug, nlx, dx, -lx0*dx, nz, dz, oz, nx, dx, ox);
            outputSmart3d("Vocig_debug", Vocig_debug, nlv, dz, -lv0*dz, nz, dz, oz, nx, dx, ox);
            outputSmart4d("Gocig_contracted_debug", Gocig_contracted_dbg, nlx, dx, -lx0*dx, nlv, dz, -lv0*dz, nz, dz, oz, nx, dx, ox);
            outputSmart4d("Gocig_contracted_ref_debug", Gocig_contracted_ref_dbg, nlx, dx, -lx0*dx, nlv, dz, -lv0*dz, nz, dz, oz, nx, dx, ox);
        }
        time_now = time(NULL);
        printf("done in %ld seconds \n", time_now-time_start);
    }
    free(Hocig_debug);
    free(Vocig_debug);
    // free(Gocig_contracted_dbg);
    // free(Gocig_contracted_ref_dbg);
    // ========= DEBUG - end ===========
    // ========= DEBUG - end ===========

    if(outputs)
    {
        time_start = time(NULL);
        printf("\n Saving H&V OCIGs ...... ");
        // outputSmart2d(image_FileName, image, nz, dz, oz, nx, dx, ox);
        outputModel2d(image_FileName, image, nz, dz, 0.0, nx, dx, 0.0, nzb, nxb);
        outputSmart3d(Hocig_FileName, Hocig, nlx, dx, -lx0*dx,  nz, dz, oz,  nx, dx, ox);
        outputSmart3d(Vocig_FileName, Vocig, nlv, dz, -lv0*dz,  nz, dz, oz,  nx, dx, ox);
        outputSmart4d(Gocig_FileName, Gocig, nlx, dx, -lx0*dx, nlv, dz, -lv0*dz, nz, dz, oz, nx, dx, ox);
        // outputSmart4d_J(Gocig_FileName, Gocig, nlx, dx, -lx0*dx, 2, nlv, dz, -lv0*dz, 2, nz, dz, oz, 2, nx, dx, ox, 2);
        time_now = time(NULL);
        printf("done in %ld seconds \n", time_now-time_start);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Transform to Lambda-Dip domain
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Transforming H&V CIGS to LambdaDip...... ");
    time_start = time(NULL);
    float *eimage_LambdaDip                = CPU_zaloc1F(nlambda * nxLD * nzLD * ndip);
    float *eimage_LambdaDip_contracted     = CPU_zaloc1F(nlambda * nxLD * nzLD * ndip);
    float *eimage_LambdaDip_contracted_ref = CPU_zaloc1F(nlambda * nxLD * nzLD * ndip);
    GPU_WRAP_lxlz2LambdaDip(0, procGPU, 0, Gocig, eimage_LambdaDip, \
                            nx, dx, nz, dz, nxb, nzb, nxLD, oxLD, nzLD, ozLD, \
                            lx0, lv0, lambdaMax, dlambda, nlambda, dipMin, dipMax, ddip, ndip);
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Saving eimage_LambdaDip ...... ");
    time_start = time(NULL);
    if(outputs)
    {
        outputSmart4d(eimage_LambdaDip_FileName, eimage_LambdaDip, nlambda, dlambda, \
                     -lambdaMax, nzLD, dz, -0.5*nzLD*dz, nxLD, dx, -0.5*nxLD*dx, ndip, ddip, dipMin);
        // outputSmart4d_J(eimage_LambdaDip_FileName, eimage_LambdaDip, nlambda, dlambda, \
        //              -lambdaMax, 2, nzLD, dz, -ozLD, 2, nxLD, dx, -oxLD, 2, ndip, ddip, dipMin, 2);
    }
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

    // Apply focusing operator
    size_t dipSZ = nlambda * nxLD * nzLD;
    int idip;
    if(iproc==0)    printf("\n\n");
    for(idip=iproc; idip<ndip; idip+=nproc)
    {
        size_t shiftMem = dipSZ * idip;
        float *eimage_tmp = eimage_LambdaDip+shiftMem;

        //  >>>>> FOCUSING IMAGE <<<<<<
        float *eimage_tmp_diff           = CPU_zaloc1F(nxLD*nzLD*nlambda);
        float *eimage_tmp_contracted     = CPU_zaloc1F(nxLD*nzLD*nlambda);
        float *eimage_tmp_contracted_ref = CPU_zaloc1F(nxLD*nzLD*nlambda);
        
        time_start = time(NULL);
        if(iproc==0)    printf("\n Applying Focusing Operator for dip %+4.1f degress ...... ", dipMin+idip*ddip);
        float zMin = 0.0;
        float zMax = nzLD*dz;
        FocusingOperator_lambda(0, eimage_tmp_diff, eimage_tmp_contracted, eimage_tmp_contracted_ref, eimage_tmp, \
                                zMin, zMax, maxDepth, nxLD, nzLD, dx, dz, ox, oz, 0, 0, lambda0, ltheta0, dtheta, \
                                otheta, iter, muFocus, epsFocus, lambdaNull, 0, iproc, nproc);
        time_now = time(NULL);
        if(iproc==0)    printf(" done in %3ld seconds", time_now-time_start);

        time_start = time(NULL);
        memcpy(eimage_LambdaDip_contracted    +shiftMem, eimage_tmp_contracted    , nxLD*nzLD*nlambda*sizeof(float));
        memcpy(eimage_LambdaDip_contracted_ref+shiftMem, eimage_tmp_contracted_ref, nxLD*nzLD*nlambda*sizeof(float));

        free(eimage_tmp_diff);
        free(eimage_tmp_contracted);
        free(eimage_tmp_contracted_ref);
        // freeMem(eimage_tmp_diff,           nxLD*nzLD*nlambda*sizeof(float));
        // freeMem(eimage_tmp_contracted,     nxLD*nzLD*nlambda*sizeof(float));
        // freeMem(eimage_tmp_contracted_ref, nxLD*nzLD*nlambda*sizeof(float));

        // if(iproc==0)    printf("\n Filtering focused image ...... ");
        // filterECIG(eimage_diff          , eiproc, eimage_msh, msh);
        // filterECIG(eimage_contracted    , eiproc, eimage_msh, msh);
        // filterECIG(eimage_contracted_ref, eiproc, eimage_msh, msh);
        // if(iproc==0)    printf(" done \n");
    }
    if(iproc==0)    printf("\n\n");
    // REDUCE of focusing results
    time_start = time(NULL);
    for(idip=0; idip<ndip; idip++)
    {
        long long int nsz = nlambda * nxLD * nzLD;
        long long int shift = ((long long int) idip) * nsz;

        float *eimage_LambdaDip_contracted_tmp     = CPU_zaloc1F(nsz);
        float *eimage_LambdaDip_contracted_ref_tmp = CPU_zaloc1F(nsz);

        float *focus = eimage_LambdaDip_contracted+shift;
        float *refer = eimage_LambdaDip_contracted_ref+shift;

        MPI_Reduce(focus,  eimage_LambdaDip_contracted_tmp    , nsz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(refer,  eimage_LambdaDip_contracted_ref_tmp, nsz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

        if(iproc==0) {
            memcpy(focus,  eimage_LambdaDip_contracted_tmp    , nsz*sizeof(float));
            memcpy(refer,  eimage_LambdaDip_contracted_ref_tmp, nsz*sizeof(float));
        }
        else {
            memset(focus, 0, nsz*sizeof(float));
            memset(refer, 0, nsz*sizeof(float));
        }
        free(eimage_LambdaDip_contracted_tmp);
        free(eimage_LambdaDip_contracted_ref_tmp);
        // freeMem(eimage_LambdaDip_contracted_tmp    , nsz*sizeof(float));
        // freeMem(eimage_LambdaDip_contracted_ref_tmp, nsz*sizeof(float));
    
        MPI_Bcast(focus, nsz, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Bcast(refer, nsz, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }
    time_now = time(NULL);
    if(iproc==0)    printf("\n\n Performed all reduce operations in %ld seconds. \n\n", time_now-time_start);



    // Saving to file the results of focusing operation
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Saving eimage_LambdaDip_contracted and eimage_LambdaDip_contracted_ref ...... ");
    time_start = time(NULL);
    if(outputs) {
        outputSmart4d(eimage_LambdaDip_contracted_FileName    , eimage_LambdaDip_contracted    , \
                      nlambda, dlambda, -lambdaMax, nzLD, dz, -ozLD, nxLD, dx, -oxLD, ndip, ddip, dipMin);
        outputSmart4d(eimage_LambdaDip_contracted_ref_FileName, eimage_LambdaDip_contracted_ref, \
                        nlambda, dlambda, -lambdaMax, nzLD, dz, -ozLD, nxLD, dx, -oxLD, ndip, ddip, dipMin);
        // outputSmart4d_J(eimage_LambdaDip_contracted_FileName    , eimage_LambdaDip_contracted    , \
        //               nlambda, dlambda, -lambdaMax, 2, nzLD, dz, -ozLD, 2, nxLD, dx, -oxLD, 2, ndip, ddip, dipMin, 1);
        // outputSmart4d_J(eimage_LambdaDip_contracted_ref_FileName, eimage_LambdaDip_contracted_ref, \
        //                 nlambda, dlambda, -lambdaMax, 2, nzLD, dz, -ozLD, 2, nxLD, dx, -oxLD, 2, ndip, ddip, dipMin, 1);
    }
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);  if(iproc==0)    printf(" done in %ld seconds \n", time_now-time_start);

    // Transforming back to H&V CIGs
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Transforming eimageLambdaDip (focused and reference) back to Gocigs (focused and reference) ...... ");
    time_start = time(NULL);
    float *Gocig_contracted     = CPU_zaloc1F(nx*nz*nlx*nlv);
    float *Gocig_contracted_ref = CPU_zaloc1F(nx*nz*nlx*nlv);
    GPU_WRAP_lxlz2LambdaDip(1, procGPU, 0, Gocig_contracted    , eimage_LambdaDip_contracted    , \
                            nx, dx, nz, dz, nxb, nzb, nxLD, oxLD, nzLD, ozLD, \
                            lx0, lv0, lambdaMax, dlambda, nlambda, dipMin, dipMax, ddip, ndip);
    GPU_WRAP_lxlz2LambdaDip(1, procGPU, 0, Gocig_contracted_ref, eimage_LambdaDip_contracted_ref, \
                            nx, dx, nz, dz, nxb, nzb, nxLD, oxLD, nzLD, ozLD, \
                            lx0, lv0, lambdaMax, dlambda, nlambda, dipMin, dipMax, ddip, ndip);
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n done %ld seconds \n", time_now-time_start);

    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Transforming Gocig (focused and reference) back to H&V CIGs ...... ");
    time_start = time(NULL);
    // Gocig2VHocig(Gocig_contracted    , Hocig_contracted    , Vocig_contracted    , nx, nz, dx, dz, lx0, lv0, lv0, iproc);
    // Gocig2VHocig(Gocig_contracted_ref, Hocig_contracted_ref, Vocig_contracted_ref, nx, nz, dx, dz, lx0, lv0, lv0, iproc);
    Gocig2VHocig(Gocig_contracted_dbg    , Hocig_contracted    , Vocig_contracted    , nx, nz, dx, dz, lx0, lv0, lv0, iproc);
    Gocig2VHocig(Gocig_contracted_ref_dbg, Hocig_contracted_ref, Vocig_contracted_ref, nx, nz, dx, dz, lx0, lv0, lv0, iproc);
    freeMem(Gocig_contracted_ref_dbg, nx*nz*nlx*nlv*sizeof(float));
    freeMem(Gocig_contracted_dbg    , nx*nz*nlx*nlv*sizeof(float));
    freeMem(Gocig_diff_dbg, nx*nz*nlx*nlv*sizeof(float));
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n done %ld seconds \n", time_now-time_start);

    MPI_Barrier(MPI_COMM_WORLD);  if(iproc==0)    printf("\n\n Saving contracted H&V OCIGs ...... ");
    time_start = time(NULL);
    if(outputs)
    {
        outputSmart3d(Hocig_contracted_FileName    , Hocig_contracted    , nlx, dx, -lx0*dx,  nz, dz,      oz, nx, dx, ox);
        outputSmart3d(Hocig_contracted_ref_FileName, Hocig_contracted_ref, nlx, dx, -lx0*dx,  nz, dz,      oz, nx, dx, ox);
        outputSmart3d(Vocig_contracted_FileName    , Vocig_contracted    , nlv, dz, -lv0*dz,  nz, dz,      oz, nx, dx, ox);
        outputSmart3d(Vocig_contracted_ref_FileName, Vocig_contracted_ref, nlv, dz, -lv0*dz,  nz, dz,      oz, nx, dx, ox);
        // outputSmart4d(Gocig_contracted_FileName    , Gocig_contracted    , nlx, dx, -lx0*dx, nlv, dz, -lv0*dz, nz, dz, oz, nx, dx, ox);
        // outputSmart4d(Gocig_contracted_ref_FileName, Gocig_contracted_ref, nlx, dx, -lx0*dx, nlv, dz, -lv0*dz, nz, dz, oz, nx, dx, ox);
        // outputSmart4d_J(Gocig_contracted_FileName    , Gocig_contracted    , nlx, dx, -lx0*dx, 2, nlv, dz, -lv0*dz, 2, nz, dz, oz, 2, nx, dx, ox, 2);
        // outputSmart4d_J(Gocig_contracted_ref_FileName, Gocig_contracted_ref, nlx, dx, -lx0*dx, 2, nlv, dz, -lv0*dz, 2, nz, dz, oz, 2, nx, dx, ox, 2);
    }
    time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);
        
    free(Gocig_contracted);  
    free(Gocig_contracted_ref);        
    
    free(eimage_LambdaDip);
    free(eimage_LambdaDip_contracted);
    free(eimage_LambdaDip_contracted_ref);
    free(Hocig);
    free(Vocig);

    MPI_Barrier(MPI_COMM_WORLD);

    // multiplyArrayByScalar(0.0, Vocig_contracted_ref, nx*nz*nlv);
    // multiplyArrayByScalar(0.0, Vocig_contracted    , nx*nz*nlv);

    if(dataset_reference != NULL)
    {
        if(iproc==0)    printf("\n Modeling reference data \n");
        extendedBornModeling(dataset_reference, NULL, Hocig_contracted_ref, Vocig_contracted_ref, vp_bg, nx, nz, lx0, lv0, lt0, \
                             dx, dz, ox, oz, jt, perc, offMax, acquisitionGeom, NULL, iproc, nproc);

        int jshot = 10;
        writeSeisData2File(dataset_reference, jshot, iproc, refDataFileName, NULL);
        if(iproc==0)    printf("\n Modeling reference data  ----->>>>> concluded\n");
    }
    if(dataset_misfit    != NULL)
    {
        if(iproc==0)    printf("\n Modeling misfit data \n");
        extendedBornModeling(dataset_misfit, NULL, Hocig_contracted, Vocig_contracted, vp_bg, nx, nz, lx0, lv0, lt0, \
                             dx, dz, ox, oz, jt, perc, offMax, acquisitionGeom, NULL, iproc, nproc);
        int jshot = 10;
        writeSeisData2File(dataset_misfit   , jshot, iproc, misDataFileName, NULL);
        if(iproc==0)    printf("\n Modeling misfit data  ----->>>>> concluded\n");
    }

    //  >>>>> DATA RESIDUALS <<<<<<  
    float *Hocig_diff = CPU_zaloc1F(nx*nz*nlx);
    float *Vocig_diff = CPU_zaloc1F(nx*nz*nlv);
    addArrays(Hocig_diff, 1.0, Hocig_contracted, -1.0, Hocig_contracted_ref, nx*nz*nlx);
    addArrays(Vocig_diff, 1.0, Vocig_contracted, -1.0, Vocig_contracted_ref, nx*nz*nlv);
    lxNullTapering(Hocig_diff, Hocig_diff, nx, nz, dx, lx0, lambdaNull);
    double objFuncH = getArrayL2Norm(Hocig_diff, nx*nz*nlx);
    double objFuncV = getArrayL2Norm(Vocig_diff, nx*nz*nlv);
    float objFunc = objFuncH + 0*objFuncV;
    printf("\n\n objFuncH=%g  objFuncV=%g  objFunc=%g \n\n", objFuncH, objFuncV, objFunc);
    multiplyArrayByScalar(0.0, Vocig_diff, nx*nz*nlv);
    if(computingGradient)
    {
        if(iproc==0)    printf("\n\n Computing data residuals. \n (Shot,time): \n");
        extendedBornModeling(dataset_residuals, NULL, Hocig_diff, Vocig_diff, vp_bg, nx, nz, lx0, lv0, lt0, \
                             dx, dz, ox, oz, jt, perc, offMax, acquisitionGeom, NULL, iproc, nproc);
    }
    // if(iproc==0)    printf("\n\n Computing data residuals. \n (Shot,time): \n");
    // objFunc = extendedBornModeling(dataset_residuals, NULL, Hocig_diff, Vocig_diff, vp_bg, nx, nz, lx0, lv0, lt0, \
    //                                dx, dz, ox, oz, jt, perc, offMax, acquisitionGeom, NULL, iproc, nproc);
    MPI_Barrier(MPI_COMM_WORLD);

    free(image);
    free(Hocig_diff);
    free(Vocig_diff);
    

return objFunc;
}


float computeODCIGs(Dataset *dataset_input, Dataset *dataset_reference, float *eimage_contracted, float *eimage_contracted_ref, \
                    float *eimage_diff, float *vp_bg, \
                    int nx, int nz, int nt, int lx0, int lz0, int lt0, int jt_lt, int ltheta0, \
                    EIProc *eiproc, Eimage_MSH *eimage_msh, float dx, float dz, float dt, float dtheta, float ox, float oz, \
                    float ot, float otheta, float offMax, int acquisitionGeom, float perc, \
                    float zMin, float zMax, float maxDepth, float muFocus, float epsFocus, float lambdaNull, int iter, int ism, int jt, int outputs, \
                    int iproc, int nproc) {


    sprintf( refDataFileName , "seismogram_referenceData_Iter%d", iter);
    sprintf( misDataFileName , "seismogram_misfitData_Iter%d", iter);
    sprintf( resDataFileName , "seismogram_residuals_Iter%d", iter);
    
    float objFunc = 0.0;

    int nlx    = 2*lx0 + 1;
    int nlz    = 2*lz0 + 1;
    int ntheta = 2*ltheta0 + 1;

    zeroArray(eimage_contracted    , nx*nz*nlx*nlz);
    zeroArray(eimage_contracted_ref, nx*nz*nlx*nlz);
    zeroArray(eimage_diff, nx*nz*nlx*nlz);

    if(iproc==0)
        printf("\n\n  nlx=%d  nlz=%d  ntheta=%d  \n\n", nlx,  nlz,  ntheta);

    float  *image, *eimage, *eimage_theta;
    image               = CPU_zaloc1F(nx*nz);
    eimage              = CPU_zaloc1F(nx*nz*nlx*nlz);
    eimage_theta        = CPU_zaloc1F(nx*nz*ntheta);

    // >>>>>>>>>>> Migrate shots <<<<<<<<<<<
    if(iproc==0)    printf("\n\n Migrating data. \n");
    migration(dataset_input, image, eimage, vp_bg, nx, nz, lx0, lz0, lt0, jt_lt, ltheta0, dx, dz, \
              dtheta, ox, oz, otheta, jt, ism, 0, iproc, nproc);
    if(iproc==0)    printf("\n\n Migrating data -->> concluded \n");
    
    if(eiproc->applyBandass2EImage)
        lxApplyBandpass(eimage, nx, nz, dz, lx0, eiproc->lambda1, eiproc->lambda2, eiproc->lambda3, eiproc->lambda4);
    MPI_Barrier(MPI_COMM_WORLD);

    if(eiproc->applyXTaper)
        lxTapering_x(eimage, eimage, nx, nz, dx, lx0, eiproc->x1, eiproc->x2, eiproc->x3, eiproc->x4);
    if(eiproc->applyZTaper) {
        lxTapering(eimage, eimage, nx, nz, dz, lx0, lz0, eiproc->z0, eiproc->z1);
        lxTapering(eimage, eimage, nx, nz, dz, lx0, lz0, eiproc->z3, eiproc->z2);
    }
    if(eiproc->applyLXTaper)
        lxTapering_lx_varDepth(eimage, eimage, nx, nz, dx, lx0, \
                               eiproc->lx1_i, eiproc->lx2_i, eiproc->lx3_i, eiproc->lx4_i, \
                               eiproc->lx1_f, eiproc->lx2_f, eiproc->lx3_f, eiproc->lx4_f, eiproc->z_i, eiproc->z_f);
    
    GPU_WRAP_lx2theta(0, procGPU, 0, eimage, eimage_theta, nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);
    // Applying cosine to ADCIG (Jacobian of the Offset to Angle transformation)
    // applyCosADCIG(eimage_theta, eimage_theta, nx, nz, ntheta, dtheta);
    
    MPI_Barrier(MPI_COMM_WORLD);

    // Save image, ODCIGs, ADCIGs to file
    if(outputs)
    {
        // outputSmart2d(image_FileName,   image                                       , nz, dz, oz, nx, dx, ox);
        outputModel2d(image_FileName,   image                                       , nz, dz, 0.0, nx, dx, 0.0, nzb, nxb);
        outputSmart3d(eimage_FileName, eimage      , nlx   , dx    , -lx0*dx        , nz, dz, oz, nx, dx, ox);
        outputSmart3d(adcigs_FileName, eimage_theta, ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(iproc==0)    printf("\n\n Applying taper \n");

    taperADCIG(eimage_theta, eimage_theta, nx, nz, ntheta, dtheta, eiproc->theta2,  eiproc->theta1, 1.0);
    taperADCIG(eimage_theta, eimage_theta, nx, nz, ntheta, dtheta, eiproc->theta3,  eiproc->theta4, 1.0);
    taperADCIG(eimage_theta, eimage_theta, nx, nz, ntheta, dtheta, eiproc->theta4, +dtheta*ltheta0, 0.0);
    taperADCIG(eimage_theta, eimage_theta, nx, nz, ntheta, dtheta, eiproc->theta1, -dtheta*ltheta0, 0.0);
    
    if(iproc==0)    printf("\n Recovering ODCIG \n");
    MPI_Barrier(MPI_COMM_WORLD);
    float *eimage_taper_angle = cg_theta2lx_gpu(procGPU, eimage_theta, nx, nz, dx, dz, lx0, ntheta, dtheta, otheta, 5);
    MPI_Barrier(MPI_COMM_WORLD); 
    
    if(outputs)
    {
        outputSmart3d(eimage_taper_angle_FileName, eimage_taper_angle, nlx, dx, -lx0*dx, nz, dz, oz, nx, dx, ox);
        outputSmart3d(adcigs_taper_angle_FileName, eimage_theta, ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz, nx, dx, ox);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(iproc==0)    printf("\n Applying Focusing Operator \n");

    //  >>>>> FOCUSING IMAGE <<<<<<
    FocusingOperator_lx_lz_lt(0, eimage_diff, eimage_contracted, eimage_contracted_ref, eimage_taper_angle, zMin, zMax, maxDepth, \
                              nx, nz, dx, dz, ox, oz, lx0, lz0, lt0, ltheta0, dtheta, otheta, \
                              eimage_msh->dipMin, eimage_msh->dipMax, eimage_msh->ddip, eimage_msh->lambdaMax, eimage_msh->dlambda, \
                              eiproc->applyADFocusing, iter, muFocus, epsFocus, lambdaNull, outputs, iproc, nproc);
    free(eimage_taper_angle);

    // Applying DSO
    // applyLinearWeight_ODCIG(eimage_diff          , eimage_diff          , dx, nx, nz, lx0);
    // applyLinearWeight_ODCIG(eimage_contracted    , eimage_contracted    , dx, nx, nz, lx0);
    // applyLinearWeight_ODCIG(eimage_contracted_ref, eimage_contracted_ref, dx, nx, nz, lx0);

    MPI_Barrier(MPI_COMM_WORLD);

    if(iproc==0)    printf("\n\n Computing data residuals. \n (Shot,time): \n");
        
    extendedBornModeling(dataset_reference, eimage_contracted_ref, NULL, NULL, vp_bg, nx, nz, lx0, lz0, lt0, \
                         dx, dz, ox, oz, jt, perc, offMax, acquisitionGeom, NULL, iproc, nproc);

    free(image);
    free(eimage);
    free(eimage_theta);

return objFunc;
}

float computeObjFunc(Dataset *dataset_misfit, Dataset *dataset_reference, float *eimage_contracted, float *vp_bg, \
                     int nx, int nz, int nt, int lx0, int lz0, int lt0, int jt_lt, float dx, float dz, float dt, float ox, float oz, \
                     float ot, float offMax, int acquisitionGeom, float perc, int jt, int imageDomain, int iproc, int nproc, int iter) {
    
    float objFunc = 0.0;

    //  >>>>> DATA RESIDUALS <<<<<<
    if(imageDomain==0)
    {
        if(iproc==0)    printf("\n\n Computing data residuals. \n (Shot,time): \n");
        extendedBornModeling(dataset_misfit, eimage_contracted, NULL, NULL, vp_bg, nx, nz, lx0, lz0, lt0, \
                             dx, dz, ox, oz, jt, perc, offMax, acquisitionGeom, NULL, iproc, nproc);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    int nshot = dataset_misfit->nshot;
    int ishot;
    for(ishot=0; ishot<nshot; ishot++) {
        removeDA(&((dataset_misfit->shot)[ishot]), &((dataset_reference->shot)[ishot]));    
        objFunc += addObjFunc(&((dataset_misfit->shot)[ishot]));
    }
    
return objFunc;
}