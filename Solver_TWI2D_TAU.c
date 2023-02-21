void Solver_TWI2D_TAU(Dataset3D *dataset, float *vel, InvParms *invparms, Parms *parms, MSH *msh, \
                      Eimage_MSH *eimsh, EIProc *eiproc, GRDProc *grdProc, int iproc, int nproc) {
    
    long long int nsz  = ((long long int) msh->nx) * ((long long int) msh->nz);
    long long int nszT = nsz * ((long long int) eimsh->ntau);
    float *image     = CPU_zaloc1F(nsz );
    float *ilumin    = CPU_zaloc1F(nsz );
    float *Tcig_dif  = CPU_zaloc1F(nszT);
    float *gradient  = CPU_zaloc1F(nsz );

    FILE *fHdr_ofunc, *fBin_ofunc;
    if(iproc==0)    initFilesSmart1d("ObjectiveFunction", &fHdr_ofunc, &fBin_ofunc, invparms->n_iter, 1, 1);

    int iter;
    for(iter=0; iter<invparms->n_iter; iter++) 
    {
        initializeFileNames(iter);

        if(iproc==0)    printf("\n\n Computing image residuals. \n");  fflush(stdout);
        int performEBM = 0;
        if(iter==0  &&  parms->performEBM)    performEBM = 1;

        double objFunc0 = ComputeResiduals_TAU(dataset, vel, Tcig_dif, image, ilumin, \
                                               msh, eimsh, eiproc, parms, performEBM, iproc, nproc);
        
        float OF = objFunc0;
        if(iproc==0)    writeSamples(&OF, fBin_ofunc, 1);
        
        if(iproc==0)  printf("\n\n Computing gradient. \n");  fflush(stdout);        
        MPI_Barrier(MPI_COMM_WORLD);
        ComputeGradient_TAU(dataset, gradient, ilumin, vel, Tcig_dif, msh, eimsh, eiproc, \
                            parms, grdProc, iter, iproc, nproc);
        
        if(iproc==0)    outputModel2d(gradient_FileName, gradient, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox, nzb, nxb);
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        if(parms->computeOnlyGrad_NoLineSearch)  break;

        if(iproc==0)    printf("\n\n Performing Line Search. \n");  fflush(stdout);
        float alpha_optm = LineSearch_TAU(dataset, vel, gradient, Tcig_dif, image, ilumin, \
                                          msh, eimsh, eiproc, parms, invparms->alpha, objFunc0, iproc, nproc);
        
        updateVel(vel, gradient, msh, alpha_optm);

        // Save updated velocity model
        if(iproc==0) {
            untransformVel(vel, msh->nx, msh->nz, msh->dx, msh->dt);
            outputModel2d(Vel_TWI_FileName, vel, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox, nzb, nxb);
            transformVel(vel, msh->nx, msh->nz, msh->dx, msh->dt);
        }

        if(iproc==0) {
            float numberIter = iter+1;
            FILE *fHdr_niter, *fBin_niter;
            initFilesSmart1d("NumberIterations",  &fHdr_niter, &fBin_niter, 1, 1, 1);
            writeSamples(&numberIter, fBin_niter, 1);
            fflush(fBin_niter);
            fclose(fBin_niter);
            fclose(fHdr_niter);
        }

    }
    free(image);
    free(ilumin);
    free(gradient);
    free(Tcig_dif);
    if(iproc==0) {
        fclose(fHdr_ofunc);
        fclose(fBin_ofunc);
    }
}

float LineSearch_TAU(Dataset3D *dataset, float *vel, float *gradient, float *Tcig_dif, float *image, float *ilumin, \
                     MSH *msh, Eimage_MSH *eimsh, EIProc *eiproc, Parms *parms, float alpha_max, double objFunc0, int iproc, int nproc)
{
    long long int nsz  = ((long long int) msh->nx) * ((long long int) msh->nz);
    float *vel_1 = CPU_zaloc1F(nsz);
    float *vel_2 = CPU_zaloc1F(nsz);

    float alpha_1 = 0.5*alpha_max;
    float alpha_2 = 1.0*alpha_max;

    memcpy(vel_1, vel, msh->nx * msh->nz * sizeof(float));
    memcpy(vel_2, vel, msh->nx * msh->nz * sizeof(float));
    updateVel(vel_1, gradient, msh, alpha_1);
    updateVel(vel_2, gradient, msh, alpha_2);

    if(iproc==0) {
        untransformVel(vel_1, msh->nx, msh->nz, msh->dx, msh->dt);
        untransformVel(vel_2, msh->nx, msh->nz, msh->dx, msh->dt);
        outputModel2d("vel_1", vel_1, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox, nzb, nxb);
        outputModel2d("vel_2", vel_2, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox, nzb, nxb);
        transformVel(vel_1, msh->nx, msh->nz, msh->dx, msh->dt);
        transformVel(vel_2, msh->nx, msh->nz, msh->dx, msh->dt);
    }

    int performEBM = 0;
    // double objFunc1 = ComputeResiduals_TAU_LS(dataset, vel_1, Tcig_dif, image, ilumin, \
    //                                           msh, eimsh, eiproc, parms, performEBM, iproc, nproc, 1);
    double objFunc1 = ComputeResiduals_TAU(dataset, vel_1, Tcig_dif, image, ilumin, \
                                              msh, eimsh, eiproc, parms, performEBM, iproc, nproc);
    
    // double objFunc2 = ComputeResiduals_TAU_LS(dataset, vel_2, Tcig_dif, image, ilumin, \
    //                                           msh, eimsh, eiproc, parms, performEBM, iproc, nproc, 2);
    double objFunc2 = ComputeResiduals_TAU(dataset, vel_2, Tcig_dif, image, ilumin, \
                                           msh, eimsh, eiproc, parms, performEBM, iproc, nproc);
    
    if(iproc==0)    printf("objFunc1=%f  objFunc2=%f\n", objFunc1, objFunc2);    fflush(stdout);

    float x0 = 0.0;
    float x1 = alpha_1;
    float x2 = alpha_2;
    float y0 = objFunc0;
    float y1 = objFunc1;
    float y2 = objFunc2;
    float ymin;
    float xmin = findParabolaMin(x0, x1, x2, y0, y1, y2, &ymin);

    if(iproc==0) {
        printf("\n\n Result from Line Search: \n");
        printf("x0   =%6.6g  y0   =%6.6g\n", x0, y0);
        printf("x1   =%6.6g  y1   =%6.6g\n", x1, y1);
        printf("x2   =%6.6g  y2   =%6.6g\n", x2, y2);
        printf("xmin =%6.6g  ymin =%6.6g  \n", xmin, ymin);
        printf("Effective value being used:    xmin = %6.6g \n", min(xmin,x2));
        printf(" <ymin is predicted value by parabolic paramters, not measured value by evaluation of objective function>  \n\n");
        fflush(stdout);
    }   

    free(vel_1);
    free(vel_2);

return min(xmin, x2);
}

void ComputeGradient_TAU(Dataset3D *dataset, float *gradient, float *ilumin, float *vel, float *Tcig_dif, \
                         MSH *msh, Eimage_MSH *eimsh, EIProc *eiproc, Parms *parms, GRDProc *grdProc, \
                         int iter, int iproc, int nproc) {
    
    MPI_Barrier(MPI_COMM_WORLD);
    memset(ilumin  , 0, msh->nx*msh->nz*sizeof(float));
    memset(gradient, 0, msh->nx*msh->nz*sizeof(float));
    int recordMovie = 0;
    int ishotMovie = 0;
    BornGradient_DataDomain_TAU(dataset, gradient, ilumin, vel, Tcig_dif, grdProc, \
                                msh, eimsh, parms, ishotMovie, recordMovie, iter, iproc, nproc);
    MPI_Barrier(MPI_COMM_WORLD);
return;
}

double ComputeResiduals_TAU(Dataset3D *dataset, float *vel, float *Tcig_dif, float *image, float *ilumin, \
                            MSH *msh, Eimage_MSH *eimsh, EIProc *eiproc, Parms *parms, int performEBM, int iproc, int nproc) {
    
    
    long long int nsz  = ((long long int) msh->nx) * ((long long int) msh->nz);
    long long int nszT = nsz * ((long long int) eimsh->ntau);
    float *Tcig_ref = CPU_zaloc1F(nszT);
    float *Tcig_foc = CPU_zaloc1F(nszT);

    memset(ilumin  , 0, nsz);
    memset(image   , 0, nsz);
    memset(Tcig_dif, 0, nszT);
    int recordMovie = 0;
    int ishotMovie = 0;
    migration_TAU(dataset, image, ilumin, Tcig_ref, vel, msh, eimsh, parms, ishotMovie, recordMovie, iproc, nproc);

    // >>> Image preconditioning - begin
    long long int rx=100, rz=5;
    modelSmooth(ilumin, ilumin, msh->nx, msh->nz, rx, rz);
    normalizeArray(ilumin,    nsz);
    normalizeArray(image,     nsz);
    // normalizeArray(Tcig_ref, nszT);
    float epsilon = 0.005;
    // applyIlum2EImage(Tcig_ref, Tcig_ref, image, image, ilumin, msh->nx, msh->nz, msh->dx, eimsh->tau0, epsilon);
    // normalizeArray( image, nx*nz);
    // normalizeArray(Tcig_ref, nszT);
    MPI_Barrier(MPI_COMM_WORLD);
    // >>> Image preconditioning - end


    if(iproc==0)    outputModel2d(image_FileName, image, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox, nzb, nxb);
    if(iproc==0)    outputSmart3d(Tcig_FileName, Tcig_ref, eimsh->ntau, eimsh->dtau, -eimsh->tau0*eimsh->dtau, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);

    if(iproc==0)    printf("\n\n Processing CIGs. \n");  fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    CIG_Preprocessing_TAU(Tcig_ref, eiproc, msh, eimsh, iproc, procGPU);
    if(iproc==0)    printf("\n\n Processing CIGs. Done \n");  fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    // Debug EBM
    if(performEBM) {
        int recordMovie = 0;
        int ism = 0;
        if(iproc==0)    printf("\n\n Performing extended Born modeling. \n");  fflush(stdout);
        ExtendedBornModeling_TAU(dataset, vel, Tcig_ref, msh, eimsh, parms, ism, recordMovie, iproc, nproc);
    }

    if(iproc==0)    outputSmart3d(Tcig_proc_FileName, Tcig_ref, eimsh->ntau, eimsh->dtau, -eimsh->tau0*eimsh->dtau, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);
    MPI_Barrier(MPI_COMM_WORLD);

    // Focusing Operator
    if(iproc==0)    printf("\n\n Applying focusing operator. \n");  fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    Apply_FocusingOperator_TAU(Tcig_ref, Tcig_foc, msh, eimsh, eiproc, iproc, nproc, procGPU);
    if(iproc==0)    printf("\n\n Applying focusing operator. Done\n");  fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    if(iproc==0)    outputSmart3d(Tcig_ref_FileName, Tcig_ref, eimsh->ntau, eimsh->dtau, -eimsh->tau0*eimsh->dtau, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);
    if(iproc==0)    outputSmart3d(Tcig_foc_FileName, Tcig_foc, eimsh->ntau, eimsh->dtau, -eimsh->tau0*eimsh->dtau, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);

    addArrays(Tcig_dif, 1.0, Tcig_foc, -1.0, Tcig_ref, nszT);
    double ObjFunc = getArrayL2Norm(Tcig_dif, nszT);
    if(iproc==0)    outputSmart3d(Tcig_dif_FileName, Tcig_dif, eimsh->ntau, eimsh->dtau, -eimsh->tau0*eimsh->dtau, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);

    MPI_Barrier(MPI_COMM_WORLD);

    free(Tcig_ref);
    free(Tcig_foc);
    

return ObjFunc;
}

double ComputeResiduals_TAU_LS(Dataset3D *dataset, float *vel, float *Tcig_dif, float *image, float *ilumin, \
                               MSH *msh, Eimage_MSH *eimsh, EIProc *eiproc, Parms *parms, int performEBM, int iproc, int nproc, int suffix) {
    
    
    long long int nsz  = ((long long int) msh->nx) * ((long long int) msh->nz);
    long long int nszT = nsz * ((long long int) eimsh->ntau);
    float *Tcig_ref = CPU_zaloc1F(nszT);
    float *Tcig_foc = CPU_zaloc1F(nszT);

    memset(ilumin  , 0, nsz);
    memset(image   , 0, nsz);
    memset(Tcig_dif, 0, nszT);
    int recordMovie = 0;
    int ishotMovie = 0;
    migration_TAU(dataset, image, ilumin, Tcig_ref, vel, msh, eimsh, parms, ishotMovie, recordMovie, iproc, nproc);

    // Saving file for debug 
    if(iproc==0) {
        char fileName[2048];
        sprintf(fileName, "Tcig_LS_%d", suffix);
        outputSmart3d(fileName, Tcig_ref, eimsh->ntau, eimsh->dtau, -eimsh->tau0*eimsh->dtau, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);
        sprintf(fileName, "image_LS_%d", suffix);
        outputModel2d(fileName, image, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox, nzb, nxb);

        if(iproc==0) {
            sprintf(fileName, "vel_LS_%d", suffix);
            untransformVel(vel, msh->nx, msh->nz, msh->dx, msh->dt);
            outputModel2d(fileName, vel, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox, nzb, nxb);
            transformVel(vel, msh->nx, msh->nz, msh->dx, msh->dt);
        }

    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(iproc==0)    printf("\n\n Processing CIGs. \n");  fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    CIG_Preprocessing_TAU(Tcig_ref, eiproc, msh, eimsh, iproc, procGPU);
    if(iproc==0)    printf("\n\n Processing CIGs. Done \n");  fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    // Saving file for debug 
    if(iproc==0) {
        char fileName[2048];
        sprintf(fileName, "Tcig_proc_LS_%d", suffix);
        outputSmart3d(fileName, Tcig_ref, eimsh->ntau, eimsh->dtau, -eimsh->tau0*eimsh->dtau, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Focusing Operator
    if(iproc==0)    printf("\n\n Applying focusing operator. \n");  fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    Apply_FocusingOperator_TAU(Tcig_ref, Tcig_foc, msh, eimsh, eiproc, iproc, nproc, procGPU);
    if(iproc==0)    printf("\n\n Applying focusing operator. Done\n");  fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    // Saving file for debug
    if(iproc==0) {
        char fileName[2048];
        sprintf(fileName, "Tcig_ref_LS_%d", suffix);
        outputSmart3d(fileName, Tcig_ref, eimsh->ntau, eimsh->dtau, -eimsh->tau0*eimsh->dtau, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);
        sprintf(fileName, "Tcig_foc_LS_%d", suffix);
        outputSmart3d(fileName, Tcig_foc, eimsh->ntau, eimsh->dtau, -eimsh->tau0*eimsh->dtau, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);
    }

    addArrays(Tcig_dif, 1.0, Tcig_foc, -1.0, Tcig_ref, nszT);
    double ObjFunc = getArrayL2Norm(Tcig_dif, nszT);
    if(iproc==0) {
        char fileName[2048];
        sprintf(fileName, "Tcig_dif_LS_%d", suffix);
        outputSmart3d(fileName, Tcig_dif, eimsh->ntau, eimsh->dtau, -eimsh->tau0*eimsh->dtau, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    free(Tcig_ref);
    free(Tcig_foc);
    

return ObjFunc;
}

void CIG_Preprocessing_TAU(float *Tcig, EIProc *eiproc, MSH *msh, Eimage_MSH *eimsh, int iproc, int procGPU)
{
    // H&V CIG filterings
    filterTCIG(Tcig, eiproc, eimsh, msh);
    if(iproc==0)    outputSmart3d("Tcig_filtered", Tcig, eimsh->ntau, eimsh->dtau, -eimsh->tau0*eimsh->dtau, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);
}

void Apply_FocusingOperator_TAU(float *Tcig_ref, float *Tcig_foc, MSH *msh, \
                                Eimage_MSH *eimsh, EIProc *eiproc, int iproc, int nproc, int procGPU)
{
    long long int nsz  = ((long long int) msh->nx) * ((long long int) msh->nz);
    long long int nszT = nsz * ((long long int) eimsh->ntau);
    float *Tcig = CPU_zaloc1F(nszT);
    
    memcpy(Tcig, Tcig_ref, nszT*sizeof(float));
    
    // Focusing Operator onto Tcig
    time_t time_start = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Applying focusing operator to Tcig ...... "); 
    int outputs = 0;
    FocusingOperator_Tcig(0, Tcig_foc, Tcig_ref, Tcig, eiproc->zMin, eiproc->zMax, \
                          msh->nx, msh->nz, msh->dx, msh->dz, msh->ox, msh->oz, msh->nxb, msh->nzb, eimsh->tau0, \
                          eiproc->muFocus, eiproc->tauNull, outputs, iproc, nproc);

    MPI_Barrier(MPI_COMM_WORLD);    
    if(iproc==0)    printf("done in %ld seconds \n", time(NULL)-time_start);
    
    
    free(Tcig);
}