void Solver_TWI2D(Dataset3D *dataset, float *vel, InvParms *invparms, Parms *parms, MSH *msh, \
                  Eimage_MSH *eimsh, EIProc *eiproc, GRDProc *grdProc, int iproc, int nproc) {
    
    long long int nsz  = ((long long int) msh->nx) * ((long long int) msh->nz);
    long long int nszH = nsz * ((long long int) eimsh->nlx);
    long long int nszV = nsz * ((long long int) eimsh->nlv);
    float *image     = CPU_zaloc1F(nsz );
    float *Hocig_dif = CPU_zaloc1F(nszH);
    float *Vocig_dif = CPU_zaloc1F(nszV);
    
    float *gradient  = CPU_zaloc1F(nsz);

    int iter;
    for(iter=0; iter<invparms->n_iter; iter++) 
    {
        initializeFileNames(iter);

        if(iproc==0)    printf("\n\n Computing image residuals. \n");
        double objFunc0 = ComputeResiduals_ORTH(dataset, vel, Hocig_dif, Vocig_dif, image, \
                                                msh, eimsh, eiproc, parms, iproc, nproc);
        
        if(iproc==0)    printf("\n\n Computing gradient. \n");  fflush(stdout);
        ComputeGradient_ORTH(dataset, gradient, vel, Hocig_dif, Vocig_dif, msh, eimsh, eiproc, parms, 
                             grdProc, iter, iproc, nproc);
        
        if(iproc==0)    outputModel2d(gradient_FileName, gradient, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox, nzb, nxb);

        MPI_Barrier(MPI_COMM_WORLD);
        if(parms->computeOnlyGrad_NoLineSearch) {
            free(Hocig_dif);
            free(Vocig_dif);
            return;
        }

        if(iproc==0)    printf("\n\n Performing Line Search. \n");  fflush(stdout);
        float alpha_optm = LineSearch_ORTH(dataset, vel, gradient, Hocig_dif, Vocig_dif, image, \
                                           msh, eimsh, eiproc, parms, invparms->alpha, objFunc0, \
                                           iproc, nproc);
        
        updateVel(vel, gradient, msh, alpha_optm);

        // Save updated velocity model
        if(iproc==0) {
            untransformVel(vel, msh->nx, msh->nz, msh->dx, msh->dt);
            outputModel2d(Vel_TWI_FileName, vel, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox, nzb, nxb);
            transformVel(vel, msh->nx, msh->nz, msh->dx, msh->dt);
        }
        
    }

    free(Hocig_dif);
    free(Vocig_dif);
    
    

}

void updateVel(float *vel, float *gradient, MSH *msh, float alpha_optm)
{
    long long int nsz  = ((long long int) msh->nx) * ((long long int) msh->nz);
    untransformVel(vel, msh->nx, msh->nz, msh->dx, msh->dt);
    addArrays(vel, 1.0, vel, -alpha_optm, gradient, nsz);
    transformVel(vel, msh->nx, msh->nz, msh->dx, msh->dt);
}

float LineSearch_ORTH(Dataset3D *dataset, float *vel, float *gradient, float *Hocig_dif, float *Vocig_dif, float *image, \
                      MSH *msh, Eimage_MSH *eimsh, EIProc *eiproc, Parms *parms, float alpha_max, double objFunc0, \
                      int iproc, int nproc)
{
    long long int nsz  = ((long long int) msh->nx) * ((long long int) msh->nz);
    float *vel_1 = CPU_zaloc1F(nsz );
    float *vel_2 = CPU_zaloc1F(nsz );

    float alpha_1 = 0.5*alpha_max;
    float alpha_2 = 1.0*alpha_max;

    untransformVel(vel, msh->nx, msh->nz, msh->dx, msh->dt);
    addArrays(vel_1, 1.0, vel, -alpha_1, gradient, nsz);
    addArrays(vel_2, 1.0, vel, -alpha_2, gradient, nsz);
    transformVel(vel  , msh->nx, msh->nz, msh->dx, msh->dt);
    transformVel(vel_1, msh->nx, msh->nz, msh->dx, msh->dt);
    transformVel(vel_2, msh->nx, msh->nz, msh->dx, msh->dt);

    double objFunc1 = ComputeResiduals_ORTH(dataset, vel_1, Hocig_dif, Vocig_dif, image, \
                                            msh, eimsh, eiproc, parms, iproc, nproc);
    
    double objFunc2 = ComputeResiduals_ORTH(dataset, vel_2, Hocig_dif, Vocig_dif, image, \
                                            msh, eimsh, eiproc, parms, iproc, nproc);
    
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

return min(xmin,x2);
}

void ComputeGradient_ORTH(Dataset3D *dataset, float *gradient, float *vel, float *Hocig_dif, float *Vocig_dif, \
                          MSH *msh, Eimage_MSH *eimsh, EIProc *eiproc, Parms *parms, GRDProc *grdProc, \
                          int iter, int iproc, int nproc) {

    long long int lv0_tmp;
    if(EIMSH_CIGS_VCIG)  lv0_tmp = eimsh->lv0;
    else                 lv0_tmp = 0;
    
    BornGradient_DataDomain(dataset, gradient, vel, NULL, Hocig_dif, Vocig_dif, grdProc, \
                            msh, eimsh, parms, 0, 0, iter, iproc, nproc);
return;
}

double ComputeResiduals_ORTH(Dataset3D *dataset, float *vel, float *Hocig_dif, float *Vocig_dif, float *image, \
                             MSH *msh, Eimage_MSH *eimsh, EIProc *eiproc, Parms *parms, int iproc, int nproc) {
    
    long long int nsz  = ((long long int) msh->nx) * ((long long int) msh->nz);
    long long int nszH = nsz * ((long long int) eimsh->nlx);
    long long int nszV = nsz * ((long long int) eimsh->nlv);
    float *Hocig_ref = CPU_zaloc1F(nszH);
    float *Hocig_foc = CPU_zaloc1F(nszH);
    float *Vocig_ref = CPU_zaloc1F(nszV);
    float *Vocig_foc = CPU_zaloc1F(nszV);

    migration_HOCIG_VOCIG(dataset, image, Hocig_ref, Vocig_ref, vel, msh, eimsh, parms, 0, 0, iproc, nproc);
    
    if(iproc==0)    outputModel2d(image_FileName, image, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox, nzb, nxb);
    if(iproc==0)    outputSmart3d(Hocig_FileName, Hocig_ref, eimsh->nlx, eimsh->dlx, -eimsh->lx0*eimsh->dlx, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);
    if(iproc==0)    outputSmart3d(Vocig_FileName, Vocig_ref, eimsh->nlv, eimsh->dlv, -eimsh->lv0*eimsh->dlv, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);

    CIG_Preprocessing_ORTH(Hocig_ref, Vocig_ref, eiproc, msh, eimsh, iproc, procGPU);

    if(iproc==0)    outputSmart3d(Hocig_proc_FileName, Hocig_ref, eimsh->nlx, eimsh->dlx, -eimsh->lx0*eimsh->dlx, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);
    if(iproc==0)    outputSmart3d(Vocig_proc_FileName, Vocig_ref, eimsh->nlv, eimsh->dlv, -eimsh->lv0*eimsh->dlv, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);
    MPI_Barrier(MPI_COMM_WORLD);

    // Focusing Operator
    Apply_FocusingOperator_ORTH(Hocig_ref, Hocig_foc, Vocig_ref, Vocig_foc, msh, eimsh, eiproc, iproc, nproc, procGPU);
    if(iproc==0)    outputSmart3d(Hocig_ref_FileName, Hocig_ref, eimsh->nlx, eimsh->dlx, -eimsh->lx0*eimsh->dlx, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);
    if(iproc==0)    outputSmart3d(Hocig_foc_FileName, Hocig_foc, eimsh->nlx, eimsh->dlx, -eimsh->lx0*eimsh->dlx, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);
    if(iproc==0)    outputSmart3d(Vocig_ref_FileName, Vocig_ref, eimsh->nlv, eimsh->dlv, -eimsh->lv0*eimsh->dlv, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);
    if(iproc==0)    outputSmart3d(Vocig_foc_FileName, Vocig_foc, eimsh->nlv, eimsh->dlv, -eimsh->lv0*eimsh->dlv, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);

    addArrays(Hocig_dif, 1.0, Hocig_foc, -1.0, Hocig_ref, nszH);
    addArrays(Vocig_dif, 1.0, Vocig_foc, -1.0, Vocig_ref, nszV);
    double ObjFunc = getArrayL2Norm(Hocig_dif, nszH) + getArrayL2Norm(Vocig_dif, nszV);
    if(iproc==0)    outputSmart3d(Hocig_dif_FileName, Hocig_dif, eimsh->nlx, eimsh->dlx, -eimsh->lx0*eimsh->dlx, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);
    if(iproc==0)    outputSmart3d(Vocig_dif_FileName, Vocig_dif, eimsh->nlv, eimsh->dlv, -eimsh->lv0*eimsh->dlv, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);

    free(Hocig_ref);
    free(Hocig_foc);
    free(Vocig_ref);
    free(Vocig_foc);

return ObjFunc;
}

void CIG_Preprocessing_ORTH(float *Hocig, float *Vocig, EIProc *eiproc, MSH *msh, Eimage_MSH *eimsh, int iproc, int procGPU)
{
    //  H&V CIG filterings  
    filterECIG(Hocig, eiproc, eimsh, msh);
    if(EIMSH_CIGS_VCIG) {
        int lx0 = eimsh->lx0;
        eimsh->lx0 = eimsh->lv0;
        filterECIG(Vocig, eiproc, eimsh, msh);
        eimsh->lx0 = lx0;
    }
    if(iproc==0)    outputSmart3d("Hocig_filtered", Hocig, eimsh->nlx, eimsh->dlx, -eimsh->lx0*msh->dx, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);

    float *Hocig_taper_angle = filterECIG_Angle(Hocig, eiproc, msh->nx, msh->nz, msh->nxb, msh->nzb, \
                                                msh->dx, msh->dz, msh->ox, msh->oz, \
                                                eimsh->ltheta0, eimsh->dtheta, eimsh->lx0, iproc, procGPU);
    float *Vocig_taper_angle = NULL;
    if(EIMSH_CIGS_VCIG) {
           Vocig_taper_angle = filterECIG_Angle(Vocig, eiproc, msh->nx, msh->nz, msh->nxb, msh->nzb, \
                                                msh->dx, msh->dz, msh->ox, msh->oz, \
                                                eimsh->ltheta0, eimsh->dtheta, eimsh->lv0, iproc, procGPU);
    }
    if(Hocig_taper_angle!=NULL) {
        memcpy(Hocig, Hocig_taper_angle, msh->nx*msh->nz*eimsh->nlx*sizeof(float));
        free(Hocig_taper_angle);
    }
    if(Vocig_taper_angle!=NULL) {
        memcpy(Vocig, Vocig_taper_angle, msh->nx*msh->nz*eimsh->nlv*sizeof(float));
        free(Vocig_taper_angle);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(iproc==0)    outputSmart3d("Hocig_filtered_1", Hocig, eimsh->nlx, eimsh->dlx, -eimsh->lx0*msh->dx, msh->nz, msh->dz, msh->oz, msh->nx, msh->dx, msh->ox);

}

void Apply_FocusingOperator_ORTH(float *Hocig_ref, float *Hocig_foc, float *Vocig_ref, float *Vocig_foc, \
                                 MSH *msh, Eimage_MSH *eimsh, EIProc *eiproc, int iproc, int nproc, int procGPU)
{
    long long int nx=msh->nx, nz=msh->nz, nxb=msh->nxb, nzb=msh->nzb;
    float         dx=msh->dx, dz=msh->dz, ox=msh->ox, oz=msh->oz;
    float         lambdaMax=eimsh->lambdaMax, dlambda=eimsh->dlambda, ddip=eimsh->ddip, dipMin=eimsh->dipMin, dipMax=eimsh->dipMax;
    long long int nlambda=eimsh->nlambda, lambda0=(nlambda-1)/2, ndip=eimsh->ndip, lx0=eimsh->lx0, lv0=eimsh->lv0;

    // Transform to DOCIG
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Transforming H&V CIGS to Docig ...... ");
    time_t time_start = time(NULL);
    float *Docig     = CPU_zaloc1F(nx*nz*nlambda*ndip);
    GPU_WRAP_VHocig2Docig(procGPU, Docig, Hocig_ref, Vocig_ref, nx, nz, dx, dz, ndip, ddip, dipMin, lx0, lv0, lambda0, iproc, EIMSH_CIGS_VCIG);
    time_t time_now = time(NULL);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);

    
    // Focusing Operator onto Docig
    float *Docig_ref  = CPU_zaloc1F(nx*nz*nlambda*ndip);
    float *Docig_foc  = CPU_zaloc1F(nx*nz*nlambda*ndip);
    float *Docig_dif = CPU_zaloc1F(nx*nz*nlambda*ndip);
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Applying focusing operator to Docig ...... "); time_start = time(NULL);
    int outputs = 0;
    FocusingOperator_Docig(0, Docig_dif, Docig_foc, Docig_ref, Docig, eiproc->zMin, eiproc->zMax, \
                           nx, nz, dx, dz, ox, oz, nxb, nzb, lambda0, ndip, dipMin, ddip, \
                           eiproc->muFocus, eiproc->epsFocus, eiproc->lambdaNull, outputs, iproc, nproc);
    MPI_Barrier(MPI_COMM_WORLD);    time_now = time(NULL); if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);
    
    
    // Transforming Docig back to H&Vocig
    MPI_Barrier(MPI_COMM_WORLD);    if(iproc==0)    printf("\n Transforming focused and reference Docigs to HVocig ...... "); time_start = time(NULL);
    int niterCG = 8;  // number of CG iterations to transform back from Docig to HVocig
    cg_VHocig2Docig(procGPU, Docig_ref, Hocig_ref, Vocig_ref, nx, nz, dx, dz, lx0, lv0, lambda0, ndip, ddip, dipMin, niterCG, iproc);
    cg_VHocig2Docig(procGPU, Docig_foc, Hocig_foc, Vocig_foc, nx, nz, dx, dz, lx0, lv0, lambda0, ndip, ddip, dipMin, niterCG, iproc);
    MPI_Barrier(MPI_COMM_WORLD);    time_now = time(NULL); if(iproc==0)    printf("done in %ld seconds \n", time_now-time_start);
    
    free(Docig_ref);
    free(Docig_foc);
    free(Docig_dif);
}