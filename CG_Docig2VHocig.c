void cg_VHocig2Docig(int procGPU, float *Docig, float *Hocig, float *Vocig, long long int nx, long long int nz, float dx, float dz, \
                     long long int lx0, long long int lv0, long long int ld0, long long int ndip, float ddip, float dipMin, int niter, int iproc)
{
    float  *A, *xh, *xv, *bh, *bv, *rh, *rv, *ph, *pv, *Aph, *Apv;
    double  beta, beta_denomin, alpha, alpha_denominator, rr;
    long long int     k, nh, nv;
    long long int     i, ix, iz;
    long long int     nlx = 2*lx0+1;
    long long int     nlv = 2*lv0+1;
    long long int     nld = 2*ld0+1;



    nh = nx*nz*nlx;
    nv = nx*nz*nlv;
    bh = CPU_zaloc1F(nh);
    bv = CPU_zaloc1F(nv);
    if(procGPU<0)
        Docig2VHocig(Docig, bh, bv, nx, nz, ndip, dx, dz, ddip, dipMin, lx0, lv0, ld0, iproc);
    else
        GPU_WRAP_Docig2VHocig(procGPU, Docig, bh, bv, nx, nz, dx, dz, ndip, ddip, dipMin, lx0, lv0, ld0, iproc, EIMSH_CIGS_VCIG);
    
    // Initial conditions
    k    = 0;
    beta = 0.0f;
    
    xh  = CPU_zaloc1F(nh);
    xv  = CPU_zaloc1F(nv);
    rh  = CPU_zaloc1F(nh);
    rv  = CPU_zaloc1F(nv);
    ph  = CPU_zaloc1F(nh);
    pv  = CPU_zaloc1F(nv);
    Aph = CPU_zaloc1F(nh);
    Apv = CPU_zaloc1F(nv);
    
    // r iteration 0
    for(i=0; i<nh; i++)    rh[i] = bh[i];
    for(i=0; i<nv; i++)    rv[i] = bv[i];
    
    // rr iteration 0
    rr = 0.0f;
    for(i=0; i<nh; i++)    rr += rh[i]*rh[i];
    for(i=0; i<nv; i++)    rr += rv[i]*rv[i];

    // Iterate
    float rr_prev;
    for(k=0; k<niter; k++)
    {
        // if(iproc==0)    printf("\n [cg_VHocig2Docig] starting iteration %lld ", k);
        rr_prev = rr;
        // Update d
        for(i=0; i<nh; i++)    ph[i] = rh[i] + beta * ph[i];
        for(i=0; i<nv; i++)    pv[i] = rv[i] + beta * pv[i];
        
        // Get Ad (action of Hessian on the image)
        getAp_VHocig2Docig(procGPU, Aph, Apv, ph, pv, nx, nz, dx, dz, lx0, lv0, ld0, ndip, ddip, dipMin, iproc);
        
        // Update alpha
        alpha_denominator = 0.0f;

        for(i=0; i<nh; i++)    {alpha_denominator += ((double) ph[i]) * ((double) Aph[i]);}
        for(i=0; i<nv; i++)    {alpha_denominator += ((double) pv[i]) * ((double) Apv[i]);}

        alpha = rr/alpha_denominator;

        // if(iproc==0)    printf("\n [cg_VHocig2Docig] alpha_denominator=%g, alpha=%g - iteration %lld ", alpha_denominator, alpha, k);
        
        // Update x
        for(i=0; i<nh; i++)    xh[i] = xh[i] + alpha * ph[i];
        for(i=0; i<nv; i++)    xv[i] = xv[i] + alpha * pv[i];

        // Update r
        for(i=0; i<nh; i++)    rh[i] = rh[i] - alpha * Aph[i];
        for(i=0; i<nv; i++)    rv[i] = rv[i] - alpha * Apv[i];
        
        // Keep denominator of beta expression
        beta_denomin = rr;

        // Update rr
        rr = 0.0f;
        for(i=0; i<nh; i++)    rr += rh[i]*rh[i];
        for(i=0; i<nv; i++)    rr += rv[i]*rv[i];
        
        beta = rr/beta_denomin;
    }

    memcpy(Hocig, xh, nh*sizeof(float));
    memcpy(Vocig, xv, nv*sizeof(float));

    free(xh);
    free(xv);    
    free(rh);
    free(rv);
    free(ph);
    free(pv);
    free(Aph);
    free(Apv);
    free(bh);
    free(bv);

return;
}

void getAp_VHocig2Docig(int procGPU, float *Aph, float *Apv, float *ph, float *pv, long long int nx, long long int nz, float dx, float dz, \
                        long long int lx0, long long int lv0, long long int ld0, long long int ndip, float ddip, float dipMin, int iproc) {
    float *pfor;
    long long int nld = 2*ld0 + 1;
    pfor = CPU_zaloc1F(nx*nz*nld*ndip);    
    if(procGPU<0) {
        VHocig2Docig(pfor,  ph,  pv, nx, nz, ndip, dx, dz, ddip, dipMin, lx0, lv0, ld0, iproc); // forward
        Docig2VHocig(pfor, Aph, Apv, nx, nz, ndip, dx, dz, ddip, dipMin, lx0, lv0, ld0, iproc); // adjoint
    }
    else {
        GPU_WRAP_VHocig2Docig(procGPU, pfor,  ph,  pv, nx, nz, dx, dz, ndip, ddip, dipMin, lx0, lv0, ld0, iproc, EIMSH_CIGS_VCIG); // forward
        GPU_WRAP_Docig2VHocig(procGPU, pfor, Aph, Apv, nx, nz, dx, dz, ndip, ddip, dipMin, lx0, lv0, ld0, iproc, EIMSH_CIGS_VCIG); // adjoint
    }
    free(pfor);
return;
}