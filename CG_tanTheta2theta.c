void getAp_tanTheta2theta(float *Ap, float *p, long long int nx, long long int nz, float dx, float dz, \
                          float xMin, float xMax, float zMin, float zMax, \
                          long long int ltheta0, float dtheta, long long int ltanth0, float dtanth, int iproc);

void cg_tanTheta2theta(float *tdcig, float *adcig, long long int nx, long long int nz, float dx, float dz, float xMin, float xMax, float zMin, float zMax, \
                       long long int ltheta0, float dtheta, long long int ltanth0, float dtanth, int niter, int iproc)
{
    float  *A, *x, *b, *r, *p, *Ap;
    double  beta, beta_denomin, alpha, alpha_denominator, rr;
    long long int     k, n;
    long long int     i, ix, iz;
    long long int     ntheta = 2*ltheta0+1;
    long long int     ntanth = 2*ltanth0+1;

    n = nx*nz*ntheta;
    b = CPU_zaloc1F(n);
    GPU_WRAP_adcig_theta2tanTheta(1, procGPU, 0, b, tdcig, nx, dx, nz, dz, ltheta0, dtheta, ltanth0, dtanth, nxb, nzb, xMin, xMax, zMin, zMax);
    
    // Initial conditions
    k    = 0;
    beta = 0.0f;
    
    x  = CPU_zaloc1F(n);
    r  = CPU_zaloc1F(n);
    p  = CPU_zaloc1F(n);
    Ap = CPU_zaloc1F(n);
    
    // r iteration 0
    for(i=0; i<n; i++)    r[i] = b[i];
    
    // rr iteration 0
    rr = 0.0f;
    for(i=0; i<n; i++)    rr += r[i]*r[i];

    // Iterate
    float rr_prev;
    for(k=0; k<niter; k++)
    {
        if(iproc==0)    printf("\n [cg_VHocig2Docig] starting iteration %lld ", k);
        rr_prev = rr;
        // Update d
        for(i=0; i<n; i++)    p[i] = r[i] + beta * p[i];
        
        // Get Ad (action of Hessian on the image)
        getAp_tanTheta2theta(Ap, p, nx, nz, dx, dz, xMin, xMax, zMin, zMax, ltheta0, dtheta, ltanth0, dtanth, iproc);
        
        // Update alpha
        alpha_denominator = 0.0f;

        for(i=0; i<n; i++)    
            alpha_denominator += ((double) p[i]) * ((double) Ap[i]);

        alpha = rr/alpha_denominator;

        if(iproc==0)    printf("\n [cg_VHocig2Docig] alpha_denominator=%g, alpha=%g - iteration %lld ", alpha_denominator, alpha, k);
        
        // Update x
        for(i=0; i<n; i++)    x[i] = x[i] + alpha * p[i];

        // Update r
        for(i=0; i<n; i++)    r[i] = r[i] - alpha * Ap[i];
        
        // Keep denominator of beta expression
        beta_denomin = rr;

        // Update rr
        rr = 0.0f;
        for(i=0; i<n; i++)    rr += r[i]*r[i];
        
        beta = rr/beta_denomin;
    }

    memcpy(adcig, x, n*sizeof(float));

    free(x);
    free(r);
    free(p);
    free(Ap);
    free(b);

return;
}

void getAp_tanTheta2theta(float *Ap, float *p, long long int nx, long long int nz, float dx, float dz, \
                          float xMin, float xMax, float zMin, float zMax, \
                          long long int ltheta0, float dtheta, long long int ltanth0, float dtanth, int iproc) {
    float *pfor;
    long long int ntanth = 2*ltanth0 + 1;
    pfor = CPU_zaloc1F(nx*nz*ntanth);    
    GPU_WRAP_adcig_theta2tanTheta(0, procGPU, 0, p , pfor, nx, dx, nz, dz, \
                                  ltheta0, dtheta, ltanth0, dtanth, \
                                  nxb, nzb, xMin, xMax, zMin, zMax); // forward
    GPU_WRAP_adcig_theta2tanTheta(1, procGPU, 0, Ap, pfor, nx, dx, nz, dz, \
                                  ltheta0, dtheta, ltanth0, dtanth, \
                                  nxb, nzb, xMin, xMax, zMin, zMax); // adjoint
    free(pfor);
return;
}