
void getAp_lx2theta_gpu(int idev, float *Ad, float *d, int nx, int nz, float dx, float dz, int lx0, int ntheta, float dtheta, float otheta);

float* cg_theta2lx_gpu(int idev, float *img_th, int nx_, int nz_, float dx_, float dz_, \
                       int lx0_, int ntheta_, float dtheta_, float otheta_, int niter)
{
    float  *A, *x, *b, *r, *p, *Ap;
    double  beta, beta_denomin, alpha, alpha_denominator, rr;
    int     k, n;
    int     i, ix, iz;
   
    int   nx     = nx_;
    int   nz     = nz_;
    float dx     = dx_; 
    float dz     = dz_;
    int   lx0    = lx0_;
    int   ntheta = ntheta_;
    float dtheta = dtheta_;
    float otheta = otheta_;
    int   nlx    = 2*lx0+1;

    // printf("\n Starting Conjugate Gradient to transform from ADCIG to ODCIG \n ");
    // printf("Residuals:  ");
    n = nx*nz*nlx;
    b = CPU_zaloc1F(n);
    GPU_WRAP_lx2theta(1, idev, 0, b, img_th, nx, dx, nz, dz, nxb, nzb, lx0, ntheta, dtheta, otheta);
    
    // Initializing algorithm
    // k      = 0
    // x0     = 0.0 (my choice)
    // p_{-1} = 0.0
    // beta_0 = 0.0;
    // r_0    = b - Ax0; >> r_0 = b; (as of my choice)
    
    // Initial conditions
    k    = 0;
    beta = 0.0f;
    
    x  = CPU_zaloc1F(n);
    r  = CPU_zaloc1F(n);
    p  = CPU_zaloc1F(n);
    Ap = CPU_zaloc1F(n);
    
    // r iteration 0
    for(i=0; i<n; i++)
        r[i] = b[i];
    
    // rr iteration 0
    rr = 0.0f;
    for(i=0; i<n; i++)
        rr += r[i]*r[i];

    // Iterate
    float rr_prev;
    for(k=0; k<niter; k++)
    {
        rr_prev = rr;
        // Update d
        for(i=0; i<n; i++)
            p[i] = r[i] + beta * p[i];
        
        // Get Ad (action of Hessian on the image)
        getAp_lx2theta_gpu(idev, Ap, p, nx, nz, dx, dz, lx0, ntheta, dtheta, otheta);
        
        // Update alpha
        alpha_denominator = 0.0f;
        for(i=0; i<n; i++)
            alpha_denominator += p[i] * Ap[i];
        alpha = rr/alpha_denominator;
        
        // Update x
        for(i=0; i<n; i++)
            x[i] = x[i] + alpha * p[i];

        // Update r
        for(i=0; i<n; i++)
            r[i] = r[i] - alpha * Ap[i];
        
        // Keep denominator of beta expression
        beta_denomin = rr;

        // Update rr
        rr = 0.0f;
        for(i=0; i<n; i++)
            rr += r[i]*r[i];
        
        // Norm of vectors
        double norm_x = 0.0;
        double norm_b = 0.0;
        for(i=0; i<n; i++) {
            norm_x += x[i]*x[i];
            norm_b += b[i]*b[i];
        }
        norm_x = sqrt(norm_x);
        norm_b = sqrt(norm_b);
        
        beta = rr/beta_denomin;        
        // printf("    iter %d, rr=%g", k, rr);        
    }
    
    free(r);
    free(p);
    free(Ap);
    free(b);

return x;
}

void getAp_lx2theta_gpu(int idev, float *Ap, float *p, int nx, int nz, float dx, float dz, int lx0, int ntheta, float dtheta, float otheta) {
    float *pfor;
    pfor = CPU_zaloc1F(nx * nz * ntheta);
    GPU_WRAP_lx2theta(0, idev, 0,  p, pfor, nx, dx, nz, dz, nxb, nzb, lx0, ntheta, dtheta, otheta);
    GPU_WRAP_lx2theta(1, idev, 0, Ap, pfor, nx, dx, nz, dz, nxb, nzb, lx0, ntheta, dtheta, otheta);
    free(pfor);
return;
}