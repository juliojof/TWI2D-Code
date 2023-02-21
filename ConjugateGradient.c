void normalizeB(float *b, int n);
void getAD(float *Ad, float *d, Dataset *dataset, float *vp_bg, int jt, int nt, int nx, int nz, float dx, float dz, \
           float dt, float ox, float oz, float ot, int lx0, int iter);

float* cg(float *image, Dataset *dataset, float *vp_bg, int jt, int nt, int nx, int nz, float dx, \
          float dz, float dt, float ox, float oz, float ot, int lx0, int niter)
{
    float  *A, *x, *b, *r, *d, *Ad;
    float   beta, beta_denomin, alpha, alpha_denominator, rr;
    int     k, n, nlx;
    int     i, ix, iz;
    
    printf("\n Starting Conjugate Gradient \n ");
    nlx = 2*lx0 + 1;
    n = nx*nz*nlx;
    b = image;

    //normalizeB(b,n);

    printf("\n\n Value of n: %d \n\n", n);

    // Save to file inverted extended reflectivity
    FILE *fImg_Hdr, *fImg_Bin;
    initFiles3d("./debug_b", &fImg_Hdr, &fImg_Bin, nlx, dx, -lx0*dx, nz, dz, oz, nx, dx, ox);
    writeSamples(b, fImg_Bin, nlx*nz*nx);
    fclose(fImg_Bin);
    fclose(fImg_Hdr);
    
    // Initializing algoritmo
    // k      = 0
    // x0     = 0.0 (my choice)
    // p_{-1} = 0.0
    // beta_0 = 0.0;
    // r_0    = b - Ax0; >> r_0 = b; (as of my choice)
    
    // Condicoes iniciais (iteracao zero)
    k    = 0;
    beta = 0.0f;
    
    x  = CPU_zaloc1F(n);
    r  = CPU_zaloc1F(n);
    d  = CPU_zaloc1F(n);
    Ad = CPU_zaloc1F(n);
    
    // r da iteracao 0
    for(i=0; i<n; i++)
        r[i] = b[i];
    
    // rr da iteracao 0
    rr = 0.0f;
    for(i=0; i<n; i++)
        rr += r[i]*r[i];
    
    
    // Save to file inverted extended reflectivity
    initFiles3d("./debug_r", &fImg_Hdr, &fImg_Bin, nlx, dx, -lx0*dx, nz, dz, oz, nx, dx, ox);
    writeSamples(r, fImg_Bin, nlx*nz*nx);
    fclose(fImg_Bin);
    fclose(fImg_Hdr);

    printf("\n Initial Values: ");
    printf("\n rr=%g   beta=%g", rr, beta);
    
    // Iteracoes
    for(k=0; k<niter; k++)
    {
        printf("\n Iteration k=%d ", k);
        printf("\n Checking 01: rr=%g \n", rr);
        // Update d
        for(i=0; i<n; i++)
            d[i] = r[i] + beta * d[i];
        
        // Save to file inverted extended reflectivity
        initFiles3d("./debug_d", &fImg_Hdr, &fImg_Bin, nlx, dx, -lx0*dx, nz, dz, oz, nx, dx, ox);
        writeSamples(d, fImg_Bin, nlx*nz*nx);
        fclose(fImg_Bin);
        
        // Get Ad (action of Hessian on the image)
        //getAd(Ad, d);
        getAD(Ad, d, dataset, vp_bg, jt, nt, nx, nz, dx, dz, dt, ox, oz, ot, lx0, k);

        // Save to file inverted extended reflectivity
        initFiles3d("./debug_Ad", &fImg_Hdr, &fImg_Bin, nlx, dx, -lx0*dx, nz, dz, oz, nx, dx, ox);
        writeSamples(Ad, fImg_Bin, nlx*nz*nx);
        fclose(fImg_Bin);
        
        // Update alpha
        alpha_denominator = 0.0f;
        for(i=0; i<n; i++)
            alpha_denominator += d[i] * Ad[i];
        
        printf("\n Checking 02: aAd=%g (alpha_denominator)\n", alpha_denominator);

        alpha = rr/alpha_denominator;
        
        printf(" Checking 03: rr=%g    alpha_denominator=%g    alpha=%g", rr, alpha, alpha_denominator);
        
        // Update x
        for(i=0; i<n; i++)
            x[i] = x[i] + alpha * d[i];
        
        // Update r
        for(i=0; i<n; i++)
            r[i] = r[i] - alpha * Ad[i];
        
        // Keep denominator of beta expression
        beta_denomin = rr;

        printf("\n Checking 04: beta_denomin=%g \n", beta_denomin);
                
        // Update rr
        rr = 0.0f;
        for(i=0; i<n; i++)
            rr += r[i]*r[i];
        
        beta = rr/beta_denomin;            
        
        printf(" Checking 05: beta_denomin=%g    beta=%g", beta_denomin, beta);
    }
    

return x;
}

void getAD(float *Ad, float *d, Dataset *dataset, float *vp_bg, int jt, int nt, int nx, int nz, float dx, float dz, \
           float dt, float ox, float oz, float ot, int lx0, int iter) {

    int nshot = ((dataset->shot)[0]).nshot;
    int ishot;
    int ntheta = 0;
    float otheta=0, dtheta=0;
    // Forward operator
    printf("\n\n Computing Ad for CG iter %d: extended modeling %d shots: \n", iter, nshot);
    for(ishot=0; ishot<nshot; ishot++) {
        printf(" Shot %d \n", ishot);
        modelShotBorn_lx(d, vp_bg, &((dataset->shot)[ishot]), jt, nt, nx, nz, dx, dz, dt, ox, oz, ot, lx0);
    }
    // Adjoint operator
    printf("\n\n Computing Ad for CG iter %d: extended migrating %d shots: \n", iter, nshot);
    for(ishot=0; ishot<nshot; ishot++) {
        printf(" Shot %d \n", ishot);
        emigrateShot_lx(NULL, NULL, NULL, Ad, NULL, vp_bg, &((dataset->shot)[ishot]), \
                        jt, nt, nx, nz, dx, dz, dt, ox, oz, ot, lx0, ntheta, dtheta, otheta);
    }
return;
}

void normalizeB(float *b, int n) {
    int i;
    float max_b = -99999.0;
    for(i=0; i<n; i++) {
        if( max_b<fabsf(b[i]) )    max_b = fabsf(b[i]);
    }
    for(i=0; i<n; i++) {
        b[i] /= max_b;
    }
return;
}