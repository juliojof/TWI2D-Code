int   writePSFs;
//int   nx;
//int   nz;
//float dx;
//float dz;
//float ox;
//float oz;
//float epsilon;

void  normalizeB(float *b, int n);

void  getAd(float *Ad, float *d, int nx, int nz, float dx, float dz, float *image_psf, int **positionPSF);

float* getInterpHessianAtPoint(int i, int nz, float dx, float dz, float *imagePSF, int **positionPSF);

float* interpK(float *A_00, float *A_01, float *A_10, float *A_11, \
               float p_00, float p_01, float p_10, float p_11, \
               int n1, int n2, float dx, float dz);

fftw_complex  *fftw2D_forw_Real2Comp(float *A, int n1, int n2);

float* fftw2D_back_Comp2Real(fftw_complex  *a, int n1, int n2);

float* cg_psf(float *image, float *image_psf, int **positionPSF, int nx_, int nz_, float dx_, \
                 float dz_, float dt_, float ox_, float oz_, int niter, float epsilon_)
{
    float  *A, *x, *b, *r, *d, *Ad;
    float   beta, beta_denomin, alpha, alpha_denominator, rr;
    int     k, n;
    int     i, ix, iz;
    int     nx = nx_;   
    int     nz = nz_;
    float   dx = dx_; 
    float   dz = dz_;
    float   ox = ox_; 
    float   oz = oz_;
    float   epsilon = epsilon_;
    
    printf("\n Starting Conjugate Gradient \n ");
    n = nx*nz;
    b = image;

    //normalizeB(b,n);

    printf("\n\n Value of n: %d \n\n", n);
    
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
    d  = CPU_zaloc1F(n);
    Ad = CPU_zaloc1F(n);
    
    // r iteration 0
    for(i=0; i<n; i++)
        r[i] = b[i];
    
    // rr iteration 0
    rr = 0.0f;
    for(i=0; i<n; i++)
        rr += r[i]*r[i];
    
    printf("\n Initial Values: ");
    printf("\n rr=%g   beta=%g", rr, beta);
    
    // Iterate
    float rr_prev;
    for(k=0; k<niter; k++)
    {
        rr_prev = rr;;
        printf("\n Iteration k=%d ", k);
        printf("\n Checking 01: rr=%g \n", rr);
        // Update d
        for(i=0; i<n; i++)
            d[i] = r[i] + beta * d[i];
        
        // Get Ad (action of Hessian on the image)
        getAd(Ad, d, nx, nz, dx, dz, image_psf, positionPSF);
        
        // Update alpha
        alpha_denominator = 0.0f;
        for(i=0; i<n; i++)
            alpha_denominator += d[i] * Ad[i];
        
        printf("\n Checking 02: aAd=%g (alpha_denominator)\n", alpha_denominator);

        alpha = rr/alpha_denominator;
        
        printf(" Checking 03: rr=%g    alpha_denominator=%g    alpha=%g", rr, alpha, alpha_denominator);
        
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
        
        if(rr>=rr_prev) {
            if(k==0)
                for(i=0; i<n; i++)    x[i] = b[i];

            break;
        }

        // Update x
        for(i=0; i<n; i++)
            x[i] = x[i] + alpha * d[i];
        
        beta = rr/beta_denomin;            
        
        printf(" Checking 05: beta_denomin=%g    beta=%g", beta_denomin, beta);
        
    }
    
    free(r);
    free(d);
    free(Ad);

return x;
}


void  getAd(float *Ad, float *d, int nx, int nz, float dx, float dz, float *image_psf, int **positionPSF)
{
    int     ix, iz, i;
    int     ix_, iz_;
    int     j, j0;
    float  *A;
    
    //A = CPU_zaloc1F( (jxPSF+1)*(jzPSF+1) );
    
    for(ix=jxPSF/2; ix<nx-jxPSF/2; ix++)
        for(iz=jzPSF/2; iz<nz-jzPSF/2; iz++)
        {
            i = iz + ix*nz;

            if(ix==nx/2  &&  iz==nz/2)    writePSFs = 1;
            else                          writePSFs = 0;
            
            Ad[i] = 0.0f;
            A = getInterpHessianAtPoint(i, nz, dx, dz, image_psf, positionPSF);
            
            for(ix_=0; ix_<jxPSF+1; ix_++)
                for(iz_=0; iz_<jzPSF+1; iz_++)
                {
                    j  = (iz+iz_-jzPSF/2) +  (ix+ix_-jxPSF/2) * nz;
                    
                    j0 = iz_ + ix_ * (jzPSF+1);
                    
                    Ad[i] += A[j0] * d[j];
                    //Ad[i] += 1.0f * d[j];
                }
            free(A);
        }

return;
}
float* getInterpHessianAtPoint(int i, int nz, float dx, float dz, float *imagePSF, int **positionPSF)
{
    int     ix0, iz0;
    int     ix1, iz1;
    int     j, j0;
    int     ix_, iz_;
    int     ix, iz;
    int     ixPSF, izPSF;
    float  *A_00, *A_01, *A_10, *A_11;
    float   d_00, d_01, d_10, d_11;
    float   p_00, p_10, p_01, p_11, p_total;
    float   A_a, A_b;
    float   filter, sigma2, d_half, d2;
    float   px_0, px_1, pz_0, pz_1;
    
    A_00 = CPU_zaloc1F( (jxPSF+1)*(jzPSF+1) );
    A_01 = CPU_zaloc1F( (jxPSF+1)*(jzPSF+1) );
    A_10 = CPU_zaloc1F( (jxPSF+1)*(jzPSF+1) );
    A_11 = CPU_zaloc1F( (jxPSF+1)*(jzPSF+1) );
    
    d_half  = 0.5*jxPSF * dx;
    sigma2  = log(2.0)/(d_half*d_half);
    
    ix  = i/nz;
    iz  = i - ix*nz;        
    ix0 = positionPSF[i][0];
    iz0 = positionPSF[i][1];
    ix1 = positionPSF[i][2];
    iz1 = positionPSF[i][3];
    
    // PSF 1    
    ixPSF = ix0;
    izPSF = iz0;
    d_00  = sqrtf( dx*(ix-ixPSF)*dx*(ix-ixPSF) + dz*(iz-izPSF)*dz*(iz-izPSF) );
    for(ix_=0; ix_<jxPSF+1; ix_++)
        for(iz_=0; iz_<jzPSF+1; iz_++)
        {
            j  = (izPSF+iz_-jzPSF/2) +  (ixPSF+ix_-jxPSF/2) * nz;
            
            j0 = iz_ + ix_* (jzPSF+1);
            
            d2 = dx*(ix_-jxPSF/2)*dx*(ix_-jxPSF/2) + dz*(iz_-jzPSF/2)*dz*(iz_-jzPSF/2);
            
            filter = 1.0;//exp(-sigma2*d2);
            
            A_00[j0] = -imagePSF[j] * filter;
        }
    // PSF 2
    ixPSF = ix0;
    izPSF = iz1;
    d_01  = sqrtf( dx*(ix-ixPSF)*dx*(ix-ixPSF) + dz*(iz-izPSF)*dz*(iz-izPSF) );
    for(ix_=0; ix_<jxPSF+1; ix_++)
        for(iz_=0; iz_<jzPSF+1; iz_++)
        {
            j  = (izPSF+iz_-jzPSF/2) +  (ixPSF+ix_-jxPSF/2) * nz;
            
            j0 = iz_ + ix_* (jzPSF+1);
            
            d2 = dx*(ix_-jxPSF/2)*dx*(ix_-jxPSF/2) + dz*(iz_-jzPSF/2)*dz*(iz_-jzPSF/2);
            
            filter = 1.0;//exp(-sigma2*d2);
            
            A_01[j0] = -imagePSF[j] * filter;
        }
    // PSF 3
    ixPSF = ix1;
    izPSF = iz0;
    d_10  = sqrtf( dx*(ix-ixPSF)*dx*(ix-ixPSF) + dz*(iz-izPSF)*dz*(iz-izPSF) );
    for(ix_=0; ix_<jxPSF+1; ix_++)
        for(iz_=0; iz_<jzPSF+1; iz_++)
        {
            j  = (izPSF+iz_-jzPSF/2) +  (ixPSF+ix_-jxPSF/2) * nz;
            
            j0 = iz_ + ix_* (jzPSF+1);
            
            d2 = dx*(ix_-jxPSF/2)*dx*(ix_-jxPSF/2) + dz*(iz_-jzPSF/2)*dz*(iz_-jzPSF/2);
            
            filter = 1.0;//exp(-sigma2*d2);
            
            A_10[j0] = -imagePSF[j] * filter;
        }
    // PSF 4
    ixPSF = ix1;
    izPSF = iz1;
    d_11  = sqrtf( dx*(ix-ixPSF)*dx*(ix-ixPSF) + dz*(iz-izPSF)*dz*(iz-izPSF) );
    for(ix_=0; ix_<jxPSF+1; ix_++)
        for(iz_=0; iz_<jzPSF+1; iz_++)
        {
            j  = (izPSF+iz_-jzPSF/2) +  (ixPSF+ix_-jxPSF/2) * nz;
            
            j0 = iz_ + ix_* (jzPSF+1);
            
            d2 = dx*(ix_-jxPSF/2)*dx*(ix_-jxPSF/2) + dz*(iz_-jzPSF/2)*dz*(iz_-jzPSF/2);
            
            filter = 1.0;//exp(-sigma2*d2);
            
            A_11[j0] = -imagePSF[j] * filter;
        }
    
    // Interpolation
    if(ix<ix0) {   
        px_0 = 1.0;
        px_1 = 0.0;
    }
    else if(ix>ix1) {
        px_0 = 0.0;         
        px_1 = 1.0;
    }
    else {              
        px_0 = 1.0f*(ix1-ix); 
        px_1 = 1.0f*(ix-ix0); 
    }
    if(iz<iz0) {   
        pz_0 = 1.0;         
        pz_1 = 0.0; 
    }
    else if(iz>iz1) {   
        pz_0 = 0.0;         
        pz_1 = 1.0;
    }
    else {              
        pz_0 = 1.0f*(iz1-iz); 
        pz_1 = 1.0f*(iz-iz0); 
    }
    p_00 = px_0*pz_0;
    p_01 = px_0*pz_1;
    p_10 = px_1*pz_0;
    p_11 = px_1*pz_1;
    /*
    p_00    = 1.0f/(d_00 + 0.003*dx);
    p_10    = 1.0f/(d_10 + 0.003*dx);
    p_01    = 1.0f/(d_01 + 0.003*dx);
    p_11    = 1.0f/(d_11 + 0.003*dx);
    */
    p_total = p_00 + p_10 + p_01 + p_11;
    p_00   = p_00/p_total;
    p_10   = p_10/p_total;
    p_01   = p_01/p_total;
    p_11   = p_11/p_total;
    float *A = interpK(A_00, A_01, A_10, A_11, p_00, p_01, p_10, p_11, (jzPSF+1), (jxPSF+1), dx, dz);

    if(writePSFs) {
    FILE *fpsf_Hdr, *fpsf_Bin;
    initFiles2d("./psf_00", &fpsf_Hdr, &fpsf_Bin, (jzPSF+1), dz, 0, (jxPSF+1), dx, 0);
    writeSamples(A_00, fpsf_Bin, (jzPSF+1)*(jxPSF+1));
    fclose(fpsf_Hdr);  fclose(fpsf_Bin);

    initFiles2d("./psf_01", &fpsf_Hdr, &fpsf_Bin, (jzPSF+1), dz, 0, (jxPSF+1), dx, 0);
    writeSamples(A_01, fpsf_Bin, (jzPSF+1)*(jxPSF+1));
    fclose(fpsf_Hdr);  fclose(fpsf_Bin);

    initFiles2d("./psf_10", &fpsf_Hdr, &fpsf_Bin, (jzPSF+1), dz, 0, (jxPSF+1), dx, 0);
    writeSamples(A_10, fpsf_Bin, (jzPSF+1)*(jxPSF+1));
    fclose(fpsf_Hdr);  fclose(fpsf_Bin);

    initFiles2d("./psf_11", &fpsf_Hdr, &fpsf_Bin, (jzPSF+1), dz, 0, (jxPSF+1), dx, 0);
    writeSamples(A_11, fpsf_Bin, (jzPSF+1)*(jxPSF+1));
    fclose(fpsf_Hdr);  fclose(fpsf_Bin);

    initFiles2d("./psf_int", &fpsf_Hdr, &fpsf_Bin, (jzPSF+1), dz, 0, (jxPSF+1), dx, 0);
    writeSamples(A, fpsf_Bin, (jzPSF+1)*(jxPSF+1));
    fclose(fpsf_Hdr);  fclose(fpsf_Bin);
    }
    
/*    
    for(ix_=0; ix_<jxPSF+1; ix_++)
        for(iz_=0; iz_<jzPSF+1; iz_++)
        {
            j0 = iz_ + ix_* (jzPSF+1);

            p_00    = 1.0f/(d_00 + 0.003*dx);
            p_10    = 1.0f/(d_10 + 0.003*dx);
            p_01    = 1.0f/(d_01 + 0.003*dx);
            p_11    = 1.0f/(d_11 + 0.003*dx);
            p_total = p_00 + p_10 + p_01 + p_11;
            p_00   /= p_total;
            p_10   /= p_total;
            p_01   /= p_total;
            p_11   /= p_total;
            
            //A[j0] = A_00[j0]*p_00  +  A_10[j0]*p_10  +  A_01[j0]*p_01  +  A_11[j0]*p_11;
            
            A_a = ( abs(ix1-ix) * A_00[j0] + abs(ix-ix0) * A_10[j0] ) / (abs(ix-ix0) + abs(ix1-ix));
            A_b = ( abs(ix1-ix) * A_01[j0] + abs(ix-ix0) * A_11[j0] ) / (abs(ix-ix0) + abs(ix1-ix));
            
            A[j0] = (abs(iz1-iz) * A_a + abs(iz-iz0) * A_b) / ( abs(iz-iz0) + abs(iz1-iz) );
                        
        }
//*/   
    
    free(A_00);
    free(A_01);
    free(A_10);
    free(A_11);
    
return A;
}


float* interpK(float *A_00, float *A_01, float *A_10, float *A_11, \
               float p_00, float p_01, float p_10, float p_11, \
               int n1, int n2, float dx, float dz) {

    fftw_complex *a_00 = fftw2D_forw_Real2Comp(A_00, n1, n2);
    fftw_complex *a_01 = fftw2D_forw_Real2Comp(A_01, n1, n2);
    fftw_complex *a_10 = fftw2D_forw_Real2Comp(A_10, n1, n2);
    fftw_complex *a_11 = fftw2D_forw_Real2Comp(A_11, n1, n2);
    
    fftw_complex *a  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n1 * n2);

    int i;
    double amp_tmp, amp;
    double maxAmp = -9999.0;
    for(i=0; i<n1*n2; i++) {
        a[i][0] = p_00*a_00[i][0] + p_01*a_01[i][0] + p_10*a_10[i][0] + p_11*a_11[i][0];
        a[i][1] = p_00*a_00[i][1] + p_01*a_01[i][1] + p_10*a_10[i][1] + p_11*a_11[i][1];
        amp_tmp = (a[i][0]*a[i][0] + a[i][1]*a[i][1]);
        if(maxAmp<amp_tmp) maxAmp = amp_tmp;
    }
    amp = 0.0;
    for(i=0; i<n1*n2; i++) {
        amp_tmp = (a[i][0]*a[i][0] + a[i][1]*a[i][1]);
        if(amp_tmp>0.05*maxAmp) amp += amp_tmp;
    }
    amp = sqrt(amp)/(n1*n2);
    //printf("\n\n rms amp = %g \n\n", amp);
    for(i=0; i<n1*n2; i++) {
        //a[i][0] += epsilon * amp;
        //a[i][1] += epsilon * amp;
        a[i][0] += 0.01e-8;
        a[i][1] += 0.01e-8;
    }
    float  *A = fftw2D_back_Comp2Real(a, n1, n2);
    
    if(writePSFs) {
    float *temp = CPU_zaloc1F(n1*n2);
    FILE *fpsf_Hdr, *fpsf_Bin;
    initFiles2d("./psf_ft_00", &fpsf_Hdr, &fpsf_Bin, (jzPSF+1), dz, 0, (jxPSF+1), dx, 0);
    for(i=0; i<n1*n2; i++) temp[i] = a_00[i][0];
    writeSamples(temp, fpsf_Bin, (jzPSF+1)*(jxPSF+1));
    fclose(fpsf_Hdr);  fclose(fpsf_Bin);

    initFiles2d("./psf_ft_01", &fpsf_Hdr, &fpsf_Bin, (jzPSF+1), dz, 0, (jxPSF+1), dx, 0);
    for(i=0; i<n1*n2; i++) temp[i] = a_01[i][0];
    writeSamples(temp, fpsf_Bin, (jzPSF+1)*(jxPSF+1));
    fclose(fpsf_Hdr);  fclose(fpsf_Bin);

    initFiles2d("./psf_ft_10", &fpsf_Hdr, &fpsf_Bin, (jzPSF+1), dz, 0, (jxPSF+1), dx, 0);
    for(i=0; i<n1*n2; i++) temp[i] = a_10[i][0];
    writeSamples(temp, fpsf_Bin, (jzPSF+1)*(jxPSF+1));
    fclose(fpsf_Hdr);  fclose(fpsf_Bin);

    initFiles2d("./psf_ft_11", &fpsf_Hdr, &fpsf_Bin, (jzPSF+1), dz, 0, (jxPSF+1), dx, 0);
    for(i=0; i<n1*n2; i++) temp[i] = a_11[i][0];
    writeSamples(temp, fpsf_Bin, (jzPSF+1)*(jxPSF+1));
    fclose(fpsf_Hdr);  fclose(fpsf_Bin);

    initFiles2d("./psf_ft", &fpsf_Hdr, &fpsf_Bin, (jzPSF+1), dz, 0, (jxPSF+1), dx, 0);
    for(i=0; i<n1*n2; i++) temp[i] = a[i][0];
    writeSamples(temp, fpsf_Bin, (jzPSF+1)*(jxPSF+1));
    fclose(fpsf_Hdr);  fclose(fpsf_Bin);

    initFiles2d("./psf_di_int", &fpsf_Hdr, &fpsf_Bin, (jzPSF+1), dz, 0, (jxPSF+1), dx, 0);
    for(i=0; i<n1*n2; i++) temp[i] = A[i];
    writeSamples(temp, fpsf_Bin, (jzPSF+1)*(jxPSF+1));
    fclose(fpsf_Hdr);  fclose(fpsf_Bin);
    free(temp);
    }

    fftw_free(a); 

return A;
}

fftw_complex* fftw2D_forw_Real2Comp(float *A, int n1, int n2) {

    fftw_complex  *orig, *tran;
    fftw_plan      plan;

    orig = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n1 * n2);
    tran = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n1 * n2);
    int i;
    for(i=0; i<n1*n2; i++)
    {
        orig[i][0] = A[i];
        orig[i][1] = 0.0f;
    }
    plan  = fftw_plan_dft_2d(n1, n2, orig, tran,  FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    fftw_free(orig);

return tran;
}
float* fftw2D_back_Comp2Real(fftw_complex  *a, int n1, int n2) {

    float *A = CPU_zaloc1F(n1*n2);
    fftw_complex  *orig;
    fftw_plan      plan;

    orig = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n1 * n2);
    plan  = fftw_plan_dft_2d(n1, n2, a, orig,  FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    int i;
    for(i=0; i<n1*n2; i++)    {
        A[i] = orig[i][0];
        //if(writePSFs)    printf("\n A[i] = %f  orig[i][0]=%f", A[i], orig[i][0]);
    }

    fftw_destroy_plan(plan);
    fftw_free(orig);

return A;
}






/*
void  getAd_new(float *Ad, float *d, float *image_psf, int **positionPSF)
{
    int     ix, iz, i;
    int     ix_, iz_;
    int     j, j0;
    float  *A;
    int     wx_psf = 2*floorf(1.2*(0.5*jxPSF))+1;
    int     wz_psf = 2*floorf(1.2*(0.5*jzPSF))+1;
    
    A = CPU_zaloc1F( wx_psf*wz_psf );
    
    for(ix=jxPSF/2; ix<nx-jxPSF/2; ix++)
        for(iz=jzPSF/2; iz<nz-jzPSF/2; iz++)
        {
            i = iz + ix*nz;
            
            Ad[i] = 0.0f;
            getInterpHessianAtPoint(i, A, image_psf, positionPSF);
            
            for(ix_=0; ix_<wx_psf; ix_++)
                for(iz_=0; iz_<wz_psf; iz_++)
                {
                    j  = ( iz + (iz_-wz_psf/2) ) +  ( ix + (ix_-wx_psf/2) ) * nz;
                    
                    j0 = iz_ + ix_ * (wz_psf+1);
                    
                    Ad[i] += A[j0] * d[j];
                }
        }
    
    free(A);

return;
}
void getInterpHessianAtPoint_new(int i, float *A, float *imagePSF, int **positionPSF)
{
    int     ix0, iz0;
    int     ix1, iz1;
    int     j, j0;
    int     ix_, iz_;
    int     ix, iz;
    int     ixPSF, izPSF;
    float  *A_00, *A_01, *A_10, *A_11;
    float   d_00, d_01, d_10, d_11;
    float   p_00, p_10, p_01, p_11, p_total;
    float   A_a, A_b;
    float   filter, sigma2, d_half, d2;
    int     wx_psf = 2*floorf(1.2*(0.5*jxPSF))+1;
    int     wz_psf = 2*floorf(1.2*(0.5*jzPSF))+1;
    
    A_00 = CPU_zaloc1F( wx_psf*wz_psf );
    A_01 = CPU_zaloc1F( wx_psf*wz_psf );
    A_10 = CPU_zaloc1F( wx_psf*wz_psf );
    A_11 = CPU_zaloc1F( wx_psf*wz_psf );
    
    d_half  = 0.5*wx_psf * dx;
    sigma2  = log(2.0)/(d_half*d_half);
    
    ix  = i/nz;
    iz  = i - ix*nz;        
    ix0 = positionPSF[i][0];
    iz0 = positionPSF[i][1];
    ix1 = positionPSF[i][2];
    iz1 = positionPSF[i][3];
    
    // PSF 1    
    ixPSF = ix0;
    izPSF = iz0;
    d_00  = sqrtf( dx*(ix-ixPSF)*dx*(ix-ixPSF) + dz*(iz-izPSF)*dz*(iz-izPSF) );
    for(ix_=0; ix_<wx_psf; ix_++)
        for(iz_=0; iz_<wz_psf; iz_++)
        {
            j  = (izPSF+iz_-wz_psf/2) +  (ixPSF+ix_-wx_psf/2) * nz;
            
            j0 = iz_ + ix_* wz_psf;
            
            d2 = dx*(ix_-wx_psf/2)*dx*(ix_-wx_psf/2) + dz*(iz_-wz_psf/2)*dz*(iz_-wz_psf/2);
            
            filter = exp(-sigma2*d2);
            
            A_00[j0] = -imagePSF[j] * filter;
        }
    // PSF 2
    ixPSF = ix0;
    izPSF = iz1;
    d_01  = sqrtf( dx*(ix-ixPSF)*dx*(ix-ixPSF) + dz*(iz-izPSF)*dz*(iz-izPSF) );
    for(ix_=0; ix_<wx_psf; ix_++)
        for(iz_=0; iz_<wz_psf; iz_++)
        {
            j  = (izPSF+iz_-wz_psf/2) +  (ixPSF+ix_-wx_psf/2) * nz;
            
            j0 = iz_ + ix_* wz_psf;
            
            d2 = dx*(ix_-wx_psf/2)*dx*(ix_-wx_psf/2) + dz*(iz_-wz_psf/2)*dz*(iz_-wz_psf/2);
            
            filter = exp(-sigma2*d2);
            
            A_01[j0] = -imagePSF[j] * filter;
        }
    // PSF 3
    ixPSF = ix1;
    izPSF = iz0;
    d_10  = sqrtf( dx*(ix-ixPSF)*dx*(ix-ixPSF) + dz*(iz-izPSF)*dz*(iz-izPSF) );
    for(ix_=0; ix_<wx_psf; ix_++)
        for(iz_=0; iz_<wz_psf; iz_++)
        {
            j  = (izPSF+iz_-wz_psf/2) +  (ixPSF+ix_-wx_psf/2) * nz;
            
            j0 = iz_ + ix_* wz_psf;
            
            d2 = dx*(ix_-wx_psf/2)*dx*(ix_-wx_psf/2) + dz*(iz_-wz_psf/2)*dz*(iz_-wz_psf/2);
            
            filter = exp(-sigma2*d2);
            
            A_10[j0] = -imagePSF[j] * filter;
        }
    // PSF 4
    ixPSF = ix1;
    izPSF = iz1;
    d_11  = sqrtf( dx*(ix-ixPSF)*dx*(ix-ixPSF) + dz*(iz-izPSF)*dz*(iz-izPSF) );
    for(ix_=0; ix_<wx_psf; ix_++)
        for(iz_=0; iz_<wz_psf; iz_++)
        {
            j  = (izPSF+iz_-wz_psf/2) +  (ixPSF+ix_-wx_psf/2) * nz;
            
            j0 = iz_ + ix_* wz_psf;
            
            d2 = dx*(ix_-wx_psf/2)*dx*(ix_-wx_psf/2) + dz*(iz_-wz_psf/2)*dz*(iz_-wz_psf/2);
            
            filter = exp(-sigma2*d2);
            
            A_11[j0] = -imagePSF[j] * filter;
        }

    for(ix_=0; ix_<wx_psf; ix_++)
        for(iz_=0; iz_<wz_psf; iz_++)
        {
            j0 = iz_ + ix_* (wz_psf+1);

            p_00    = 1.0f/(d_00 + 0.003*dx);
            p_10    = 1.0f/(d_10 + 0.003*dx);
            p_01    = 1.0f/(d_01 + 0.003*dx);
            p_11    = 1.0f/(d_11 + 0.003*dx);
            p_total = p_00 + p_10 + p_01 + p_11;
            p_00   /= p_total;
            p_10   /= p_total;
            p_01   /= p_total;
            p_11   /= p_total;
            
            //A[j0] = A_00[j0]*p_00  +  A_10[j0]*p_10  +  A_01[j0]*p_01  +  A_11[j0]*p_11;
            
            A_a = ( abs(ix1-ix) * A_00[j0] + abs(ix-ix0) * A_10[j0] ) / (abs(ix-ix0) + abs(ix1-ix));
            A_b = ( abs(ix1-ix) * A_01[j0] + abs(ix-ix0) * A_11[j0] ) / (abs(ix-ix0) + abs(ix1-ix));
            
            A[j0] = (abs(iz1-iz) * A_a + abs(iz-iz0) * A_b) / ( abs(iz-iz0) + abs(iz1-iz) );
        }

    free(A_00);
    free(A_01);
    free(A_10);
    free(A_11);

return;
}
//*/