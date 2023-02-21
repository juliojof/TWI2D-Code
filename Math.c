float* transform2Envelope(float *traces, int nt, float dt, int ntrace) {
    int itr, it;
    float *envelopeGather = CPU_zaloc1F(ntrace*nt);
    for(itr=0; itr<ntrace; itr++) {
        
        for(it=0; it<nt; it++)    
            envelopeGather[it+itr*nt] = traces[it+itr*nt];        
            
        envelope(envelopeGather+itr*nt, nt, dt);

    }
return envelopeGather;
}

void envelope(float *trace, int nt, float dt)
{
    int it;
    float *traceHilbert = CPU_zaloc1F(nt);
    float *instPhase    = CPU_zaloc1F(nt);

    for(it=0; it<nt; it++) {
        traceHilbert[it] = trace[it];
    }
    hilbert(traceHilbert, nt, dt);

    for(it=0; it<nt; it++) {
        trace[it] = trace[it]*trace[it] + traceHilbert[it]*traceHilbert[it];
        trace[it] = sqrtf(trace[it]);
    }

    free(traceHilbert);
    free(instPhase);

return;
}

float* transform2InstantaneousPhase(float *traces, int nt, float dt, int ntrace) {
    int itr, it;
    float *instPhaseGather = CPU_zaloc1F(ntrace*nt);
    for(itr=0; itr<ntrace; itr++) {
        
        for(it=0; it<nt; it++)    
            instPhaseGather[it+itr*nt] = traces[it+itr*nt];        
            
        instantaneousPhase(instPhaseGather+itr*nt, nt, dt);

    }
return instPhaseGather;
}

void instantaneousPhase(float *trace, int nt, float dt)
{
    int it;
    float *traceHilbert = CPU_zaloc1F(nt);
    float *instPhase    = CPU_zaloc1F(nt);

    for(it=0; it<nt; it++) {
        traceHilbert[it] = trace[it];
    }
    hilbert(traceHilbert, nt, dt);
/*
    for(it=0; it<nt; it++) {
        if(trace[it]!=0.0 || traceHilbert[it]!=0.0)
            trace[it] = atan2f(trace[it], traceHilbert[it]) - 0.0*M_PI;
        else
            trace[it] = 0.0f;
    }
*/
    for(it=0; it<nt; it++)
        trace[it] = atan2f(trace[it], traceHilbert[it]) - 0.0*M_PI;
    for(it=0; it<nt; it++)
        trace[it+1] += trace[it];  
    
    free(traceHilbert);
    free(instPhase);

return;
}

void hilbert(float *trace, int nt, float dt)
{

    fftw_complex  *orig, *tran;
    fftw_plan      planF, planB;
    int            nf, ifreq, it;
    float          dfreq, omega;
    float          a, b, alfa, beta, mod, freq;
    float          ap, bp, am, bm;
    
    nf   = nt;
    dfreq   = 1.0/(nt*dt);

    orig  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt);
    tran  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt);

    planF = fftw_plan_dft_1d(nt, orig, tran,  FFTW_FORWARD, FFTW_ESTIMATE);
    planB = fftw_plan_dft_1d(nt, tran, orig, FFTW_BACKWARD, FFTW_ESTIMATE);

    ap = cos(+M_PI/2.0f);
    bp = sin(+M_PI/2.0f);
    am = cos(-M_PI/2.0f);
    bm = sin(-M_PI/2.0f);

    for(it=0; it<nt; it++)
    {
        orig[it][0] = trace[it]/sqrtf((1.0f*nt));
        orig[it][1] = 0.0;
    }

    fftw_execute(planF);        
        
    tran[0][0] = 0.0f;
    tran[0][1] = 0.0f;
    if( (nf%2)==0 )
    {
        tran[nf/2][0] = 0.0f;
        tran[nf/2][1] = 0.0f;
    }   
    
    for(ifreq=1; ifreq<(nf-1)/2; ifreq++)
    {
        freq  = ifreq*dfreq;
        omega = 2.0f*M_PI*freq;

        alfa  = tran[ifreq][0];
        beta  = tran[ifreq][1];

        tran[ifreq   ][0] = (am*alfa - bm*beta);
        tran[ifreq   ][1] = (am*beta + bm*alfa);

        alfa  = tran[nf-ifreq][0];
        beta  = tran[nf-ifreq][1];
       
        //tran[nf-ifreq][0] = (ap*alfa - bp*beta)
        //tran[nf-ifreq][1] = (ap*beta + bp*alfa);   
    }
        
    fftw_execute(planB);
        
    for(it=0; it<nt; it++)
    {
        trace[it] = orig[it][0]/sqrtf((1.0f*nt));
    }
    
    fftw_destroy_plan(planF);
    fftw_destroy_plan(planB);

    fftw_free(orig);
    fftw_free(tran);

return;
}


// Sinc function
float sinc(float x) {
    if(x==0.0f)
        return 1.0f;
    else
        return (sinf(M_PI*x)/(M_PI*x));
}
// Kaiser window
float kaiser(float x, float r, float b) {
    float window;
    if(fabsf(x)<=r)
        window = mBFZ(b*sqrtf(1.0 - (x/r)*(x/r))) / mBFZ(b);
    else 
        window = 0.0;
return window;
}
// Modified Bessel First Kind Zero Order
float mBFZ(float x) {
    double add;
    double out  = 0.0;
    double samp = 1.0;
    double arg = 0.25f * x*x;
    long long int kfat  = 1;
    long long int kfat1 = 1;
    long long int k;
    k = 0;
    do {
        // out += ((samp/(kfat))/kfat1);
        add = samp/(kfat*kfat1);
        out += add;
        samp *= arg;
        k++;
        kfat  *= k;
        kfat1 *= (k+1);
    } while(add>0.0001*arg);

return ((float) out);
}


void applyHalfDerivative(float *inp, int nt, float dt, float sign)
{

    fftw_complex  *trcT, *trcF;
    fftw_plan      planForw, planBack;
    int            nfreq, ifreq, it;
    float          dfreq, omega;
    float          a, b, RE, IM, mod, freq;
    float          ap, bp, am, bm;

    
    nfreq   = nt;
    dfreq   = 1.0/(nt*dt);

    trcT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt);
    trcF = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt);

    planForw = fftw_plan_dft_1d(nt, trcT, trcF, FFTW_FORWARD,  FFTW_ESTIMATE);
    planBack = fftw_plan_dft_1d(nt, trcF, trcT, FFTW_BACKWARD, FFTW_ESTIMATE);

    ap = cos(sign*M_PI/4.0f);
    bp = sin(sign*M_PI/4.0f);
    am = cos(sign*3.0f*M_PI/4.0f);
    bm = sin(sign*3.0f*M_PI/4.0f);

    for(it=0; it<nt; it++)
    {
        trcT[it][0] = inp[it]/sqrtf((1.0f*nt));
        trcT[it][1] = 0.0;
    }

    fftw_execute(planForw);        

    trcF[0][0] = 0.0f;
    trcF[0][1] = 0.0f;
    
    if( (nfreq%2)==0 )
    {
        trcF[nfreq/2][0] = 0.0f;
        trcF[nfreq/2][1] = 0.0f;
    }   
    
    for(ifreq=1; ifreq<(nfreq-1)/2; ifreq++)
    {
        freq  = ifreq*dfreq;
        omega = 2.0f*M_PI*freq;

        RE  = trcF[ifreq][0];
        IM  = trcF[ifreq][1];

        trcF[ifreq][0] = sqrtf(omega) * (ap*RE - bp*IM);
        trcF[ifreq][1] = sqrtf(omega) * (ap*IM + bp*RE);
       
        trcF[nfreq-ifreq][0] =  trcF[ifreq][0];
        trcF[nfreq-ifreq][1] = -trcF[ifreq][1];
        
    }
        
    fftw_execute(planBack);
        
    for(it=0; it<nt; it++)
    {
        inp[it] = trcT[it][0]/sqrtf((1.f*nt));
    }

    fftw_destroy_plan(planForw);
    fftw_destroy_plan(planBack);

    fftw_free(trcT);
    fftw_free(trcF);

return;
}


void bandpass(float *inp, int nt, float dt, float f1, float f2, float f3, float f4)
{

    fftw_complex  *trcT, *trcF;
    fftw_plan      planForw, planBack;
    int            nfreq, ifreq, it;
    float          dfreq, omega;
    float          filter, freq;

    nfreq   = nt;
    dfreq   = 1.0/(nt*dt);

    trcT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt);
    trcF = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt);

    planForw = fftw_plan_dft_1d(nt, trcT, trcF, FFTW_FORWARD,  FFTW_ESTIMATE);
    planBack = fftw_plan_dft_1d(nt, trcF, trcT, FFTW_BACKWARD, FFTW_ESTIMATE);


    for(it=0; it<nt; it++)
    {
        trcT[it][0] = inp[it]/sqrtf((1.0f*nt));
        trcT[it][1] = 0.0;
    }

    fftw_execute(planForw);        
   

    for(ifreq=0; ifreq<=nfreq/2; ifreq++)
    {
        freq  = ifreq*dfreq;
        
        if(freq<f1)
            filter = 0.0f;
        
        if(freq>=f1  &&  freq<f2)
            filter = (freq-f1)/(f2-f1);
        
        if(freq>=f2  &&  freq<=f3)
            filter = 1.0f;
        
        if(freq>f3  &&  freq<=f4)
            filter = (f4-freq)/(f4-f3);
        
        if(freq>f4)
            filter = 0.0f;
                
            trcF[ifreq][0] *= filter;
            trcF[ifreq][1] *= filter;
        
        if(ifreq>0 && ifreq<nfreq-ifreq)
        {
            trcF[nfreq-ifreq][0] *= filter;
            trcF[nfreq-ifreq][1] *= filter;
        }
        
    }
        
    fftw_execute(planBack);
        
    for(it=0; it<nt; it++)
    {
        inp[it] = trcT[it][0]/sqrtf((1.0f*nt));
    }
    
    fftw_destroy_plan(planForw);
    fftw_destroy_plan(planBack);

    fftw_free(trcT);
    fftw_free(trcF);

return;
}

float getVolumeRMS(float *array, int n) {
    int i;
    double rms = 0;
    for(i=0; i<n; i++)
        rms += (array[i]*array[i])/n;
    rms = sqrt(rms);
return ((float) rms);
}
void applyLinearWeight_TPCIG(float *tpcig_dso, float *tpcig, float dp, int nx, int ntau, int lp0) {
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int ix, itau, lp, np=2*lp0+1;
    for(ix=0;ix<nx;ix++)
        for(itau=0;itau<ntau;itau++) {
            for(lp=-lp0; lp<=lp0; lp++)
            {
                int itp = itau + (lp+lp0)*ntau + ix*np*ntau;
                tpcig_dso[itp] = fabsf(lp*dp) * fabsf(lp*dp) * tpcig[itp];
            }
        }
return;
}
void applyLinearWeight_ODCIG(float *eimage_dso, float *eimage, float dx, int nx, int nz, int lx0) {
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;
    int ix, iz, lx, nlx=2*lx0+1;
    for(ix=lx0;ix<nx-lx0;ix++)
        for(iz=0;iz<nz;iz++) {
            for(lx=-lx0; lx<=lx0; lx++)
            {
                int i  = id3_eic_lx(lx,iz,ix);
                eimage_dso[i] = fabsf(lx*dx) * eimage[i];
            }
        }
return;
}
void applyLinearDSO_ODCIG(float *eimage_dso, float *eimage, float dx, int nx, int nz, int lx0, float linear_coefficient) {
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;
    int ix, iz, lx, nlx=2*lx0+1;
    for(ix=lx0;ix<nx-lx0;ix++)
        for(iz=0;iz<nz;iz++) {
            for(lx=-lx0; lx<=lx0; lx++)
            {
                int i  = id3_eic_lx(lx,iz,ix);
                eimage_dso[i] = linear_coefficient * fabsf(lx*dx) * eimage[i];
            }
        }
return;
}
void applyAngleDerivative(float *out, float *eimage_theta, int nx, int nz, int ntheta, float dtheta) {

    int ix, iz, itheta, i;
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;
    int i_th0, i_th1, i_tha;
    int ltheta0 = ntheta/2;
    
    float *tmp = CPU_zaloc1F(nx*nz*ntheta);

    
    for(ix=0; ix<nx; ix++)
        for(iz=0; iz<nz; iz++) {
            i_th0  = id3_eic_th(ltheta0-1,iz,ix);
            i_tha  = id3_eic_th(ltheta0  ,iz,ix);
            i_th1  = id3_eic_th(ltheta0+1,iz,ix);
            tmp[i_tha] = (eimage_theta[i_th1] - eimage_theta[i_th0])/(2.0*dtheta);
            for(itheta=1; itheta<ltheta0-1; itheta++)
            {
                i_th0  = id3_eic_th(itheta-1+ltheta0,iz,ix);
                i_th1  = id3_eic_th(itheta+1+ltheta0,iz,ix);
                tmp[i_th0] = (eimage_theta[i_th1] - eimage_theta[i_th0])/(2.0f*dtheta);

                i_th0  = id3_eic_th(ltheta0-(itheta-1),iz,ix);
                i_th1  = id3_eic_th(ltheta0-(itheta+1),iz,ix);
                tmp[i_th0] = (eimage_theta[i_th1] - eimage_theta[i_th0])/(2.0f*dtheta);
            }
        }
    for(i=0; i<nx*nz*ntheta; i++) {
        if(out==NULL)    eimage_theta[i] = tmp[i];
        else             out[i]          = tmp[i];
    }
    free(tmp);
return;
}

void getFlatADCIG(float *out, float *eimage_theta, int nx, int nz, int ntheta, float dz, float dtheta, float oz, float theta_max) {

    int ix, iz, itheta, i;
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;
    int i_th0, i_th;
    int ltheta0 = ntheta/2;

    int ntheta_max = floorf(theta_max/dtheta);
    
    // Stack
    for(ix=0; ix<nx; ix++)
        for(iz=0; iz<nz; iz++) {
            i_th0  = id3_eic_th(0,iz,ix);
            for(itheta=0; itheta<ntheta_max; itheta++)
            {
                i_th = id3_eic_th(+itheta,iz,ix);
                out[i_th0] += eimage_theta[i_th]/(2*ntheta_max);
                i_th = id3_eic_th(-itheta,iz,ix);
                out[i_th0] += eimage_theta[i_th]/(2*ntheta_max);
            }
        }
    outputSmart2d("stack", out+((nx/2)*ntheta*nz), ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz);
    float *out_envelope = transform2Envelope(out, nz, dz, ntheta*nx);
    divideArrays(out, out, out_envelope, nx*nz*ntheta);
    outputSmart2d("stack_norm", out+((nx/2)*ntheta*nz), ntheta, dtheta, -ltheta0*dtheta, nz, dz, oz);
    free(out_envelope);
    // Spray
    for(ix=0; ix<nx; ix++)
        for(iz=0; iz<nz; iz++) {
            i_th0  = id3_eic_th(0,iz,ix);
            for(itheta=1; itheta<=ltheta0; itheta++)
            {
                i_th = id3_eic_th(+itheta,iz,ix);
                out[i_th] = out[i_th0];
                i_th = id3_eic_th(-itheta,iz,ix);
                out[i_th] = out[i_th0];
            }
        }
return;
}

void taperADCIG(float *out, float *eimage_theta, int nx, int nz, int ntheta, float dtheta, float thetaIni, float thetaFin, float power) {

    int ix, iz, itheta, i;
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;
    int i_th0, i_th;
    int ltheta0 = ntheta/2;
    int ithetaIni = thetaIni/dtheta;
    int ithetaFin = thetaFin/dtheta;

    int ithetaMin, ithetaMax;
    
    float taper, theta;
    int step;
    if(ithetaIni<ithetaFin) {
        ithetaMin = ithetaIni;
        ithetaMax = ithetaFin;
    }
    else {
        ithetaMax = ithetaIni;
        ithetaMin = ithetaFin;
    }
    for(ix=0; ix<nx; ix++)
        for(iz=0; iz<nz; iz++)
            for(itheta=ithetaMin; itheta<=ithetaMax; itheta++)
            {
                theta = itheta*dtheta;
                i_th = id3_eic_th(itheta,iz,ix);
                
                taper = (thetaFin-theta)/(thetaFin-thetaIni);
                if(power==0.0)    taper = 0.0;
                else              taper = pow(taper,power);
                out[i_th] = taper * eimage_theta[i_th];
            }
return;
}
void makeReciprocalADCIG(float *out, float *inp, int nx, int nz, int nxb, int nzb, int ntheta, int acqGeom) {

    int ix, iz, itheta;
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;
    int ltheta0 = ntheta/2;

    for(ix=ixMin; ix<ixMax; ix++)
    {
        for(iz=izMin; iz<izMax; iz++)
            for(itheta=1; itheta<=ltheta0; itheta++)
            {
                int i_inp = id3_eic_th(-acqGeom*itheta, iz, ix);
                int i_out = id3_eic_th(+acqGeom*itheta, iz, ix);
                out[i_out] = inp[i_inp];
            }
    }
return;
}

void applyCosADCIG(float *out, float *eimage_theta, int nx, int nz, int ntheta, float dtheta) {

    int ix, iz, itheta, i;
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;
    int i_th0, i_th;
    int ltheta0 = ntheta/2;

    float theta;
    for(ix=0; ix<nx; ix++)
        for(iz=0; iz<nz; iz++)
            for(itheta=-ltheta0; itheta<=ltheta0; itheta++)
            {
                theta     = itheta*dtheta * M_PI/180.0;
                i_th      = id3_eic_th(itheta,iz,ix);
                out[i_th] = cosf(theta) * eimage_theta[i_th];
            }
return;
}

void flatSpikeADCIG(float *out, int nx, int nz, float dz, int ntheta, float dtheta, float z) {

    int ix, iz, itheta, i;
    int i_th;
    int ltheta0 = ntheta/2;

    iz = floorf(0.5+z/dz) + nzb;

    for(ix=0; ix<nx; ix++)
        for(itheta=-ltheta0; itheta<=ltheta0; itheta++)
        {
            i_th = id3_eic_th(itheta,iz,ix);
            out[i_th] = 1.0;
        }
    
return;
}
double getArrayL2Norm(float *array, long long int n) {
    long long int i;
    double L2Norm = 0.0;
    for(i=0; i<n; i++) {
        L2Norm += array[i]*array[i];
    }
    L2Norm = sqrt(L2Norm);
return L2Norm;
}
void normalizeArray(float *array, long long int n) {
    long long int i;
    float max = -9999.0;
    for(i=0; i<n; i++) {
        if(max<fabsf(array[i]))    max = fabsf(array[i]);
    }
    if(max>0.0) {
        for(i=0; i<n; i++)    array[i] /= max;
    }
return;
}
void normalizeArraysJointly(float *array1, long long int n1, float *array2, long long int n2) {
    long long int i;
    float max = -9999.0;
    for(i=0; i<n1; i++) {
        if(max<fabsf(array1[i]))    max = fabsf(array1[i]);
    }
    for(i=0; i<n2; i++) {
        if(max<fabsf(array2[i]))    max = fabsf(array2[i]);
    }
    if(max>0.0) {
        for(i=0; i<n1; i++)    array1[i] /= max;
        for(i=0; i<n2; i++)    array2[i] /= max;
    }
return;
}
void takeAbsluteValueArrays(float *out, float *inp, long long int n) {
    long long int i;
    for(i=0; i<n; i++) {
        out[i] = fabsf(inp[i]);
    }
return;
}
void takeSquareRootArrays(float *out, float *inp, long long int n) {
    long long int i;
    for(i=0; i<n; i++) {
        out[i] = sqrtf(fabsf(inp[i]));
    }
return;
}
void multiplyArrays(float *out, float *inp1, float *inp2, long long int n) {
    long long int i;
    for(i=0; i<n; i++) {
        out[i] = inp1[i]*inp2[i];
    }
return;
}
void addArrays(float *out, float alpha1, float *inp1, float alpha2, float *inp2, long long int n) {
    long long int i;
    for(i=0; i<n; i++) {
        out[i] = alpha1*inp1[i] + alpha2*inp2[i];
    }
return;
}
void divideArraysPow(float *out, float pow1, float *inp1, float pow2, float *inp2, long long int n) {
    long long int i;
    for(i=0; i<n; i++) {
        out[i] = powf(inp1[i],pow1) / powf(inp2[i],pow2);
    }
return;
}

void multiplyArrayByScalar(float alpha, float *array, long long int n) {
    long long int i;
    for(i=0; i<n; i++) {
        array[i] *= alpha;
    }
return;
}

double computeL2NormArray(float *array, long long int n) {
    long long int i;
    double L2Norm = 0.0;
    for(i=0; i<n; i++) {
        L2Norm += array[i]*array[i];
    }
    L2Norm /= n;
return L2Norm;
}

void divideArrays(float *out, float *inp1, float *inp2, long long int n) {
    long long int i;
    double avg = 0;
    double max = -9999999;
    long long int count = 0;
    for(i=0; i<n; i++) {
        if(inp2[i]!=0.0)   {
            avg += inp2[i]*inp2[i];
            if(max < inp2[i])    max = inp2[i];
            count++;
        }
    }
    avg = sqrt(avg)/count;
    for(i=0; i<n; i++) {
        // out[i] = inp1[i]/(inp2[i]+3000.0*avg);
        out[i] = inp1[i]/(inp2[i]+0.4*max);
    }
    printf("\n\n avg=%g  max=%g\n\n", avg, max);
return;
}
void divideArraysWhiteNoise(float *out, float *num, float *den, long long int n, float noiseLevel) {
    long long int i;
    for(i=0; i<n; i++)    out[i] = num[i]/(den[i]+noiseLevel);
return;
}
//*
float* residualMoveOut(float *eimage_theta, int nx, int nz, int ntheta, float dx, float dz, float lagMax) {

    int ix, itheta;
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;

    int ltheta0 = ntheta/2;

    float window = 2*lagMax;
    float dz_ = dz/10.0f;
    int   nz_ = floorf(nz*dz/dz_);
    int   nwin = window/dz_;
    if(nwin%2==0)    nwin++;  // Make it odd
    
    float *resMO_lags   = CPU_zaloc1F(nx*nz*ntheta);
    float *resMO_corr   = CPU_zaloc1F(nx*nz*ntheta);

    for(ix=ixMin; ix<ixMax; ix++)
        for(itheta=0; itheta<(ntheta/2); itheta++)
        {
            float *thetaTrace       = CPU_zaloc1F(nz );
            float *thetaTraceInterp = CPU_zaloc1F(nz_);
            getThetaTrace(eimage_theta, thetaTrace, nx, nz, ntheta, itheta, ix);
            interpTraceLinear(thetaTrace, thetaTraceInterp, nz, dz, nz_, dz_);
            free(thetaTrace);
            free(thetaTraceInterp);
        }
    free(resMO_corr);

return resMO_lags;
}
void getXCorr(float *trcA, float *trcB, int *lags, float *corr, int n, int shiftMax) {
    int i, lag, is;
    for(i=0; i<n; i++) {
        lags[i] = 0;
        corr[i] = -999999.0f;
    }
    for(lag=-shiftMax; lag<=shiftMax; lag++) {         // loops over lags
        for(i=0; i<n; i++) {                           // loops over trace samples
            float corr_ = 0;
            float normA = 0;
            float normB = 0;
            for(is=-shiftMax; is<=shiftMax; is++) {    // sums over a window
                corr_ += trcA[i+is] * trcB[i+lag+is];
                normA += trcA[i+is] * trcA[i+lag+is];
                normB += trcA[i+is] * trcA[i+lag+is];
            }
            if(corr[i]<corr_) {                        // Updats corr and lags
                corr[i] = corr_;
                lags[i] = lag;
            }

        }
    }
}


void getThetaTrace(float *eimage_theta, float *thetaTrace, int nx, int nz, int ntheta, int itheta, int ix) {
    int iz;
    int ltheta0 = ntheta/2;
    for(iz=0; iz<nz; iz++)
    {
        int i_th = id3_eic_th(itheta,iz,ix);
        thetaTrace[iz] = eimage_theta[i_th];
    }
return;
}
void transp(float *ou, float *in, int n1, int n2) {
    int i2, i1;
    for(i2=0; i2<n2; i2++) {
        for(i1=0; i1<n1; i1++) {
            ou[i2+i1*n2] = in[i1+i2*n1];
        }
    }
return;
}
void smoothTrace(float *gather, int ntrc, int nsmp, int nsmt) {
    int itrc, ismp, ismt;
    for(itrc=0; itrc<ntrc; itrc++)
    {
        float *tmp = CPU_zaloc1F(nsmp);
            for(ismp=nsmt; ismp<nsmp-nsmt; ismp++)
            {
                for(ismt=-nsmt; ismt<=nsmt; ismt++)
                    tmp[ismp] += gather[ismp+ismt+itrc*nsmp]/(2*nsmt);
            }
            for(ismp=nsmt; ismp<nsmp-nsmt; ismp++)
                gather[ismp+itrc*nsmp] = tmp[ismp];
        free(tmp);
    }

return;
}
//*/
float* getADCIGDipGrad(float *eimage_theta, int nx, int nz, int ntheta, float dx, float dtheta) {

    int ix, iz, itheta, i;

    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;

    int ltheta0 = ntheta/2;
    
    float *tmp_x   = CPU_zaloc1F(nx*nz*ntheta);
    float *tmp_z   = CPU_zaloc1F(nx*nz*ntheta);
    float *dipgrad = CPU_zaloc1F(nx*nz*ntheta);

    for(ix=0; ix<nx; ix++)
        for(iz=0; iz<nz-1; iz++)
            for(itheta=0; itheta<ntheta-1; itheta++) {
                int i_th0  = id3_eic_th(itheta+0,iz+0,ix);
                int i_th1  = id3_eic_th(itheta+1,iz+0,ix);
                int i_th2  = id3_eic_th(itheta+0,iz+1,ix);
                tmp_x[i_th0] = (eimage_theta[i_th1] - eimage_theta[i_th0])/dtheta;
                tmp_z[i_th0] = (eimage_theta[i_th2] - eimage_theta[i_th0])/dx;
            }
    for(i=0; i<nx*nz*ntheta; i++) {
        dipgrad[i] = atan2(tmp_z[i],tmp_x[i]);
    }
    free(tmp_x);
    free(tmp_z);
return dipgrad;
}
void dsoTauP(float *adcigtp, float *adcigtp_dso, int nx, int np, int ntau, float dp, float dtau) {

    int ix, ip, itau;

    int ixMin = nxb;
    int ixMax = nx-nxb;
    int np0 = np/2;

    for(ix=ixMin; ix<ixMax; ix++)
        for(itau=0; itau<ntau; itau++)
            for(ip=1; ip<np-1; ip++)
            {
                int itp   = itau +  ip   *ntau + ix*np*ntau;
                int itp_m = itau + (ip-1)*ntau + ix*np*ntau;
                int itp_p = itau + (ip+1)*ntau + ix*np*ntau;
                //adcigtp_dso[itp] = (adcigtp[itp_p] - adcigtp[itp_m])/(2.0*dp);
                adcigtp_dso[itp] = (adcigtp[itp_p] + adcigtp[itp_m] - 2.0f*adcigtp[itp])/(dp*dp);
            }
return;
}
void dsoTauP_env(float *adcigtp, float *adcigtp_dso, int nx, int np, int ntau, float dp, float dtau) {

    int ix, ip, itau;

    int ixMin = nxb;
    int ixMax = nx-nxb;
    int np0 = np/2;

    for(ix=ixMin; ix<ixMax; ix++)
        for(itau=0; itau<ntau; itau++)
            for(ip=1; ip<np-1; ip++)
            {
                int itp   = itau +  ip   *ntau + ix*np*ntau;
                int itp_m = itau + (ip-1)*ntau + ix*np*ntau;
                int itp_p = itau + (ip+1)*ntau + ix*np*ntau;
                adcigtp_dso[itp] = (adcigtp[itp_p] - adcigtp[itp_m])/(2.0*dp);
                //adcigtp_dso[itp] = (adcigtp[itp_p] + adcigtp[itp_m] - 2.0f*adcigtp[itp])/(dp*dp);
            }
return;
}
//*
void ODCIG_localSlantStack(int adj, float *eimage_localSlantStack, float *eimage, \
                           int nx, float dx, int nz, float dz, int lx0, int ntheta, float dtheta) {

    int ix, iz, lx, iz_, lx_, il, itheta;
    int nlx=2*lx0+1;
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;

    int ltheta0 = ntheta/2;
    float theta;
    int   nl = 11;
    float L  = (nl-1)*dx;
    int   nl0 = (nl-1)/2;
    float dl = L/(nl-1);
    float l, Delta_lx, Delta_z;
    float stack_value, plx_0, plx_1, pz_0, pz_1;
    float lx_pos, z_pos;

    printf("\n\n ntheta=%d  ltheta0=%d  nlx=%d  lx0=%d \n\n", ntheta, ltheta0, nlx, lx0);

    // ODCIG local slant-stack
    int   i, i_z0_lx0, i_z0_lx1, i_z1_lx0, i_z1_lx1;
    int  cond;
    for(ix=ixMin;ix<ixMax;ix++)
        for(iz=izMin;iz<izMax;iz++) 
            for(lx=-lx0; lx<=lx0; lx++) 
            {
                // DeMoveout
                lx_pos = lx*dx;
                z_pos  = iz*dz;

                for(itheta=-ltheta0; itheta<=ltheta0; itheta++)
                {
                    printf("\n\n  ix=%d  iz=%d  lx=%d  itheta=%d \n\n", ix, iz, lx, itheta);
                    // Slant-Stack
                    if(adj==0)
                        stack_value = 0.0f;
                    else if(adj==1)
                        stack_value = eimage_localSlantStack[id3_eic_lx_lst(itheta,lx,iz,ix)];

                    theta = itheta * dtheta * M_PI/180.0;

                    for(il=0; il<nl; il++)
                    {
                        l = (il-nl0)*dl;
                        Delta_lx = + l * cosf(theta);
                        Delta_z  = - l * sinf(theta);

                        lx_   = floorf((lx_pos+Delta_lx)/dx);
                        plx_1 = ((lx_pos+Delta_lx) - lx_*dx)/dx;
                        plx_0 = 1.0 - plx_1;

                        iz_   = floorf((z_pos+Delta_z)/dz);
                        pz_1  = ((z_pos+Delta_z) - iz_*dz)/dz;
                        pz_0  = 1.0 - pz_1;

                        i_z0_lx0 = id3_eic_lx(lx_  ,iz_  ,ix);
                        i_z1_lx0 = id3_eic_lx(lx_  ,iz_+1,ix);
                        i_z0_lx1 = id3_eic_lx(lx_+1,iz_  ,ix);
                        i_z1_lx1 = id3_eic_lx(lx_+1,iz_+1,ix);

                        cond = ( (-lx0)<=lx_ && lx_<=(lx0-1) )  &&  ( 0<=iz_  &&  iz_<nz-1 );
                        if( !cond )    continue;

                        if(adj==0) {
                            stack_value += plx_0*pz_0 * eimage[i_z0_lx0]
                                         + plx_0*pz_1 * eimage[i_z1_lx0]
                                         + plx_1*pz_0 * eimage[i_z0_lx1]
                                         + plx_1*pz_1 * eimage[i_z1_lx1];
                        }
                        else if(adj==1) {
                            eimage[i_z0_lx0] += stack_value * plx_0*pz_0;
                            eimage[i_z1_lx0] += stack_value * plx_0*pz_1; 
                            eimage[i_z0_lx1] += stack_value * plx_1*pz_0;
                            eimage[i_z1_lx1] += stack_value * plx_1*pz_1;
                        }
                    } // loop along slant-stack line segment
                    if(adj==0)
                        eimage_localSlantStack[id3_eic_lx_lst(itheta,lx,iz,ix)] += stack_value;
                } // loop over theta
            } // loop over lx
}
void ODCIG_MoveIn_LocalSlantStack(float *eimage_SlantStack, float *eimage_SlantStack_deMO, \
                                int nx, float dx, int nz, float dz, int lx0, \
                                float contraction_factor, float lx_null, int ntheta, float dtheta) {

    int ix, iz, lx, iz_, lx_, il, itheta;
    int nlx=2*lx0+1;
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;

    int ltheta0 = ntheta/2;
    float theta;

    float stack_value, plx_0, plx_1, pz_0, pz_1;
    float lx_pos, z_pos, lx_move, z_move, move, sign_move;

    // ODCIG contraction
    int   i, i_z0_lx0, i_z0_lx1, i_z1_lx0, i_z1_lx1;
    int  cond;
    for(ix=0;ix<nx;ix++)
        for(iz=0;iz<nz;iz++) 
            for(lx=-lx0; lx<=lx0; lx++) 
            {
                // Original position
                lx_pos = lx*dx; 
                z_pos  = iz*dz;
                for(itheta=-ltheta0; theta<=ltheta0; itheta++)
                {
                    theta = itheta * dtheta * M_PI/180.0;

                    stack_value = eimage_SlantStack[id3_eic_lx_lst(itheta,lx,iz,ix)];

                    // Compute movein to be applied
                    if ( abs(lx)*dx>=lx_null )
                        move = (1.0-contraction_factor) * (abs(lx)*dx-lx_null);
                    else
                        move = 0.0;                    
                    if(fsignf(theta)==fsignf(lx)) move *= -1.0;
                    if(fsignf(theta)!=fsignf(lx)) move *= +1.0;

                    // New position after move in
                    lx_move = lx_pos + move * sinf(theta);
                     z_move =  z_pos + move * cosf(theta);

                    lx_   = floorf(lx_move/dx);
                    plx_1 = (lx_move - lx_*dx)/dx;
                    plx_0 = 1.0 - plx_1;

                    iz_   = floorf(z_move/dz);
                    pz_1  = (z_move - iz_*dz)/dz;
                    pz_0  = 1.0 - pz_1;

                    i_z0_lx0 = id3_eic_lx_lst(itheta, lx_  ,iz_  ,ix);
                    i_z1_lx0 = id3_eic_lx_lst(itheta, lx_  ,iz_+1,ix);
                    i_z0_lx1 = id3_eic_lx_lst(itheta, lx_+1,iz_  ,ix);
                    i_z1_lx1 = id3_eic_lx_lst(itheta, lx_+1,iz_+1,ix);

                    cond = ( (-lx0)<=lx_ && lx_<=(lx0-1) )  &&  ( 0<=iz_  &&  iz_<nz-1 );
                    if( !cond )    continue;

                    eimage_SlantStack_deMO[i_z0_lx0] += plx_0*pz_0 * stack_value;
                    eimage_SlantStack_deMO[i_z1_lx0] += plx_0*pz_1 * stack_value;
                    eimage_SlantStack_deMO[i_z0_lx1] += plx_1*pz_0 * stack_value;
                    eimage_SlantStack_deMO[i_z1_lx1] += plx_1*pz_1 * stack_value;
                    
                } // loop over theta
            } // loop over lx
return;
}
void lxMoveIn(int add, float *eimage, float *eimage_movein, \
              int nx, float dx, int nz, float dz, int lx0, \
              int ntheta, float dtheta, \
              float contraction_factor, float lx_null) {

    int ix, iz, lx;
    int nlx=2*lx0+1;
    int adj;

    float *eimage_localSlantStack        = CPU_zaloc1F(nx*nz*nlx*ntheta);
    float *eimage_localSlantStack_moveIn = CPU_zaloc1F(nx*nz*nlx*ntheta);

    if(add == 0)
    {
        for(ix=0;ix<nx;ix++)
            for(iz=0;iz<nz;iz++) 
                for(lx=-lx0; lx<=lx0; lx++) {
                    int   i  = id3_eic_lx(lx,iz,ix);
                    eimage_movein[i] = 0.0f;
                }
    }

    adj = 0;
    ODCIG_localSlantStack(adj, eimage_localSlantStack, eimage, nx, dx, nz, dz, lx0, ntheta, dtheta);
    // ODCIG_MoveIn_LocalSlantStack(eimage_localSlantStack, eimage_localSlantStack_moveIn, nx, dx, nz, dz, lx0, contraction_factor, lx_null, ntheta, dtheta);
    adj = 1;
    ODCIG_localSlantStack(adj, eimage_localSlantStack, eimage_movein, nx, dx, nz, dz, lx0, ntheta, dtheta);
    // ODCIG_localSlantStack(adj, eimage_localSlantStack_moveIn, eimage_movein, nx, dx, nz, dz, lx0, ntheta, dtheta);

    free(eimage_localSlantStack);
    free(eimage_localSlantStack_moveIn);

    
return;
}
//*/
void lx2DeMoveout_forward3(int add, float *eimage, float *eimage_deMO, \
                          int nx, float dx, int nz, float dz, int lx0, \
                          float contraction_factor, float lx1, float lx2, float lx3, float lx4) {

    int ix, iz, lx, iz_, lx_, il;
    int nlx=2*lx0+1;
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;

    float thetaMin = -75.0 * M_PI/180.0;
    float thetaMax = +75.0 * M_PI/180.0;
    float dtheta = 1.0;
    float theta;
    int   nl = 11;
    float L  = (nl-1)*dx;
    int   nl0 = (nl-1)/2;
    float dl = L/(nl-1);
    float l, Delta_lx, Delta_z;
    float stack_value, plx_0, plx_1, pz_0, pz_1;
    float lx_pos, z_pos, lx_move, z_move, move, sign_move;

    if(add == 0)
    {
        for(ix=ixMin;ix<ixMax;ix++)
            for(iz=izMin;iz<izMax;iz++) 
                for(lx=-lx0; lx<=lx0; lx++) {
                    int   i_th  = id3_eic_lx(lx,iz,ix);
                    eimage_deMO[i_th] = 0.0f;
                }
    }


    int   i, i_z0_lx0, i_z0_lx1, i_z1_lx0, i_z1_lx1;
    int  cond;
    for(ix=0;ix<nx;ix++)
        for(iz=0;iz<nz;iz++) 
            for(lx=-lx0; lx<=lx0; lx++) 
            {
                lx_pos = lx*dx; 
                z_pos  = iz*dz;
                for(theta=thetaMin; theta<=thetaMax; theta+=dtheta)
                {
                    // Slant-Stack
                    stack_value = 0.0f;
                    for(il=0; il<nl; il++) {
                        l = (il-nl0)*dl;
                        Delta_lx = + l * cosf(theta);
                        Delta_z  = - l * sinf(theta);

                        lx_   = floorf((lx_pos+Delta_lx)/dx);
                        plx_1 = ((lx_pos+Delta_lx) - lx_*dx)/dx;
                        plx_0 = 1.0 - plx_1;

                        iz_   = floorf((z_pos+Delta_z)/dz);
                        pz_1  = ((z_pos+Delta_z) - iz_*dz)/dz;
                        pz_0  = 1.0 - pz_1;

                        i_z0_lx0 = id3_eic_lx(lx_  ,iz_  ,ix);
                        i_z1_lx0 = id3_eic_lx(lx_  ,iz_+1,ix);
                        i_z0_lx1 = id3_eic_lx(lx_+1,iz_  ,ix);
                        i_z1_lx1 = id3_eic_lx(lx_+1,iz_+1,ix);

                        // cond = (<=lx_ && lx_<lx0-1) && (0<=iz_ && iz_<nz-1);
                        cond = ( (-lx0)<=lx_ && lx_<=(lx0-1) )  &&  ( 0<=iz_  &&  iz_<nz-1 );
                        if( !cond )    continue;

                        stack_value += plx_0*pz_0 * eimage[i_z0_lx0]
                                     + plx_0*pz_1 * eimage[i_z1_lx0]
                                     + plx_1*pz_0 * eimage[i_z0_lx1]
                                     + plx_1*pz_1 * eimage[i_z1_lx1];
                    }

                    // Demoveout                    
                    float ll = (abs(lx)*dx);
                    if( ll<lx1 )
                        move = 0.0;
                    else if ( lx1<=ll  &&  ll<lx2 )
                        move = (1.0-contraction_factor) * (ll-lx1);
                    else if ( lx2<=ll  &&  ll<lx3 )
                        move = (1.0-contraction_factor) * (lx2-lx1);
                    else if ( lx3<=ll  &&  ll<lx4 )
                        move = (1.0-contraction_factor) * (lx2-lx1) * (lx4-ll)/(lx4-lx3);
                    else if ( lx4<=ll )
                        move = 0.0f;

                    if(fsignf(theta)==fsignf(lx)) move *= -1.0;
                    if(fsignf(theta)!=fsignf(lx)) move *= +1.0;
                    lx_move = lx_pos + move * sinf(theta);
                     z_move =  z_pos + move * cosf(theta);

                    // Slant-Spread
                    for(il=0; il<nl; il++) {
                        l = (il-nl0)*dl;
                        Delta_lx = + l * cosf(theta);
                        Delta_z  = - l * sinf(theta);

                        lx_   = floorf((lx_move+Delta_lx)/dx);
                        plx_1 = ((lx_move+Delta_lx) - lx_*dx)/dx;
                        plx_0 = 1.0 - plx_1;

                        iz_   = floorf((z_move+Delta_z)/dz);
                        pz_1  = ((z_move+Delta_z) - iz_*dz)/dz;
                        pz_0  = 1.0 - pz_1;

                        i_z0_lx0 = id3_eic_lx(lx_  ,iz_  ,ix);
                        i_z1_lx0 = id3_eic_lx(lx_  ,iz_+1,ix);
                        i_z0_lx1 = id3_eic_lx(lx_+1,iz_  ,ix);
                        i_z1_lx1 = id3_eic_lx(lx_+1,iz_+1,ix);

                        cond = ( (-lx0)<=lx_ && lx_<=(lx0-1) )  &&  ( 0<=iz_  &&  iz_<nz-1 );
                        if( !cond )    continue;

                        eimage_deMO[i_z0_lx0] += plx_0*pz_0 * stack_value;
                        eimage_deMO[i_z1_lx0] += plx_0*pz_1 * stack_value;
                        eimage_deMO[i_z0_lx1] += plx_1*pz_0 * stack_value;
                        eimage_deMO[i_z1_lx1] += plx_1*pz_1 * stack_value;
                    }
                } // loop over theta
            } // loop over lx

return;
}
void lx2DeMoveout_forward(int add, float *eimage, float *eimage_deMO, \
                          int nx, float dx, int nz, float dz, int lx0, \
                          float xMin, float xMax, float zMin, float zMax,
                          float contraction_factor, float lx_null) {

    int ix, iz, lx, iz_, lx_, il;
    int nlx=2*lx0+1;
    int ixMin = nxb + xMin/dx;
    int ixMax = nxb + xMax/dx;
    int izMin = nzb + zMin/dz;
    int izMax = nzb + zMax/dz;

    float thetaMin = -55.0 * M_PI/180.0;
    float thetaMax = +55.0 * M_PI/180.0;
    float dtheta = 1.0 * M_PI/180.0;
    float theta;
    int   nl = 11;
    float L  = (nl-1)*dx;
    int   nl0 = (nl-1)/2;
    float dl = L/(nl-1);
    float l, Delta_lx, Delta_z;
    float stack_value, plx_0, plx_1, pz_0, pz_1;
    float lx_pos, z_pos, lx_move, z_move, move, sign_move;

    if(add == 0)
    {
        for(ix=ixMin;ix<ixMax;ix++)
            for(iz=izMin;iz<izMax;iz++) 
                for(lx=-lx0; lx<=lx0; lx++) {
                    int   i_th  = id3_eic_lx(lx,iz,ix);
                    eimage_deMO[i_th] = 0.0f;
                }
    }

    // ODCIG contraction
    int   i, i_z0_lx0, i_z0_lx1, i_z1_lx0, i_z1_lx1;
    int  cond;
    for(ix=ixMin;ix<ixMax;ix++)
        for(iz=izMin;iz<izMax;iz++)
            for(lx=-lx0; lx<=lx0; lx++)
            {
                // DeMoveout
                lx_pos = lx*dx; 
                z_pos  = iz*dz;
                for(theta=thetaMin; theta<=thetaMax; theta+=dtheta)
                {
                    // Slant-Stack
                    stack_value = 0.0f;
                    for(il=0; il<nl; il++) {
                        l = (il-nl0)*dl;
                        Delta_lx = + l * cosf(theta);
                        Delta_z  = - l * sinf(theta);

                        lx_   = floorf((lx_pos+Delta_lx)/dx);
                        plx_1 = ((lx_pos+Delta_lx) - lx_*dx)/dx;
                        plx_0 = 1.0 - plx_1;

                        iz_   = floorf((z_pos+Delta_z)/dz);
                        pz_1  = ((z_pos+Delta_z) - iz_*dz)/dz;
                        pz_0  = 1.0 - pz_1;

                        i_z0_lx0 = id3_eic_lx(lx_  ,iz_  ,ix);
                        i_z1_lx0 = id3_eic_lx(lx_  ,iz_+1,ix);
                        i_z0_lx1 = id3_eic_lx(lx_+1,iz_  ,ix);
                        i_z1_lx1 = id3_eic_lx(lx_+1,iz_+1,ix);

                        // cond = (<=lx_ && lx_<lx0-1) && (0<=iz_ && iz_<nz-1);
                        cond = ( (-lx0)<=lx_ && lx_<=(lx0-1) )  &&  ( 0<=iz_  &&  iz_<nz-1 );
                        if( !cond )    continue;

                        stack_value += plx_0*pz_0 * eimage[i_z0_lx0]
                                     + plx_0*pz_1 * eimage[i_z1_lx0]
                                     + plx_1*pz_0 * eimage[i_z0_lx1]
                                     + plx_1*pz_1 * eimage[i_z1_lx1];
                    }
                    // Slant-Spread
                    /*
                    if ( (abs(lx)*dx)>=lx_null )
                        move = (1.0-contraction_factor) * (abs(lx)*dx-lx_null);// * cosf(theta);
                    else
                        move = 0.0;
                    //*/

                    //*
                    if ( (abs(lx)*dx)>=lx_null )
                        move = (1.0-contraction_factor) * lx_null;
                    else
                        move = (1.0-contraction_factor) * abs(lx)*dx / cosf(theta);
                    //*/

                    // sign_move = - fsignf(lx) * fsignf(theta);
                    // move *= sign_move;
                    // lx_move = lx_pos + move*sinf(fabsf(theta));
                    //  z_move =  z_pos + move*cosf(fabsf(theta));
                    
                    if(fsignf(theta)==fsignf(lx)) move *= -1.0;
                    if(fsignf(theta)!=fsignf(lx)) move *= +1.0;
                    lx_move = lx_pos + move * sinf(theta);
                     z_move =  z_pos + move * cosf(theta);
                    //  z_move =  z_pos + move;

                    for(il=0; il<nl; il++)
                    {
                        l = (il-nl0)*dl;
                        Delta_lx = + l * cosf(theta);
                        Delta_z  = - l * sinf(theta);

                        lx_   = floorf((lx_move+Delta_lx)/dx);
                        plx_1 = ((lx_move+Delta_lx) - lx_*dx)/dx;
                        plx_0 = 1.0 - plx_1;

                        iz_   = floorf((z_move+Delta_z)/dz);
                        pz_1  = ((z_move+Delta_z) - iz_*dz)/dz;
                        pz_0  = 1.0 - pz_1;

                        i_z0_lx0 = id3_eic_lx(lx_  ,iz_  ,ix);
                        i_z1_lx0 = id3_eic_lx(lx_  ,iz_+1,ix);
                        i_z0_lx1 = id3_eic_lx(lx_+1,iz_  ,ix);
                        i_z1_lx1 = id3_eic_lx(lx_+1,iz_+1,ix);

                        cond = ( (-lx0)<=lx_ && lx_<=(lx0-1) )  &&  ( 0<=iz_  &&  iz_<nz-1 );
                        if( !cond )    continue;

                        eimage_deMO[i_z0_lx0] += plx_0*pz_0 * stack_value;
                        eimage_deMO[i_z1_lx0] += plx_0*pz_1 * stack_value;
                        eimage_deMO[i_z0_lx1] += plx_1*pz_0 * stack_value;
                        eimage_deMO[i_z1_lx1] += plx_1*pz_1 * stack_value;
                    }
                } // loop over theta
            } // loop over lx
return;
}
void lx2DeMoveout_forward2(int add, float *eimage, float *eimage_deMO, \
                          int nx, float dx, int nz, float dz, int lx0, \
                          float contraction_factor, float lx_null) {

    int ix, iz, lx, iz_, lx_, il;
    int nlx=2*lx0+1;
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;

    float thetaMin = -75.0 * M_PI/180.0;
    float thetaMax = +75.0 * M_PI/180.0;
    float dtheta = 1.0;
    float theta;
    int   nl = 11;
    float L  = (nl-1)*dx;
    int   nl0 = (nl-1)/2;
    float dl = L/(nl-1);
    float l, Delta_lx, Delta_z;
    float stack_value, plx_0, plx_1, pz_0, pz_1;
    float lx_pos, z_pos, lx_move, z_move, move, sign_move;

    int   nl_ = nx;
    int   nl0_ = (nl_-1)/2;

    if(add == 0)
    {
        for(ix=ixMin;ix<ixMax;ix++)
            for(iz=izMin;iz<izMax;iz++) 
                for(lx=-lx0; lx<=lx0; lx++) {
                    int   i_th  = id3_eic_lx(lx,iz,ix);
                    eimage_deMO[i_th] = 0.0f;
                }
    }

    // ODCIG contraction
    int   i, i_z0_lx0, i_z0_lx1, i_z1_lx0, i_z1_lx1;
    int  cond;
    for(ix=0;ix<nx;ix++)
        for(iz=0;iz<nz;iz++) 
            for(lx=-lx0; lx<=lx0; lx++) 
            {
                // DeMoveout
                lx_pos = lx*dx; 
                z_pos  = iz*dz;
                for(theta=thetaMin; theta<=thetaMax; theta+=dtheta) {
                    // Slant-Stack
                    stack_value = 0.0f;
                    for(il=0; il<nl; il++) {
                        l = (il-nl0)*dl;
                        Delta_lx = + l * cosf(theta);
                        Delta_z  = - l * sinf(theta);

                        lx_   = floorf((lx_pos+Delta_lx)/dx);
                        plx_1 = ((lx_pos+Delta_lx) - lx_*dx)/dx;
                        plx_0 = 1.0 - plx_1;

                        iz_   = floorf((z_pos+Delta_z)/dz);
                        pz_1  = ((z_pos+Delta_z) - iz_*dz)/dz;
                        pz_0  = 1.0 - pz_1;

                        i_z0_lx0 = id3_eic_lx(lx_  ,iz_  ,ix);
                        i_z1_lx0 = id3_eic_lx(lx_  ,iz_+1,ix);
                        i_z0_lx1 = id3_eic_lx(lx_+1,iz_  ,ix);
                        i_z1_lx1 = id3_eic_lx(lx_+1,iz_+1,ix);

                        // cond = (<=lx_ && lx_<lx0-1) && (0<=iz_ && iz_<nz-1);
                        cond = ( (-lx0)<=lx_ && lx_<=(lx0-1) )  &&  ( 0<=iz_  &&  iz_<nz-1 );
                        if( !cond )    continue;

                        stack_value += plx_0*pz_0 * eimage[i_z0_lx0]
                                     + plx_0*pz_1 * eimage[i_z1_lx0]
                                     + plx_1*pz_0 * eimage[i_z0_lx1]
                                     + plx_1*pz_1 * eimage[i_z1_lx1];
                    }
                    // Slant-Spread
                    if ( (abs(lx)*dx)>=lx_null )
                        move = (1.0-contraction_factor) * (abs(lx)*dx-lx_null);
                    else
                        move = 0.0;

                    // sign_move = - fsignf(lx) * fsignf(theta);
                    // move *= sign_move;
                    // lx_move = lx_pos + move*sinf(fabsf(theta));
                    //  z_move =  z_pos + move*cosf(fabsf(theta));
                    
                    if(fsignf(theta)==fsignf(lx)) move *= -1.0;
                    if(fsignf(theta)!=fsignf(lx)) move *= +1.0;
                    lx_move = lx_pos + move * sinf(theta);
                     z_move =  z_pos + move * cosf(theta);


                    for(il=0; il<nl_; il++) {
                        l = (il-nl0_)*dl;
                        Delta_lx = + l * cosf(theta);
                        Delta_z  = - l * sinf(theta);

                        lx_   = floorf((lx_move+Delta_lx)/dx);
                        plx_1 = ((lx_move+Delta_lx) - lx_*dx)/dx;
                        plx_0 = 1.0 - plx_1;

                        iz_   = floorf((z_move+Delta_z)/dz);
                        pz_1  = ((z_move+Delta_z) - iz_*dz)/dz;
                        pz_0  = 1.0 - pz_1;

                        i_z0_lx0 = id3_eic_lx(lx_  ,iz_  ,ix);
                        i_z1_lx0 = id3_eic_lx(lx_  ,iz_+1,ix);
                        i_z0_lx1 = id3_eic_lx(lx_+1,iz_  ,ix);
                        i_z1_lx1 = id3_eic_lx(lx_+1,iz_+1,ix);

                        cond = ( (-lx0)<=lx_ && lx_<=(lx0-1) )  &&  ( 0<=iz_  &&  iz_<nz-1 );
                        if( !cond )    continue;

                        eimage_deMO[i_z0_lx0] += plx_0*pz_0 * stack_value;
                        eimage_deMO[i_z1_lx0] += plx_0*pz_1 * stack_value;
                        eimage_deMO[i_z0_lx1] += plx_1*pz_0 * stack_value;
                        eimage_deMO[i_z1_lx1] += plx_1*pz_1 * stack_value;
                    }
                } // loop over theta
            } // loop over lx
return;
}
void lx2DeMoveout_new(int add, float *eimage, float *eimage_deMO, \
                      int nx, float dx, int nz, float dz, int lx0, \
                      float contraction_factor, float lx_null) {

    int ix, iz, lx, iz_, lx_, il;
    int nlx=2*lx0+1;
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;

    float thetaMin = -75.0 * M_PI/180.0;
    float thetaMax = +75.0 * M_PI/180.0;
    float dtheta = 1.0;
    float theta;
    int   nl = 11;
    float L  = (nl-1)*dx;
    int   nl0 = (nl-1)/2;
    float dl = L/(nl-1);
    float l, Delta_lx, Delta_z;
    float stack_value, plx_0, plx_1, pz_0, pz_1;
    float lx_pos, z_pos, lx_move, z_move, move, sign_move;

    if(add == 0)
    {
        for(ix=ixMin;ix<ixMax;ix++)
            for(iz=izMin;iz<izMax;iz++) 
                for(lx=-lx0; lx<=lx0; lx++) {
                    int   i_th  = id3_eic_lx(lx,iz,ix);
                    eimage_deMO[i_th] = 0.0f;
                }
    }

    // ODCIG contraction
    int   i, i_z0_lx0, i_z0_lx1, i_z1_lx0, i_z1_lx1;
    int  cond;
    for(ix=0;ix<nx;ix++)
        for(iz=0;iz<nz;iz++) 
            for(lx=-lx0; lx<=lx0; lx++) 
            {
                // DeMoveout
                lx_pos = lx*dx; 
                z_pos  = iz*dz;
                for(theta=thetaMin; theta<=thetaMax; theta+=dtheta) {
                    // Slant-Stack
                    stack_value = 0.0f;
                    for(il=0; il<nl; il++) {
                        l = (il-nl0)*dl;
                        Delta_lx = + l * cosf(theta);
                        Delta_z  = - l * sinf(theta);

                        lx_   = floorf((lx_pos+Delta_lx)/dx);
                        plx_1 = ((lx_pos+Delta_lx) - lx_*dx)/dx;
                        plx_0 = 1.0 - plx_1;

                        iz_   = floorf((z_pos+Delta_z)/dz);
                        pz_1  = ((z_pos+Delta_z) - iz_*dz)/dz;
                        pz_0  = 1.0 - pz_1;

                        i_z0_lx0 = id3_eic_lx(lx_  ,iz_  ,ix);
                        i_z1_lx0 = id3_eic_lx(lx_  ,iz_+1,ix);
                        i_z0_lx1 = id3_eic_lx(lx_+1,iz_  ,ix);
                        i_z1_lx1 = id3_eic_lx(lx_+1,iz_+1,ix);

                        // cond = (<=lx_ && lx_<lx0-1) && (0<=iz_ && iz_<nz-1);
                        cond = ( (-lx0)<=lx_ && lx_<=(lx0-1) )  &&  ( 0<=iz_  &&  iz_<nz-1 );
                        if( !cond )    continue;

                        stack_value += plx_0*pz_0 * eimage[i_z0_lx0]
                                     + plx_0*pz_1 * eimage[i_z1_lx0]
                                     + plx_1*pz_0 * eimage[i_z0_lx1]
                                     + plx_1*pz_1 * eimage[i_z1_lx1];
                    }
                    // Slant-Spread
                    if ( (abs(lx)*dx)>=lx_null )
                        move = (1.0-contraction_factor) * (abs(lx)*dx-lx_null);
                    else
                        move = 0.0;

                    // sign_move = - fsignf(lx) * fsignf(theta);
                    // move *= sign_move;
                    // lx_move = lx_pos + move*sinf(fabsf(theta));
                    //  z_move =  z_pos + move*cosf(fabsf(theta));
                    
                    if(fsignf(theta)==fsignf(lx)) move *= -1.0;
                    if(fsignf(theta)!=fsignf(lx)) move *= +1.0;
                    lx_move = lx_pos + move * sinf(theta);
                     z_move =  z_pos + move * cosf(theta);


                    for(il=0; il<nl; il++) {
                        l = (il-nl0)*dl;
                        Delta_lx = + l * cosf(theta);
                        Delta_z  = - l * sinf(theta);

                        lx_   = floorf((lx_move+Delta_lx)/dx);
                        plx_1 = ((lx_move+Delta_lx) - lx_*dx)/dx;
                        plx_0 = 1.0 - plx_1;

                        iz_   = floorf((z_move+Delta_z)/dz);
                        pz_1  = ((z_move+Delta_z) - iz_*dz)/dz;
                        pz_0  = 1.0 - pz_1;

                        i_z0_lx0 = id3_eic_lx(lx_  ,iz_  ,ix);
                        i_z1_lx0 = id3_eic_lx(lx_  ,iz_+1,ix);
                        i_z0_lx1 = id3_eic_lx(lx_+1,iz_  ,ix);
                        i_z1_lx1 = id3_eic_lx(lx_+1,iz_+1,ix);

                        cond = ( (-lx0)<=lx_ && lx_<=(lx0-1) )  &&  ( 0<=iz_  &&  iz_<nz-1 );
                        if( !cond )    continue;

                        eimage_deMO[i_z0_lx0] += plx_0*pz_0 * stack_value;
                        eimage_deMO[i_z1_lx0] += plx_0*pz_1 * stack_value;
                        eimage_deMO[i_z0_lx1] += plx_1*pz_0 * stack_value;
                        eimage_deMO[i_z1_lx1] += plx_1*pz_1 * stack_value;
                    }
                } // loop over theta
            } // loop over lx
return;
}
//*/
void lx2Contraction_forward(int add, float *eimage, float *eimage_contracted, \
                            int nx, float dx, int nz, int lx0, float contraction_factor) {

    int ix, iz, lx;
    int nlx=2*lx0+1;
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;

    if(add == 0)
    {
        for(ix=ixMin;ix<ixMax;ix++)
            for(iz=izMin;iz<izMax;iz++) 
                for(lx=-lx0; lx<=lx0; lx++) {
                    int   i_th  = id3_eic_lx(lx,iz,ix);
                    eimage_contracted[i_th] = 0.0f;
                }
    }

    // ODCIG contraction
    float lambda, lambda_;
    float plx_0, plx_1, plx_tot;
    int   lx_0, lx_1;
    int   i, i_0, i_1;
/*
    for(ix=0;ix<nx;ix++)
        for(iz=0;iz<nz;iz++) 
            for(lx=-lx0; lx<=lx0; lx++) 
            {
                lambda  = lx * dx;
                lambda_ = lambda * contraction_factor;
                lx_0    = floorf(lambda_/dx);
                lx_1    = lx_0 + 1;

                plx_0 = (dx*lx_1 - lambda_);
                plx_1 = (lambda_ - dx*lx_0);
                plx_tot = plx_0 + plx_1;
                plx_0 /= plx_tot;
                plx_1 /= plx_tot;

                // printf("\n plx_0=%g    plx_1=%g", plx_0, plx_1);

                i      = id3_eic_lx(lx  ,iz,ix);
                i_0    = id3_eic_lx(lx_0,iz,ix);
                i_1    = id3_eic_lx(lx_1,iz,ix);
                eimage_contracted[i_0] += contraction_factor * plx_0 * eimage[i];
                eimage_contracted[i_1] += contraction_factor * plx_1 * eimage[i];
            }
//*/
    for(ix=0;ix<nx;ix++)
        for(iz=0;iz<nz;iz++) 
            for(lx=-lx0; lx<=lx0; lx++) 
            {
                lambda  = lx * dx;
                lambda_ = lambda * contraction_factor;
                lx_0    = floorf(lambda_/dx);

                float alpha = (lambda_ - dx*lx_0)/dx;
                float *coef = interSinc1D_8P(alpha);

                i = id3_eic_lx(lx, iz, ix);
                int i_;
                for(i_=0; i_<8; i_++)
                {
                    int ll = lx_0-3+i_;
                    int ii =  id3_eic_lx(ll, iz, ix);
                    if(-lx0<=ll  &&  ll<=lx0)
                        eimage_contracted[ii] += contraction_factor * coef[i_] * eimage[i];
                }
                free(coef);
            }
return;
}
/*
void taperGrad(float *grad, int nx, int nz, float dx, float dz, float z1, float z2, float z3, float z4) {
    int ix, iz;
    float z, x, taper;
    for(ix=nxb; ix<nx-nxb; ix++)
        for(iz=nzb; iz<nz-nzb; iz++) {
            z = (iz-nzb) * dz;
            if(z<z1)    taper = 0.0;
            else if(z1<=z  &&  z<z2)    taper = (z-z1)/(z2-z1);
            else if(z2<=z1  )

        }
}
//*/

void lxApplyBandpass(float *eimage, int nx, int nz, float dz, int lx0, float lambda1, float lambda2, float lambda3, float lambda4)
{
    int ix, iz, lx;
    int nlx=2*lx0+1;

    float k1, k2, k3, k4;

    // printf("\n\n lambda1=%f  lambda2=%f  lambda3=%f  lambda4=%f  nz=%d  dz=%f \n\n", lambda1, lambda2, lambda3, lambda4, nz, dz);

    if(lambda1==0.0)    k1 = 1.0/(2.0*dz);
    else                k1 = 1.0/lambda1;

    if(lambda2==0.0)    k2 = 1.0/(2.0*dz);
    else                k2 = 1.0/lambda2;
    
    if(lambda3==0.0)    k3 = 1.0/(2.0*dz);
    else                k3 = 1.0/lambda3;
    
    if(lambda4==0.0)    k4 = 1.0/(2.0*dz);
    else                k4 = 1.0/lambda4;

    // printf("\n\n k1=%f  k2=%f  k3=%f  k4=%f \n\n", k1, k2, k3, k4);

    // ODCIG transpose
    float *eimage_transp = CPU_zaloc1F(nx*nz*nlx);
    for(ix=0; ix<nx; ix++)
        for(iz=0; iz<nz; iz++)
            for(lx=-lx0; lx<=lx0; lx++) 
            {
                int i  = id3_eic_lx(lx,iz,ix);
                int i_ = iz + (lx+lx0)*nz + ix*nz*nlx;
                eimage_transp[i_] = eimage[i]; 
            }
    // ODCIG bandpass
    for(ix=0; ix<nx; ix++)
        for(lx=-lx0; lx<=lx0; lx++) 
        {
            size_t shift = (lx+lx0)*nz + ix*nz*nlx;
            bandpass(eimage_transp+shift, nz, dz, k1, k2, k3, k4);
        }
    // ODCIG transpose again
    for(ix=0; ix<nx; ix++)
        for(iz=0; iz<nz; iz++)
            for(lx=-lx0; lx<=lx0; lx++) 
            {
                int i  = id3_eic_lx(lx,iz,ix);
                int i_ = iz + (lx+lx0)*nz + ix*nz*nlx;
                eimage[i] = eimage_transp[i_]; 
            }
    free(eimage_transp);

return;
}
void modelBandpass(float *model, int nx, int nz, float dz, float k1, float k2, float k3, float k4) {
    int ix;
    for(ix=0; ix<nx; ix++)    bandpass(model+(ix*nz), nz, dz, k1, k2, k3, k4);
}
void modelBandpass_Lambdas(float *model, int nx, int nz, float dz, float lambda1, float lambda2, float lambda3, float lambda4) {
    int ix;
    
    float k1, k2, k3, k4;

    if(lambda1==0.0)    k1 = 1.0/(2.0*dz);
    else                k1 = 1.0/lambda1;

    if(lambda2==0.0)    k2 = 1.0/(2.0*dz);
    else                k2 = 1.0/lambda2;
    
    if(lambda3==0.0)    k3 = 1.0/(2.0*dz);
    else                k3 = 1.0/lambda3;
    
    if(lambda4==0.0)    k4 = 1.0/(2.0*dz);
    else                k4 = 1.0/lambda4;
    
    for(ix=0; ix<nx; ix++)    bandpass(model+(ix*nz), nz, dz, k1, k2, k3, k4);
}
void modelSmooth_Target(float *out, float *inp, int nx, int nz, int rx, int rz, int min_z, int max_z) {
    int ix, iz, ix_, iz_;
    int nn, i, i_;
    float *tmp = CPU_zaloc1F(nx*nz);

    for(ix=0; ix<nx; ix++) 
    {
        for(iz=min_z; iz<max_z; iz++)
        {
            nn = 0;
            i = iz + ix*nz;
            int ix_min = max(0   ,ix-rx);
            int ix_max = min(nx-1,ix+rx);
            int iz_min = max(max(0   ,iz-rz), min_z);
            int iz_max = min(min(nz-1,iz+rz), max_z);
            for(ix_=ix_min; ix_<=ix_max; ix_++)
            {
                for(iz_=iz_min; iz_<=iz_max; iz_++)
                {
                    nn++;
                    i_ = iz_ + ix_*nz;
                    tmp[i] += inp[i_];
                }
            }
            if(nn>0)
                tmp[i] /= nn;
            else
                tmp[i] = 0;
        }
    }

    for(ix=0; ix<nx; ix++) 
        for(iz=0; iz<nz; iz++)
        {
            i = iz + ix*nz;
            if(min_z<=iz &&  iz<max_z)
                out[i] = tmp[i];
            else 
                out[i] = inp[i];
        }
    
    free(tmp);
}
void modelSmooth(float *out, float *inp, int nx, int nz, int rx, int rz) {
    int ix, iz, ix_, iz_;
    int nn, i, i_;
    float *tmp = CPU_zaloc1F(nx*nz);

    for(ix=0; ix<nx; ix++) 
    {
        for(iz=0; iz<nz; iz++)
        {
            nn = 0;
            i = iz + ix*nz;
            int ix_min = max(0   ,ix-rx);
            int ix_max = min(nx-1,ix+rx);
            int iz_min = max(0   ,iz-rz);
            int iz_max = min(nz-1,iz+rz);
            for(ix_=ix_min; ix_<=ix_max; ix_++)
            {
                for(iz_=iz_min; iz_<=iz_max; iz_++)
                {
                    if(iz_>=nz  ||  iz_<0  ||  ix_>=nx  ||  ix_<0)
                    {
                        printf("\n\n WEIRD ix_=%d  iz_=%d  ix=%d  iz=%d  rx=%d  rz=%d\n", ix_, iz_, ix, iz, rx, rz);
                        printf("ix_min=%d  ix_max=%d  iz_min=%d  iz_max=%d", ix_min,  ix_max,  iz_min,  iz_max);
                        printf("max(0,ix-rx)=%d    min(nx-1,ix+rx)=%d\n", max(0,ix-rx), min(nx-1,ix+rx));
                        printf("max(0,iz-rz)=%d    min(nz-1,iz+rz)=%d\n\n", max(0,iz-rz), min(nz-1,iz+rz));
                    }
                    nn++;
                    i_ = iz_ + ix_*nz;
                    tmp[i] += inp[i_];
                }
            }
            if(nn>0)
                tmp[i] /= nn;
            else
                tmp[i] = 0;
        }
    }

    for(ix=0; ix<nx; ix++) 
        for(iz=0; iz<nz; iz++)
        {
            i = iz + ix*nz;
            out[i] = tmp[i];
        }
    
    free(tmp);
}

void modelTaperingZ(float *model, float *model_tapered, \
                    int nx, int nz, float dz, float z0, float z1) {

    int ix, iz;
    float taper, z;

    // Tapering
    for(ix=0;ix<nx;ix++)
        for(iz=0;iz<nz;iz++)  
        {
            if(z0<z1)
            {
                z = dz*(iz-nzb);
                if(z0<=z && z<=z1)    taper = (z-z0)/(z1-z0);
                else if(z0>z)         taper = 0.0;
                else                  taper = 1.0;
            }
            else if(z0>z1)
            {
                z = dz*(iz-nzb);
                if(z1<=z && z<=z0)    taper = (z0-z)/(z0-z1);
                else if(z0<z)         taper = 0.0;
                else                  taper = 1.0;
            }
            else
            {
                taper = 1.0;
            }
            int i = iz + ix*nz;
            model_tapered[i] = taper * model[i];
        }
return;
}
void modelTaperingX(float *model, float *model_tapered, \
                    int nx, int nz, float dx, float x1, \
                    float x2, float x3, float x4) {

    int ix, iz;
    float taper, x, ox;

    // ODCIG contraction
    ox = -nxb*dx;
    for(ix=0;ix<nx;ix++)
    {
        x = ix*dx + ox;
        if(x<x1)                     taper = 0.0;
        else if(x1<=x  &&  x<x2 )    taper = (x-x1)/(x2-x1);
        else if(x2<=x  &&  x<=x3)    taper = 1.0;
        else if(x3<x   &&  x<=x4)    taper = (x4-x)/(x4-x3);
        else if(x>x4)                taper = 0.0;

        for(iz=0;iz<nz;iz++)  
        {
            int i = iz + ix*nz;
            model_tapered[i] = taper * model[i];
        }
    }
return;
}
void modelTaperingX_old(float *model, float *model_tapered, \
                        int nx, int nz, float dx, float x0, float x1) {

    int ix, iz;
    float taper, x;

    // ODCIG contraction
    for(ix=0;ix<nx;ix++)
        for(iz=0;iz<nz;iz++)  
        {
            if(x0<x1)
            {
                x = dx*(ix-nxb);
                if(x0<=x && x<=x1)    taper = (x-x0)/(x1-x0);
                else if(x0>x)         taper = 0.0;
                else                  taper = 1.0;
            }
            else
            {
                x = dx*(ix-nxb);
                if(x1<=x && x<=x0)    taper = (x0-x)/(x0-x1);
                else if(x0<x)         taper = 0.0;
                else                  taper = 1.0;
            }
            int i = iz + ix*nz;
            model_tapered[i] = taper * model[i];
        }
return;
}
void lxTapering(float *eimage, float *eimage_tapered, \
               int nx, int nz, float dz, int lx0, int lz0, float z0, float z1) {

    int ix, iz, lx, lz;
    int nlx=2*lx0+1;
    int nlz=2*lz0+1;

    // ODCIG contraction
    for(ix=0;ix<nx;ix++)
        for(iz=0;iz<nz;iz++)  
        {
            float taper;
            float z = dz*(iz-nzb);
            
            if(z0<=z1)
            {
                if(z0<=z && z<=z1)    taper = (z-z0)/(z1-z0);
                else if(z0>z)         taper = 0.0;
                else                  taper = 1.0;
            }
            else
            {
                if(z1<=z && z<=z0)    taper = (z-z1)/(z0-z1);
                else if(z0<z)         taper = 0.0;
                else                  taper = 1.0;
            }

            // printf("\n Weird stuff: z=%f  z0=%f  z1=%f  taper=%f  dz=%f  iz=%d  nzb=%d", z, z0, z1, taper, dz, iz, nzb);

            // if(taper<1.0 && ix==nx/2)  printf("\n Will it perform the filtering Z taper=%f ?\n", taper);
            
            // for(lx=-lx0; lx<=lx0; lx++) 
            //     for(lz=-lz0; lz<=lz0; lz++) 
            //     {
            //         int i = id3_eic_lx_lz(lx,lz,iz,ix);
            //         eimage_tapered[i] = taper * eimage[i];
            //     }
            for(lx=-lx0; lx<=lx0; lx++) 
            {
                int i = (lx+lx0) + iz*nlx + ix*nz*nlx; 
                eimage_tapered[i] = taper * eimage[i];
            }
        }
return;
}
void lxTapering_MultiDim(float *eimage, float *eimage_tapered, \
                         int nx, int nz, float dx, float dz, int lx0, \
                         float z0, float z1, \
                         float x0, float x1) {

    int ix, iz, lx;
    int nlx=2*lx0+1;

    // ODCIG contraction
    for(ix=0;ix<nx;ix++)
    {
        float taper_x;
        float x = dx*(ix-nxb);
        if(x0<=x && x<=x1)    taper_x = (x-x0)/(x1-x0);
        else if(x0>x)         taper_x = 0.0;
        else                  taper_x = 1.0;
        for(iz=0;iz<nz;iz++)
        {
            float taper_z;
            float z = dz*(iz-nzb);
            if(z0<=z && z<=z1)    taper_z = (z-z0)/(z1-z0);
            else if(z0>z)         taper_z = 0.0;
            else                  taper_z = 1.0;
            for(lx=-lx0; lx<=lx0; lx++) 
            {
                float taper = taper_x * taper_z;
                int i = id3_eic_lx(lx,iz,ix);
                eimage_tapered[i] = taper * eimage[i];
            }
        }
    }
return;
}
void lxTapering_x(float *eimage, float *eimage_tapered, \
                  int nx, int nz, float dx, int lx0, \
                  float x1, float x2, float x3, float x4) {

    int ix, iz, lx;
    int nlx=2*lx0+1;

    // ODCIG contraction
    for(ix=0;ix<nx;ix++)
    {
        float taper_x;
        float x = dx*(ix-nxb);
        if(x1<=x && x<=x2)         taper_x = (x-x1)/(x2-x1);
        else if(x1>x)              taper_x = 0.0;
        else if(x3<=x && x<=x4)    taper_x = (x4-x)/(x4-x3);
        else if(x4<x)              taper_x = 0.0;
        else                       taper_x = 1.0;
        for(iz=0;iz<nz;iz++)
        {
            for(lx=-lx0; lx<=lx0; lx++) 
            {
                int i = id3_eic_lx(lx,iz,ix);
                eimage_tapered[i] = taper_x * eimage[i];
            }
        }
    }
return;
}
void applyIlum2EImage(float *eimage, float *eimage_ilumin, \
                      float *image, float *image_ilumin, float *ilumin, \
                      int nx, int nz, float dx, int lx0, float epsilon) {

    int ix, iz, lx, i;
    int nlx=2*lx0+1;
    double ilm, ilMax=-1e20, ilMin=1e+20, ilAvg=0.0;
    double eps;

    int n = nx*nz;
    for(i=0; i<n; i++) {
        if(ilMax<ilumin[i])    ilMax = ilumin[i];
    }
    for(i=0; i<n; i++) {
        ilumin[i] /= ilMax;
    }
    int nn=0;
    for(i=0; i<n; i++)
    {
        if(ilMin>ilumin[i]  &&  ilumin[i]>0.0)
            ilMin = ilumin[i];

        if(ilumin[i]>0.0)
        {
            ilAvg += ilumin[i];
            nn++;
        }
    }
    ilAvg /= nn;

    if(epsilon == 0.0)
        eps = max(0.001*ilAvg, ilMin);
    else
        eps = epsilon;

    // Divide by illumination
    for(ix=0; ix<nx; ix++)
        for(iz=0; iz<nz; iz++)
        {
            if(nxb<=ix  &&  ix<nx-nxb  &&  nzb<=iz  &&  iz<nz-nzb)
                ilm = 1.0/(ilMax*(ilumin[iz+ix*nz]+eps));
            else
                ilm = 0.0;

            if(image!=NULL  &&  image_ilumin!=NULL)
                image_ilumin[iz+ix*nz] = ilm * image[iz+ix*nz];
            
            if(eimage!=NULL  &&  eimage_ilumin!=NULL)
            {
                for(lx=-lx0; lx<=lx0; lx++)
                {
                    int i = id3_eic_lx(lx,iz,ix);
                    eimage_ilumin[i] = ilm * eimage[i];
                }
            }
        }
return;
}
void lxTapering_lx(float *eimage, float *eimage_tapered, \
                   int nx, int nz, float dx, int lx0, \
                   float lx1, float lx2, float lx3, float lx4) {

    int ix, iz, lx;
    int nlx=2*lx0+1;

    // ODCIG contraction
    for(ix=0;ix<nx;ix++)
    {
        for(iz=0;iz<nz;iz++)
        {
            for(lx=-lx0; lx<=lx0; lx++) 
            {
                float taper_lx;
                float lxf = dx*lx;
                if(lx1<=lxf && lxf<=lx2)        taper_lx = (lxf-lx1)/(lx2-lx1);
                else if(lx1>lxf)                taper_lx = 0.0;
                else if(lx3<=lxf && lxf<=lx4)   taper_lx = (lx4-lxf)/(lx4-lx3);
                else if(lx4<lxf)                taper_lx = 0.0;
                else                            taper_lx = 1.0;

                int i = id3_eic_lx(lx,iz,ix);
                eimage_tapered[i] = taper_lx * eimage[i];
            }
        }
    }
return;
}
void lxNullTapering(float *eimage, float *eimage_tapered, int nx, int nz, float dx, int lx0, float lxNull) {

    int ix, iz, lx;
    int nlx=2*lx0+1;

    if(lxNull<dx)    return;

    // ODCIG contraction
    float lxNull_a = 0.5*lxNull;
    float lxNull_b = 1.0*lxNull;
    for(ix=0;ix<nx;ix++)
    {
        for(iz=0;iz<nz;iz++)
        {
            for(lx=-lx0; lx<=lx0; lx++) 
            {
                float taper_lx;
                float lxf = dx*lx;
                if(fabsf(lxf) <= lxNull_a)                     taper_lx = 0.0;
                else if(lxNull_a <= fabsf(lxf) <= lxNull_b)    taper_lx = (fabsf(lxf)-lxNull_a) / (lxNull_b-lxNull_a);
                else                                           taper_lx = 1.0;

                int i = id3_eic_lx(lx,iz,ix);
                eimage_tapered[i] = taper_lx * eimage[i];
                // eimage_tapered[i] *= 1.5f;//*eimage[i];
            }
        }
    }
return;
}
//*
void lxTapering_lx_varXZ(float *eimage, float *eimage_tapered, \
                         int nx, int nz, float dx, int lx0, \
                         TapeLX *tapelx, int ntape) {

    int ix, iz, lx;
    int nlx=2*lx0+1;
    float lx1, lx2, lx3, lx4, z, deltaZ;

    float lx1_i, lx2_i, lx3_i, lx4_i;
    float lx1_f, lx2_f, lx3_f, lx4_f;
    float z_i, z_f;

    // Applying taper
    for(ix=0; ix<nx;ix++)
    {
        // Finding adjacent taper point to the left
        float x = (ix-nxb)*dx;
        int itape = 0;
        while( x>tapelx[itape].xpos  &&  itape<ntape)    itape++;

        if(itape>0) itape--;

        printf("\n x=%f  itape=%d  xpos=%f  nxb=%d", x, itape, tapelx[itape].xpos, nxb);
        
        // Finding weighs
        if( (itape==ntape-1)  ||  (itape==0  &&  x<tapelx[itape].xpos))
        {
            lx1_i = tapelx[itape].lx1_i;
            lx2_i = tapelx[itape].lx2_i;
            lx3_i = tapelx[itape].lx3_i;
            lx4_i = tapelx[itape].lx4_i;
            lx1_f = tapelx[itape].lx1_f;
            lx2_f = tapelx[itape].lx2_f;
            lx3_f = tapelx[itape].lx3_f;
            lx4_f = tapelx[itape].lx4_f;
            z_i = tapelx[itape].z_i;
            z_f = tapelx[itape].z_f;
        }
        else
        {
            float p0 = tapelx[itape+1].xpos - x;
            float p1 = x - tapelx[itape].xpos;
            float ptot = p0 + p1;

            p0 /= ptot;
            p1 /= ptot;

            lx1_i = p0*tapelx[itape].lx1_i + p1*tapelx[itape+1].lx1_i;
            lx2_i = p0*tapelx[itape].lx2_i + p1*tapelx[itape+1].lx2_i;
            lx3_i = p0*tapelx[itape].lx3_i + p1*tapelx[itape+1].lx3_i;
            lx4_i = p0*tapelx[itape].lx4_i + p1*tapelx[itape+1].lx4_i;

            lx1_f = p0*tapelx[itape].lx1_f + p1*tapelx[itape+1].lx1_f;
            lx2_f = p0*tapelx[itape].lx2_f + p1*tapelx[itape+1].lx2_f;
            lx3_f = p0*tapelx[itape].lx3_f + p1*tapelx[itape+1].lx3_f;
            lx4_f = p0*tapelx[itape].lx4_f + p1*tapelx[itape+1].lx4_f;

            z_i   = p0*tapelx[itape].z_i   + p1*tapelx[itape+1].z_i;
            z_f   = p0*tapelx[itape].z_f   + p1*tapelx[itape+1].z_f;
        }

        for(iz=0; iz<nz;iz++)
        {
            z = (iz-nzb)*dx;
            if(z<z_i) {
                lx1 = lx1_i;
                lx2 = lx2_i;
                lx3 = lx3_i;
                lx4 = lx4_i;
            }
            else if(z>z_f) {
                lx1 = lx1_f;
                lx2 = lx2_f;
                lx3 = lx3_f;
                lx4 = lx4_f;
            }
            else {
                deltaZ = (z-z_i)/(z_f-z_i);
                lx1 = lx1_i + deltaZ * (lx1_f-lx1_i);
                lx2 = lx2_i + deltaZ * (lx2_f-lx2_i);
                lx3 = lx3_i + deltaZ * (lx3_f-lx3_i);
                lx4 = lx4_i + deltaZ * (lx4_f-lx4_i);
            }

            for(lx=-lx0; lx<=lx0; lx++) 
            {
                float taper_lx;
                float lxp = dx*lx;
                if(lx1<=lxp && lxp<=lx2)        taper_lx = (lxp-lx1)/(lx2-lx1);
                else if(lx1>lxp)                taper_lx = 0.0;
                else if(lx3<=lxp && lxp<=lx4)   taper_lx = (lx4-lxp)/(lx4-lx3);
                else if(lx4<lxp)                taper_lx = 0.0;
                else                            taper_lx = 1.0;

                int i = id3_eic_lx(lx,iz,ix);
                eimage_tapered[i] = taper_lx * eimage[i];
            }
        }
    }
return;
}
//*/
void lxTapering_lx_varDepth(float *eimage, float *eimage_tapered, \
                            int nx, int nz, float dx, int lx0, \
                            float lx1_i, float lx2_i, float lx3_i, float lx4_i, \
                            float lx1_f, float lx2_f, float lx3_f, float lx4_f, \
                            float z_i, float z_f) {

    int ix, iz, lx;
    int nlx=2*lx0+1;
    float lx1, lx2, lx3, lx4, z, deltaZ;

    // ODCIG contraction
    for(ix=0; ix<nx;ix++)
    {
        for(iz=0; iz<nz;iz++)
        {
            z = (iz-nzb)*dx;
            if(z<z_i) {
                lx1 = lx1_i;
                lx2 = lx2_i;
                lx3 = lx3_i;
                lx4 = lx4_i;
            }
            else if(z>z_f) {
                lx1 = lx1_f;
                lx2 = lx2_f;
                lx3 = lx3_f;
                lx4 = lx4_f;
            }
            else {
                deltaZ = (z-z_i)/(z_f-z_i);
                lx1 = lx1_i + deltaZ * (lx1_f-lx1_i);
                lx2 = lx2_i + deltaZ * (lx2_f-lx2_i);
                lx3 = lx3_i + deltaZ * (lx3_f-lx3_i);
                lx4 = lx4_i + deltaZ * (lx4_f-lx4_i);
            }

            for(lx=-lx0; lx<=lx0; lx++) 
            {
                float taper_lx;
                float lxp = dx*lx;
                if(lx1<=lxp && lxp<=lx2)        taper_lx = (lxp-lx1)/(lx2-lx1);
                else if(lx1>lxp)                taper_lx = 0.0;
                else if(lx3<=lxp && lxp<=lx4)   taper_lx = (lx4-lxp)/(lx4-lx3);
                else if(lx4<lxp)                taper_lx = 0.0;
                else                            taper_lx = 1.0;

                int i = id3_eic_lx(lx,iz,ix);
                eimage_tapered[i] = taper_lx * eimage[i];
            }
        }
    }
return;
}
void stacking_eimage_threads(int nthreads, float *eimage_input, float *eimage_reduced, int nx, int nz, int lx0) {
    int ix, iz, lx, ithread;
    int nlx=2*lx0+1;
    // ODCIG "reduction"
    printf("\n\n  Reducing eimage.  nthreads=%d  \n\n", nthreads);
    for(ithread=0; ithread<nthreads; ithread++)
    {
        size_t shift = nlx*nx*nz * ithread; 
        for(ix=0;ix<nx;ix++)
            for(iz=0;iz<nz;iz++)
                for(lx=-lx0; lx<=lx0; lx++) 
                {
                    int i = id3_eic_lx(lx,iz,ix);
                    eimage_reduced[i] += eimage_input[i + shift];
                }
    }
return;
}
void cpyAngleODCIG(float *eimage_temp, float *eimage_ampAngle, int nx, int nz, int lx0, int itheta) { 

    int ix, iz, lx;
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;
    int nlx=2*lx0+1;

    int shiftMem = nx*nz*nlx*itheta;

    for(ix=ixMin;ix<ixMax;ix++)
        for(iz=izMin;iz<izMax;iz++) {
            for(lx=-lx0; lx<=lx0; lx++)
            {
                int i = id3_eic_lx(lx,iz,ix);
                eimage_ampAngle[i] = eimage_temp[i];
            }
        }
return;
}
void cpyAngle(float *eimage_theta_temp, float *eimage_theta, int nx, int nz, int ntheta, int itheta) {

    int ix, iz;
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;

    int ltheta0 = ntheta/2;

    for(ix=ixMin; ix<ixMax; ix++)
        for(iz=izMin; iz<izMax; iz++)
        {
            int i_th = id3_eic_th(itheta,iz,ix);
            eimage_theta_temp[i_th] = eimage_theta[i_th];
        }
return;
}
void ADCIG2TauP_forward(float *adcigtp, float *eimage_theta, int nx, int nz, int ntheta, int np, int ntau, float dx, float dtheta, float dp, float dtau) {

    int ix, itau, ip, itheta;
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;

    int np0 = np/2;

    int ltheta0 = ntheta/2;
    float dz = dx;

    printf("\n\n  tau p: ntheta=%d np=%d ntau=%d dtheta=%f dp=%f dtau=%f \n\n", ntheta, np, ntau, dtheta, dp, dtau);

    for(ix=ixMin; ix<ixMax; ix++)
        for(itau=0; itau<ntau; itau++)
            for(ip=0; ip<np; ip++)
            {
                float p = (ip-np0)*dp;
                for(itheta=-ltheta0; itheta<=ltheta0; itheta++)
                {
                    float theta  = itheta * dtheta;
                    float tau    = itau*dtau;
                    float z      = tau + p*theta;
                    int iz0      = floor(z/dz);
                    int iz1      = iz0 + 1;
                    float w0     = (iz1*dz-z)/dz;
                    float w1     = 1.0f-w0;
                    iz0 += nzb;
                    iz1 += nzb;
                    if(iz1<=izMax && iz0>=izMin)
                    {
                        int i_th0    = id3_eic_th(itheta,iz0,ix);
                        int i_th1    = id3_eic_th(itheta,iz1,ix);
                        int itp      = itau + ip*ntau + ix*np*ntau;
                        adcigtp[itp] += w0 * eimage_theta[i_th0] + w1 * eimage_theta[i_th1];
                    }
                }
            }
return;
}

void ADCIG2TauP_adjoint(float *adcigtp, float *eimage_theta, int nx, int nz, int ntheta, int np, int ntau, float dx, float dtheta, float dp, float dtau) {

    int ix, itau, ip, itheta, i;
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;

    int np0 = np/2;

    int ltheta0 = ntheta/2;

    float dz = dx;

    for(i=0; i<nx*nz*ntheta; i++)
        eimage_theta[i] = 0.0;
                
    for(ix=ixMin; ix<ixMax; ix++)
        for(itau=0; itau<ntau; itau++)
            for(ip=0; ip<np; ip++)
            {
                float p = (ip-np0)*dp;
                for(itheta=0; itheta<ntheta; itheta++)
                {
                    float theta  = (itheta-ltheta0) * dtheta;
                    float tau    = itau*dtau;
                    float z      = tau + p*theta;
                    int iz0      = floor(z/dz);
                    int iz1      = iz0 + 1;
                    float w0     = (iz1*dz-z)/dz;
                    float w1     = 1.0f-w0;
                    iz0 += nzb;
                    iz1 += nzb;
                    if(iz1<=izMax && iz0>=izMin)
                    {
                        int i_th0    = id3_eic_th(itheta,iz0,ix);
                        int i_th1    = id3_eic_th(itheta,iz1,ix);
                        int itp      = itau + ip*ntau + ix*np*ntau;
                        eimage_theta[i_th0] += w0 * adcigtp[itp];
                        eimage_theta[i_th1] += w1 * adcigtp[itp];
                    }
                }
            }
return;
}

float fsignf(float x) {
    return (x/fabsf(x));
}
double fsign(double x) {
    return (x/fabs(x));
}
int sign(int x) {
    return (x/abs(x));
}

void transposeArray(float *array, int n1, int n2) {
    int i1, i2, i;
    float *transp = (float *) calloc(n2*n1, sizeof(float));
    for(i1=0; i1<n1; i1++)
        for(i2=0; i2<n2; i2++)    transp[i2+i1*n2] = array[i1+i2*n1];
    for(i=0; i<n1*n2; i++)   array[i] = transp[i];
    free(transp);
}

float findParabolaMin(float x0, float x1, float x2, float y0, float y1, float y2, float *ymin) {
    float a, b, c, xmin;
    getParabolaParams(&a, &b, &c, x0, x1, x2, y0, y1, y2);
    if(a<=0.0) {
        printf("\n\n  Warning: parabola has negative curvature! a=%g  \
        \n -->> Will return either of x0, x1, x2 corresponding to the least objective Function value value  y0, y1, y2 \n\n", a);
        if     (y0<y1 && y0<y2)    return x0;
        else if(y1<y2)             return x1;
        else                       return x2;
    } 
    printf("\n Parabola parameters: a=%g  b=%g  c=%g \n", a, b, c);
    xmin = -b/(2.0*a);
    *ymin = a*a*xmin + b*xmin + c;
return xmin;    
}

void getParabolaParams(float *a, float *b, float *c, double x0, double x1, double x2, double y0, double y1, double y2) {
    double d10 = y1-y0;
    double d20 = y2-y0;
    *a = ( d20 + d10*(x0-x2)/(x1-x0) ) / ( (x2*x2-x0*x0) + (x0*x0 - x1*x1) * (x2-x0)/(x1-x0) );
    *b = d10/(x1-x0) - (*a) * (x1*x1-x0*x0)/(x1-x0);
    *c = y0 - (*a)*x0*x0 - (*b)*x0;
return;
}

float* correlation(float *a, float *b, int na, int nb) {
    
    int t, tau, n;
    n = na + nb - 1;
    float *c = CPU_zaloc1F(n);

    int tmin = -na+1;
    int tmax =  na-1;

    for(t=tmin; t<=tmax; t++)
    {
        int tau_max = min(na-1-t,nb-1);
        int tau_min = max(-t,0);
        for(tau=tau_min; tau<=tau_max; tau++)
        {
             c[t-tmin] += b[tau] * a[t+tau];
        }
    }

return c;
}

float* convolution(float *a, float *b, int na, int nb) {
    
    int t, tau, n;
    n = na + nb - 1;
    float *c = CPU_zaloc1F(n);

    for(t=0; t<n; t++)
    {
        int tau_max = min(t,nb-1);
        for(tau=0; tau<=tau_max; tau++)
        {
             if( (t-tau)<na )    c[t] += b[tau] * a[t-tau];
        }
    }

return c;
}

float* codeGen(int n1, int j1, int ncodes, int seed)
{
    int i1, icode, jj;
    float *codes = CPU_zaloc1F(n1*ncodes);

    srand(seed);

    for(icode=0; icode<ncodes; icode++)
    {
        jj = j1;
        for(i1=0; i1<n1; i1+=jj)
        {
            float sample = ((1.0*rand())/(1.0f*RAND_MAX)) - 0.5;
            float spike;
            if(sample >= 0.0)    spike = +1.0;
            else                 spike = -1.0;
            if(i1>0)
                codes[i1+icode*n1] = spike;
            // codes[i1+icode*n1] = sample;
            jj = j1 * (1.0+2.0*sample);
            // printf("\n sample=%f  spike=%f  jj=%d  j1=%d", sample, spike, jj, j1);
        }
    }
return codes;
}

void resamp(float *trace, float *trace_resamp, int nt, float dt, int nt_resamp, float dt_resamp)
{
    int   i, it, it_resamp, r_=3;
    float t, alpha;
    for(it_resamp=0; it_resamp<nt_resamp; it_resamp++)
    {
        t = it_resamp * dt_resamp;
        it = floorf(t/dt);

        if(dt_resamp>dt)
        {
            alpha = (t-it*dt)/dt;
            float *coef = interSinc1D_8P(alpha);
            for(i=0; i<8; i++)
                if( (it-r_+i)>=0  &&  (it-r_+i)<nt )  trace_resamp[it_resamp] += trace[it-r_+i];
            free(coef);
        }
        else
        {
            float p0   = (it+1) * dt - t;
            float p1   = t - it*dt;
            float ptot = p0 + p1;
            p0 /= ptot;
            p1 /= ptot;
            if(it<nt-1) trace_resamp[it_resamp] = p1*trace[it+1] + p0*trace[it];
        }
    }
return;
}

void eimage_NonGeologicDipFiltering(float *eimage, int nx, int nz, float dx, float dz, int lx0, int lz0)
{
    int ix, iz, lx, lz;
    int nlx, nlz;
    int ikx, ikz;

    nlx = 2*lx0 + 1;
    nlz = 2*lz0 + 1;

    float dkx  = 1.0f/(nx * dx);
    float dkz  = 1.0f/(nz * dz);

    // Transposing to make offset the slowest dimension
    // float *eimage_transp_temp = transp_dim4(eimage            , nlx, nlz,  nz, nx, 13);
    // float *eimage_transp      = transp_dim4(eimage_transp_temp,  nz, nlz, nlx, nx, 24);
    float *eimage_transp_temp = transp_dim3(eimage            , nlx*nlz,      nz, nx, 12);
    float *eimage_transp      = transp_dim3(eimage_transp_temp,      nz, nlx*nlz, nx, 23);
    free(eimage_transp_temp);

    // Allocating arrays that will transformed CIGs
    fftw_complex *eimage_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*nz*nlx*nlz);

    // Allocating arrays that will receive corrected CIGs
    fftw_complex *eimage_filtered_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*nz*nlx*nlz);

    float *eimage_filtered_transp = CPU_zaloc1F(nx*nz*nlx*nlz);

    // Transforming to kx-kz
    for(lz=-lz0; lz<=lz0; lz++)
        for(lx=-lx0; lx<=lx0; lx++)
        {
            size_t shift = (lx+lx0)*nx*nz + (lz+lz0) * nz*nx*nlx;
            applyFFTW2D_forwardR2I(eimage_fft+shift, eimage_transp+shift, nx, nz);
       }
    // Applying dip correction to HOCIG
    int foundNaN = 0;
    for(lz=-lz0; lz<=lz0; lz++)
        for(lx=-lx0; lx<=lx0; lx++)
        {
            // float alpha_apparent = atan2f(-dz*lz,dx*lx); 
            double hz = dz*lz;
            double hx = dx*lx;
            double hh = sqrt(hx*hx + hz*hz);
            double alpha_lambda;
            if(hh>0.0)  alpha_lambda = acosf(hx/hh);
            else        alpha_lambda = 0.0;
            if( hx*hz > 0 )    alpha_lambda *= -1;
            
            for(ikz=0; ikz<nz; ikz++)
                for(ikx=0; ikx<nx; ikx++)
                {
                    double kx, kz;

                    if(ikx<=nx/2)    kx =      ikx  * dkx;
                    else             kx = -(nx-ikx) * dkx;
                    if(ikz<=nz/2)    kz =      ikz  * dkz;
                    else             kz = -(nz-ikz) * dkz;

                    // double k = sqrt(kx*kx + kz*kz) + 0.01 * 1.0/(dx*nx);  // 1.0/(dx*nx) i s dkx
                    // double alpha_dip = asinf(kx/k);

                    double k = sqrt(kx*kx + kz*kz);  // 1.0/(dx*nx) i s dkx
                    double alpha_dip;
                    if(k>0.0)      alpha_dip = asinf(kx/k);
                    else           alpha_dip = 0.0;

                    // double sigma = 10.0 * M_PI/180.0; // equivalent to 3 degrees
                    // double sigma2 = sigma*sigma;
                    // double filter = exp( - (alpha_lambda - alpha_dip)*(alpha_lambda - alpha_dip) / (2.0*sigma2) );

                    double sigma = cosf(88.0 * M_PI/180.0); // equivalent to 3 degrees
                    double sigma2 = sigma*sigma;
                    double projection;
                    double kx_normalized, kz_normalized, hx_normalized, hz_normalized;
                    if(k>0.0)  {
                        kx_normalized = kx / k;
                        kz_normalized = kz / k;
                    }
                    else {
                        kx_normalized = 0.0;
                        kz_normalized = 0.0;
                    }
                    if(hh>0.0)  {
                        hx_normalized = hx / hh;
                        hz_normalized = hz / hh;
                    }
                    else {
                        hx_normalized = 0.0;
                        hz_normalized = 0.0;
                    }
                    // if(k>0.0  &&  hh>0.0)  projection = (hx*kx + hz*kz) / (k*hh);
                    // else                   projection = 0.0;
                    projection = hx_normalized * kx_normalized + hz_normalized * kz_normalized;
                    double filter = exp( - projection*projection / (2.0*sigma2) );

                    if(isnan(filter))
                    {
                        foundNaN++;
                        printf("\n  sigma=%g  sigma2=%g  alpha_lambda=%g  alpha_dip=%g  filter=%g   ", sigma, sigma2, alpha_lambda, alpha_dip, filter); 
                    }

                    int i   = ikz + ikx*nz + (lx+lx0)*nz*nx + (lz+lz0)*nz*nx*nlx;
                    eimage_filtered_fft[i][0] = filter * eimage_fft[i][0];
                    eimage_filtered_fft[i][1] = filter * eimage_fft[i][1];
                }
        }
    printf("\n\n foundNaN=%d \n\n", foundNaN);
    // Transform back to space coordinates
    for(lz=-lz0; lz<=lz0; lz++)
        for(lx=-lx0; lx<=lx0; lx++) 
        {
            size_t shift = (lx+lx0)*nx*nz + (lz+lz0) * nz*nx*nlx;
            applyFFTW2D_backwardI2R(eimage_filtered_transp+shift, eimage_filtered_fft+shift, nx, nz);
            // applyFFTW2D_backwardI2R(eimage_filtered_transp+shift, eimage_fft+shift, nx, nz);
        }

    // float *eimage_filtered_transp_temp = transp_dim4(eimage_filtered_transp     , nz , nx, nlx, nlz, 13);
    // float *eimage_filtered             = transp_dim4(eimage_filtered_transp_temp, nlx, nx,  nz, nlz, 24);
    float *eimage_filtered_transp_temp = transp_dim3(eimage_filtered_transp     , nz,      nx, nlx*nlz, 23);
    float *eimage_filtered             = transp_dim3(eimage_filtered_transp_temp, nz, nlx*nlz,      nx, 12);
    free(eimage_filtered_transp_temp); 

    free(eimage_transp);
    free(eimage_filtered_transp);
    memcpy(eimage, eimage_filtered, nx*nz*nlx*nlz*sizeof(float));
    free(eimage_filtered);

return;
}

void VHocig_dispersalCorrection(float *Hocig, float *Vocig, int nx, int nz, float dx, float dz, int lx0, int lv0)
{
    int ix, iz, lx, lv;
    int nlx, nlv;
    int ikx, ikz;

    nlx = 2*lx0 + 1;
    nlv = 2*lv0 + 1;

    float dkx  = 1.0f/(nx * dx);
    float dkz  = 1.0f/(nz * dz);

    // Transposing to make offset the slowest dimension
    float *Hocig_transp = transp_dim3(Hocig, nlx, nz, nx, 13);
    float *Vocig_transp = transp_dim3(Vocig, nlv, nz, nx, 13);

    // Allocating arrays that will transformed CIGs
    fftw_complex *Hocig_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*nz*nlx);
    fftw_complex *Vocig_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*nz*nlv);

    fftw_complex *Hocig_corrected_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*nz*nlx);
    fftw_complex *Vocig_corrected_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*nz*nlv);

    float *Hocig_corrected_transp = CPU_zaloc1F(nx*nz*nlx);
    float *Vocig_corrected_transp = CPU_zaloc1F(nx*nz*nlv);

    // Transforming to kx-kz
    for(lx=-lx0; lx<=lx0; lx++) {
        size_t shift = (lx+lx0)*nx*nz;
        applyFFTW2D_forwardR2I(Hocig_fft+shift, Hocig_transp+shift, nx, nz);
    }
    for(lv=-lv0; lv<=lv0; lv++) {
        size_t shift = (lv+lv0)*nx*nz;
        applyFFTW2D_forwardR2I(Vocig_fft+shift, Vocig_transp+shift, nx, nz);
    }

    // Applying dip correction to HOCIG
    for(lx=-lx0; lx<=lx0; lx++)
    {
        for(ikz=0; ikz<nz; ikz++)
            for(ikx=0; ikx<nx; ikx++)
            {
                float kx = ikx * dkx;
                float kz = ikz * dkz;
                float alpha = - atan2f(kx,kz);
                float llx  = lx*dx;
                float llx_ = llx * cosf(alpha);
                float lx_a = floorf(llx_/dx);
                float lx_b = lx_a+1;
                float pa = (lx_b*dx - llx_)/dx;
                float pb = (llx_ - lx_a*dx)/dx;
                float ptot = pa + pb;
                pa /= ptot;
                pb /= ptot;

                int i   = ikx + ikz*nx + (lx  +lx0)*nz*nx;
                int i_a = ikx + ikz*nx + (lx_a+lx0)*nz*nx;
                int i_b = ikx + ikz*nx + (lx_b+lx0)*nz*nx;
                
                if(-lx0<=lx_a  &&  lx_a<=+lx0)
                {
                    Hocig_corrected_fft[i_a][0] += pa*Hocig_fft[i][0];
                    Hocig_corrected_fft[i_a][1] += pa*Hocig_fft[i][1];
                }
                if(-lx0<=lx_b  &&  lx_b<=+lx0)
                {
                    Hocig_corrected_fft[i_b][0] += pb*Hocig_fft[i][0];
                    Hocig_corrected_fft[i_b][1] += pb*Hocig_fft[i][1];
                }
            }
    }
    // Applying dip correction to VOCIG
    //*
    for(lv=-lv0; lv<=lv0; lv++)
    {
        for(ikz=0; ikz<nz; ikz++)
            for(ikx=0; ikx<nx; ikx++)
            {
                float kx = ikx * dkx;
                float kz = ikz * dkz;
                float alpha = - atan2f(kx,kz);
                float llv  = lv*dz;
                float llv_ = llv * sinf(alpha);
                float lv_a = floorf(llv_/dz);
                float lv_b = lv_a+1;
                float pa = (lv_b*dz - llv_)/dz;
                float pb = (llv_ - lv_a*dz)/dz;
                float ptot = pa + pb;
                pa /= ptot;
                pb /= ptot;

                int i   = ikx + ikz*nx + (lv  +lv0) * nz*nx;
                int i_a = ikx + ikz*nx + (lv_a+lv0) * nz*nx;
                int i_b = ikx + ikz*nx + (lv_b+lv0) * nz*nx;

                if(-lv0<=lv_a  &&  lv_a<=+lv0)
                {
                    Vocig_corrected_fft[i_a][0] += pa*Vocig_fft[i][0];
                    Vocig_corrected_fft[i_a][1] += pa*Vocig_fft[i][1];
                }
                if(-lv0<=lv_b  &&  lv_b<=+lv0)
                {
                    Vocig_corrected_fft[i_b][0] += pb*Vocig_fft[i][0];
                    Vocig_corrected_fft[i_b][1] += pb*Vocig_fft[i][1];
                }
            }
    }
    //*/
    

    // Transform back to space coordinates
    for(lx=-lx0; lx<=lx0; lx++) 
    {
        size_t shift = (lx+lx0)*nx*nz;
        applyFFTW2D_backwardI2R(Hocig_corrected_transp+shift, Hocig_corrected_fft+shift, nx, nz);
    }
    for(lv=-lv0; lv<=lv0; lv++) 
    {
        size_t shift = (lv+lv0)*nx*nz;
        applyFFTW2D_backwardI2R(Vocig_corrected_transp+shift, Vocig_corrected_fft+shift, nx, nz);
    }

    float *Hocig_corrected = transp_dim3(Hocig_corrected_transp, nx, nz, nlx, 13);
    float *Vocig_corrected = transp_dim3(Vocig_corrected_transp, nx, nz, nlv, 13);

    memcpy(Hocig, Hocig_corrected, nx*nz*nlx*sizeof(float));
    memcpy(Vocig, Vocig_corrected, nx*nz*nlv*sizeof(float));

    free(Hocig_transp);
    free(Vocig_transp);

    fftw_free(Hocig_fft);
    fftw_free(Vocig_fft);

    fftw_free(Hocig_corrected_fft);
    fftw_free(Vocig_corrected_fft);

    free(Hocig_corrected_transp);
    free(Vocig_corrected_transp);

    free(Hocig_corrected);
    free(Vocig_corrected);

return;
}
float* VHocig2Gocig(float *Hocig, float *Vocig, long long int nx, long long int nz, float dx, float dz, \
                    long long int lx0, long long int lv0, long long int lz0, int iproc)
{
    long long int i, ix, iz, lx, lv, lz;
    long long int nlx, nlv, nlz;
    long long int ikx, ikz;

    nlx = 2*lx0 + 1;
    nlv = 2*lv0 + 1;
    nlz = 2*lz0 + 1;

    float dkx  = 1.0f/(nx * dx);
    float dkz  = 1.0f/(nz * dz);

    // Allocating arrays that will transformed CIGs
    // fftw_complex *Hocig_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*nz*nlx);
    fftw_complex *Hocig_fft = CPU_zalocFFTW1(nx*nz*nlx);
    // fftw_complex *Vocig_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*nz*nlv);
    fftw_complex *Vocig_fft = CPU_zalocFFTW1(nx*nz*nlv);
    // fftw_complex *Gocig_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*nz*nlx*nlz);
    fftw_complex *Gocig_fft = CPU_zalocFFTW1(nx*nz*nlx*nlz);

    float *Gocig_transp     = CPU_zaloc1F(nx*nz*nlx*nlz);

    // Transposing to make offset the slowest dimension
    float *Hocig_transp = transp_dim3(Hocig, nlx, nz, nx, 13);
    float *Vocig_transp = transp_dim3(Vocig, nlv, nz, nx, 13);

    // Transforming to kx-kz
    for(lx=-lx0; lx<=lx0; lx++) {
        size_t shift = (lx+lx0)*nx*nz;
        applyFFTW2D_forwardR2I(Hocig_fft+shift, Hocig_transp+shift, nx, nz);
    }
    for(lv=-lv0; lv<=lv0; lv++) {
        size_t shift = (lv+lv0)*nx*nz;
        applyFFTW2D_forwardR2I(Vocig_fft+shift, Vocig_transp+shift, nx, nz);
    }
    /*
    // Computing amplitude spectrum for debuggin purposes
    float *ampSpec = CPU_zaloc1F(nx*nz);
    for(lx=-lx0; lx<=lx0; lx++) {
        size_t shift = (lx+lx0)*nx*nz;
        for(i=0; i<nx*nz; i++)
        {
            fftw_complex *spec = Hocig_fft+shift;
            ampSpec[i] = sqrtf(spec[i][0]*spec[i][0] + spec[i][1]*spec[i][1]);
        }
    }
    if(iproc==0)    outputSmart2d("AmplitudeSpectrum", ampSpec, nx, dx, 0, nz, dz, 0);
    free(ampSpec);
    if(iproc==0)
    {
        float *Hocig_fft_real = CPU_zaloc1F(nx*nz*nlx);
        float *Hocig_fft_imag = CPU_zaloc1F(nx*nz*nlx);
        long long int i;
        for(i=0; i<nx*nz*nlx; i++)      Hocig_fft_real[i] = Hocig_fft[i][0];
        for(i=0; i<nx*nz*nlx; i++)      Hocig_fft_imag[i] = Hocig_fft[i][1];
        outputSmart3d("Hocig_transp",   Hocig_transp, nx, dx, 0, nz, dz, 0, nlx   , dx    , -lx0*dx);
        outputSmart3d("Hocig_fft_real", Hocig_fft_real, nx, dx, 0, nz, dz, 0, nlx   , dx    , -lx0*dx);
        outputSmart3d("Hocig_fft_imag", Hocig_fft_imag, nx, dx, 0, nz, dz, 0, nlx   , dx    , -lx0*dx);
        free(Hocig_fft_real);
        free(Hocig_fft_imag);
    }
    //*/

    // Applying dip correction to HOCIG
    for(lx=-lx0; lx<=lx0; lx++)
    {
        for(ikz=0; ikz<nz; ikz++)
            for(ikx=0; ikx<nx; ikx++)
            {
                float kx, kz;
                if(ikx<nx/2.0)    kx =       ikx  * dkx;
                else              kx = -((nx-ikx) * dkx);
                if(ikz<nz/2.0)    kz =       ikz  * dkz;
                else              kz = -((nz-ikz) * dkz);
                float alpha    = atan2f(kx,kz);

                float llx      =   lx*dx;
                float lambda   =   llx * cosf(alpha);
                float lambda_x =   lambda * cosf(alpha);
                float lambda_z = - lambda * sinf(alpha);

                long long int lx_a     = floorf(lambda_x/dx);
                long long int lx_b     = lx_a+1;
                long long int lz_a     = floorf(lambda_z/dz);
                long long int lz_b     = lz_a+1;
                
                float plx_a     = (lx_b*dx  - lambda_x)/dx;
                float plx_b     = (lambda_x - lx_a*dx )/dx;
                float plz_a     = (lz_b*dz  - lambda_z)/dz;
                float plz_b     = (lambda_z - lz_a*dz )/dz;
                
                float ptot_x   = plx_a + plx_b;
                plx_a /= ptot_x;
                plx_b /= ptot_x;

                float ptot_z   = plz_a + plz_b;
                plz_a /= ptot_z;
                plz_b /= ptot_z;

                float p1   = plx_a * plz_a;
                float p2   = plx_b * plz_a;
                float p3   = plx_a * plz_b;
                float p4   = plx_b * plz_b;
                float ptot = p1 + p2 + p3 + p4;
                float cos2 = cosf(alpha)*cosf(alpha);
                p1 *= cos2/ptot;  
                p2 *= cos2/ptot;  
                p3 *= cos2/ptot;  
                p4 *= cos2/ptot;

                long long int i   = ikx + ikz*nx + (lx  +lx0)*nz*nx;
                long long int i_1 = ikx + ikz*nx + (lx_a+lx0)*nz*nx + (lz_a+lz0)*nlx*nz*nx;
                long long int i_2 = ikx + ikz*nx + (lx_b+lx0)*nz*nx + (lz_a+lz0)*nlx*nz*nx;
                long long int i_3 = ikx + ikz*nx + (lx_a+lx0)*nz*nx + (lz_b+lz0)*nlx*nz*nx;
                long long int i_4 = ikx + ikz*nx + (lx_b+lx0)*nz*nx + (lz_b+lz0)*nlx*nz*nx;

                if(-lx0<=lx_a  &&  lx_a<=+lx0  &&  -lz0<=lz_a  &&  lz_a<=+lz0)
                {
                    Gocig_fft[i_1][0] += p1*Hocig_fft[i][0];
                    Gocig_fft[i_1][1] += p1*Hocig_fft[i][1];
                }
                if(-lx0<=lx_b  &&  lx_b<=+lx0  &&  -lz0<=lz_a  &&  lz_a<=+lz0)
                {
                    Gocig_fft[i_2][0] += p2*Hocig_fft[i][0];
                    Gocig_fft[i_2][1] += p2*Hocig_fft[i][1];
                }
                if(-lx0<=lx_a  &&  lx_a<=+lx0  &&  -lz0<=lz_b  &&  lz_b<=+lz0)
                {
                    Gocig_fft[i_3][0] += p3*Hocig_fft[i][0];
                    Gocig_fft[i_3][1] += p3*Hocig_fft[i][1];
                }
                if(-lx0<=lx_b  &&  lx_b<=+lx0  &&  -lz0<=lz_b  &&  lz_b<=+lz0)
                {
                    Gocig_fft[i_4][0] += p4*Hocig_fft[i][0];
                    Gocig_fft[i_4][1] += p4*Hocig_fft[i][1];
                }
            }
    }
    // Applying dip correction to VOCIG
    /*
    for(lv=-lv0; lv<=lv0; lv++)
    {
        for(ikz=0; ikz<nz; ikz++)
            for(ikx=0; ikx<nx; ikx++)
            {
                float kx, kz;
                if(ikx<nx/2.0)    kx =       ikx  * dkx;
                else              kx = -((nx-ikx) * dkx);
                if(ikz<nz/2.0)    kz =       ikz  * dkz;
                else              kz = -((nz-ikz) * dkz);
                float alpha    = atan2f(kx,kz);

                float llv      =   lv*dz;
                float lambda   =   llv * sinf(alpha);
                float lambda_x =   lambda * cosf(alpha);
                float lambda_z = - lambda * sinf(alpha);

                long long int lx_a     = floorf(lambda_x/dx);
                long long int lx_b     = lx_a+1;
                long long int lz_a     = floorf(lambda_z/dz);
                long long int lz_b     = lz_a+1;
                
                float plx_a     = (lx_b*dx  - lambda_x)/dx;
                float plx_b     = (lambda_x - lx_a*dx )/dx;
                float plz_a     = (lz_b*dz  - lambda_z)/dz;
                float plz_b     = (lambda_z - lz_a*dz )/dz;

                float ptot_x   = plx_a + plx_b;
                plx_a /= ptot_x;
                plx_b /= ptot_x;

                float ptot_z   = plz_a + plz_b;
                plz_a /= ptot_z;
                plz_b /= ptot_z;

                float p1   = plx_a * plz_a;
                float p2   = plx_b * plz_a;
                float p3   = plx_a * plz_b;
                float p4   = plx_b * plz_b;
                float ptot = p1 + p2 + p3 + p4;
                float sin2 = sinf(alpha)*sinf(alpha);
                p1 *= sin2/ptot;  
                p2 *= sin2/ptot;  
                p3 *= sin2/ptot;  
                p4 *= sin2/ptot;

                long long int i   = ikx + ikz*nx + (lv  +lv0)*nz*nx;
                long long int i_1 = ikx + ikz*nx + (lx_a+lx0)*nz*nx + (lz_a+lz0)*nlx*nz*nx;
                long long int i_2 = ikx + ikz*nx + (lx_b+lx0)*nz*nx + (lz_a+lz0)*nlx*nz*nx;
                long long int i_3 = ikx + ikz*nx + (lx_a+lx0)*nz*nx + (lz_b+lz0)*nlx*nz*nx;
                long long int i_4 = ikx + ikz*nx + (lx_b+lx0)*nz*nx + (lz_b+lz0)*nlx*nz*nx;

                if(-lx0<=lx_a  &&  lx_a<=+lx0  &&  -lz0<=lz_a  &&  lz_a<=+lz0)
                {
                    Gocig_fft[i_1][0] += p1*Vocig_fft[i][0];
                    Gocig_fft[i_1][1] += p1*Vocig_fft[i][1];
                }
                if(-lx0<=lx_b  &&  lx_b<=+lx0  &&  -lz0<=lz_a  &&  lz_a<=+lz0)
                {
                    Gocig_fft[i_2][0] += p2*Vocig_fft[i][0];
                    Gocig_fft[i_2][1] += p2*Vocig_fft[i][1];
                }
                if(-lx0<=lx_a  &&  lx_a<=+lx0  &&  -lz0<=lz_b  &&  lz_b<=+lz0)
                {
                    Gocig_fft[i_3][0] += p3*Vocig_fft[i][0];
                    Gocig_fft[i_3][1] += p3*Vocig_fft[i][1];
                }
                if(-lx0<=lx_b  &&  lx_b<=+lx0  &&  -lz0<=lz_b  &&  lz_b<=+lz0)
                {
                    Gocig_fft[i_4][0] += p4*Vocig_fft[i][0];
                    Gocig_fft[i_4][1] += p4*Vocig_fft[i][1];
                }                
            }
    }*/

    // Transform back to space coordinates
    for(lz=-lz0; lz<=lz0; lz++)
        for(lx=-lx0; lx<=lx0; lx++)
        {
            size_t shift = (lx+lx0)*nx*nz + (lz+lz0)*nx*nz*nlx;
            applyFFTW2D_backwardI2R(Gocig_transp+shift, Gocig_fft+shift, nx, nz);
        }
    // float *Gocig_tmp = transp_dim4(Gocig_transp,  nx, nz, nlx, nlz, 13);
    // float *Gocig     = transp_dim4(Gocig_tmp,    nlx, nz,  nx, nlz, 24);
    // float *Gocig = transp_dim4(Gocig_transp,  nx, nz, nlx, nlz, 13);
    float *Gocig = transp_dim3(Gocig_transp,  nx, nz, nlx*nlz, 13);
    free(Gocig_transp);
    free(Hocig_transp);
    free(Vocig_transp);
    // freeMem(Gocig_transp, nx*nz*nlx*nlz*sizeof(float));
    // freeMem(Hocig_transp, nx*nz*nlx*sizeof(float));
    // freeMem(Vocig_transp, nx*nz*nlv*sizeof(float));
    // free(Gocig_tmp);
    fftw_free(Hocig_fft);
    fftw_free(Vocig_fft);
    fftw_free(Gocig_fft);
    // fftw_freeMem(Gocig_fft, nx*nz*nlx*nlz);
    // fftw_freeMem(Hocig_fft, nx*nz*nlx    );
    // fftw_freeMem(Vocig_fft, nx*nz*nlv    );

return Gocig;
}

void Gocig2VHocig(float *Gocig, float *Hocig, float *Vocig, long long int nx, long long int nz, float dx, float dz, \
                  long long int lx0, long long int lv0, long long int lz0, int iproc)
{
    long long int i, ix, iz, lx, lv, lz;
    long long int nlx, nlv, nlz;
    long long int ikx, ikz;

    nlx = 2*lx0 + 1;
    nlv = 2*lv0 + 1;
    nlz = 2*lz0 + 1;

    float dkx  = 1.0f/(nx * dx);
    float dkz  = 1.0f/(nz * dz);

    // Allocating arrays that will transformed CIGs
    float *Hocig_transp = CPU_zaloc1F(nx*nz*nlx);
    float *Vocig_transp = CPU_zaloc1F(nx*nz*nlv);
    fftw_complex *Hocig_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*nz*nlx);
    fftw_complex *Vocig_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*nz*nlv);
    fftw_complex *Gocig_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*nz*nlx*nlz);

    float *Gocig_transp = transp_dim3(Gocig, nlx*nlz, nz, nx, 13);

    // Transform back to space coordinates
    for(lz=-lz0; lz<=lz0; lz++)
        for(lx=-lx0; lx<=lx0; lx++)
        {
            size_t shift = (lx+lx0)*nx*nz + (lz+lz0)*nx*nz*nlx;
            applyFFTW2D_forwardR2I(Gocig_fft+shift, Gocig_transp+shift, nx, nz);
        }
    free(Gocig_transp);
    // freeMem(Gocig_transp, nlx*nlz*nz*nx);    

    // Applying dip correction to recover HOCIG
    for(lx=-lx0; lx<=lx0; lx++)
    {
        for(ikz=0; ikz<nz; ikz++)
            for(ikx=0; ikx<nx; ikx++)
            {
                float kx, kz;
                if(ikx<nx/2.0)    kx =       ikx  * dkx;
                else              kx = -((nx-ikx) * dkx);
                if(ikz<nz/2.0)    kz =       ikz  * dkz;
                else              kz = -((nz-ikz) * dkz);

                float alpha    = atan2f(kx,kz);

                float llx      =   lx*dx;
                float lambda   =   llx * cosf(alpha);
                float lambda_x =   lambda * cosf(alpha);
                float lambda_z = - lambda * sinf(alpha);

                long long int lx_a     = floorf(lambda_x/dx);
                long long int lx_b     = lx_a+1;
                long long int lz_a     = floorf(lambda_z/dz);
                long long int lz_b     = lz_a+1;
                
                float plx_a     = (lx_b*dx  - lambda_x)/dx;
                float plx_b     = (lambda_x - lx_a*dx )/dx;
                float plz_a     = (lz_b*dz  - lambda_z)/dz;
                float plz_b     = (lambda_z - lz_a*dz )/dz;
                
                float ptot_x   = plx_a + plx_b;
                plx_a /= ptot_x;
                plx_b /= ptot_x;

                float ptot_z   = plz_a + plz_b;
                plz_a /= ptot_z;
                plz_b /= ptot_z;

                float p1   = plx_a * plz_a;
                float p2   = plx_b * plz_a;
                float p3   = plx_a * plz_b;
                float p4   = plx_b * plz_b;
                float ptot = p1 + p2 + p3 + p4;
                float cos2 = cosf(alpha)*cosf(alpha);
                p1 *= cos2/ptot;  
                p2 *= cos2/ptot;  
                p3 *= cos2/ptot;  
                p4 *= cos2/ptot;
                // p1 /= ptot;  
                // p2 /= ptot;  
                // p3 /= ptot;  
                // p4 /= ptot;

                if(isfinite(ptot)==0  ||  ptot_x==0  ||  ptot_z==0  ||  isfinite(p1)==0  ||  isfinite(p2)==0  ||  isfinite(p3)==0  ||  isfinite(p4)==0)
                {
                    printf("\n\n  ATENTION (HOCIG): ptot=%f  ptot_x=%f  ptot_z=%f  p1=%f  p2=%f  p3=%f  p4=%f \n\n", ptot, ptot_x, ptot_z, p1, p2, p3, p4);
                }

                long long int i   = ikx + ikz*nx + (lx  +lx0)*nz*nx;
                long long int i_1 = ikx + ikz*nx + (lx_a+lx0)*nz*nx + (lz_a+lz0)*nlx*nz*nx;
                long long int i_2 = ikx + ikz*nx + (lx_b+lx0)*nz*nx + (lz_a+lz0)*nlx*nz*nx;
                long long int i_3 = ikx + ikz*nx + (lx_a+lx0)*nz*nx + (lz_b+lz0)*nlx*nz*nx;
                long long int i_4 = ikx + ikz*nx + (lx_b+lx0)*nz*nx + (lz_b+lz0)*nlx*nz*nx;

                if(-lx0<=lx_a  &&  lx_a<=+lx0  &&  -lz0<=lz_a  &&  lz_a<=+lz0)
                {
                    Hocig_fft[i][0] += p1*Gocig_fft[i_1][0];
                    Hocig_fft[i][1] += p1*Gocig_fft[i_1][1];
                }
                if(-lx0<=lx_b  &&  lx_b<=+lx0  &&  -lz0<=lz_a  &&  lz_a<=+lz0)
                {
                    Hocig_fft[i][0] += p2*Gocig_fft[i_2][0];
                    Hocig_fft[i][1] += p2*Gocig_fft[i_2][1];
                }
                if(-lx0<=lx_a  &&  lx_a<=+lx0  &&  -lz0<=lz_b  &&  lz_b<=+lz0)
                {
                    Hocig_fft[i][0] += p3*Gocig_fft[i_3][0];
                    Hocig_fft[i][1] += p3*Gocig_fft[i_3][1];
                }
                if(-lx0<=lx_b  &&  lx_b<=+lx0  &&  -lz0<=lz_b  &&  lz_b<=+lz0)
                {
                    Hocig_fft[i][0] += p4*Gocig_fft[i_4][0];
                    Hocig_fft[i][1] += p4*Gocig_fft[i_4][1];
                }
            }
    }
    //*
    // Applying dip correction to VOCIG
    for(lv=-lv0; lv<=lv0; lv++)
    {
        for(ikz=0; ikz<nz; ikz++)
            for(ikx=0; ikx<nx; ikx++)
            {
                float kx, kz;
                if(ikx<nx/2.0)    kx =       ikx  * dkx;
                else              kx = -((nx-ikx) * dkx);
                if(ikz<nz/2.0)    kz =       ikz  * dkz;
                else              kz = -((nz-ikz) * dkz);
                float alpha    = atan2f(kx,kz);

                float llv      =   lv*dz;
                float lambda   =   llv * sinf(alpha);
                float lambda_x =   lambda * cosf(alpha);
                float lambda_z = - lambda * sinf(alpha);

                long long int lx_a = floorf(lambda_x/dx);
                long long int lx_b = lx_a+1;
                long long int lz_a = floorf(lambda_z/dz);
                long long int lz_b = lz_a+1;
                
                float plx_a = (lx_b*dx  - lambda_x)/dx;
                float plx_b = (lambda_x - lx_a*dx )/dx;
                float plz_a = (lz_b*dz  - lambda_z)/dz;
                float plz_b = (lambda_z - lz_a*dz )/dz;
                

                float ptot_x = plx_a + plx_b;
                plx_a /= ptot_x;
                plx_b /= ptot_x;

                float ptot_z = plz_a + plz_b;
                plz_a /= ptot_z;
                plz_b /= ptot_z;

                float p1   = plx_a * plz_a;
                float p2   = plx_b * plz_a;
                float p3   = plx_a * plz_b;
                float p4   = plx_b * plz_b;
                float ptot = p1 + p2 + p3 + p4;
                float sin2 = sinf(alpha)*sinf(alpha);
                p1 *= sin2/ptot;
                p2 *= sin2/ptot;
                p3 *= sin2/ptot;
                p4 *= sin2/ptot;
                // p1 /= ptot;  
                // p2 /= ptot;  
                // p3 /= ptot;  
                // p4 /= ptot;

                if(isfinite(ptot)==0  ||  ptot_x==0  ||  ptot_z==0  ||  isfinite(p1)==0  ||  isfinite(p2)==0  ||  isfinite(p3)==0  ||  isfinite(p4)==0)
                {
                    printf("\n\n  ATENTION (VOCIG): ptot=%f  ptot_x=%f  ptot_z=%f  p1=%f  p2=%f  p3=%f  p4=%f \n\n", ptot, ptot_x, ptot_z, p1, p2, p3, p4);
                }

                long long int i   = ikx + ikz*nx + (lv  +lv0)*nz*nx;
                long long int i_1 = ikx + ikz*nx + (lx_a+lx0)*nz*nx + (lz_a+lz0)*nlx*nz*nx;
                long long int i_2 = ikx + ikz*nx + (lx_b+lx0)*nz*nx + (lz_a+lz0)*nlx*nz*nx;
                long long int i_3 = ikx + ikz*nx + (lx_a+lx0)*nz*nx + (lz_b+lz0)*nlx*nz*nx;
                long long int i_4 = ikx + ikz*nx + (lx_b+lx0)*nz*nx + (lz_b+lz0)*nlx*nz*nx;

                if(-lx0<=lx_a  &&  lx_a<=+lx0  &&  -lz0<=lz_a  &&  lz_a<=+lz0)
                {
                    Vocig_fft[i][0] += 0*p1*Gocig_fft[i_1][0];
                    Vocig_fft[i][1] += 0*p1*Gocig_fft[i_1][1];
                }
                if(-lx0<=lx_b  &&  lx_b<=+lx0  &&  -lz0<=lz_a  &&  lz_a<=+lz0)
                {
                    Vocig_fft[i][0] += 0*p2*Gocig_fft[i_2][0];
                    Vocig_fft[i][1] += 0*p2*Gocig_fft[i_2][1];
                }
                if(-lx0<=lx_a  &&  lx_a<=+lx0  &&  -lz0<=lz_b  &&  lz_b<=+lz0)
                {
                    Vocig_fft[i][0] += 0*p3*Gocig_fft[i_3][0];
                    Vocig_fft[i][1] += 0*p3*Gocig_fft[i_3][1];
                }
                if(-lx0<=lx_b  &&  lx_b<=+lx0  &&  -lz0<=lz_b  &&  lz_b<=+lz0)
                {
                    Vocig_fft[i][0] += 0*p4*Gocig_fft[i_4][0];
                    Vocig_fft[i][1] += 0*p4*Gocig_fft[i_4][1];
                }

            }
    }
    //*/
    fftw_free(Gocig_fft);
    // fftw_freeMem(Gocig_fft, nlx*nlz*nz*nx); printf("\n Gocig2VHocig 12 iproc=%d", iproc); printf("\n");

    // Transforming to kx-kz
    for(lx=-lx0; lx<=lx0; lx++) {
        size_t shift = (lx+lx0)*nx*nz;
        applyFFTW2D_backwardI2R(Hocig_transp+shift, Hocig_fft+shift, nx, nz);
    }
    for(lv=-lv0; lv<=lv0; lv++) {
        size_t shift = (lv+lv0)*nx*nz;
        applyFFTW2D_backwardI2R(Vocig_transp+shift, Vocig_fft+shift, nx, nz);
    }
    fftw_free(Hocig_fft);
    fftw_free(Vocig_fft);
    
    // Transposing to make offset the slowest dimension
    float *Hocig_tmp = transp_dim3(Hocig_transp, nx, nz, nlx, 13);
    float *Vocig_tmp = transp_dim3(Vocig_transp, nx, nz, nlv, 13);
    free(Hocig_transp);
    free(Vocig_transp);

    memcpy(Hocig, Hocig_tmp, nx*nz*nlx*sizeof(float));
    memcpy(Vocig, Vocig_tmp, nx*nz*nlv*sizeof(float));
    free(Hocig_tmp);
    free(Vocig_tmp);


return;
}
void VHocig2Docig_testField(float *Docig, float *Hocig, float *Vocig, \
                  long long int nx, long long int nz, long long int ndip, \
                  float dx, float dz, float ddip, float dipMin, \
                  long long int lx0, long long int lv0, long long int ld0, int iproc)
{
 
    dipMin *= M_PI/180.0;
    ddip   *= M_PI/180.0;

    size_t szw = sizeof(fftw_complex);

    long long int i, i_a, i_b, ix, iz, lx, lv, ld, idip;
    long long int nlx, nlv, nld;
    long long int ikx, ikz;
    float p_dip_a, p_dip_b, p_dip_tot;

    nlx = 2*lx0 + 1;
    nlv = 2*lv0 + 1;
    nld = 2*ld0 + 1;

    float dkx  = 1.0f/(nx * dx);
    float dkz  = 1.0f/(nz * dz);

    // Allocating arrays that will transformed CIGs
    fftw_complex *Hocig_fft = CPU_zalocFFTW1(nx*nz*nlx);
    fftw_complex *Vocig_fft = CPU_zalocFFTW1(nx*nz*nlv);
    fftw_complex *Docig_fft = CPU_zalocFFTW1(nx*nz*nld*ndip);
    float *Docig_transp     = CPU_zaloc1F(nx*nz*nld*ndip);

    // Transposing to make offset the slowest dimension
    float *Hocig_transp = transp_dim3(Hocig, nlx, nz, nx, 13);
    float *Vocig_transp = transp_dim3(Vocig, nlv, nz, nx, 13);

    // Transforming to kx-kz
    for(lx=-lx0; lx<=lx0; lx++)
    {
        size_t shift = (lx+lx0)*nx*nz;
        fftw_complex *Hocig_fft_tmp = CPU_zalocFFTW1(nx*nz);
        applyFFTW2D_forwardR2I(Hocig_fft_tmp, Hocig_transp+shift, nx, nz); //wayhome
        for(i=0; i<nx*nz; i++)
        {
            Hocig_fft[shift+i][0] = Hocig_fft_tmp[i][0];
            Hocig_fft[shift+i][1] = Hocig_fft_tmp[i][1];
        }
        free(Hocig_fft_tmp);
    }
    for(lv=-lv0; lv<=lv0; lv++)
    {
        size_t shift = (lv+lv0)*nx*nz;
        fftw_complex *Vocig_fft_tmp = CPU_zalocFFTW1(nx*nz);
        applyFFTW2D_forwardR2I(Vocig_fft_tmp, Vocig_transp+shift, nx, nz);
        for(i=0; i<nx*nz; i++)
        {
            Vocig_fft[shift+i][0] = Vocig_fft_tmp[i][0];
            Vocig_fft[shift+i][1] = Vocig_fft_tmp[i][1];
        }
        free(Vocig_fft_tmp);
    }

    // Applying dip correction to HOCIG
    for(lx=-lx0; lx<=lx0; lx++)
    {
        for(ikz=0; ikz<nz; ikz++)
            for(ikx=0; ikx<nx; ikx++)
            {
                i   = ikx + ikz*nx + (lx+lx0)*nz*nx;
                Docig_fft[i][0] =  Hocig_fft[i][0];
                Docig_fft[i][1] =  Hocig_fft[i][1];
            }
    }

    // Transform back to space coordinates
    for(lx=-lx0; lx<=lx0; lx++)
    {
        size_t shift = (lx+lx0)*nx*nz;
        fftw_complex *Docig_fft_tmp = CPU_zalocFFTW1(nx*nz);
        for(i=0; i<nx*nz; i++)
        {
            Docig_fft_tmp[i][0] = Docig_fft[shift+i][0];
            Docig_fft_tmp[i][1] = Docig_fft[shift+i][1];
        }
        applyFFTW2D_backwardI2R(Docig_transp+shift, Docig_fft_tmp, nx, nz); // gethere
        free(Docig_fft_tmp);
    }
    float *Docig_tmp = transp_dim3(Docig_transp, nx, nz, nlx, 13);
    memcpy(Hocig, Docig_tmp, nx*nz*nlx*sizeof(float));
    free(Docig_tmp);
    // for(idip=0; idip<ndip; idip++)
    //     for(ld=-ld0; ld<=ld0; ld++)
    //     {
    //         size_t shift = idip*nld*nx*nz + (ld+ld0)*nx*nz;
    //         fftw_complex *Docig_fft_tmp = CPU_zalocFFTW1(nx*nz);
    //         for(i=0; i<nx*nz; i++)
    //         {
    //             Docig_fft_tmp[i][0] = Docig_fft[shift+i][0];
    //             Docig_fft_tmp[i][1] = Docig_fft[shift+i][1];
    //         }
    //         applyFFTW2D_backwardI2R(Docig_transp+shift, Docig_fft_tmp, nx, nz); // gethere
    //         free(Docig_fft_tmp);
    //     }

    // for(idip=0; idip<ndip; idip++)
    // {
    //     size_t shift = idip*nld*nx*nz;
    //     float *Docig_tmp = transp_dim3(Docig_transp+shift, nx, nz, nld, 13);
    //     memcpy(Docig+shift, Docig_tmp, nx*nz*nld*sizeof(float));
    //     free(Docig_tmp);
    // }

    free(Docig_transp);
    free(Hocig_transp);
    free(Vocig_transp);

    fftw_free(Hocig_fft);
    fftw_free(Vocig_fft);
    fftw_free(Docig_fft);

return;
}
void VHocig2Docig(float *Docig, float *Hocig, float *Vocig, \
                  long long int nx, long long int nz, long long int ndip, \
                  float dx, float dz, float ddip, float dipMin, \
                  long long int lx0, long long int lv0, long long int ld0, int iproc)
{
 
    dipMin *= M_PI/180.0;
    ddip   *= M_PI/180.0;

    size_t szw = sizeof(fftw_complex);

    long long int i, i_a, i_b, ix, iz, lx, lv, ld, idip;
    long long int nlx, nlv, nld;
    long long int ikx, ikz;
    float p_dip_a, p_dip_b, p_dip_tot;

    nlx = 2*lx0 + 1;
    nlv = 2*lv0 + 1;
    nld = 2*ld0 + 1;

    float dkx  = 1.0f/(nx * dx);
    float dkz  = 1.0f/(nz * dz);

    // Allocating arrays that will transformed CIGs
    fftw_complex *Hocig_fft = CPU_zalocFFTW1(nx*nz*nlx);
    fftw_complex *Vocig_fft = CPU_zalocFFTW1(nx*nz*nlv);
    fftw_complex *Docig_fft = CPU_zalocFFTW1(nx*nz*nld*ndip);
    float *Docig_transp     = CPU_zaloc1F(nx*nz*nld*ndip);

    // Transposing to make offset the slowest dimension
    float *Hocig_transp = transp_dim3(Hocig, nlx, nz, nx, 13);
    float *Vocig_transp = transp_dim3(Vocig, nlv, nz, nx, 13);

    // Transforming to kx-kz
    for(lx=-lx0; lx<=lx0; lx++)
    {
        size_t shift = (lx+lx0)*nx*nz;
        fftw_complex *Hocig_fft_tmp = CPU_zalocFFTW1(nx*nz);
        applyFFTW2D_forwardR2I(Hocig_fft_tmp, Hocig_transp+shift, nx, nz); //homeishere
        for(i=0; i<nx*nz; i++)
        {
            Hocig_fft[shift+i][0] = Hocig_fft_tmp[i][0];
            Hocig_fft[shift+i][1] = Hocig_fft_tmp[i][1];
        }
        free(Hocig_fft_tmp);
    }
    for(lv=-lv0; lv<=lv0; lv++)
    {
        size_t shift = (lv+lv0)*nx*nz;
        fftw_complex *Vocig_fft_tmp = CPU_zalocFFTW1(nx*nz);
        applyFFTW2D_forwardR2I(Vocig_fft_tmp, Vocig_transp+shift, nx, nz);
        for(i=0; i<nx*nz; i++)
        {
            Vocig_fft[shift+i][0] = Vocig_fft_tmp[i][0];
            Vocig_fft[shift+i][1] = Vocig_fft_tmp[i][1];
        }
        free(Vocig_fft_tmp);
    }

    // Applying dip correction to HOCIG
    //*
    for(lx=-lx0; lx<=lx0; lx++)
    {
        for(ikz=0; ikz<nz; ikz++)
            for(ikx=0; ikx<nx; ikx++)
            {
                float kx, kz;
                if(ikx<nx/2.0)    kx =       ikx  * dkx;
                else              kx = -((nx-ikx) * dkx);
                if(ikz<nz/2.0)    kz =       ikz  * dkz;
                else              kz = -((nz-ikz) * dkz);
                float alpha = atan2f(kx,kz);

                // idip = (long long int) (0.5f + (alpha-dipMin)/ddip);
                idip = (long long int) floorf((alpha-dipMin)/ddip);

                if(idip>=ndip)
                {
                    idip = ndip-1;
                    float dipMax = dipMin + ddip * (ndip-1);
                    float sigma = logf(5.0)/ddip;
                    p_dip_a = exp(-sigma*(alpha-dipMax));
                    p_dip_b = 0.0;
                }
                else if(idip<0)
                {
                    idip = 0;
                    float sigma = logf(5.0)/ddip;
                    p_dip_a = exp(-sigma*(dipMin-alpha));
                    p_dip_b = 0.0;
                }
                else
                {
                    float dip_a = dipMin +  idip   *ddip;
                    float dip_b = dipMin + (idip+1)*ddip;
                    p_dip_a = (dip_b - alpha)/ddip;
                    p_dip_b = (alpha - dip_a)/ddip;
                    p_dip_tot = p_dip_a + p_dip_b;
                    p_dip_a /= p_dip_tot;
                    p_dip_b /= p_dip_tot;
                }

                float llx       =   lx*dx;
                float lambda    =   llx * cosf(alpha);
                // float lambda    =   llx * cosf(alpha)/(0.001+cosf(alpha-dip));

                long long int ld_a = floorf(lambda/dx);
                long long int ld_b = ld_a+1;
                
                float pld_a  = (ld_b*dx - lambda)/dx;
                float pld_b  = (lambda  - ld_a*dx)/dx;
                float ptot_d = pld_a + pld_b;
                float cos2   = cosf(alpha)*cosf(alpha);
                pld_a *= cos2/ptot_d;
                pld_b *= cos2/ptot_d;

                i   = ikx + ikz*nx + (lx  +lx0)*nz*nx;
                i_a = ikx + ikz*nx + (ld_a+ld0)*nz*nx + idip*nld*nz*nx;
                i_b = ikx + ikz*nx + (ld_b+ld0)*nz*nx + idip*nld*nz*nx;
                if(-ld0<=ld_a  &&  ld_a<=+ld0)
                {
                    Docig_fft[i_a][0] += pld_a * p_dip_a * Hocig_fft[i][0];
                    Docig_fft[i_a][1] += pld_a * p_dip_a * Hocig_fft[i][1];
                }
                if(-ld0<=ld_b  &&  ld_b<=+ld0)
                {
                    Docig_fft[i_b][0] += pld_b * p_dip_a * Hocig_fft[i][0];
                    Docig_fft[i_b][1] += pld_b * p_dip_a * Hocig_fft[i][1];
                }
                if(idip<ndip-1)
                {
                i   = ikx + ikz*nx + (lx  +lx0)*nz*nx;
                i_a = ikx + ikz*nx + (ld_a+ld0)*nz*nx + (idip+1)*nld*nz*nx;
                i_b = ikx + ikz*nx + (ld_b+ld0)*nz*nx + (idip+1)*nld*nz*nx;
                if(-ld0<=ld_a  &&  ld_a<=+ld0)
                {
                    Docig_fft[i_a][0] += pld_a * p_dip_b * Hocig_fft[i][0];
                    Docig_fft[i_a][1] += pld_a * p_dip_b * Hocig_fft[i][1];
                }
                if(-ld0<=ld_b  &&  ld_b<=+ld0)
                {
                    Docig_fft[i_b][0] += pld_b * p_dip_b * Hocig_fft[i][0];
                    Docig_fft[i_b][1] += pld_b * p_dip_b * Hocig_fft[i][1];
                }
                }

            }
    }
    //*/
    // Applying dip correction to VOCIG
    //*
    for(lv=-lv0; lv<=lv0; lv++)
    {
        for(ikz=0; ikz<nz; ikz++)
            for(ikx=0; ikx<nx; ikx++)
            {
                float kx, kz;
                if(ikx<nx/2.0)    kx =       ikx  * dkx;
                else              kx = -((nx-ikx) * dkx);
                if(ikz<nz/2.0)    kz =       ikz  * dkz;
                else              kz = -((nz-ikz) * dkz);
                float alpha    = atan2f(kx,kz);

                idip = (long long int) (0.5f + (alpha-dipMin)/ddip);

                if(idip>=ndip)
                {
                    idip = ndip-1;
                    float dipMax = dipMin + ddip * (ndip-1);
                    float sigma = logf(5.0)/ddip;
                    p_dip_a = exp(-sigma*(alpha-dipMax));
                    p_dip_b = 0.0;
                }
                else if(idip<0)
                {
                    idip = 0;
                    float sigma = logf(5.0)/ddip;
                    p_dip_a = exp(-sigma*(dipMin-alpha));
                    p_dip_b = 0.0;
                }
                else
                {
                    float dip_a = dipMin +  idip   *ddip;
                    float dip_b = dipMin + (idip+1)*ddip;
                    p_dip_a = (dip_b - alpha)/ddip;
                    p_dip_b = (alpha - dip_a)/ddip;
                    p_dip_tot = p_dip_a + p_dip_b;
                    p_dip_a /= p_dip_tot;
                    p_dip_b /= p_dip_tot;
                }

                float llv      =   lv*dz;
                float lambda   =   llv * sinf(alpha);
                // float lambda   =   llv * sinf(alpha)/(0.001+sinf(alpha-dip));

                long long int ld_a = floorf(lambda/dx);
                long long int ld_b = ld_a+1;
                
                float pld_a  = (ld_b*dx - lambda)/dx;
                float pld_b  = (lambda  - ld_a*dx)/dx;
                float ptot_d = pld_a + pld_b;
                float sin2   = sinf(alpha)*sinf(alpha);
                pld_a *= sin2/ptot_d;
                pld_b *= sin2/ptot_d;

                i   = ikx + ikz*nx + (lv  +lv0)*nz*nx;
                i_a = ikx + ikz*nx + (ld_a+ld0)*nz*nx + idip*nld*nz*nx;
                i_b = ikx + ikz*nx + (ld_b+ld0)*nz*nx + idip*nld*nz*nx;
                if(-ld0<=ld_a  &&  ld_a<=+ld0)
                {
                    Docig_fft[i_a][0] += pld_a * p_dip_a * Vocig_fft[i][0];
                    Docig_fft[i_a][1] += pld_a * p_dip_a * Vocig_fft[i][1];
                }
                if(-ld0<=ld_b  &&  ld_b<=+ld0)
                {
                    Docig_fft[i_b][0] += pld_b * p_dip_a * Vocig_fft[i][0];
                    Docig_fft[i_b][1] += pld_b * p_dip_a * Vocig_fft[i][1];
                }
                
                if(idip<ndip-1)
                {
                i   = ikx + ikz*nx + (lv  +lv0)*nz*nx;
                i_a = ikx + ikz*nx + (ld_a+ld0)*nz*nx + (idip+1)*nld*nz*nx;
                i_b = ikx + ikz*nx + (ld_b+ld0)*nz*nx + (idip+1)*nld*nz*nx;
                if(-ld0<=ld_a  &&  ld_a<=+ld0)
                {
                    Docig_fft[i_a][0] += pld_a * p_dip_b * Vocig_fft[i][0];
                    Docig_fft[i_a][1] += pld_a * p_dip_b * Vocig_fft[i][1];
                }
                if(-ld0<=ld_b  &&  ld_b<=+ld0)
                {
                    Docig_fft[i_b][0] += pld_b * p_dip_b * Vocig_fft[i][0];
                    Docig_fft[i_b][1] += pld_b * p_dip_b * Vocig_fft[i][1];
                }
                }

                
            }
    }//*/

    // Transform back to space coordinates
    for(idip=0; idip<ndip; idip++)
        for(ld=-ld0; ld<=ld0; ld++)
        {
            size_t shift = idip*nld*nx*nz + (ld+ld0)*nx*nz;
            fftw_complex *Docig_fft_tmp = CPU_zalocFFTW1(nx*nz);
            for(i=0; i<nx*nz; i++)
            {
                Docig_fft_tmp[i][0] = Docig_fft[shift+i][0];
                Docig_fft_tmp[i][1] = Docig_fft[shift+i][1];
            }
            applyFFTW2D_backwardI2R(Docig_transp+shift, Docig_fft_tmp, nx, nz); // gethere
            free(Docig_fft_tmp);
        }

    for(idip=0; idip<ndip; idip++)
    {
        size_t shift = idip*nld*nx*nz;
        float *Docig_tmp = transp_dim3(Docig_transp+shift, nx, nz, nld, 13);
        memcpy(Docig+shift, Docig_tmp, nx*nz*nld*sizeof(float));
        free(Docig_tmp);
    }

    free(Docig_transp);
    free(Hocig_transp);
    free(Vocig_transp);

    fftw_free(Hocig_fft);
    fftw_free(Vocig_fft);
    fftw_free(Docig_fft);

return;
}
void Docig2VHocig(float *Docig, float *Hocig, float *Vocig, \
                  long long int nx, long long int nz, long long int ndip, \
                  float dx, float dz, float ddip, float dipMin, \
                  long long int lx0, long long int lv0, long long int ld0, int iproc)
{
    dipMin *= M_PI/180.0;
    ddip   *= M_PI/180.0;

    size_t szw = sizeof(fftw_complex);

    long long int i, i_a, i_b, ix, iz, lx, lv, ld, idip;
    long long int nlx, nlv, nld;
    long long int ikx, ikz;
    float p_dip_a, p_dip_b, p_dip_tot;

    nlx = 2*lx0 + 1;
    nlv = 2*lv0 + 1;
    nld = 2*ld0 + 1;

    float dkx  = 1.0f/(nx * dx);
    float dkz  = 1.0f/(nz * dz);

    // Allocating arrays that will transformed CIGs
    fftw_complex *Hocig_fft = CPU_zalocFFTW1(nx*nz*nlx);
    fftw_complex *Vocig_fft = CPU_zalocFFTW1(nx*nz*nlv);
    fftw_complex *Docig_fft = CPU_zalocFFTW1(nx*nz*nld*ndip);
    
    float *Docig_transp     = CPU_zaloc1F(nx*nz*nld*ndip);
    float *Hocig_transp     = CPU_zaloc1F(nx*nz*nlx);
    float *Vocig_transp     = CPU_zaloc1F(nx*nz*nlv);


    for(idip=0; idip<ndip; idip++)
    {
        size_t shift = idip*nld*nx*nz;
        float *Docig_tmp = transp_dim3(Docig+shift, nld, nz, nx, 13);
        memcpy(Docig_transp+shift, Docig_tmp, nx*nz*nld*sizeof(float));
        free(Docig_tmp);
    }
    for(idip=0; idip<ndip; idip++)
        for(ld=-ld0; ld<=ld0; ld++)
        {
            size_t shift = idip*nld*nx*nz + (ld+ld0)*nx*nz;
            fftw_complex *Docig_fft_tmp = CPU_zalocFFTW1(nx*nz);
            applyFFTW2D_forwardR2I(Docig_fft_tmp, Docig_transp+shift, nx, nz); 
            for(i=0; i<nx*nz; i++) {
                Docig_fft[shift+i][0] = Docig_fft_tmp[i][0];
                Docig_fft[shift+i][1] = Docig_fft_tmp[i][1];
            }
            free(Docig_fft_tmp);
        }

    // Applying dip correction to HOCIG
    for(lx=-lx0; lx<=lx0; lx++)
    {
        for(ikz=0; ikz<nz; ikz++)
            for(ikx=0; ikx<nx; ikx++)
            {
                float kx, kz;
                if(ikx<nx/2.0)    kx =       ikx  * dkx;
                else              kx = -((nx-ikx) * dkx);
                if(ikz<nz/2.0)    kz =       ikz  * dkz;
                else              kz = -((nz-ikz) * dkz);
                float alpha = atan2f(kx,kz);

                // idip = (long long int) (0.5f + (alpha-dipMin)/ddip);
                idip = (long long int) floorf((alpha-dipMin)/ddip);

                if(idip>=ndip)
                {
                    idip = ndip-1;
                    float dipMax = dipMin + ddip * idip;
                    float sigma = logf(5.0)/ddip;
                    p_dip_a = exp(-sigma*(alpha-dipMax));
                    p_dip_b = 0.0;
                }
                else if(idip<0)
                {
                    idip = 0;
                    float sigma = logf(5.0)/ddip;
                    p_dip_a = exp(-sigma*(dipMin-alpha));
                    p_dip_b = 0.0;
                }
                else
                {
                    float dip_a = dipMin +  idip   *ddip;
                    float dip_b = dipMin + (idip+1)*ddip;
                    p_dip_a = (dip_b - alpha)/ddip;
                    p_dip_b = (alpha - dip_a)/ddip;
                    p_dip_tot = p_dip_a + p_dip_b;
                    p_dip_a /= p_dip_tot;
                    p_dip_b /= p_dip_tot;
                }

                float llx       =   lx*dx;
                float lambda    =   llx * cosf(alpha);
                // float lambda    =   llx * cosf(alpha)/(0.001+cosf(alpha-dip));

                long long int ld_a = floorf(lambda/dx);
                long long int ld_b = ld_a+1;
                
                float pld_a  = (ld_b*dx - lambda)/dx;
                float pld_b  = (lambda  - ld_a*dx)/dx;
                float ptot_d = pld_a + pld_b;
                float cos2   = cosf(alpha)*cosf(alpha);
                pld_a *= cos2/ptot_d;
                pld_b *= cos2/ptot_d;

                i   = ikx + ikz*nx + (lx  +lx0)*nz*nx;
                i_a = ikx + ikz*nx + (ld_a+ld0)*nz*nx + idip*nld*nz*nx;
                i_b = ikx + ikz*nx + (ld_b+ld0)*nz*nx + idip*nld*nz*nx;
                if(-ld0<=ld_a  &&  ld_a<=+ld0)
                {
                    Hocig_fft[i][0] += pld_a * p_dip_a * Docig_fft[i_a][0];
                    Hocig_fft[i][1] += pld_a * p_dip_a * Docig_fft[i_a][1];
                }
                if(-ld0<=ld_b  &&  ld_b<=+ld0)
                {
                    Hocig_fft[i][0] += pld_b * p_dip_a * Docig_fft[i_b][0];
                    Hocig_fft[i][1] += pld_b * p_dip_a * Docig_fft[i_b][1];
                }
                
                if(idip<ndip-1)
                {
                i   = ikx + ikz*nx + (lx  +lx0)*nz*nx;
                i_a = ikx + ikz*nx + (ld_a+ld0)*nz*nx + (idip+1)*nld*nz*nx;
                i_b = ikx + ikz*nx + (ld_b+ld0)*nz*nx + (idip+1)*nld*nz*nx;
                if(-ld0<=ld_a  &&  ld_a<=+ld0)
                {
                    Hocig_fft[i][0] += pld_a * p_dip_b * Docig_fft[i_a][0];
                    Hocig_fft[i][1] += pld_a * p_dip_b * Docig_fft[i_a][1];
                }
                if(-ld0<=ld_b  &&  ld_b<=+ld0)
                {
                    Hocig_fft[i][0] += pld_b * p_dip_b * Docig_fft[i_b][0];
                    Hocig_fft[i][1] += pld_b * p_dip_b * Docig_fft[i_b][1];
                }
                }
            }
    }
    // Applying dip correction to VOCIG
    //*
    for(lv=-lv0; lv<=lv0; lv++)
    {
        for(ikz=0; ikz<nz; ikz++)
            for(ikx=0; ikx<nx; ikx++)
            {
                float kx, kz;
                if(ikx<nx/2.0)    kx =       ikx  * dkx;
                else              kx = -((nx-ikx) * dkx);
                if(ikz<nz/2.0)    kz =       ikz  * dkz;
                else              kz = -((nz-ikz) * dkz);
                float alpha    = atan2f(kx,kz);

                idip = (long long int) (0.5f + (alpha-dipMin)/ddip);

                if(idip>=ndip)
                {
                    idip = ndip-1;
                    float dipMax = dipMin + ddip * (ndip-1);
                    float sigma = logf(5.0)/ddip;
                    p_dip_a = exp(-sigma*(alpha-dipMax));
                    p_dip_b = 0.0;
                }
                else if(idip<0)
                {
                    idip = 0;
                    float sigma = logf(5.0)/ddip;
                    p_dip_a = exp(-sigma*(dipMin-alpha));
                    p_dip_b = 0.0;
                }
                else
                {
                    float dip_a = dipMin +  idip   *ddip;
                    float dip_b = dipMin + (idip+1)*ddip;
                    p_dip_a = (dip_b - alpha)/ddip;
                    p_dip_b = (alpha - dip_a)/ddip;
                    p_dip_tot = p_dip_a + p_dip_b;
                    p_dip_a /= p_dip_tot;
                    p_dip_b /= p_dip_tot;
                }

                float llv      =   lv*dz;
                float lambda   =   llv * sinf(alpha);
                // float lambda   =   llv * sinf(alpha)/(0.001+sinf(alpha-dip));

                long long int ld_a = floorf(lambda/dx);
                long long int ld_b = ld_a+1;
                
                float pld_a  = (ld_b*dx - lambda)/dx;
                float pld_b  = (lambda  - ld_a*dx)/dx;
                float ptot_d = pld_a + pld_b;
                float sin2   = sinf(alpha)*sinf(alpha);
                pld_a *= sin2/ptot_d;
                pld_b *= sin2/ptot_d;

                //*
                i   = ikx + ikz*nx + (lv  +lv0)*nz*nx;
                i_a = ikx + ikz*nx + (ld_a+ld0)*nz*nx + idip*nld*nz*nx;
                i_b = ikx + ikz*nx + (ld_b+ld0)*nz*nx + idip*nld*nz*nx;
                if(-ld0<=ld_a  &&  ld_a<=+ld0)
                {
                    Vocig_fft[i][0] += pld_a * p_dip_a * Docig_fft[i_a][0];
                    Vocig_fft[i][1] += pld_a * p_dip_a * Docig_fft[i_a][1];
                }
                if(-ld0<=ld_b  &&  ld_b<=+ld0)
                {
                    Vocig_fft[i][0] += pld_b * p_dip_a * Docig_fft[i_b][0];
                    Vocig_fft[i][1] += pld_b * p_dip_a * Docig_fft[i_b][1];
                }
                
                if(idip<ndip-1)
                {
                i   = ikx + ikz*nx + (lv  +lv0)*nz*nx;
                i_a = ikx + ikz*nx + (ld_a+ld0)*nz*nx + (idip+1)*nld*nz*nx;
                i_b = ikx + ikz*nx + (ld_b+ld0)*nz*nx + (idip+1)*nld*nz*nx;
                if(-ld0<=ld_a  &&  ld_a<=+ld0)
                {
                    Vocig_fft[i][0] += pld_a * p_dip_b * Docig_fft[i_a][0];
                    Vocig_fft[i][1] += pld_a * p_dip_b * Docig_fft[i_a][1];
                }
                if(-ld0<=ld_b  &&  ld_b<=+ld0)
                {
                    Vocig_fft[i][0] += pld_b * p_dip_b * Docig_fft[i_b][0];
                    Vocig_fft[i][1] += pld_b * p_dip_b * Docig_fft[i_b][1];
                }
                }
                //*/
            }
    }

    // Transforming: from kx-kz to x-z
    for(lx=-lx0; lx<=lx0; lx++) {
        size_t shift = (lx+lx0)*nx*nz;
        fftw_complex *Hocig_fft_tmp = CPU_zalocFFTW1(nx*nz);
        for(i=0; i<nx*nz; i++) {
            Hocig_fft_tmp[i][0] = Hocig_fft[shift+i][0];
            Hocig_fft_tmp[i][1] = Hocig_fft[shift+i][1];
        }
        applyFFTW2D_backwardI2R(Hocig_transp+shift, Hocig_fft_tmp, nx, nz);
        free(Hocig_fft_tmp);
    }
    for(lv=-lv0; lv<=lv0; lv++) {
        size_t shift = (lv+lv0)*nx*nz;
        fftw_complex *Vocig_fft_tmp = CPU_zalocFFTW1(nx*nz);
        for(i=0; i<nx*nz; i++) {
            Vocig_fft_tmp[i][0] = Vocig_fft[shift+i][0];
            Vocig_fft_tmp[i][1] = Vocig_fft[shift+i][1];
        }
        applyFFTW2D_backwardI2R(Vocig_transp+shift, Vocig_fft_tmp, nx, nz); //comehere
        free(Vocig_fft_tmp);
    }

    // Transposing to make offset the slowest dimension
    float *Hocig_tmp = transp_dim3(Hocig_transp, nx, nz, nlx, 13);
    float *Vocig_tmp = transp_dim3(Vocig_transp, nx, nz, nlv, 13);

    memcpy(Hocig, Hocig_tmp, nx*nz*nlx*sizeof(float));
    memcpy(Vocig, Vocig_tmp, nx*nz*nlv*sizeof(float));

    free(Hocig_tmp);
    free(Vocig_tmp);

    free(Docig_transp);
    free(Hocig_transp);
    free(Vocig_transp);

    fftw_free(Hocig_fft);
    fftw_free(Vocig_fft);
    fftw_free(Docig_fft);

return;
}

float* VHocig2Gocig_simple(float **Hocig_new, float *Hocig, float *Vocig, int nx, int nz, float dx, float dz, int lx0, int lv0, int lz0, int iproc)
{
    int ix, iz, lx, lv, lz;
    int nlx, nlv, nlz;
    int ikx, ikz;

    nlx = 2*lx0 + 1;
    nlv = 2*lv0 + 1;
    lz0 = 0;
    nlz = 2*lz0 + 1;

    float dkx  = 1.0f/(nx * dx);
    float dkz  = 1.0f/(nz * dz);

    // Allocating arrays that will transformed CIGs
    fftw_complex *Hocig_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*nz*nlx);
    fftw_complex *Gocig_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*nz*nlx*nlz);

    float *Gocig_transp     = CPU_zaloc1F(nx*nz*nlx*nlz);

    // Transposing to make offset the slowest dimension
    // float *Hocig_transp = transp_dim3(Hocig, nlx, nz, nx, 13);
    float *Hocig_transp = transp_dim4(Hocig, nlx, nz, nx, 1, 13);

    // Transforming to kx-kz
    for(lx=-lx0; lx<=lx0; lx++) {
        size_t shift = (lx+lx0)*nx*nz;
        applyFFTW2D_forwardR2I(Hocig_fft+shift, Hocig_transp+shift, nx, nz);
    }

    if(iproc==0)
    {
        float *Hocig_fft_real = CPU_zaloc1F(nx*nz*nlx);
        float *Hocig_fft_imag = CPU_zaloc1F(nx*nz*nlx);
        int i;
        for(i=0; i<nx*nz*nlx; i++)    Hocig_fft_real[i] = Hocig_fft[i][0];
        for(i=0; i<nx*nz*nlx; i++)    Hocig_fft_imag[i] = Hocig_fft[i][1];
        outputSmart3d("Hocig_transp", Hocig_transp, nx, dx, 0, nz, dz, 0, nlx   , dx    , -lx0*dx);
        outputSmart3d("Hocig_fft_real", Hocig_fft_real, nx, dx, 0, nz, dz, 0, nlx   , dx    , -lx0*dx);
        outputSmart3d("Hocig_fft_imag", Hocig_fft_imag, nx, dx, 0, nz, dz, 0, nlx   , dx    , -lx0*dx);
        free(Hocig_fft_real);
        free(Hocig_fft_imag);
    }

    // Applying dip correction to HOCIG
    for(lx=-lx0; lx<=lx0; lx++)
    {
        for(ikz=0; ikz<nz; ikz++)
            for(ikx=0; ikx<nx; ikx++)
            {
                int i   = ikx + ikz*nx + (lx  +lx0)*nz*nx;

                Gocig_fft[i][0] = Hocig_fft[i][0];
                Gocig_fft[i][1] = Hocig_fft[i][1];
                
            }
    }    

    // Transform back to space coordinates
    for(lx=-lx0; lx<=lx0; lx++)
    {
        size_t shift = (lx+lx0)*nx*nz;
        applyFFTW2D_backwardI2R(Hocig_transp+shift, Hocig_fft+shift, nx, nz);
    }
    // Transposing to make offset the slowest dimension
    *Hocig_new = transp_dim4(Hocig_transp, nx, nz, nlx, 1, 13);

    for(lz=-lz0; lz<=lz0; lz++)
        for(lx=-lx0; lx<=lx0; lx++)
        {
            size_t shift = (lx+lx0)*nx*nz + (lz+lz0)*nx*nz*nlx;
            applyFFTW2D_backwardI2R(Gocig_transp+shift, Gocig_fft+shift, nx, nz);
        }

    // float *Gocig_tmp = transp_dim4(Gocig_transp,  nx, nz, nlx, nlz, 13);
    // float *Gocig     = transp_dim4(Gocig_tmp,    nlx, nz,  nx, nlz, 24);
    float *Gocig = transp_dim4(Gocig_transp,  nx, nz, nlx, nlz, 13);
    // float *Gocig = transp_dim3(Gocig_transp,  nx, nz, nlx, 13);

    free(Gocig_transp);
    free(Hocig_transp);
    // free(Gocig_tmp);
    fftw_free(Hocig_fft);
    fftw_free(Gocig_fft);

return Gocig;
}
/*
void VHocig_dispersalCorrection_Dip(float *Hocig, float *Vocig, int nx, int nz, float dx, float dz, \
                                    int lx0, int lv0, \
                                    int nlambda, float lambdaMax, float dlambda, \
                                    float dipMin, float dipMax, float ddip)
{
    int ix, iz, lx, lv;
    int nlx, nlv, nk;
    int ikx, ikz, ik;

    nlx = 2*lx0 + 1;
    nlv = 2*lv0 + 1;

    float dkx  = 1.0f/(nx * dx);
    float dkz  = 1.0f/(nz * dz);

    nk = 

    // Transposing to make offset the slowest dimension
    float *Hocig_transp = transp_dim3(Hocig, nlx, nz, nx, 13);
    float *Vocig_transp = transp_dim3(Vocig, nlv, nz, nx, 13);

    // Allocating arrays that will transformed CIGs
    fftw_complex *Hocig_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*nz*nlx);
    fftw_complex *Vocig_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*nz*nlv);

    fftw_complex *Hocig_corrected_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*nz*nlx);
    fftw_complex *Vocig_corrected_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*nz*nlv);

    float *Hocig_corrected_transp = CPU_zaloc1F(nx*nz*nlx);
    float *Vocig_corrected_transp = CPU_zaloc1F(nx*nz*nlv);

    // Transforming to kx-kz
    for(lx=-lx0; lx<=lx0; lx++) {
        size_t shift = (lx+lx0)*nx*nz;
        applyFFTW2D_forwardR2I(Hocig_fft+shift, Hocig_transp+shift, nx, nz);
    }
    for(lv=-lv0; lv<=lv0; lv++) {
        size_t shift = (lv+lv0)*nx*nz;
        applyFFTW2D_forwardR2I(Vocig_fft+shift, Vocig_transp+shift, nx, nz);
    }

    // Applying dip correction to HOCIG
    for(lx=-lx0; lx<=lx0; lx++)
    {
        for(dip=dipMin; dip<dipMax; dip+=ddip)
        {
            int idip = (int) ( 0.5f + (dip-dipMin)/ddip );
            float alpha = dip * M_PI/180.0f;
            for(ik=0; ik<nk; ik++)
            {
                float kx = ikx * dkx;
                float kz = ikz * dkz;
                float alpha = - atan2f(kx,kz);
                float llx  = lx*dx;
                float llx_ = llx * cosf(alpha);
                float lx_a = floorf(llx_/dx);
                float lx_b = lx_a+1;
                float pa = (lx_b*dx - llx_)/dx;
                float pb = (llx_ - lx_a*dx)/dx;
                float ptot = pa + pb;
                pa /= ptot;
                pb /= ptot;

                int i   = ikx + ikz*nx + (lx  +lx0)*nz*nx;
                int i_a = ikx + ikz*nx + (lx_a+lx0)*nz*nx;
                int i_b = ikx + ikz*nx + (lx_b+lx0)*nz*nx;
                
                if(-lx0<=lx_a  &&  lx_a<=+lx0)
                {
                    Hocig_corrected_fft[i_a][0] += pa*Hocig_fft[i][0];
                    Hocig_corrected_fft[i_a][1] += pa*Hocig_fft[i][1];
                }
                if(-lx0<=lx_b  &&  lx_b<=+lx0)
                {
                    Hocig_corrected_fft[i_b][0] += pb*Hocig_fft[i][0];
                    Hocig_corrected_fft[i_b][1] += pb*Hocig_fft[i][1];
                }
            }
    }
    // Applying dip correction to VOCIG
    //*
    for(lv=-lv0; lv<=lv0; lv++)
    {
        for(ikz=0; ikz<nz; ikz++)
            for(ikx=0; ikx<nx; ikx++)
            {
                float kx = ikx * dkx;
                float kz = ikz * dkz;
                float alpha = - atan2f(kx,kz);
                float llv  = lv*dz;
                float llv_ = llv * sinf(alpha);
                float lv_a = floorf(llv_/dz);
                float lv_b = lv_a+1;
                float pa = (lv_b*dz - llv_)/dz;
                float pb = (llv_ - lv_a*dz)/dz;
                float ptot = pa + pb;
                pa /= ptot;
                pb /= ptot;

                int i   = ikx + ikz*nx + (lv  +lv0) * nz*nx;
                int i_a = ikx + ikz*nx + (lv_a+lv0) * nz*nx;
                int i_b = ikx + ikz*nx + (lv_b+lv0) * nz*nx;

                if(-lv0<=lv_a  &&  lv_a<=+lv0)
                {
                    Vocig_corrected_fft[i_a][0] += pa*Vocig_fft[i][0];
                    Vocig_corrected_fft[i_a][1] += pa*Vocig_fft[i][1];
                }
                if(-lv0<=lv_b  &&  lv_b<=+lv0)
                {
                    Vocig_corrected_fft[i_b][0] += pb*Vocig_fft[i][0];
                    Vocig_corrected_fft[i_b][1] += pb*Vocig_fft[i][1];
                }
            }
    }
    // Transform back to space coordinates
    for(lx=-lx0; lx<=lx0; lx++) 
    {
        size_t shift = (lx+lx0)*nx*nz;
        applyFFTW2D_backwardI2R(Hocig_corrected_transp+shift, Hocig_corrected_fft+shift, nx, nz);
    }
    for(lv=-lv0; lv<=lv0; lv++) 
    {
        size_t shift = (lv+lv0)*nx*nz;
        applyFFTW2D_backwardI2R(Vocig_corrected_transp+shift, Vocig_corrected_fft+shift, nx, nz);
    }

    float *Hocig_corrected = transp_dim3(Hocig_corrected_transp, nx, nz, nlx, 13);
    float *Vocig_corrected = transp_dim3(Vocig_corrected_transp, nx, nz, nlv, 13);

    memcpy(Hocig, Hocig_corrected, nx*nz*nlx*sizeof(float));
    memcpy(Vocig, Vocig_corrected, nx*nz*nlv*sizeof(float));

    free(Hocig_transp);
    free(Vocig_transp);

    fftw_free(Hocig_fft);
    fftw_free(Vocig_fft);

    fftw_free(Hocig_corrected_fft);
    fftw_free(Vocig_corrected_fft);

    free(Hocig_corrected_transp);
    free(Vocig_corrected_transp);

    free(Hocig_corrected);
    free(Vocig_corrected);

return;
}
*/
float* transp_dim3(float *m_input, long long int n1, long long int n2, long long int n3, int plane)
{
    float *m_transp = CPU_zaloc1F(n1*n2*n3);
    long long int i1, i2, i3;

    if(plane==12)
    {
        for(i3=0; i3<n3; i3++)
          for(i2=0; i2<n2; i2++)
            for(i1=0; i1<n1; i1++)
            {
                long long int i_in = i1 + i2*n1 + i3*n2*n1;
                long long int i_tr = i2 + i1*n2 + i3*n2*n1;
                m_transp[i_tr] = m_input[i_in];
            }
    }
    if(plane==13)
    {
        for(i3=0; i3<n3; i3++)
          for(i2=0; i2<n2; i2++)
            for(i1=0; i1<n1; i1++)
            {
                long long int i_in = i1 + i2*n1 + i3*n2*n1;
                long long int i_tr = i3 + i2*n3 + i1*n2*n3;
                m_transp[i_tr] = m_input[i_in];
            }
    }
    if(plane==23)
    {
        for(i3=0; i3<n3; i3++)
          for(i2=0; i2<n2; i2++)
            for(i1=0; i1<n1; i1++)
            {
                long long int i_in = i1 + i2*n1 + i3*n2*n1;
                long long int i_tr = i1 + i3*n1 + i2*n3*n1;
                m_transp[i_tr] = m_input[i_in];
            }
    }

return m_transp;
}
float* transp_dim4(float *m_input, long long int n1, long long int n2, long long int n3, long long int n4, int plane)
{
    float *m_transp = CPU_zaloc1F(n1*n2*n3*n4);
    long long int i1, i2, i3, i4, i_in, i_tr;

    if(plane==12)
    {
        for(i4=0; i4<n4; i4++)
          for(i3=0; i3<n3; i3++)
            for(i2=0; i2<n2; i2++)
              for(i1=0; i1<n1; i1++)
              {
                  i_in = i1 + i2*n1 + i3*n2*n1 + i4*n3*n2*n1;
                  i_tr = i2 + i1*n2 + i3*n2*n1 + i4*n3*n2*n1;
                  m_transp[i_tr] = m_input[i_in];
              }
    }
    if(plane==13)
    {
        for(i4=0; i4<n4; i4++)
          for(i3=0; i3<n3; i3++)
            for(i2=0; i2<n2; i2++)
              for(i1=0; i1<n1; i1++)
              {
                  i_in = i1 + i2*n1 + i3*n2*n1 + i4*n3*n2*n1;
                  i_tr = i3 + i2*n3 + i1*n2*n3 + i4*n3*n2*n1;
                  m_transp[i_tr] = m_input[i_in];
              }
    }
    if(plane==23)
    {
        for(i4=0; i4<n4; i4++)
          for(i3=0; i3<n3; i3++)
            for(i2=0; i2<n2; i2++)
              for(i1=0; i1<n1; i1++)
              {
                  i_in = i1 + i2*n1 + i3*n2*n1 + i4*n3*n2*n1;
                  i_tr = i1 + i3*n1 + i2*n3*n1 + i4*n3*n2*n1;
                  m_transp[i_tr] = m_input[i_in];
              }
    }
    if(plane==14)
    {
        for(i4=0; i4<n4; i4++)
          for(i3=0; i3<n3; i3++)
            for(i2=0; i2<n2; i2++)
              for(i1=0; i1<n1; i1++)
              {
                  i_in = i1 + i2*n1 + i3*n2*n1 + i4*n3*n2*n1;
                  i_tr = i4 + i2*n4 + i3*n2*n4 + i1*n3*n2*n4;
                  m_transp[i_tr] = m_input[i_in];
              }
    }
    if(plane==24)
    {
        for(i4=0; i4<n4; i4++)
          for(i3=0; i3<n3; i3++)
            for(i2=0; i2<n2; i2++)
              for(i1=0; i1<n1; i1++)
              {
                  i_in = i1 + i2*n1 + i3*n2*n1 + i4*n3*n2*n1;
                  i_tr = i1 + i4*n1 + i3*n4*n1 + i2*n3*n4*n1;
                  m_transp[i_tr] = m_input[i_in];
              }
    }
    if(plane==34)
    {
        for(i4=0; i4<n4; i4++)
          for(i3=0; i3<n3; i3++)
            for(i2=0; i2<n2; i2++)
              for(i1=0; i1<n1; i1++)
              {
                  i_in = i1 + i2*n1 + i3*n2*n1 + i4*n3*n2*n1;
                  i_tr = i1 + i2*n1 + i4*n2*n1 + i3*n4*n2*n1;
                  m_transp[i_tr] = m_input[i_in];
              }
    }
    

return m_transp;
}
void applyFFTW2D_forwardR2I(fftw_complex *out, float *inp, int n1, int n2)
{
    // n1: fastest, n2:slowest
    fftw_complex   *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n1*n2);
    fftw_plan planForw = fftw_plan_dft_2d(n2, n1, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    int i;
    for(i=0; i<n1*n2; i++)
    {
        in[i][0] = inp[i]/sqrtf((1.0f*n1*n2));
        in[i][1] = 0.0;
    }
    fftw_execute(planForw);
    fftw_destroy_plan(planForw);
    fftw_free(in);
return;
}
void applyFFTW2D_backwardI2R(float *out, fftw_complex *inp, long long int n1, long long int n2)
{
    long long int i;

    fftw_complex  *in = inp;
    fftw_complex  *ou = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n1*n2);

    fftw_plan planBack = fftw_plan_dft_2d(n2, n1, in, ou, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(planBack);

    for(i=0; i<n1*n2; i++)    out[i] = ou[i][0]/sqrtf((1.0f*n1*n2));
        
    
    fftw_destroy_plan(planBack);

return;
}
void applyFFTW2D_backwardI2R_debug(float *out, fftw_complex *inp, int n1, int n2)
{
    int i1, i2, i3, i_p, i_n, i;

    fftw_complex  *in = inp;
    fftw_complex  *ou = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n1*n2);

    fftw_plan planBack = fftw_plan_dft_2d(n1, n2, in, ou, FFTW_BACKWARD, FFTW_ESTIMATE);

    /*for(i2=1; i2<n2/2; i2++)
        for(i1=1; i1<n1/2; i1++)
        {
            int i_p     = i1 + i2*n1;

            int i_n_1   = (n1-i1) +     i2 *n1,  s_1   = -1;
            int i_n_2   =     i1  + (n2-i2)*n1,  s_2   = -1;
            int i_n_12  = (n1-i1) + (n2-i2)*n1,  s_12  = +1;

            in[i_n_1][0] =         in[i_p][0];
            in[i_n_1][1] = s_1   * in[i_p][1];

            in[i_n_2][0] =         in[i_p][0];
            in[i_n_2][1] = s_2   * in[i_p][1];

            in[i_n_12][0] =        in[i_p][0];
            in[i_n_12][1] = s_12 * in[i_p][1];
        }*/

    fftw_execute(planBack);

    // printf("\n sqrtf((1.0f*n1*n2))=%g   1/sqrtf((1.0f*n1*n2))=%g", sqrtf((1.0f*n1*n2)), 1.0f/sqrtf((1.0f*n1*n2)));
    int foundNaN=0, inpNaN_0=0, inpNaN_1=0, ouNaN_0=0, ouNaN_1=0;
    for(i=0; i<n1*n2; i++) {
        out[i] = ou[i][0]/sqrtf((1.0f*n1*n2));
        if(isnan(out[i])) {
            // printf("\n n1=%d  n2=%d  sqrtf((1.0f*n1*n2))=%g   1/sqrtf((1.0f*n1*n2))=%g   ou[i][0]=%g    ou[i][0]/sqrtf((1.0f*n1*n2))=%g", \
            n1, n2, sqrtf((1.0f*n1*n2)), 1.0f/sqrtf((1.0f*n1*n2)), ou[i][0], ou[i][0]/sqrtf((1.0f*n1*n2))); // getthere
            foundNaN = 1;
        }
        if(isnan(ou[i][0]))   ouNaN_0  = 1;
        if(isnan(ou[i][1]))   ouNaN_1  = 1;
        if(isnan(inp[i][0])  ||  !isfinite(inp[i][0]))  inpNaN_0 = 1;
        if(isnan(inp[i][1])  ||  !isfinite(inp[i][1]))  inpNaN_1 = 1;
    }

    double max = 0.0;
    double min = 1.0e+20;
    for(i=0; i<n1*n2; i++) {
        if(max<fabs(in[i][0]))    max = fabs(in[i][0]);
        if(max<fabs(in[i][1]))    max = fabs(in[i][1]);
        if(min>fabs(in[i][0]))    min = fabs(in[i][0]);
        if(min>fabs(in[i][1]))    min = fabs(in[i][1]);
    }

    if(inpNaN_0)    printf("  == inp0 has Nan  ==  ");
    if(inpNaN_1)    printf("  == inp1 has Nan  ==  ");
    if(ouNaN_0)     printf("  == ou0  has Nan  ==  ");
    if(ouNaN_1)     printf("  == ou1  has Nan  ==  ");
    if(foundNaN)    printf("  !! NaN was found !!  sqrtf((1.0f*n1*n2)=%g    max=%g  min=%g", sqrtf((1.0f*n1*n2)), min, max);
        
    
    fftw_destroy_plan(planBack);

return;
}
/*
fftw_complex* applyFFTW2D_forwardR2I(float *inp, int n1, int n2)
{
    // n1: fastest, n2:slowest
    fftw_complex   *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n1*n2);
    fftw_complex   *ou = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n1*n2);
    fftw_plan planForw = fftw_plan_dft_2d(n2, n1, in, ou, FFTW_FORWARD,  FFTW_ESTIMATE);
    int i;
    for(i=0; i<n1*n2; i++)
    {
        in[i][0] = inp[i]/sqrtf((1.0f*n1*n2));
        in[i][1] = 0.0;
    }
    fftw_execute(planForw);
    fftw_destroy_plan(planForw);
    fftw_free(in);
return ou;
}
float* applyFFTW2D_backwardI2R(fftw_complex *inp, int n1, int n2)
{
    int            i1, i2, i3, i_p, i_n, i;
    float         *out = CPU_zaloc1F(n1*n2);

    fftw_complex  *in = inp;
    fftw_complex  *ou = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n1*n2);

    fftw_plan planBack = fftw_plan_dft_3d(n2, n1, in, ou, FFTW_BACKWARD, FFTW_ESTIMATE);

    for(i2=1; i2<n2/2; i2++)
        for(i1=1; i1<n1/2; i1++)
        {
            int i_p     = i1 + i2*n1;

            int i_n_1   = (n1-i1) +     i2 *n1,  s_1   = -1;
            int i_n_2   =     i1  + (n2-i2)*n1,  s_2   = -1;
            int i_n_12  = (n1-i1) + (n2-i2)*n1,  s_12  = +1;

            in[i_n_1][0] =           in[i_p][0];
            in[i_n_1][1] = s_1     * in[i_p][1];

            in[i_n_2][0] =           in[i_p][0];
            in[i_n_2][1] = s_2     * in[i_p][1];

            in[i_n_12][0] =          in[i_p][0];
            in[i_n_12][1] = s_12   * in[i_p][1];
        }

    fftw_execute(planBack);

    for(i=0; i<n1*n2; i++)
        out[i] = ou[i][0]/sqrtf((1.0f*n1*n2));
    
    fftw_destroy_plan(planBack);

return out;
}
//*/

fftw_complex* applyFFTW3D_forwardR2I(float *inp, int n1, int n2, int n3)
{
    // n1: fastest, n3:slowest
    fftw_complex   *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n1*n2*n3);
    fftw_complex   *ou = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n1*n2*n3);
    fftw_plan planForw = fftw_plan_dft_3d(n3, n2, n1, in, ou, FFTW_FORWARD,  FFTW_ESTIMATE);
    int i;
    for(i=0; i<n1*n2*n3; i++)
    {
        in[i][0] = inp[i]/sqrtf((1.0f*n1*n2*n3));
        in[i][1] = 0.0;
    }
    fftw_execute(planForw);
    fftw_destroy_plan(planForw);
    fftw_free(in);
return ou;
}
float* applyFFTW3D_backwardI2R(fftw_complex *inp, int n1, int n2, int n3)
{
    int            i1, i2, i3, i_p, i_n, i;
    float         *out = CPU_zaloc1F(n1*n2*n3);

    fftw_complex  *in = inp;
    fftw_complex  *ou = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n1*n2*n3);

    fftw_plan planBack = fftw_plan_dft_3d(n3, n2, n1, in, ou, FFTW_BACKWARD, FFTW_ESTIMATE);

    for(i3=1; i3<n3/2; i3++)
        for(i2=1; i2<n2/2; i2++)
            for(i1=1; i1<n1/2; i1++)
            {
                int i_p     = i1 + i2*n1 + i3*n1*n2;

                int i_n_1   = (n1-i1) +     i2 *n1 +      i3*n1*n2,  s_1   = -1;
                int i_n_2   =     i1  + (n2-i2)*n1 +      i3*n1*n2,  s_2   = -1;
                int i_n_12  = (n1-i1) + (n2-i2)*n1 +      i3*n1*n2,  s_12  = +1;

                int i_n_3   =     i1  +     i2 *n1 + (n3-i3)*n1*n2,  s_3   = -1;
                int i_n_13  = (n1-i1) +     i2 *n1 + (n3-i3)*n1*n2,  s_13  = +1;
                int i_n_23  =     i1  + (n2-i2)*n1 + (n3-i3)*n1*n2,  s_23  = +1;
                int i_n_123 = (n1-i1) + (n2-i2)*n1 + (n3-i3)*n1*n2,  s_123 = -1;

                in[i_n_1][0] =           in[i_p][0];
                in[i_n_1][1] = s_1     * in[i_p][1];

                in[i_n_2][0] =           in[i_p][0];
                in[i_n_2][1] = s_2     * in[i_p][1];

                in[i_n_12][0] =          in[i_p][0];
                in[i_n_12][1] = s_12   * in[i_p][1];

                in[i_n_3][0] =           in[i_p][0];
                in[i_n_3][1] = s_3     * in[i_p][1];

                in[i_n_13][0] =          in[i_p][0];
                in[i_n_13][1] = s_13   * in[i_p][1];

                in[i_n_23][0] =          in[i_p][0];
                in[i_n_23][1] = s_23   * in[i_p][1];

                in[i_n_123][0] =         in[i_p][0];
                in[i_n_123][1] = s_123 * in[i_p][1];
                
                // in[i_n][0] = + in[i_p][0];
                // in[i_n][1] = - in[i_p][1];
            }

    fftw_execute(planBack);

    for(i=0; i<n1*n2*n3; i++)
        out[i] = ou[i][0]/sqrtf((1.0f*n1*n2*n3));
    
    fftw_destroy_plan(planBack);

return out;
}

/*
fftw_complex *applyFFTW_forwardR2I(float *inp, int nt, float dt)
{
    fftw_complex  *trcT, *trcF;
    fftw_plan      planForw;
    int            nfreq, ifreq, it;
    nfreq   = nt;
    trcT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt);
    trcF = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt);
    planForw = fftw_plan_dft_1d(nt, trcT, trcF, FFTW_FORWARD,  FFTW_ESTIMATE);
    for(it=0; it<nt; it++)
    {
        trcT[it][0] = inp[it]/sqrtf((1.0f*nt));
        trcT[it][1] = 0.0;
    }
    fftw_execute(planForw);
    fftw_destroy_plan(planForw);
    fftw_free(trcT);
return trcF;
}
float *applyFFTW_backwardI2R(fftw_complex *inp, int nt, float dt)
{
    fftw_complex  *trcT, *trcF;
    fftw_plan      planBack;
    int            nfreq, ifreq, it;
    float         *out = CPU_zaloc1F(nt);
    nfreq   = nt;
    trcT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt);
    trcF = inp;
    planBack = fftw_plan_dft_1d(nt, trcF, trcT, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(planBack);
    for(it=0; it<nt; it++)
        out[it] = trcT[it][0]/sqrtf((1.0f*nt));
    
    fftw_destroy_plan(planBack);

return out;
}
*/

/*
void applyFFTW_template(float *inp, int nt, float dt)
{

    fftw_complex  *trcT, *trcF;
    fftw_plan      planForw, planBack;
    int            nfreq, ifreq, it;
    float          dfreq, omega;
    float          freq;

    nfreq   = nt;
    dfreq   = 1.0/(nt*dt);

    trcT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt);
    trcF = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt);

    planForw = fftw_plan_dft_1d(nt, trcT, trcF, FFTW_FORWARD,  FFTW_ESTIMATE);
    planBack = fftw_plan_dft_1d(nt, trcF, trcT, FFTW_BACKWARD, FFTW_ESTIMATE);


    for(it=0; it<nt; it++)
    {
        trcT[it][0] = inp[it]/sqrtf((1.0f*nt));
        trcT[it][1] = 0.0;
    }

    fftw_execute(planForw);
    fftw_execute(planBack);
        
    for(it=0; it<nt; it++)
    {
        inp[it] = trcT[it][0]/sqrtf((1.0f*nt));
    }
    
    fftw_destroy_plan(planForw);
    fftw_destroy_plan(planBack);

    fftw_free(trcT);
    fftw_free(trcF);

return;
}
//*/