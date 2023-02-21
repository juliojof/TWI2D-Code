void MT_aplicaFiltroPassaBanda(float *pulso, int nt, float dt, float f1, float f2, float f3, float f4)
{

    fftw_complex  *tracoTempo, *tracoFreq;
    fftw_plan      planoForward, planoBackward;
    int            nf, ifreq, it;
    float          df, omega;
    float          filtro, freq;

    
    nf   = nt;
    df   = 1.0/(nt*dt);

    tracoTempo = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt);
    tracoFreq  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt);

    planoForward  = fftw_plan_dft_1d(nt, tracoTempo, tracoFreq,  FFTW_FORWARD,  FFTW_ESTIMATE);
    planoBackward = fftw_plan_dft_1d(nt, tracoFreq,  tracoTempo, FFTW_BACKWARD, FFTW_ESTIMATE);


    for(it=0; it<nt; it++)
    {
        tracoTempo[it][0] = pulso[it]/sqrtf((1.0f*nt));
        tracoTempo[it][1] = 0.0;
    }

    fftw_execute(planoForward);        
   

    for(ifreq=0; ifreq<=nf/2; ifreq++)
    {
        freq  = ifreq*df;
        
        if(freq<f1)
            filtro = 0.0f;
        
        if(freq>=f1  &&  freq<f2)
            filtro = (freq-f1)/(f2-f1);
        
        if(freq>=f2  &&  freq<=f3)
            filtro = 1.0f;
        
        if(freq>f3  &&  freq<=f4)
            filtro = (f4-freq)/(f4-f3);
        
        if(freq>f4)
            filtro = 0.0f;
        
        
            tracoFreq[ifreq   ][0] *= filtro;
            tracoFreq[ifreq   ][1] *= filtro;
        
        if(ifreq>0 && ifreq<nf-ifreq)
        {
            tracoFreq[nf-ifreq][0] *= filtro;
            tracoFreq[nf-ifreq][1] *= filtro;
        }
        
    }
        
    fftw_execute(planoBackward);
        
    for(it=0; it<nt; it++)
    {
        pulso[it] = tracoTempo[it][0]/sqrtf((1.0f*nt));
    }
    

    fftw_destroy_plan(planoForward);
    fftw_destroy_plan(planoBackward);

    fftw_free(tracoTempo);
    fftw_free(tracoFreq);

return;
}