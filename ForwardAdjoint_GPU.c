void modelShot_GPU(float *vp, Shot *shot, int jt, int nt, int nx, int nz, float dx, float dz, float dt, \
                   float ox, float oz, float ot, int recordMovie, int ishotMovie) {
    /*
    int it;
    float *wav1, *wav2, *lap, *sigX, *sigZ, *sigma, *sigmaInv;
    wav1 = CPU_zaloc1F(nx*nz);
    wav2 = CPU_zaloc1F(nx*nz);
    lap  = CPU_zaloc1F(nx*nz);
    sigX = CPU_zaloc1F(nx);
    sigZ = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);

    // applyHalfDerivative2Source(shot);

    FILE *fWavMovie_Hdr, *fWavMovie_Bin;
    if(shot->ishot==ishotMovie && recordMovie) {
        initFiles3d("./movie_modeling", &fWavMovie_Hdr, &fWavMovie_Bin, nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
    }
    for(it=0; it<nt; it++)
    {
        recordData(shot, wav2, nx, nz, it);
        if(it%jt==0  &&  shot->ishot==ishotMovie && recordMovie)
            writeSamples(wav2, fWavMovie_Bin, nx*nz);
        
        injectSource_Acoustic(wav1, shot->source, it, shot->ixs, shot->izs, nx, nz, dx, dz);
        //propagate_Acoustic(wav1, wav2, lap, vp, nx, nz);
        // propagate_Acoustic_abs(wav1, wav2, lap, vp, sigma, nx, nz);
        propagate_Acoustic_abs2(wav1, wav2, lap, vp, sigma, sigmaInv, nx, nz);
        float *tmp = wav2;
        wav2=wav1;
        wav1=tmp;
    }
    if(shot->ishot==ishotMovie && recordMovie) {
        fclose(fWavMovie_Bin);
        fclose(fWavMovie_Hdr);
    }
    free(sigma);
    free(sigX);
    free(sigZ);
    free(wav1);
    free(wav2);
    free(lap);
*/

return;
}