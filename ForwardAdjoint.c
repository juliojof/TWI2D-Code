void modelShot(float *vp, Shot *shot, int jt, int nt, int nx, int nz, float dx, float dz, float dt, \
               float ox, float oz, float ot, int recordMovie, int ishotMovie) {
    int it;
    float *wav1, *wav2, *lap, *sigX, *sigZ, *sigma, *sigmaInv, *tmp;
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
        initFilesSmart3d("movie_modeling", &fWavMovie_Hdr, &fWavMovie_Bin, nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
    }
    for(it=0; it<nt; it++)
    {
        recordData(shot, wav2, nx, nz, it);
        if(it%jt==0  &&  shot->ishot==ishotMovie && recordMovie)
            writeSamples(wav2, fWavMovie_Bin, nx*nz);
        
        injectSource_Acoustic(wav1, shot->source, it, (shot->ixs)[0], (shot->izs)[0], nx, nz, dx, dz);
        //propagate_Acoustic(wav1, wav2, lap, vp, nx, nz);
        // propagate_Acoustic_abs(wav1, wav2, lap, vp, sigma, nx, nz);
        propagate_Acoustic_abs2(wav1, wav2, lap, vp, sigma, sigmaInv, nx, nz);
        tmp = wav2;
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
    free(sigmaInv);

return;
}

void modelShot_GPU(int gpu_device, float *vp, Shot3D *shot, int jt, int nt, int nx, int nz, int nxb, int nzb, \
                   float dx, float dz, float dt, float ox, float oz, float ot, int recordMovie, int ishotMovie) {
    
    int it;
    float *wav1, *wav2, *lap, *sigX, *sigZ, *sigma, *sigmaInv;
    wav1 = CPU_zaloc1F(nx*nz);
    wav2 = CPU_zaloc1F(nx*nz);
    lap  = CPU_zaloc1F(nx*nz);
    sigX = CPU_zaloc1F(nx);
    sigZ = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);


    // Initializing GPU environment for current shot
    int lx0=0, lz0=0, lt0=0;
    GPU_OPER_initEnvironment(gpu_device, nx, nz, nt, dx, dz, dt, shot->nsou, shot->nrec, lx0, lz0, lt0, nxb, nzb, MODE_CIGS_FULL);
    GPU_OPER_alocArrays_acoustic();
    GPU_OPER_getVelP_CPU2GPU(vp);
    GPU_OPER_getSigmas_CPU2GPU(sigma, sigmaInv);
    GPU_OPER_getSourceWavelet_CPU2GPU(shot->source);
    GPU_OPER_getRecordedData_CPU2GPU(shot->seismogram);
    GPU_OPER_getSourceReceiverIndexes_CPU2GPU(shot->ixs, shot->izs, shot->ixr, shot->izr);

    FILE *fWavMovie_Hdr, *fWavMovie_Bin;
    if(shot->ishot==ishotMovie && recordMovie) {
        initFilesSmart3d("movie_modeling", &fWavMovie_Hdr, &fWavMovie_Bin, nx, dx, ox, nz, dz, oz, 1 + nt/jt, dt*jt, ot);
    }
    for(it=0; it<nt; it++)
    {
        GPU_OPER_setTimeStep(it);
        GPU_WRAP_recordInjectData(0);

        if(it%jt==0  &&  shot->ishot==ishotMovie  &&  recordMovie) {
            GPU_OPER_copyPresentWavefield_GPU2CPU(wav2);
            writeSamples(wav2, fWavMovie_Bin, nx*nz);
        }
        
        GPU_WRAP_injectSource_Acoustic();
        GPU_WRAP_propagate_Acoustic_abs2();
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
    free(sigmaInv);

    GPU_OPER_getRecordedData_GPU2CPU(shot->seismogram);
    // demultiplexSeismogram(shot->seismogram, shot->nt, shot->nrec);
    GPU_OPER_freeArrays_acoustic();

return;
}

void modelShotBorn_GPU(int gpu_device, float *vp, float *ref, Shot *shot, int jt, \
                       int nt, int nx, int nz, float dx, float dz, float dt, \
                       float ox, float oz, float ot, int recordMovie, int ishotMovie) {
    int it;
    float *wav1, *wav2, *lap, *sigX, *sigZ, *sigma, *sigmaInv;
    wav1 = CPU_zaloc1F(nx*nz);
    wav2 = CPU_zaloc1F(nx*nz);
    lap  = CPU_zaloc1F(nx*nz);
    sigX = CPU_zaloc1F(nx);
    sigZ = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);


    // Initializing GPU environment for current shot
    int lx0=0, lz0=0, lt0=0;
    GPU_OPER_initEnvironment(gpu_device, nx, nz, nt, dx, dz, dt, shot->nsou, shot->nrec, lx0, lz0, lt0, nxb, nzb, MODE_CIGS_FULL);
    GPU_OPER_alocArrays_acoustic();
    GPU_OPER_alocArrays_acoustic_BornSct();
    GPU_OPER_getVelP_CPU2GPU(vp);
    GPU_OPER_copyRef_CPU2GPU(ref);
    GPU_OPER_getSigmas_CPU2GPU(sigma, sigmaInv);
    GPU_OPER_getSourceWavelet_CPU2GPU(shot->source);
    GPU_OPER_getSourceReceiverIndexes_CPU2GPU(shot->ixs, shot->izs, shot->ixr, shot->izr);

    FILE *fp = fopen("source_debug", "w");
    fwrite(shot->source, sizeof(float), nt, fp);
    fclose(fp);

    if(shot->ishot==ishotMovie && recordMovie) {
        printf("\n Source location  ixs=%d  izs=%d  \n", shot->ixs[0], shot->izs[0]);
        printf("\n dx=%f  dz=%f  dt=%f  \n", dx, dz, dt);
    }

    FILE *fMovieInc_Hdr, *fMovieInc_Bin;
    FILE *fMovieSct_Hdr, *fMovieSct_Bin;
    if(shot->ishot==ishotMovie && recordMovie) {
        initFilesSmart3d("movie_inc_wave", &fMovieInc_Hdr, &fMovieInc_Bin, nx, dx, ox, nz, dz, oz, 1 + nt/jt, dt*jt, ot);
        initFilesSmart3d("movie_sct_wave", &fMovieSct_Hdr, &fMovieSct_Bin, nx, dx, ox, nz, dz, oz, 1 + nt/jt, dt*jt, ot);
    }
    for(it=0; it<nt; it++)
    {
        GPU_OPER_setTimeStep(it);
        GPU_WRAP_recordInjectData(2);

        if(it%jt==0  &&  shot->ishot==ishotMovie  &&  recordMovie) {
            GPU_OPER_copyPresentWavefield_GPU2CPU(wav1);
            GPU_OPER_copyPresentScatteredWavefield_GPU2CPU(wav2);
            writeSamples(wav1, fMovieInc_Bin, nx*nz);
            writeSamples(wav2, fMovieSct_Bin, nx*nz);
        }
        
        GPU_OPER_wavefieldTimeSecondDerivative_part1(0);
        GPU_WRAP_injectSource_Acoustic();
        GPU_WRAP_propagate_Acoustic_abs2();
        GPU_OPER_wavefieldTimeSecondDerivative_part2(0);
        GPU_WRAP_scatter_Acoustic();
        GPU_WRAP_propagate_Acoustic_Scattered_abs2();
    }
    if(shot->ishot==ishotMovie && recordMovie) {
        fclose(fMovieInc_Bin);
        fclose(fMovieInc_Hdr);
        fclose(fMovieSct_Bin);
        fclose(fMovieSct_Hdr);
    }
    free(sigma);
    free(sigX);
    free(sigZ);
    free(wav1);
    free(wav2);
    free(lap);
    free(sigmaInv);

    GPU_OPER_getRecordedData_GPU2CPU(shot->seismogram);
    GPU_OPER_freeArrays_acoustic();
    GPU_OPER_freeArrays_acoustic_BornSct();

return;
}

void migrateShot(float *image, float *vp, Shot *shot, int jt, int nt, int nx, int nz, \
                 float dx, float dz, float dt, float ox, float oz, float ot, int recordMovie, int ishotMovie) {

    int it;
    float *wav1, *wav2, *lap, *sigX, *sigZ, *sigma, *sigmaInv, *ws;
    wav1 = CPU_zaloc1F(nx*nz);
    wav2 = CPU_zaloc1F(nx*nz);
    lap  = CPU_zaloc1F(nx*nz);
    ws   = CPU_zaloc1F(nx*nz);
    sigX = CPU_zaloc1F(nx);
    sigZ = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);

    //applyHalfDerivative2Source(shot);
    // applyHalfDerivative2Data(shot);

    // Source wavefield
    FILE *fWavMovie_Hdr, *fWavMovie_Bin, *fWS_Bin, *fWS_Hdr;
    if(shot->ishot==ishotMovie) {
        initFilesSmart3d("movie_WS", &fWavMovie_Hdr, &fWavMovie_Bin, nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
        initFilesSmart3d("WS_imaging", &fWS_Hdr, &fWS_Bin, nz, dz, oz, nx, dx, ox, nt, dt, ot);
    }
    for(it=0; it<nt; it++)
    {
        fwrite(wav2, sizeof(float), nx*nz, fWS_Bin);
        if(it%jt==0)
            writeSamples(wav2, fWavMovie_Bin, nx*nz);
        
        injectSource_Acoustic(wav1, shot->source, it, (shot->ixs)[0], (shot->izs)[0], nx, nz, dx, dz);
        propagate_Acoustic_abs(wav1, wav2, lap, vp, sigma, nx, nz);
        float *tmp = wav2;
        wav2=wav1;
        wav1=tmp;
    }
    if(shot->ishot==ishotMovie) {
        fclose(fWavMovie_Bin);
        fclose(fWavMovie_Hdr);
    }

    zeroArray(wav2,nx*nz);
    zeroArray(wav1,nx*nz);

    // Receiver wavefield
    if(shot->ishot==ishotMovie) {
        initFilesSmart3d("movie_WR", &fWavMovie_Hdr, &fWavMovie_Bin, nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, 0.0f);
    }
    for(it=nt-1; it>=0; it--)
    {
        injectData(shot, wav1, nx, nz, dx, dz, it);
        //injectDataTaper(shot, wav1, nx, nz, dx, dz, it);
        fseek(fWS_Bin, it*nx*nz*sizeof(float), 0);
        size_t err = fread(ws, sizeof(float), nx*nz, fWS_Bin);
        imagingCondition(image, ws, wav2, nx, nz);
        if(it%jt==0 && shot->ishot==ishotMovie)
            writeSamples(ws, fWavMovie_Bin, nx*nz);
            //writeSamples(wav2, fWavMovie_Bin, nx*nz);
        propagate_Acoustic_abs(wav1, wav2, lap, vp, sigma, nx, nz);
        float *tmp = wav2;
        wav2=wav1;
        wav1=tmp;
    }
    fclose(fWS_Bin);
    if(shot->ishot==ishotMovie) {
        fclose(fWavMovie_Bin);
        fclose(fWavMovie_Hdr);
    }

    free(sigma);
    free(sigX);
    free(sigZ);
    free(wav1);
    free(wav2);
    free(lap);
    free(ws);
    free(sigmaInv);

return;
}
void emigrateShot_lx(float *image, float *eimage_shot, float *eimage_shot_theta, \
                     float *eimage, float *eimage_theta, float *vp, Shot *shot, \
                     int jt, int nt, int nx, int nz, float dx, float dz, float dt, \
                     float ox, float oz, float ot, int lx0, int ntheta, float dtheta, float otheta, \
                     int recordMovie, int ishotMovie) {

    int it;
    float *wav1, *wav2, *lap, *sigX, *sigZ, *sigma, *sigmaInv, *ws;
    wav1 = CPU_zaloc1F(nx*nz);
    wav2 = CPU_zaloc1F(nx*nz);
    lap  = CPU_zaloc1F(nx*nz);
    ws   = CPU_zaloc1F(nx*nz);
    sigX = CPU_zaloc1F(nx);
    sigZ = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);

    // Source wavefield
    FILE *fWavMovie_Hdr, *fWavMovie_Bin;
    FILE *fWS_Bin, *fWS_Hdr;
    FILE *fWR_Bin, *fWR_Hdr;
    FILE *fIM_Bin, *fIM_Hdr;
    if(shot->ishot==ishotMovie && recordMovie) {
        initFilesSmart3d("movie_WS", &fWavMovie_Hdr, &fWavMovie_Bin, nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
    }
    initFilesSmart3d("WS_imaging",     &fWS_Hdr,       &fWS_Bin, nz, dz, oz, nx, dx, ox, nt   , dt   , ot);

    //applyHalfDerivative2Source(shot);
    // applyHalfDerivative2Data(shot);

    for(it=0; it<nt; it++)
    {
        if(it%jt==0) fwrite(wav2, sizeof(float), nx*nz, fWS_Bin);
        if(shot->ishot==ishotMovie && recordMovie) {
            if(it%jt==0)    writeSamples(wav2, fWavMovie_Bin, nx*nz);
        }
        
        injectSource_Acoustic(wav1, shot->source, it, (shot->ixs)[0], (shot->izs)[0], nx, nz, dx, dz);
        propagate_Acoustic_abs(wav1, wav2, lap, vp, sigma, nx, nz);
        float *tmp = wav2;
        wav2=wav1;
        wav1=tmp;
    }
    if(shot->ishot==ishotMovie && recordMovie) {
        fclose(fWavMovie_Bin);
        fclose(fWavMovie_Hdr);
    }

    zeroArray(wav2,nx*nz);
    zeroArray(wav1,nx*nz);

    // Receiver wavefield
    if(shot->ishot==ishotMovie && recordMovie) {
        initFilesSmart3d("movie_WR", &fWavMovie_Hdr, &fWavMovie_Bin, nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, 0.0f);
        initFilesSmart3d("movie_IM",       &fIM_Hdr,       &fIM_Bin, nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, 0.0f);
    }
    for(it=nt-1; it>=0; it--)
    {
        injectData(shot, wav1, nx, nz, dx, dz, it);
        //injectDataTaper(shot, wav1, nx, nz, dx, dz, it);
        if(it%jt==0) {
            fseek(fWS_Bin, it/jt*nx*nz*sizeof(float), 0);
            size_t err = fread(ws, sizeof(float), nx*nz, fWS_Bin);
        }
        if(it%jt==0) {
            //if(eimage_theta!=NULL)  eimagingCondition_lx_theta(eimage, eimage_shot, eimage_theta, eimage_shot_theta, ws, wav2, nx, dx, nz, dz, lx0, ntheta, dtheta, otheta);
            //else                    eimagingCondition_lx(eimage, ws, wav2, nx, nz, lx0);
            imagingCondition(image, ws, wav2, nx, nz);
            eimagingCondition_lx(eimage_shot, ws, wav2, nx, nz, lx0);

        }
        if(shot->ishot==ishotMovie && recordMovie) {
            if(it%jt==0) {
                writeSamples(wav2, fWavMovie_Bin, nx*nz);
                writeSamples(image,      fIM_Bin, nx*nz);
            }
        }
        propagate_Acoustic_abs(wav1, wav2, lap, vp, sigma, nx, nz);
        float *tmp = wav2;
        wav2=wav1;
        wav1=tmp;
    }
    lx2theta_forward(1, eimage_shot, eimage_shot_theta, nx, dx, nz,  dz, lx0, ntheta, dtheta, otheta);
    accumulate(eimage, eimage_shot, nx*nz*(2*lx0+1));
    accumulate(eimage_theta, eimage_shot_theta, nx*nz*ntheta);
    fclose(fWS_Bin);
    fclose(fWS_Hdr);
    if(shot->ishot==ishotMovie && recordMovie) {
        fclose(fIM_Bin);
        fclose(fIM_Hdr);
        fclose(fWavMovie_Bin);
        fclose(fWavMovie_Hdr);
    }

    free(sigma);
    free(sigX);
    free(sigZ);
    free(wav1);
    free(wav2);
    free(lap);
    free(ws);
    free(sigmaInv);

return;
}

void stripBoundary(float *wholeArray, float *strippedArray, int nx, int nz, int mode) {
    int ix, iz, i, i_;
    for(ix=nxb; ix<nx-nxb; ix++) 
        for(iz=nzb; iz<nz-nzb; iz++)
        {
            i  = iz + ix*nz;
            i_ = (iz-nzb) + (ix-nxb) * (nz-2*nzb);
            if(mode==0)
                strippedArray[i_] = wholeArray[i];
            if(mode==1)
                wholeArray[i] = strippedArray[i_]; 
        }
}
void stripBoundary2(float *wholeArray, float *strippedArray, int nx, int nz, int mode) {
    int iz, shift_whl, shift_str;
    size_t nn = nx-2*nxb;
    for(iz=nzb; iz<nz-nzb; iz++) 
    {
        // shift_whl = ix*nz + nzb;
        // shift_str = (ix-nxb)*(nz-2*nzb);
        shift_whl = iz*nx + nxb;
        shift_str = (iz-nzb)*(nx-2*nxb);
        if(mode==0)
            memcpy(strippedArray+shift_str, wholeArray+shift_whl   , nn);
        if(mode==1)
            memcpy(wholeArray+shift_whl   , strippedArray+shift_str, nn);
    }
}
void stripProp(float *wholeArray, float *strippedArray, int nx, int nz, int ixMin, int ixMax, int mode) {
    int iz, shift_whl, shift_str;
    size_t nxProp = ixMax-ixMin+1;
    for(iz=nzb; iz<nz-nzb; iz++) 
    {
        shift_whl = iz*nx + ixMin;
        shift_str = (iz-nzb) * nxProp;
        if(mode==0)    memcpy(strippedArray+shift_str, wholeArray+shift_whl   , nxProp);
        if(mode==1)    memcpy(wholeArray+shift_whl   , strippedArray+shift_str, nxProp);
    }
}

void emigrateShot_lx_lz_lt_GPU(int gpu_device, float *image, float *eimage_shot, float *eimage_shot_theta, \
                            float *eimage, float *vp, Shot *shot, \
                            int jt, int nt, int nx, int nz, float dx, float dz, float dt, \
                            float ox, float oz, float ot, int lx0, int lz0, int lt0, int ntheta, float dtheta, float otheta, \
                            int recordMovie, int ishotMovie) {

    int it;
    float *wav1, *wav2, *lap, *sigX, *sigZ, *sigma, *sigmaInv, *ws, *wstp;
    wav1 = CPU_zaloc1F(nx*nz);
    wav2 = CPU_zaloc1F(nx*nz);
    lap  = CPU_zaloc1F(nx*nz);
    ws   = CPU_zaloc1F(nx*nz*(nt/jt));
    wstp = CPU_zaloc1F((nx-2*nxb)*(nz-2*nzb));
    sigX = CPU_zaloc1F(nx);
    sigZ = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);

    int ixMin = shot->ixMin;
    int ixMax = shot->ixMax;
    int nxProp = ixMax - ixMin + 1;
    
    // Initializing GPU environment for current shot
    GPU_OPER_initEnvironment(gpu_device, nx, nz, nt, dx, dz, dt, shot->nsou, shot->nrec, lx0, lz0, lt0, nxb, nzb, MODE_CIGS_FULL);
    GPU_OPER_alocArrays_acoustic();
    GPU_OPER_getVelP_CPU2GPU(vp);
    GPU_OPER_getSigmas_CPU2GPU(sigma, sigmaInv);
    GPU_OPER_getSourceWavelet_CPU2GPU(shot->source);
    GPU_OPER_getRecordedData_CPU2GPU(shot->seismogram);
    GPU_OPER_getSourceReceiverIndexes_CPU2GPU(shot->ixs, shot->izs, shot->ixr, shot->izr);

    // Source wavefield
    FILE *fWavMovie_Hdr, *fWavMovie_Bin;
    FILE *fWS_Bin, *fWS_Hdr;
    FILE *fWR_Bin, *fWR_Hdr;
    FILE *fIM_Bin, *fIM_Hdr;
    char movieSouFileName[1024];
    char WSFileName[1024];
    if(shot->ishot==ishotMovie && recordMovie) {
        sprintf(movieSouFileName, "movie_WS_ishot%d", ishotMovie);
        initFilesSmart3d(movieSouFileName, &fWavMovie_Hdr, &fWavMovie_Bin, nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
        // initFilesSmart3d("./movie_WS", &fWavMovie_Hdr, &fWavMovie_Bin, nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
    }
    sprintf(WSFileName, "WS_imaging");
    initFilesSmart3d(WSFileName, &fWS_Hdr, &fWS_Bin, nz, dz, oz, nx, dx, ox, nt, dt, ot);
    // initFilesSmart3d("WS_imaging", &fWS_Hdr, &fWS_Bin, nz, dz, oz, nx, dx, ox, nt, dt, ot);

    // float *wstp = CPU_zaloc1F((nx-2*nxb)*(nz-2*nzb));
    for(it=0; it<nt; it++)
    {
        GPU_OPER_setTimeStep(it);
        if(it%jt==0) {
            // GPU_OPER_copyPresentWavefield_GPU2CPU(wav2);
            GPU_OPER_copyPresentWavefield_secDer_GPU2CPU(wav2);
            // fwrite(wav2, sizeof(float), nx*nz, fWS_Bin);
            // fwrite(wav2+(nxb*nz), sizeof(float), nz*(nx-2*nxb), fWS_Bin);

            // This is the good one
            fwrite(wav2+(nzb*nx), sizeof(float), nx*(nz-2*nzb), fWS_Bin);

            // stripProp(wav2, wstp, nx, nz, ixMin, ixMax, 0);
            // fwrite(wstp, sizeof(float), nxProp*(nz-2*nzb), fWS_Bin);

            // stripBoundary2(wav2, wstp, nx, nz, 0);
            // fwrite(wstp, sizeof(float), (nx-2*nxb)*(nz-2*nzb), fWS_Bin);
        }
        if(shot->ishot==ishotMovie && recordMovie) {
            if(it%jt==0)    writeSamples(wav2, fWavMovie_Bin, nx*nz);
        }
        GPU_OPER_wavefieldTimeSecondDerivative_part1(0);
        GPU_WRAP_injectSource_Acoustic();
        GPU_WRAP_propagate_Acoustic_abs2();
        GPU_OPER_wavefieldTimeSecondDerivative_part2(0);
    }
    if(shot->ishot==ishotMovie && recordMovie) {
        fclose(fWavMovie_Bin);
        fclose(fWavMovie_Hdr);
    }

    GPU_OPER_zeroWavefieldArrays();

    // Receiver wavefield
    if(shot->ishot==ishotMovie && recordMovie) {
        initFilesSmart3d("movie_WR", &fWavMovie_Hdr, &fWavMovie_Bin, nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, 0.0f);
        initFilesSmart3d("movie_IM",       &fIM_Hdr,       &fIM_Bin, nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, 0.0f);
    }
    for(it=nt-1; it>=0; it--)
    {
        GPU_OPER_setTimeStep(it);
        GPU_WRAP_recordInjectData(1);
        if(it%jt==0) {
            // fseek(fWS_Bin, it/jt*nx*nz*sizeof(float), 0);
            // size_t err = fread(ws, sizeof(float), nx*nz, fWS_Bin);

            // fseek(fWS_Bin, it/jt*nz*(nx-2*nxb)*sizeof(float), 0);
            // fread(ws+(nxb*nz), sizeof(float), nz*(nx-2*nxb), fWS_Bin);

            // This is the good one
            fseek(fWS_Bin, (it/jt)*nx*(nz-2*nzb)*sizeof(float), 0);
            fread(ws+(nzb*nx), sizeof(float), nx*(nz-2*nzb), fWS_Bin);

            // fseek(fWS_Bin, (it/jt)*(nz-2*nzb)*nxProp*sizeof(float), 0);
            // fread(wstp, sizeof(float), nxProp*(nz-2*nzb), fWS_Bin);
            // stripProp(ws, wstp, nx, nz, ixMin, ixMax, 1);

            // fseek(fWS_Bin, (it/jt)*(nz-2*nzb)*(nx-2*nxb)*sizeof(float), 0);
            // fread(wstp, sizeof(float), (nx-2*nxb)*(nz-2*nzb), fWS_Bin);
            // stripBoundary2(ws, wstp, nx, nz, 1);

            GPU_OPER_copyPresentWavefieldWS_CPU2GPU(ws);
            GPU_WRAP_eimagingCondition_lx();

        }
        if(shot->ishot==ishotMovie && recordMovie) {
            if(it%jt==0) {
                writeSamples(wav2, fWavMovie_Bin, nx*nz);
                writeSamples(image,      fIM_Bin, nx*nz);
            }
        }
        GPU_WRAP_propagate_Acoustic_abs2();
    }
    GPU_OPER_copyEImage_GPU2CPU(eimage_shot);
    GPU_OPER_freeArrays_acoustic();

    if(eimage_shot_theta!=NULL) {
        // lx2theta_forward(1, eimage_shot, eimage_shot_theta, nx, dx, nz,  dz, lx0, ntheta, dtheta, otheta);
        GPU_WRAP_lx2theta(0, gpu_device, 0, eimage_shot, eimage_shot_theta, nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);
    }

    accumulate(eimage, eimage_shot, nx*nz*(2*lx0+1));

    fclose(fWS_Bin);
    fclose(fWS_Hdr);
    if(shot->ishot==ishotMovie && recordMovie) {
        fclose(fIM_Bin);
        fclose(fIM_Hdr);
        fclose(fWavMovie_Bin);
        fclose(fWavMovie_Hdr);
        remove(movieSouFileName);
    }
    remove(WSFileName);

    transposeArray(image, nx, nz);

    free(sigma);
    free(sigX);
    free(sigZ);
    free(wav1);
    free(wav2);
    free(wstp);
    free(lap);
    free(ws);
    free(sigmaInv);

return;
}

//*
void emigrateShot_TLCIG_Mem_GPU(int gpu_device, float *image, float *ilumin, float *TLcig, float *vp, Shot *shot, 
                                long long int jt, long long int nt, long long int nx, long long int nz, float dx, float dz, float dt, \
                                float ox, float oz, float ot, long long int lt0, int recordMovie, int ishotMovie) {
    long long int it;
    long long int jt_movie = jt/5;
    float *wav2, *wav3, *sigX, *sigZ, *sigma, *sigmaInv, *image_shot, *ilumin_shot;
    image_shot  = CPU_zaloc1F(nx*nz);
    ilumin_shot = CPU_zaloc1F(nx*nz);
    wav2 = CPU_zaloc1F(nx*nz*(1+(nt/jt_movie)));
    wav3 = CPU_zaloc1F(nx*nz*(1+(nt/jt_movie)));
    sigX = CPU_zaloc1F(nx);
    sigZ = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);

    // Initializing GPU environment for current shot
    
    GPU_OPER_zeroArrays_acoustic();
    GPU_OPER_getVelP_CPU2GPU(vp);
    GPU_OPER_getSigmas_CPU2GPU(sigma, sigmaInv);
    GPU_OPER_getSourceWavelet_CPU2GPU(shot->source);
    GPU_OPER_getRecordedData_CPU2GPU(shot->seismogram);
    GPU_OPER_getSourceReceiverIndexes_CPU2GPU(shot->ixs, shot->izs, shot->ixr, shot->izr);

    int mode, side;
    int store_inc = 1, store_sct = 0;
    GPU_OPER_alocArrays_EIC(jt, store_inc, store_sct);

    // Source wavefield
    for(it=0; it<nt; it++)
    {
        GPU_OPER_setTimeStep(it);
        if(it%jt_movie==0)    GPU_WRAP_imagingCondition(IC_REGULAR_ILUMIN);
        // IC_BORNSCT_IMAGE2 stores time second derivative of wavefield
        // IC_BORNSCT_IMAGE  stores wavefield itself
        // IC_SOU_WAVE stores as source side wavefield
        // IC_REC_WAVE stores as receiver side wavefield
        GPU_OPER_EIC_storeWavefields_GPU2GPU(it, IC_SOU_WAVE, IC_BORNSCT_IMAGE2);

        GPU_OPER_wavefieldTimeSecondDerivative_part1(0);
        GPU_WRAP_injectSource_Acoustic();
        GPU_WRAP_propagate_Acoustic_abs2();
        GPU_OPER_wavefieldTimeSecondDerivative_part2(0);
    }
    GPU_OPER_zeroWavefieldArrays();
    
    char SourceWavefield_FileName[1024];
    sprintf(SourceWavefield_FileName, "Source_Wavefield_ishot%d", ishotMovie);
    if(recordMovie)    outputSmart3d(SourceWavefield_FileName, wav2, nz, dz, oz, nx, dx, ox, 1+(nt/jt_movie), dt*jt_movie, 0);

    // Receiver wavefield
    for(it=nt-1; it>=0; it--)
    {        
        GPU_OPER_setTimeStep(it);
        GPU_WRAP_recordInjectData(1);
        if(it%jt_movie==0)
        {
            long long int shift = nx*nz*(it/jt_movie);
            if(recordMovie)    GPU_OPER_copyPresentWavefield_GPU2CPU(wav3+shift);

            GPU_OPER_EIC_storeWavefields_GPU2GPU(it, IC_REC_WAVE, IC_BORNSCT_IMAGE2);
            GPU_WRAP_imagingCondition_smart(IC_REGULAR_IMAGE, it);
        }
        GPU_WRAP_propagate_Acoustic_abs2();
    }

    char ReceiverWavefield_FileName[1024];
    sprintf(ReceiverWavefield_FileName, "Receiver_Wavefield_ishot%d", ishotMovie);
    if(recordMovie)    outputSmart3d(ReceiverWavefield_FileName, wav3, nz, dz, oz, nx, dx, ox, 1+(nt/jt_movie), dt*jt_movie, 0);
    
    GPU_WRAP_eimagingCondition_lt();

    GPU_OPER_freeArrays_EIC();

    free(sigma);
    free(sigX);
    free(sigZ);
    free(wav2);
    free(wav3);
    free(sigmaInv);
    free(image_shot);
    free(ilumin_shot);

return;
}
//*/
     
void emigrateShot_TAU_Mem_GPU(int gpu_device, float *image, float *ilumin, float *Tcig, float *vp, \
                              Shot3D *shot, MSH *msh, Eimage_MSH *eimsh, int recordMovie, int ishotMovie) {
    
    long long int nxb  = msh->nxb;
    long long int nyb  = msh->nyb;
    long long int nzb  = msh->nzb;
    long long int nx   = msh->nx;
    long long int ny   = msh->ny;
    long long int nz   = msh->nz;
    long long int nt   = msh->nt;
    float         dx   = msh->dx;
    float         dy   = msh->dy;
    float         dz   = msh->dz;
    float         dt   = msh->dt;
    float         ox   = msh->ox;
    float         oy   = msh->oy;
    float         oz   = msh->oz;
    float         ot   = msh->ot;
    long long int jt   = msh->jt;
    long long int jtmv = msh->jt * 10;
    
    long long int it;
    long long int lt0 = 0;

    float *wav2, *wav3, *sigX, *sigZ, *sigma, *sigmaInv;
    wav2 = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
    wav3 = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
    sigX = CPU_zaloc1F(nx);
    sigZ = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);

    // Initializing GPU environment for current shot
    GPU_OPER_initEnvironment(gpu_device, nx, nz, nt, dx, dz, dt, shot->nsou, shot->nrec, eimsh->lx0, eimsh->lz0, lt0, nxb, nzb, MODE_CIGS_TIME);
    GPU_OPER_initTAU(eimsh->tau0, eimsh->dtau);
    GPU_OPER_alocArrays_acoustic();
    GPU_OPER_getVelP_CPU2GPU(vp);
    GPU_OPER_getSigmas_CPU2GPU(sigma, sigmaInv);
    GPU_OPER_getSourceWavelet_CPU2GPU(shot->source);
    GPU_OPER_getRecordedData_CPU2GPU(shot->seismogram);
    GPU_OPER_getSourceReceiverIndexes_CPU2GPU(shot->ixs, shot->izs, shot->ixr, shot->izr);

    int mode, side;
    int store_inc = 1, store_sct = 0;
    GPU_OPER_alocArrays_EIC(jt, store_inc, store_sct);

    // Source wavefield
    for(it=0; it<nt; it++)
    {
        GPU_OPER_setTimeStep(it);
        if(it%jt==0)
        {
            GPU_WRAP_imagingCondition(IC_REGULAR_ILUMIN);
            if(recordMovie)
            {
                long long int shift = nx*nz*(it/jtmv);
                GPU_OPER_copyPresentWavefield_secDer_GPU2CPU(wav2+shift);
                // GPU_OPER_copyPresentWavefield_GPU2CPU(wav2+shift);
            }
        }
        // IC_BORNSCT_IMAGE2 stores time second derivative of wavefield
        // IC_BORNSCT_IMAGE  stores wavefield itself
        // IC_SOU_WAVE stores as source side wavefield
        // IC_REC_WAVE stores as receiver side wavefield
        GPU_OPER_EIC_storeWavefields_GPU2GPU(it, IC_SOU_WAVE, IC_BORNSCT_IMAGE2);

        GPU_OPER_wavefieldTimeSecondDerivative_part1(0);
        GPU_WRAP_injectSource_Acoustic();
        GPU_WRAP_propagate_Acoustic_abs2();
        GPU_OPER_wavefieldTimeSecondDerivative_part2(0);
    }
    GPU_OPER_zeroWavefieldArrays();
    char SourceWavefield_FileName[1024];
    sprintf(SourceWavefield_FileName, "Source_Wavefield_ishot%d", ishotMovie);
    if(recordMovie)    outputSmart3d(SourceWavefield_FileName, wav2, nx, dx, ox, nz, dz, oz, (nt/jtmv), dt*jtmv, 0);

    // Receiver wavefield
    for(it=nt-1; it>=0; it--)
    {        
        GPU_OPER_setTimeStep(it);
        if(it%jt==0)
        {
            if(recordMovie &&  it%jtmv==0) {
                long long int shift = nx*nz*(it/jtmv);
                GPU_OPER_copyPresentWavefield_GPU2CPU(wav3+shift);
            }
            // GPU_WRAP_eimagingCondition_OCIGS_smart(EIMSH_CIGS_VCIG, it);
            GPU_WRAP_imagingCondition_smart(IC_REGULAR_IMAGE, it);
        }
        // IC_BORNSCT_IMAGE2 stores time second derivative of wavefield
        // IC_BORNSCT_IMAGE  stores wavefield itself
        // IC_SOU_WAVE stores as source side wavefield
        // IC_REC_WAVE stores as receiver side wavefield
        GPU_OPER_EIC_storeWavefields_GPU2GPU(it, IC_REC_WAVE, IC_BORNSCT_IMAGE);
        GPU_WRAP_recordInjectData(1);
        GPU_WRAP_propagate_Acoustic_abs2();
    }

    // Apply extended imaging condition
    GPU_WRAP_eimagingCondition_TAU();

    char ReceiverWavefield_FileName[1024];
    sprintf(ReceiverWavefield_FileName, "Receiver_Wavefield_ishot%d", ishotMovie);
    if(recordMovie)    outputSmart3d(ReceiverWavefield_FileName, wav3, nx, dx, ox, nz, dz, oz, (nt/jtmv), dt*jtmv, 0);
    
    GPU_OPER_copyTCIG_GPU2CPU(Tcig);
    GPU_OPER_copyIlumin_GPU2CPU(ilumin);
    GPU_OPER_copyImage_GPU2CPU(image);
    GPU_OPER_freeArrays_acoustic();
    GPU_OPER_freeArrays_EIC();

    free(sigma);
    free(sigX);
    free(sigZ);
    free(wav2);
    free(wav3);
    free(sigmaInv);

return;
}

void emigrateShot_HOCIG_VOCIG_Mem_GPU(int gpu_device, float *image, float *ilumin, float *Hocig, float *Vocig, \
                                      float *vp, Shot3D *shot, MSH *msh, Eimage_MSH *eimsh, int recordMovie, int ishotMovie) {
    
    long long int nxb  = msh->nxb;
    long long int nyb  = msh->nyb;
    long long int nzb  = msh->nzb;
    long long int nx   = msh->nx;
    long long int ny   = msh->ny;
    long long int nz   = msh->nz;
    long long int nt   = msh->nt;
    float         dx   = msh->dx;
    float         dy   = msh->dy;
    float         dz   = msh->dz;
    float         dt   = msh->dt;
    float         ox   = msh->ox;
    float         oy   = msh->oy;
    float         oz   = msh->oz;
    float         ot   = msh->ot;
    long long int jt   = msh->jt;
    long long int jtmv = msh->jt * 10;
    
    long long int it;
    long long int lt0 = 0;

    float *wav2, *wav3, *sigX, *sigZ, *sigma, *sigmaInv;
    wav2 = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
    wav3 = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
    sigX = CPU_zaloc1F(nx);
    sigZ = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);

    // Initializing GPU environment for current shot
    GPU_OPER_initEnvironment(gpu_device, nx, nz, nt, dx, dz, dt, shot->nsou, shot->nrec, eimsh->lx0, eimsh->lz0, lt0, nxb, nzb, MODE_CIGS_ORTH);
    GPU_OPER_alocArrays_acoustic();
    GPU_OPER_getVelP_CPU2GPU(vp);
    GPU_OPER_getSigmas_CPU2GPU(sigma, sigmaInv);
    GPU_OPER_getSourceWavelet_CPU2GPU(shot->source);
    GPU_OPER_getRecordedData_CPU2GPU(shot->seismogram);
    GPU_OPER_getSourceReceiverIndexes_CPU2GPU(shot->ixs, shot->izs, shot->ixr, shot->izr);

    int mode, side;
    int store_inc = 1, store_sct = 0;
    GPU_OPER_alocArrays_EIC(jt, store_inc, store_sct);

    // Source wavefield
    for(it=0; it<nt; it++)
    {
        GPU_OPER_setTimeStep(it);
        if(it%jt==0)
        {
            GPU_WRAP_imagingCondition(IC_REGULAR_ILUMIN);
            if(recordMovie)
            {
                long long int shift = nx*nz*(it/jtmv);
                GPU_OPER_copyPresentWavefield_secDer_GPU2CPU(wav2+shift);
                // GPU_OPER_copyPresentWavefield_GPU2CPU(wav2+shift);
            }
        }
        // IC_BORNSCT_IMAGE2 stores time second derivative of wavefield
        // IC_BORNSCT_IMAGE  stores wavefield itself
        // IC_SOU_WAVE stores as source side wavefield
        // IC_REC_WAVE stores as receiver side wavefield
        GPU_OPER_EIC_storeWavefields_GPU2GPU(it, IC_SOU_WAVE, IC_BORNSCT_IMAGE2);

        GPU_OPER_wavefieldTimeSecondDerivative_part1(0);
        GPU_WRAP_injectSource_Acoustic();
        GPU_WRAP_propagate_Acoustic_abs2();
        GPU_OPER_wavefieldTimeSecondDerivative_part2(0);
    }
    GPU_OPER_zeroWavefieldArrays();
    char SourceWavefield_FileName[1024];
    sprintf(SourceWavefield_FileName, "Source_Wavefield_ishot%d", ishotMovie);
    if(recordMovie)    outputSmart3d(SourceWavefield_FileName, wav2, nx, dx, ox, nz, dz, oz, (nt/jtmv), dt*jtmv, 0);

    // Receiver wavefield
    for(it=nt-1; it>=0; it--)
    {        
        GPU_OPER_setTimeStep(it);
        GPU_WRAP_recordInjectData(1);
        if(it%jt==0)
        {
            if(recordMovie &&  it%jtmv==0) {
                long long int shift = nx*nz*(it/jtmv);
                GPU_OPER_copyPresentWavefield_GPU2CPU(wav3+shift);
            }
            GPU_WRAP_eimagingCondition_OCIGS_smart(EIMSH_CIGS_VCIG, it);
            GPU_WRAP_imagingCondition_smart(IC_REGULAR_IMAGE, it);
        }
        GPU_WRAP_propagate_Acoustic_abs2();
    }

    char ReceiverWavefield_FileName[1024];
    sprintf(ReceiverWavefield_FileName, "Receiver_Wavefield_ishot%d", ishotMovie);
    if(recordMovie)    outputSmart3d(ReceiverWavefield_FileName, wav3, nx, dx, ox, nz, dz, oz, (nt/jtmv), dt*jtmv, 0);
    
    GPU_OPER_copyOCIGS_GPU2CPU(Hocig, Vocig);
    GPU_OPER_copyIlumin_GPU2CPU(ilumin);
    GPU_OPER_copyImage_GPU2CPU(image);
    GPU_OPER_freeArrays_acoustic();
    GPU_OPER_freeArrays_EIC();

    free(sigma);
    free(sigX);
    free(sigZ);
    free(wav2);
    free(wav3);
    free(sigmaInv);
return;
}

void emigrateShot_lx_lz_lt_Mem_GPU(int gpu_device, float *image, float *ilumin, float *eimage_shot, \
                                   float *eimage_shot_theta, float *eimage, float *vp, Shot *shot, int jt, \
                                   int nt, int nx, int nz, float dx, float dz, float dt, float ox, float oz, \
                                   float ot, int lx0, int lz0, int lt0, int jt_lt, int ntheta, float dtheta, float otheta, \
                                   int recordMovie, int ishotMovie) {
    int it;
    float *wav2, *sigX, *sigZ, *sigma, *sigmaInv, *image_shot, *ilumin_shot;
    image_shot  = CPU_zaloc1F(nx*nz);
    ilumin_shot = CPU_zaloc1F(nx*nz);
    wav2 = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
    sigX = CPU_zaloc1F(nx);
    sigZ = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);

    // Initializing GPU environment for current shot
    GPU_OPER_initEnvironment(gpu_device, nx, nz, nt, dx, dz, dt, shot->nsou, shot->nrec, lx0, lz0, lt0, nxb, nzb, MODE_CIGS_FULL);
    GPU_OPER_alocArrays_acoustic();
    GPU_OPER_getVelP_CPU2GPU(vp);
    GPU_OPER_getSigmas_CPU2GPU(sigma, sigmaInv);
    GPU_OPER_getSourceWavelet_CPU2GPU(shot->source);
    GPU_OPER_getRecordedData_CPU2GPU(shot->seismogram);
    GPU_OPER_getSourceReceiverIndexes_CPU2GPU(shot->ixs, shot->izs, shot->ixr, shot->izr);

    int mode, side;
    int store_inc=1, store_sct=0;
    mode = IC_BORNSCT_IMAGE;
    GPU_OPER_alocArrays_EIC(jt_lt, store_inc, store_sct);

    // Source wavefield
    for(it=0; it<nt; it++)
    {
        GPU_OPER_setTimeStep(it);
        if(it%jt==0)
        {
            long long int shift = nx*nz*(it/jt);
            // GPU_OPER_copyPresentWavefield_GPU2CPU(wav2+shift);
            GPU_OPER_copyPresentWavefield_secDer_GPU2CPU(wav2+shift);
            GPU_WRAP_imagingCondition(IC_REGULAR_ILUMIN);
        }
        side = IC_SOU_WAVE;
        GPU_OPER_EIC_storeWavefields_GPU2GPU(it, side, mode);

        GPU_OPER_wavefieldTimeSecondDerivative_part1(0);
        GPU_WRAP_injectSource_Acoustic();
        GPU_WRAP_propagate_Acoustic_abs2();
        GPU_OPER_wavefieldTimeSecondDerivative_part2(0);
    }

    GPU_OPER_zeroWavefieldArrays();

    // Receiver wavefield
    for(it=nt-1; it>=0; it--)
    {
        GPU_OPER_setTimeStep(it);
        GPU_WRAP_recordInjectData(1);
        if(it%jt==0  &&  jt_lt<=0)
        {
            long long int shift = nx*nz*(it/jt);
            GPU_OPER_copyPresentWavefieldWS_CPU2GPU(wav2+shift);
            if(lz0==0)    GPU_WRAP_eimagingCondition_lx();
            else          GPU_WRAP_eimagingCondition_lx_lz();
            GPU_WRAP_imagingCondition(IC_REGULAR_IMAGE);
        }
        side = IC_REC_WAVE;
        GPU_OPER_EIC_storeWavefields_GPU2GPU(it, side, mode);
        GPU_WRAP_propagate_Acoustic_abs2();
    }
    GPU_WRAP_eimagingCondition_lx_lz_lt();

    GPU_OPER_copyEImage_GPU2CPU(eimage_shot);
    GPU_OPER_copyIlumin_GPU2CPU(ilumin_shot);
    GPU_OPER_copyImage_GPU2CPU(image_shot);
    GPU_OPER_freeArrays_acoustic();
    GPU_OPER_freeArrays_EIC();

    if(eimage_shot_theta!=NULL)
    {
        // lx2theta_forward(1, eimage_shot, eimage_shot_theta, nx, dx, nz,  dz, lx0, ntheta, dtheta, otheta);
        GPU_WRAP_lx2theta(0, gpu_device, 0, eimage_shot, eimage_shot_theta, nx, dx, nz,  dz, nxb, nzb, lx0, ntheta, dtheta, otheta);
    }

    accumulate(eimage, eimage_shot, nx*nz*(2*lx0+1)*(2*lz0+1)*(2*lt0+1));
    accumulate(ilumin, ilumin_shot, nx*nz);
    accumulate(image, image_shot, nx*nz);

    // transposeArray(image, nx, nz);

    free(sigma);
    free(sigX);
    free(sigZ);
    free(wav2);
    free(sigmaInv);
    free(image_shot);
    free(ilumin_shot);

return;
}

void modelShotBorn(float *ref, float *vp, Shot *shot, int jt, int nt, int nx, int nz, float dx, float dz, \
                   float dt, float ox, float oz, float ot) {
    int it;
    float *wav1, *wav2;
    float *sct1, *sct2;
    float *lap, *sigX, *sigZ, *sigma, *sigmaInv, *tmp;
    wav1 = CPU_zaloc1F(nx*nz);
    wav2 = CPU_zaloc1F(nx*nz);
    sct1 = CPU_zaloc1F(nx*nz);
    sct2 = CPU_zaloc1F(nx*nz);
    lap  = CPU_zaloc1F(nx*nz);
    sigX = CPU_zaloc1F(nx);
    sigZ = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);

    //applyHalfDerivative2Source(shot);

    // Source wavefield
    //int recordBornModeling = (shot->ishot==(shot->nshot/2));
    int recordBornModeling = (shot->ishot==0);
    FILE *fwav_Hdr, *fwav_Bin;
    FILE *fsct_Hdr, *fsct_Bin;
    if(recordBornModeling) {
        initFilesSmart3d("movie_incident" , &fwav_Hdr, &fwav_Bin, nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
        initFilesSmart3d("movie_scattered", &fsct_Hdr, &fsct_Bin, nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
    }
    for(it=0; it<nt; it++)
    {
        recordData(shot, sct2, nx, nz, it);

        if(it%jt==0 && recordBornModeling) {
            writeSamples(wav2, fwav_Bin, nx*nz);
            writeSamples(sct2, fsct_Bin, nx*nz);
        }

        injectSource_Acoustic(wav1, shot->source, it, (shot->ixs)[0], (shot->izs)[0], nx, nz, dx, dz);
        scatter_Acoustic(ref, sct1, 2.0f, wav2, -1.0f, wav1, nx, nz, dt);
        propagate_Acoustic_abs(wav1, wav2, lap, vp, sigma, nx, nz);
        scatter_Acoustic(ref, sct1, -1.0f, wav2, 0.0f, wav1, nx, nz, dt);;
        propagate_Acoustic_abs(sct1, sct2, lap, vp, sigma, nx, nz);
        tmp = wav2;
        wav2=wav1;
        wav1=tmp;
        tmp = sct2;
        sct2=sct1;
        sct1=tmp;
    }
    if(recordBornModeling) {
        fclose(fwav_Bin);
        fclose(fwav_Hdr);
        fclose(fsct_Bin);
        fclose(fsct_Hdr);
    }
    free(sigma);
    free(sigX);
    free(sigZ);
    free(wav1);
    free(wav2);
    free(sct1);
    free(sct2);
    free(lap);
    free(sigmaInv);

return;
}
void migrateExtBorn_lx(float *image1, float *image2, float *ilumin1, float *ilumin2, float *eref, float *vp, Shot *shot, int jt, int nt, int nx, int nz, float dx, float dz, \
                       float dt, float ox, float oz, float ot, int lx0, int ishotMovie, int recordMovie) {
    int it;
    size_t err;
    float *wav1, *wav2;
    float *sct1, *sct2;
    float *lap, *sigX, *sigZ, *sigma, *sigmaInv, *tmp, *ws, *wc, *secDer;
    secDer = CPU_zaloc1F(nx*nz);
    wav1   = CPU_zaloc1F(nx*nz);
    wav2   = CPU_zaloc1F(nx*nz);
    ws     = CPU_zaloc1F(nx*nz);
    wc     = CPU_zaloc1F(nx*nz);
    sct1   = CPU_zaloc1F(nx*nz);
    sct2   = CPU_zaloc1F(nx*nz);
    lap    = CPU_zaloc1F(nx*nz);
    sigX   = CPU_zaloc1F(nx);
    sigZ   = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);

    // applyHalfDerivative2Source(shot);
    // applyHalfDerivative2Data(shot);

    // Source wavefield
    //int recordBornModeling = (shot->ishot==(shot->nshot/2));
    int recordBornModeling = (shot->ishot==0);
    FILE *fwav_Hdr, *fwav_Bin;
    FILE *fWavRecInc_Bin, *fWavRecInc_Hdr;
    FILE *fWavRecSct_Bin, *fWavRecSct_Hdr;
    FILE *fsct_Hdr, *fsct_Bin;
    FILE *fWS_Bin, *fWS_Hdr;
    FILE *fWC_Bin, *fWC_Hdr;
    if(recordMovie && shot->ishot==ishotMovie) {
        initFilesSmart3d("movie_incident_extendedBorn" , &fwav_Hdr, &fwav_Bin, nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
        initFilesSmart3d("movie_scattered_extendedBorn", &fsct_Hdr, &fsct_Bin, nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
        initFilesSmart3d("movie_incident_extendedBorn_receiver" , &fWavRecInc_Hdr, &fWavRecInc_Bin, nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
        initFilesSmart3d("movie_scattered_extendedBorn_receiver", &fWavRecSct_Hdr, &fWavRecSct_Bin, nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
    }
    initFilesSmart3d("WS_imaging", &fWS_Hdr, &fWS_Bin, nz, dz, oz, nx, dx, ox, nt   , dt   , ot);
    initFilesSmart3d("WC_imaging", &fWC_Hdr, &fWC_Bin, nz, dz, oz, nx, dx, ox, nt   , dt   , ot);
    
    // applyHalfDerivative2Data(shot);

    // Forward
    for(it=0; it<nt; it++)
    {
        fwrite(wav2, sizeof(float), nx*nz, fWS_Bin);
        fwrite(sct2, sizeof(float), nx*nz, fWC_Bin);

        if(it%jt==0 && recordMovie && shot->ishot==ishotMovie)
        {
            writeSamples(wav2, fwav_Bin, nx*nz);
            writeSamples(sct2, fsct_Bin, nx*nz);
        }
        /*
        injectSource_Acoustic(wav1, shot->source, it, (shot->ixs)[0], (shot->izs)[0], nx, nz, dx, dz);
        scatter_Acoustic_lx(eref, sct1, 2.0f, wav2, -1.0f, wav1, nx, nz, dt, lx0);
        propagate_Acoustic_abs2(wav1, wav2, lap, vp, sigma, sigmaInv, nx, nz);
        scatter_Acoustic_lx(eref, sct1, -1.0f, wav2, 0.0f, wav1, nx, nz, dt, lx0);
        propagate_Acoustic_abs2(sct1, sct2, lap, vp, sigma, sigmaInv, nx, nz);
        */
        zeroArray(secDer,nx*nz);
        injectSource_Acoustic(wav1, shot->source, it, (shot->ixs)[0], (shot->izs)[0], nx, nz, dx, dz);
        wavefieldTimeSecondDerivative(secDer, -2.0f, wav2, +1.0f, wav1, nx, nz, dt);
        propagate_Acoustic_abs2(wav1, wav2, lap, vp, sigma, sigmaInv, nx, nz);
        wavefieldTimeSecondDerivative(secDer, +1.0f, wav2,  0.0f, wav1, nx, nz, dt);

        scatter_Acoustic_lx_opt(eref, sct1, secDer, nx, nz, dt, lx0);
        propagate_Acoustic_abs2(sct1, sct2, lap, vp, sigma, sigmaInv, nx, nz);

        if(it%jt==0) {
            imagingCondition(ilumin1, wav2, wav2, nx, nz);
            imagingCondition(ilumin2, sct2, sct2, nx, nz);
        }

        tmp  = wav2;
        wav2 = wav1;
        wav1 = tmp;
        tmp  = sct2;
        sct2 = sct1;
        sct1 = tmp;
    }

    zeroArray(wav2,nx*nz);
    zeroArray(wav1,nx*nz);
    zeroArray(sct2,nx*nz);
    zeroArray(sct1,nx*nz);

    // Backward
    for(it=nt-1; it>=0; it--)
    {
        injectData(shot, wav1, nx, nz, dx, dz, it);
        fseek(fWS_Bin, it*nx*nz*sizeof(float), 0);
        fseek(fWC_Bin, it*nx*nz*sizeof(float), 0);
        err = fread(ws, sizeof(float), nx*nz, fWS_Bin);
        err = fread(wc, sizeof(float), nx*nz, fWC_Bin);
        if(it%jt==0) {
            imagingCondition(image1, ws, sct2, nx, nz);
            imagingCondition(image2, wc, wav2, nx, nz);
        }
        if(shot->ishot==ishotMovie && recordMovie) {
            if(it%jt==0) {
                writeSamples(wav2, fWavRecInc_Bin, nx*nz);
                writeSamples(sct2, fWavRecSct_Bin, nx*nz);
                // writeSamples(image,      fIM_Bin, nx*nz);
            }
        }
        /*
        scatter_Acoustic_lx(eref, sct1,  2.0f, wav2, -1.0f, wav1, nx, nz, dt, lx0);
        propagate_Acoustic_abs2(wav1, wav2, lap, vp, sigma, sigmaInv, nx, nz);
        scatter_Acoustic_lx(eref, sct1, -1.0f, wav2,  0.0f, wav1, nx, nz, dt, lx0);
        propagate_Acoustic_abs2(sct1, sct2, lap, vp, sigma, sigmaInv, nx, nz);
        */
        zeroArray(secDer,nx*nz);
        wavefieldTimeSecondDerivative(secDer, -2.0f, wav2, +1.0f, wav1, nx, nz, dt);
        propagate_Acoustic_abs2(wav1, wav2, lap, vp, sigma, sigmaInv, nx, nz);
        wavefieldTimeSecondDerivative(secDer, +1.0f, wav2,  0.0f, wav1, nx, nz, dt);

        scatter_Acoustic_lx_opt(eref, sct1, secDer, nx, nz, dt, lx0);
        propagate_Acoustic_abs2(sct1, sct2, lap, vp, sigma, sigmaInv, nx, nz);

        tmp  = wav2;
        wav2 = wav1;
        wav1 = tmp;
        tmp  = sct2;
        sct2 = sct1;
        sct1 = tmp;
    }

    fclose(fWS_Bin);
    fclose(fWS_Hdr);
    if(recordMovie && shot->ishot==ishotMovie) {
        fclose(fwav_Bin);
        fclose(fwav_Hdr);
        fclose(fsct_Bin);
        fclose(fsct_Hdr);
    }

    free(sigma);
    free(sigX);
    free(sigZ);
    free(ws);
    free(wc);
    free(wav1);
    free(wav2);
    free(sct1);
    free(sct2);
    free(lap);
    free(sigmaInv);

return;
}
void migrateExtBorn_lx_lz_lt_GPU(float *imageSouSide, float *imageRecSide, \
                           float *imageSouSide_shot, float *imageRecSide_shot, \
                           float *iluminSouSide, float *iluminRecSide, \
                           float *eref, float *vp, Shot *shot, int jt, int nt, \
                           int nx, int nz, float dx, float dz, \
                           float dt, float ox, float oz, float ot, int lx0, int lz0, int lt0, \
                           int ishotMovie, int recordMovie) {
    int it;
    size_t err;
    float *wav1, *wav2;
    float *sct1, *sct2;
    float *lap, *sigX, *sigZ, *sigma, *sigmaInv, *tmp, *ws, *wc, *secDer, *wstp;
    secDer = CPU_zaloc1F(nx*nz);
    wav1   = CPU_zaloc1F(nx*nz);
    wav2   = CPU_zaloc1F(nx*nz);
    ws     = CPU_zaloc1F(nx*nz);
    wstp = CPU_zaloc1F((nx-2*nxb)*(nz-2*nzb));
    wc     = CPU_zaloc1F(nx*nz);
    sct1   = CPU_zaloc1F(nx*nz);
    sct2   = CPU_zaloc1F(nx*nz);
    lap    = CPU_zaloc1F(nx*nz);
    sigX   = CPU_zaloc1F(nx);
    sigZ   = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);

    int ixMin = shot->ixMin;
    int ixMax = shot->ixMax;
    int nxProp = ixMax - ixMin + 1;

    // Initializing GPU environment for current shot
    GPU_OPER_initEnvironment(gpuDev, nx, nz, nt, dx, dz, dt, shot->nsou, shot->nrec, lx0, lz0, lt0, nxb, nzb, MODE_CIGS_FULL);
    GPU_OPER_alocArrays_acoustic();
    GPU_OPER_alocArrays_acoustic_ExtBornSct();
    GPU_OPER_getVelP_CPU2GPU(vp);
    GPU_OPER_copyEImage_CPU2GPU(eref);
    GPU_OPER_getSigmas_CPU2GPU(sigma, sigmaInv);
    GPU_OPER_getSourceWavelet_CPU2GPU(shot->source);
    GPU_OPER_getRecordedData_CPU2GPU(shot->seismogram);
    GPU_OPER_getSourceReceiverIndexes_CPU2GPU(shot->ixs, shot->izs, shot->ixr, shot->izr);

    // Source wavefield
    int recordBornModeling = ((shot->ishot==ishotMovie) && recordMovie);
    FILE *fwav_Hdr, *fwav_Bin;
    FILE *fWavRecInc_Bin, *fWavRecInc_Hdr;
    FILE *fWavRecSct_Bin, *fWavRecSct_Hdr;
    FILE *fsct_Hdr, *fsct_Bin;
    FILE *fWS_Bin, *fWS_Hdr;
    FILE *fWC_Bin, *fWC_Hdr;
    char movieIncSouFileName[1024];
    char movieSctSouFileName[1024];
    char movieIncRecFileName[1024];
    char movieSctRecFileName[1024];
    char WSFileName[1024];
    char WCFileName[1024];
    if(recordBornModeling)
    {
        sprintf(movieIncSouFileName, "movie_inc_extBorn_sou_ishot%d", ishotMovie);
        sprintf(movieSctSouFileName, "movie_sct_extBorn_sou_ishot%d", ishotMovie);
        sprintf(movieIncRecFileName, "movie_inc_extBorn_ishot%d", ishotMovie);
        sprintf(movieSctRecFileName, "movie_sct_extBorn_ishot%d", ishotMovie);
        initFilesSmart3d(movieIncSouFileName, &fwav_Hdr, &fwav_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieSctSouFileName, &fsct_Hdr, &fsct_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieIncRecFileName, &fWavRecInc_Hdr, &fWavRecInc_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieSctRecFileName, &fWavRecSct_Hdr, &fWavRecSct_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        // initFilesSmart3d("./movie_incident_extendedBorn" , &fwav_Hdr, &fwav_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        // initFilesSmart3d("./movie_scattered_extendedBorn", &fsct_Hdr, &fsct_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        // initFilesSmart3d("./movie_incident_extendedBorn_receiver" , &fWavRecInc_Hdr, &fWavRecInc_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        // initFilesSmart3d("./movie_scattered_extendedBorn_receiver", &fWavRecSct_Hdr, &fWavRecSct_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
    }
    sprintf(WSFileName, "WS_imaging");
    sprintf(WCFileName, "WC_imaging");
    initFilesSmart3d(WSFileName, &fWS_Hdr, &fWS_Bin, nx, dx, ox, nz, dz, oz, nt, dt, ot);
    initFilesSmart3d(WCFileName, &fWC_Hdr, &fWC_Bin, nx, dx, ox, nz, dz, oz, nt, dt, ot);
    // initFilesSmart3d("WS_imaging", &fWS_Hdr, &fWS_Bin, nx, dx, ox, nz, dz, oz, nt, dt, ot);
    // initFilesSmart3d("WC_imaging", &fWC_Hdr, &fWC_Bin, nx, dx, ox, nz, dz, oz, nt, dt, ot);

    // Forward
    for(it=0; it<nt; it++)
    {
        GPU_OPER_setTimeStep(it);

        if(it%jt==0)
        {
            // GPU_OPER_copyPresentWavefield_secDer_GPU2CPU(wav2);
            // GPU_OPER_copyPresentScatteredWavefield_secDer_GPU2CPU(sct2);
            GPU_OPER_copyPresentWavefield_GPU2CPU(wav2);
            GPU_OPER_copyPresentScatteredWavefield_GPU2CPU(sct2);
            // fwrite(wav2, sizeof(float), nx*nz, fWS_Bin);
            // fwrite(sct2, sizeof(float), nx*nz, fWC_Bin);
            // fwrite(wav2+(nxb*nz), sizeof(float), nz*(nx-2*nxb), fWS_Bin);
            // fwrite(sct2+(nxb*nz), sizeof(float), nz*(nx-2*nxb), fWC_Bin);

            // This is the good one
            fwrite(wav2+(nzb*nx), sizeof(float), nx*(nz-2*nzb), fWS_Bin);
            fwrite(sct2+(nzb*nx), sizeof(float), nx*(nz-2*nzb), fWC_Bin);

            // stripProp(wav2, wstp, nx, nz, ixMin, ixMax, 0);
            // fwrite(wstp, sizeof(float), nxProp*(nz-2*nzb), fWS_Bin);
            // stripProp(sct2, wstp, nx, nz, ixMin, ixMax, 0);
            // fwrite(wstp, sizeof(float), nxProp*(nz-2*nzb), fWC_Bin);
            
            // stripBoundary2(wav2, wstp, nx, nz, 0);
            // fwrite(wstp, sizeof(float), (nx-2*nxb)*(nz-2*nzb), fWS_Bin);
            // stripBoundary2(sct2, wstp, nx, nz, 0);
            // fwrite(wstp, sizeof(float), (nx-2*nxb)*(nz-2*nzb), fWC_Bin);

            GPU_WRAP_imagingCondition(IC_BORNSCT_ILUMIN);
        }
        if(it%jt==0 && recordBornModeling)
        {
            writeSamples(wav2, fwav_Bin, nx*nz);
            writeSamples(sct2, fsct_Bin, nx*nz);
        }
        GPU_OPER_wavefieldTimeSecondDerivative_part1(0);
        GPU_OPER_wavefieldTimeSecondDerivative_part1(1);
        GPU_WRAP_injectSource_Acoustic();
        GPU_WRAP_propagate_Acoustic_abs2();
        GPU_OPER_wavefieldTimeSecondDerivative_part2(0);
        GPU_WRAP_scatter_Acoustic_lx_opt();
        GPU_WRAP_propagate_Acoustic_Scattered_abs2();
        GPU_OPER_wavefieldTimeSecondDerivative_part2(1);
    }

    GPU_OPER_zeroWavefieldArrays();
    GPU_OPER_zeroScatteredWavefieldArrays();

    // Backward
    for(it=nt-1; it>=0; it--)
    {
        GPU_OPER_setTimeStep(it);

        if(it%jt==0) {
            // fseek(fWS_Bin, (it/jt)*nx*nz*sizeof(float), 0);
            // fseek(fWC_Bin, (it/jt)*nx*nz*sizeof(float), 0);
            // err = fread(ws, sizeof(float), nx*nz, fWS_Bin);
            // err = fread(wc, sizeof(float), nx*nz, fWC_Bin);

            // fseek(fWS_Bin, (it/jt)*nz*(nx-2*nxb)*sizeof(float), 0);
            // fseek(fWC_Bin, (it/jt)*nz*(nx-2*nxb)*sizeof(float), 0);
            // err = fread(ws+(nxb*nz), sizeof(float), nz*(nx-2*nxb), fWS_Bin);
            // err = fread(wc+(nxb*nz), sizeof(float), nz*(nx-2*nxb), fWC_Bin);

            // This is the good one
            fseek(fWS_Bin, (it/jt)*nx*(nz-2*nzb)*sizeof(float), 0);
            fseek(fWC_Bin, (it/jt)*nx*(nz-2*nzb)*sizeof(float), 0);
            err = fread(ws+(nzb*nx), sizeof(float), nx*(nz-2*nzb), fWS_Bin);
            err = fread(wc+(nzb*nx), sizeof(float), nx*(nz-2*nzb), fWC_Bin);


            // fseek(fWS_Bin, (it/jt)*(nz-2*nzb)*nxProp*sizeof(float), 0);
            // fread(wstp, sizeof(float), nxProp*(nz-2*nzb), fWS_Bin);
            // stripProp(ws, wstp, nx, nz, ixMin, ixMax, 1);
            // fseek(fWC_Bin, (it/jt)*(nz-2*nzb)*nxProp*sizeof(float), 0);
            // fread(wstp, sizeof(float), nxProp*(nz-2*nzb), fWC_Bin);
            // stripProp(wc, wstp, nx, nz, ixMin, ixMax, 1);


            // fseek(fWS_Bin, (it/jt)*(nz-2*nzb)*(nx-2*nxb)*sizeof(float), 0);
            // fread(wstp, sizeof(float), (nx-2*nxb)*(nz-2*nzb), fWS_Bin);
            // stripBoundary2(ws, wstp, nx, nz, 1);
            // fseek(fWC_Bin, (it/jt)*(nz-2*nzb)*(nx-2*nxb)*sizeof(float), 0);
            // fread(wstp, sizeof(float), (nx-2*nxb)*(nz-2*nzb), fWC_Bin);
            // stripBoundary2(wc, wstp, nx, nz, 1);

            GPU_OPER_copyPresentWavefieldWS_CPU2GPU(ws);
            GPU_OPER_copyPresentWavefieldWC_CPU2GPU(wc);
            GPU_WRAP_imagingCondition(IC_BORNSCT_IMAGE);
        }
        if(recordBornModeling && it%jt==0)
        {
            GPU_OPER_copyPresentWavefield_GPU2CPU(wav2);
            GPU_OPER_copyPresentScatteredWavefield_GPU2CPU(sct2);
            // GPU_OPER_copyPresentWavefieldWS_GPU2CPU(wav2);
            // GPU_OPER_copyPresentWavefieldWC_GPU2CPU(sct2);
            // writeSamples(ws, fWavRecInc_Bin, nx*nz);
            // writeSamples(wc, fWavRecSct_Bin, nx*nz);
            writeSamples(wav2, fWavRecInc_Bin, nx*nz);
            writeSamples(sct2, fWavRecSct_Bin, nx*nz);            
            
            
        }
        GPU_OPER_wavefieldTimeSecondDerivative_part1(0);
        GPU_WRAP_recordInjectData(1);
        GPU_WRAP_propagate_Acoustic_abs2();
        GPU_OPER_wavefieldTimeSecondDerivative_part2(0);
        GPU_WRAP_scatter_Acoustic_lx_opt();
        GPU_WRAP_propagate_Acoustic_Scattered_abs2();
    }

    GPU_OPER_copyImage_BRNSCT_GPU2CPU(imageSouSide_shot,  imageRecSide_shot, 0);
    GPU_OPER_copyImage_BRNSCT_GPU2CPU(imageSouSide,  imageRecSide, 1);
    GPU_OPER_copyIlumin_BRNSCT_GPU2CPU(iluminSouSide, iluminRecSide, 1);

    GPU_OPER_freeArrays_acoustic();
    GPU_OPER_freeArrays_acoustic_ExtBornSct();

    fclose(fWS_Bin);
    fclose(fWS_Hdr);
    fclose(fWC_Bin);
    fclose(fWC_Hdr);
    if(recordMovie && shot->ishot==ishotMovie) {
        fclose(fwav_Bin);
        fclose(fwav_Hdr);
        fclose(fsct_Bin);
        fclose(fsct_Hdr);
    }
    free(sigma);
    free(sigX);
    free(sigZ);
    free(ws);
    free(wc);
    free(wstp);
    free(wav1);
    free(wav2);
    free(sct1);
    free(sct2);
    free(lap);
    free(sigmaInv);


    if(recordBornModeling)
    {
        remove(movieIncSouFileName);
        remove(movieSctSouFileName);
        remove(movieIncRecFileName);
        remove(movieSctRecFileName);
    }
    remove(WSFileName);
    remove(WCFileName);

return;
}

void emigrateExplodingExtendedImages_HOCIG_VOCIG_Mem_GPU(int gpu_device, float *image,     float *ilumin, \
                                                         float *Hocig_background_original, float *Hocig_background_phaseShift, \
                                                         float *Hocig_residual_original,   float *Hocig_residual_phaseShift, \
                                                         float *Vocig_background_original, float *Vocig_background_phaseShift, \
                                                         float *Vocig_residual_original,   float *Vocig_residual_phaseShift, \
                                                         float *vp, long long int jt, int ntMax, \
                                                         long long int nt, long long int nx, long long int nz, \
                                                         float dx, float dz, float dt, float ox, float oz, float ot, \
                                                         long long int lx0, long long int lz0, int recordMovie) {
    long long int it;
    int it_min;
    float *sigX, *sigZ, *sigma, *sigmaInv, *ilumin_shot;
    sigX = CPU_zaloc1F(nx);
    sigZ = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);

    int mode_cigs = MODE_CIGS_ORTH;

    // Initializing GPU environment
    GPU_OPER_initEnvironment(gpu_device, nx, nz, nt, dx, dz, dt, 0, 0, lx0, lz0, 0, nxb, nzb, mode_cigs);
    GPU_OPER_alocArrays_ExtExpRef(jt);
    
    // Copying H&V OCIGS (background & residual, original & phaseShift)
    GPU_OPER_copyOCIGS_CPU2GPU(Hocig_background_original, Vocig_background_original);
    GPU_OPER_copyResidualOCIGS_CPU2GPU(Hocig_residual_original, Vocig_residual_original);
    GPU_OPER_copyOCIGS_phaseShift_CPU2GPU(Hocig_background_phaseShift, Vocig_background_phaseShift);
    GPU_OPER_copyResidualOCIGS_phaseShift_CPU2GPU(Hocig_residual_phaseShift, Vocig_residual_phaseShift);

    GPU_OPER_getVelP_CPU2GPU(vp);
    GPU_OPER_getSigmas_CPU2GPU(sigma, sigmaInv);

    FILE *f_ExpWavFor_Hdr, *f_ExpWavFor_Bin; 
    FILE *f_ExpSctFor_Hdr, *f_ExpSctFor_Bin;
    char movieExpWavForFileName[1024];
    char movieExpSctForFileName[1024];
    sprintf(movieExpWavForFileName, "movie_exp_wav_for");
    sprintf(movieExpSctForFileName, "movie_exp_sct_for");
    
    if(recordMovie)
    {
        initFilesSmart3d(movieExpWavForFileName, &f_ExpWavFor_Hdr, &f_ExpWavFor_Bin, nx, dx, ox, nz, dz, oz, ntMax/jt, dt*jt, ot);
        initFilesSmart3d(movieExpSctForFileName, &f_ExpSctFor_Hdr, &f_ExpSctFor_Bin, nx, dx, ox, nz, dz, oz, ntMax/jt, dt*jt, ot);
    }
    
    float  *wav2 = CPU_zaloc1F(nx*nz);
    float  *sct2 = CPU_zaloc1F(nx*nz);

    int stopInjectingCIGs = 0;
    int ntTotal = ntMax;
    for(it=0; it<ntTotal; it++)
    {
        GPU_OPER_setTimeStep(it);
        
        // Injecting OCIGs
        if(stopInjectingCIGs==0)
        {
            stopInjectingCIGs = GPU_WRAP_injectExtRef_Acoustic(EIMSH_CIGS_VCIG);
            if(stopInjectingCIGs)
            {
                ntTotal = it + nt;
                if(ntTotal>ntMax)    ntTotal = ntMax;
                printf("\n Stopping injection of OCIGs at it=%lld. Setting ntTotal=%d (ntMax=%d and nt=%lld). \n", \
                         it, ntTotal, ntMax, nt);
            }
        }

        // if(it%200==0)    printf("\n In ExtExpRef propagating it=%lld.", it);

        // GPU_WRAP_injectResidualExtRef_Acoustic(EIMSH_CIGS_VCIG);
        // int it_tmp = GPU_WRAP_injectBackgroundExtRef_Acoustic(EIMSH_CIGS_VCIG);
        // if(it_tmp>0)    it_min = it_tmp;

        GPU_WRAP_propagate_Acoustic_abs2();
        GPU_WRAP_propagate_Acoustic_Scattered_abs2();

        GPU_WRAP_imagingCondition(IC_REGULAR_EXPREF);

        if(it%jt==0  &&  it>=0)
        {
            GPU_OPER_copyPresentWavefield_GPU2CPU(wav2);
            GPU_OPER_copyPresentScatteredWavefield_GPU2CPU(sct2);
            if(recordMovie)
            {
                writeSamples(wav2, f_ExpWavFor_Bin, nx*nz);
                writeSamples(sct2, f_ExpSctFor_Bin, nx*nz);
            }
        }
    }

    GPU_OPER_zeroWavefieldArrays();
    GPU_OPER_zeroScatteredWavefieldArrays();

    printf("\n\n In new gradient computation: nt=%lld  jt=%lld  it_min=%d\n\n", nt, jt, it_min);

    GPU_OPER_copyImage_GPU2CPU(image);
    GPU_OPER_copyIlumin_GPU2CPU(ilumin);
    GPU_OPER_freeArrays_ExtExpRef();

    free(wav2);
    free(sct2);
    if(recordMovie)
    {
        fclose(f_ExpWavFor_Hdr);
        fclose(f_ExpWavFor_Bin);
        fclose(f_ExpSctFor_Hdr);
        fclose(f_ExpSctFor_Bin);
    }
    free(sigX);
    free(sigZ);
    free(sigma);
    free(sigmaInv);

return;
}

void emigrateExplodingExtendedImages_DOCIG_Mem_GPU(int gpu_device, float *image, float *ilumin, \
                                                   float *Docig_background_original, float *Docig_background_phaseShift, \
                                                   float *Docig_residual_original,   float *Docig_residual_phaseShift, \
                                                   float *vp, long long int jt, int ntMax, \
                                                   long long int nt, long long int nx, long long int nz, \
                                                   float dx, float dz, float dt, float ox, float oz, float ot, \
                                                   long long int ndip, float dipMin, float ddip, long long int lambda0, int recordMovie) {
    long long int it;
    int jt_movie = 4*jt;
    // int jt_movie = 1;
    int it_min;
    float *sigX, *sigZ, *sigma, *sigmaInv, *ilumin_shot;
    sigX = CPU_zaloc1F(nx);
    sigZ = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);

    int mode_cigs = MODE_CIGS_ORTH;

    // Initializing GPU environment
    int lx0 = 0;
    int lz0 = 0;
    GPU_OPER_initEnvironment(gpu_device, nx, nz, nt, dx, dz, dt, 0, 0, lx0, lz0, 0, nxb, nzb, mode_cigs);
    GPU_OPER_init_ExtExpRef_Docig(jt, lambda0, ndip, dipMin, ddip);
    printf("\n  [emigrateExplodingExtendedImages_DOCIG_Mem_GPU]:  ndip=%lld  dipMin=%f  ddip=%f  \n", ndip, dipMin, ddip);
    
    // Copying H&V OCIGS (background & residual, original & phaseShift)
    GPU_OPER_copyDOCIG_ExtExpRef_CPU2GPU(Docig_background_original, Docig_residual_original, Docig_background_phaseShift, Docig_residual_phaseShift);

    GPU_OPER_getVelP_CPU2GPU(vp);
    GPU_OPER_getSigmas_CPU2GPU(sigma, sigmaInv);

    FILE *f_ExpWavFor_Hdr, *f_ExpWavFor_Bin;
    FILE *f_ExpSctFor_Hdr, *f_ExpSctFor_Bin;
    FILE *f_ExpGrdFor_Hdr, *f_ExpGrdFor_Bin;
    char movieExpWavForFileName[1024];
    char movieExpSctForFileName[1024];
    char movieExpGrdForFileName[1024];
    sprintf(movieExpWavForFileName, "movie_exp_wav_for");
    sprintf(movieExpSctForFileName, "movie_exp_sct_for");
    sprintf(movieExpGrdForFileName, "movie_exp_grd_for");
    
    if(recordMovie)
    {
        initFilesSmart3d(movieExpWavForFileName, &f_ExpWavFor_Hdr, &f_ExpWavFor_Bin, nx, dx, ox, nz, dz, oz, ntMax/jt_movie, dt*jt_movie, ot);
        initFilesSmart3d(movieExpSctForFileName, &f_ExpSctFor_Hdr, &f_ExpSctFor_Bin, nx, dx, ox, nz, dz, oz, ntMax/jt_movie, dt*jt_movie, ot);
        initFilesSmart3d(movieExpGrdForFileName, &f_ExpGrdFor_Hdr, &f_ExpGrdFor_Bin, nz, dz, oz, nx, dx, ox, ntMax/jt_movie, dt*jt_movie, ot);
    }
    
    float  *wav2 = CPU_zaloc1F(nx*nz);
    float  *sct2 = CPU_zaloc1F(nx*nz);
    float  *grad = CPU_zaloc1F(nx*nz);

    int stopInjectingCIGs = 0;
    int ntTotal = ntMax;
    for(it=0; it<ntTotal; it++)
    {
        GPU_OPER_setTimeStep(it);
        
        // Injecting OCIGs
        if(stopInjectingCIGs==0)
        {
            stopInjectingCIGs = GPU_WRAP_injectExtRef_Acoustic_Docig();
            if(stopInjectingCIGs)
            {
                ntTotal = it + nt;
                if(ntTotal>ntMax)    ntTotal = ntMax;
                printf("\n Stopping injection of OCIGs at it=%lld. Setting ntTotal=%d (ntMax=%d and nt=%lld). \n", \
                         it, ntTotal, ntMax, nt);
            }
        }

        float lambda1=2.0*dx, lambda2=4.0*dx, lambda3=1000*dx, lambda4=2000*dx;
        if(it%(5*jt)==0) GPU_WRAP_wavefield_bandpass_kx_kz_iso_GPU2GPU(gpu_device, nx, nz, dx, dz, lambda1, lambda2, lambda3, lambda4);

        GPU_OPER_wavefieldTimeSecondDerivative_part1(0);
        GPU_WRAP_propagate_Acoustic_abs2();
        GPU_WRAP_propagate_Acoustic_Scattered_abs2();
        GPU_OPER_wavefieldTimeSecondDerivative_part2(0);

        if((it+jt)%jt==0)    GPU_WRAP_imagingCondition_ExtExpRef(IC_REGULAR_EXPREF_SECDER);
        // if(it%jt==0)    GPU_WRAP_imagingCondition_ExtExpRef(IC_REGULAR_EXPREF_SECDER);

        if(it%jt_movie==0  &&  it>=0  &&  recordMovie)
        {
            GPU_OPER_copyPresentWavefield_GPU2CPU(wav2);
            GPU_OPER_copyPresentScatteredWavefield_GPU2CPU(sct2);
            GPU_OPER_copyImage_GPU2CPU(grad);
            writeSamples(wav2, f_ExpWavFor_Bin, nx*nz);
            writeSamples(sct2, f_ExpSctFor_Bin, nx*nz);
            writeSamples(grad, f_ExpGrdFor_Bin, nx*nz);
        }
    }

    printf("\n\n In new gradient computation: nt=%lld  jt=%lld  it_min=%d\n\n", nt, jt, it_min);

    GPU_OPER_copyImage_GPU2CPU(image);
    GPU_OPER_copyIlumin_GPU2CPU(ilumin);
    GPU_OPER_freeArrays_ExtExpRef_Docig();

    free(wav2);
    free(sct2);
    free(grad);
    if(recordMovie)
    {
        fclose(f_ExpWavFor_Hdr);
        fclose(f_ExpWavFor_Bin);
        fclose(f_ExpSctFor_Hdr);
        fclose(f_ExpSctFor_Bin);
        fclose(f_ExpGrdFor_Hdr);
        fclose(f_ExpGrdFor_Bin);
    }

    free(sigX);
    free(sigZ);
    free(sigma);
    free(sigmaInv);

return;
}

void migrateExtBorn_lx_lz_lt_Mem_GPU_old(int idev, float *imageSouSide, float *imageRecSide, \
                                     float *imageSouSide_shot, float *imageRecSide_shot, \
                                     float *iluminSouSide, float *iluminRecSide, \
                                     float *eref, float *Heref, float *Veref, float *vp, Shot *shot, \
                                     long long int jt, long long int jt_EIC, long long int nt, \
                                     long long int nx, long long int nz, float dx, float dz, \
                                     float dt, float ox, float oz, float ot, long long int lx0, long long int lz0, long long int lt0, \
                                     int ishotMovie, int recordMovie, int imageDomain, int OPTIM) {

    int mode_cigs;
    if(Heref!=NULL  &&  Veref!=NULL)
        mode_cigs = MODE_CIGS_ORTH;
    else 
        mode_cigs = MODE_CIGS_FULL;
    
    long long int     it;
    size_t  err;
    float  *wav1, *wav2;
    float  *sct1, *sct2;
    float  *lap, *sigX, *sigZ, *sigma, *sigmaInv, *tmp, *ws, *wc, *secDer, *wstp;
    secDer = CPU_zaloc1F(nx*nz);
    wav1   = CPU_zaloc1F(nx*nz);
    wav2   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
    ws     = CPU_zaloc1F(nx*nz);
    wstp   = CPU_zaloc1F((nx-2*nxb)*(nz-2*nzb));
    wc     = CPU_zaloc1F(nx*nz);
    sct1   = CPU_zaloc1F(nx*nz);
    sct2   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
    lap    = CPU_zaloc1F(nx*nz);
    sigX   = CPU_zaloc1F(nx);
    sigZ   = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);

    int mode, side;
    int store_inc = 1, store_sct = 1;
    if(imageDomain==0)    mode=IC_BORNSCT_IMAGE2;
    else                  mode=IC_BORNSCT_IMAGE;

    long long int ixMin = shot->ixMin;
    long long int ixMax = shot->ixMax;
    long long int nxProp = ixMax - ixMin + 1;

    // Initializing GPU environment for current shot
    if(!OPTIM)
    {
        GPU_OPER_initEnvironment(idev, nx, nz, nt, dx, dz, dt, shot->nsou, shot->nrec, lx0, lz0, lt0, nxb, nzb, mode_cigs);
        GPU_OPER_alocArrays_acoustic();
        GPU_OPER_alocArrays_acoustic_ExtBornSct();
        if(mode_cigs == MODE_CIGS_FULL) {
            GPU_OPER_copyEImage_CPU2GPU(eref);
        }
        else if(mode_cigs == MODE_CIGS_ORTH) {
            GPU_OPER_copyOCIGS_CPU2GPU(Heref, Veref);
        }
    }
    else
    {
        GPU_OPER_zeroArrays_acoustic();
        GPU_OPER_zeroArrays_ExtBornSct();
    }
    GPU_OPER_getVelP_CPU2GPU(vp);
    GPU_OPER_getSigmas_CPU2GPU(sigma, sigmaInv);
    GPU_OPER_getSourceWavelet_CPU2GPU(shot->source);
    GPU_OPER_getRecordedData_CPU2GPU(shot->seismogram);
    GPU_OPER_getSourceReceiverIndexes_CPU2GPU(shot->ixs, shot->izs, shot->ixr, shot->izr);

    if(jt_EIC)   GPU_OPER_alocArrays_EIC(jt_EIC, store_inc, store_sct);
    if(OPTIM)    GPU_OPER_alocArrays_EIC(jt, store_inc, store_sct);

    // Source wavefield
    int recordBornModeling = ((shot->ishot==ishotMovie) && recordMovie);

    float  *inc_sou, *inc_rec, *sct_sou, *sct_rec, *img_sou, *img_rec;
    if(recordBornModeling)
    {
        inc_sou   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
        inc_rec   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
        sct_sou   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
        sct_rec   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
        img_sou   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
        img_rec   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
    }

    FILE *fimg_sou_Hdr, *fimg_sou_Bin;
    FILE *fimg_rec_Hdr, *fimg_rec_Bin;
    FILE *fWavSouInc_Bin, *fWavSouInc_Hdr;
    FILE *fWavSouSct_Bin, *fWavSouSct_Hdr;
    FILE *fWavRecInc_Bin, *fWavRecInc_Hdr;
    FILE *fWavRecSct_Bin, *fWavRecSct_Hdr;
    char movieIncSouFileName[1024];
    char movieSctSouFileName[1024];
    char movieIncRecFileName[1024];
    char movieSctRecFileName[1024];
    char movieImgSouFileName[1024];
    char movieImgRecFileName[1024];
    if(recordBornModeling)
    {
        sprintf(movieIncSouFileName, "movie_inc_extBorn_sou_ishot%d", ishotMovie);
        sprintf(movieSctSouFileName, "movie_sct_extBorn_sou_ishot%d", ishotMovie);
        sprintf(movieImgSouFileName, "movie_img_extBorn_sou_ishot%d", ishotMovie);
        sprintf(movieIncRecFileName, "movie_inc_extBorn_rec_ishot%d", ishotMovie);
        sprintf(movieSctSouFileName, "movie_sct_extBorn_rec_ishot%d", ishotMovie);
        sprintf(movieImgRecFileName, "movie_img_extBorn_rec_ishot%d", ishotMovie);
        initFilesSmart3d(movieIncSouFileName, &fWavSouInc_Hdr, &fWavSouInc_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieSctSouFileName, &fWavSouSct_Hdr, &fWavSouSct_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieIncRecFileName, &fWavRecInc_Hdr, &fWavRecInc_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieSctRecFileName, &fWavRecSct_Hdr, &fWavRecSct_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieImgSouFileName, &fimg_sou_Hdr  , &fimg_sou_Bin  , nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieImgRecFileName, &fimg_rec_Hdr  , &fimg_rec_Bin  , nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
    }

    // Forward
    for(it=0; it<nt; it++)
    {
        GPU_OPER_setTimeStep(it);

        if(it%jt==0)
        {
            long long int shift = nx*nz * (it/jt);
            if(!OPTIM) // IF COMMENTED, IS DEBUG
            {
                // GPU_OPER_copyPresentWavefield_GPU2CPU(wav2+shift);
                // GPU_OPER_copyPresentScatteredWavefield_GPU2CPU(sct2+shift);
                GPU_OPER_copyPresentWavefield_secDer_GPU2CPU(wav2+shift);
                if(imageDomain==0)
                    GPU_OPER_copyPresentScatteredWavefield_secDer_GPU2CPU(sct2+shift);
                else
                    GPU_OPER_copyPresentScatteredWavefield_GPU2CPU(sct2+shift);
            }
            GPU_WRAP_imagingCondition(IC_BORNSCT_ILUMIN);

            if(recordBornModeling)
            {
                GPU_OPER_copyPresentWavefield_secDer_GPU2CPU(inc_sou+shift);
                GPU_OPER_copyPresentScatteredWavefield_secDer_GPU2CPU(sct_sou+shift);
                // GPU_OPER_copyPresentWavefield_secDer_GPU2CPU(wav2+shift);
                // GPU_OPER_copyPresentScatteredWavefield_GPU2CPU(sct2+shift);
            }
        }
        
        // printf("\n idev=%d  it=%d  forward ", idev, it);
        side = IC_SOU_WAVE;
        GPU_OPER_EIC_storeWavefields_GPU2GPU(it, side, mode);
        
        GPU_OPER_wavefieldTimeSecondDerivative_part1(0);
        GPU_OPER_wavefieldTimeSecondDerivative_part1(1);
        GPU_WRAP_injectSource_Acoustic();
        GPU_WRAP_propagate_Acoustic_abs2();
        GPU_OPER_wavefieldTimeSecondDerivative_part2(0);
        if(mode_cigs == MODE_CIGS_FULL)
        {
            if(lz0==0)    GPU_WRAP_scatter_Acoustic_lx_opt();
            else          GPU_WRAP_scatter_Acoustic_lx_lz_opt();
        }
        else if(mode_cigs == MODE_CIGS_ORTH)
        {
            if(OPTIM)    GPU_WRAP_scatter_Acoustic_OCIGS_opt(EIMSH_CIGS_VCIG);
            else         GPU_WRAP_scatter_Acoustic_OCIGS_opt(EIMSH_CIGS_VCIG);
        }
        GPU_WRAP_propagate_Acoustic_Scattered_abs2();
        GPU_OPER_wavefieldTimeSecondDerivative_part2(1);
    }

    GPU_OPER_zeroWavefieldArrays();
    GPU_OPER_zeroScatteredWavefieldArrays();
    // Backward
    for(it=nt-1; it>=0; it--)
    {
        GPU_OPER_setTimeStep(it);

        long long int shift = nx*nz * (it/jt);

        if(it%jt==0  &&  jt_EIC<=0)
        {
            if(!OPTIM)
            {
                GPU_OPER_copyPresentWavefieldWS_CPU2GPU(wav2+shift);
                GPU_OPER_copyPresentWavefieldWC_CPU2GPU(sct2+shift);
                if(imageDomain==0)
                    GPU_WRAP_imagingCondition(IC_BORNSCT_IMAGE);
                else
                    GPU_WRAP_imagingCondition(IC_BORNSCT_IMAGE2);
            }
            // Debug - begin
            /*
            GPU_OPER_copyPresentWavefieldWS_CPU2GPU(wav2+shift);
            GPU_OPER_copyPresentWavefieldWC_CPU2GPU(sct2+shift);
            if(imageDomain==0)
                GPU_WRAP_imagingCondition(IC_BORNSCT_IMAGE);
            else
                GPU_WRAP_imagingCondition(IC_BORNSCT_IMAGE2);
            //*/
            // Debug - end
        }
        //*
        if(OPTIM)
        {
            if(imageDomain==0)
                GPU_WRAP_imagingCondition_smart(IC_BORNSCT_IMAGE, it);
            else
                GPU_WRAP_imagingCondition_smart(IC_BORNSCT_IMAGE2, it);
        }
        //*/

        if(jt_EIC>0) 
        {
            side = IC_REC_WAVE;
            GPU_OPER_EIC_storeWavefields_GPU2GPU(it, side, mode);
        }

        // DEBUG - begin ////////////////////
        if(it%jt==0)
        {
            if(recordBornModeling)
            {
                GPU_OPER_copyPresentWavefield_GPU2CPU(inc_rec+shift);
                GPU_OPER_copyPresentScatteredWavefield_GPU2CPU(sct_rec+shift);
                GPU_OPER_copyImage_BRNSCT_GPU2CPU(img_sou+shift, img_rec+shift, 0);
                // GPU_OPER_copyPresentWavefield_secDer_GPU2CPU(wav1);
                // GPU_OPER_copyPresentScatteredWavefield_GPU2CPU(sct1);
                // writeSamples(wav1, fWavRecInc_Bin, nx*nz);
                // writeSamples(sct1, fWavRecSct_Bin, nx*nz);
            }
        }
        // DEBUG - end ////////////////////

        GPU_OPER_wavefieldTimeSecondDerivative_part1(0);
        GPU_WRAP_recordInjectData(1);
        GPU_WRAP_propagate_Acoustic_abs2();
        GPU_OPER_wavefieldTimeSecondDerivative_part2(0);
        if(mode_cigs == MODE_CIGS_FULL)
        {
            if(lz0==0)    GPU_WRAP_scatter_Acoustic_lx_opt();
            else          GPU_WRAP_scatter_Acoustic_lx_lz_opt();
        }
        else if(mode_cigs == MODE_CIGS_ORTH)
        {
            if(OPTIM)    GPU_WRAP_scatter_Acoustic_OCIGS_opt(EIMSH_CIGS_VCIG);
            else         GPU_WRAP_scatter_Acoustic_OCIGS_opt(EIMSH_CIGS_VCIG);
        }
        GPU_WRAP_propagate_Acoustic_Scattered_abs2();
    }
    if(jt_EIC>0)    GPU_WRAP_imagingCondition_transmission();

    if(!OPTIM)
    {
        GPU_OPER_copyImage_BRNSCT_GPU2CPU(  imageSouSide_shot,   imageRecSide_shot, 0);
        GPU_OPER_copyImage_BRNSCT_GPU2CPU(  imageSouSide     ,   imageRecSide     , 1);
        GPU_OPER_copyIlumin_BRNSCT_GPU2CPU(iluminSouSide     ,  iluminRecSide     , 1);
        GPU_OPER_freeArrays_acoustic_ExtBornSct();
        GPU_OPER_freeArrays_acoustic();
    }
    if(OPTIM  ||  jt_EIC>0)  GPU_OPER_freeArrays_EIC();


    if(recordBornModeling)
    {
        for(it=0; it<nt/jt; it++)
        {
            long long int shift = it*nx*nz;
            writeSamples(inc_sou+shift, fWavSouInc_Bin, nx*nz);
            writeSamples(sct_sou+shift, fWavSouSct_Bin, nx*nz);
            writeSamples(inc_rec+shift, fWavRecInc_Bin, nx*nz);
            writeSamples(sct_rec+shift, fWavRecSct_Bin, nx*nz);
            writeSamples(img_sou+shift,   fimg_sou_Bin, nx*nz);
            writeSamples(img_rec+shift,   fimg_rec_Bin, nx*nz);
        }
        fclose(fWavSouInc_Bin);
        fclose(fWavSouInc_Hdr);

        fclose(fWavSouSct_Bin);
        fclose(fWavSouSct_Hdr);

        fclose(fWavRecInc_Bin);
        fclose(fWavRecInc_Hdr);

        fclose(fWavRecSct_Hdr);
        fclose(fWavRecSct_Bin);

        fclose(fimg_sou_Bin);
        fclose(fimg_sou_Hdr);

        fclose(fimg_rec_Bin);
        fclose(fimg_rec_Hdr);

        free(inc_sou);
        free(inc_rec);
        free(sct_sou);
        free(sct_rec);
        free(img_sou);
        free(img_rec);
    }

    
    free(sigma);
    free(sigX);
    free(sigZ);
    free(ws);
    free(wc);
    free(wstp);
    free(wav1);
    free(wav2);
    free(sct1);
    free(sct2);
    free(lap);
    free(sigmaInv);

    // if(recordBornModeling)
    // {
    //     outputSmart2d("grad_souSide_debug", imageSouSide, nz, dz, oz, nx, dx, ox);
    //     outputSmart2d("grad_recSide_debug", imageRecSide, nz, dz, oz, nx, dx, ox);
    // }

    /*
    if(recordBornModeling)
    {
        remove(movieIncSouFileName);
        remove(movieSctSouFileName);
        remove(movieIncRecFileName);
        remove(movieSctRecFileName);
    }
    */
return;
}



void ExtendedBornModeling_TAU_Mem_GPU(int idev, float *Tcig, float *vp, Shot3D *shot, \
                                      MSH *msh, Eimage_MSH *eimsh, int ishotMovie, int recordMovie) {

    long long int nxb    = msh->nxb;
    long long int nyb    = msh->nyb;
    long long int nzb    = msh->nzb;
    long long int nx     = msh->nx;
    long long int ny     = msh->ny;
    long long int nz     = msh->nz;
    long long int nt     = msh->nt;
    float         dx     = msh->dx;
    float         dy     = msh->dy;
    float         dz     = msh->dz;
    float         dt     = msh->dt;
    float         ox     = msh->ox;
    float         oy     = msh->oy;
    float         oz     = msh->oz;
    float         ot     = msh->ot;
    long long int jt     = msh->jt;
    long long int jt_EIC = msh->jt_EIC;
    long long int lt0    = 0;
    long long int it;
    
    int mode_cigs = MODE_CIGS_TIME;

    size_t  err;
    float  *wav1, *wav2;
    float  *sct1, *sct2;
    float  *lap, *sigX, *sigZ, *sigma, *sigmaInv, *tmp, *ws, *wc, *secDer, *wstp;
    secDer = CPU_zaloc1F(nx*nz);
    wav1   = CPU_zaloc1F(nx*nz);
    wav2   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
    ws     = CPU_zaloc1F(nx*nz);
    wstp   = CPU_zaloc1F((nx-2*nxb)*(nz-2*nzb));
    wc     = CPU_zaloc1F(nx*nz);
    sct1   = CPU_zaloc1F(nx*nz);
    sct2   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
    lap    = CPU_zaloc1F(nx*nz);
    sigX   = CPU_zaloc1F(nx);
    sigZ   = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);

    int store_inc = 1, store_sct = 0;

    long long int ixMin = shot->ixMin;
    long long int ixMax = shot->ixMax;
    long long int nxProp = ixMax - ixMin + 1;

    // Initializing GPU environment for current shot
    GPU_OPER_initEnvironment(idev, nx, nz, nt, dx, dz, dt, shot->nsou, shot->nrec, eimsh->lx0, eimsh->lz0, 0, nxb, nzb, MODE_CIGS_TIME);
    GPU_OPER_initTAU(eimsh->tau0, eimsh->dtau);
    GPU_OPER_alocArrays_acoustic();
    GPU_OPER_alocArrays_acoustic_ExtBornSct();
    GPU_OPER_copyTCIG_CPU2GPU(Tcig);

    GPU_OPER_zeroArrays_acoustic();
    GPU_OPER_zeroArrays_ExtBornSct();
    
    GPU_OPER_getVelP_CPU2GPU(vp);
    GPU_OPER_getSigmas_CPU2GPU(sigma, sigmaInv);
    GPU_OPER_getSourceWavelet_CPU2GPU(shot->source);
    GPU_OPER_getRecordedData_CPU2GPU(shot->seismogram);
    GPU_OPER_getSourceReceiverIndexes_CPU2GPU(shot->ixs, shot->izs, shot->ixr, shot->izr);

    GPU_OPER_alocArrays_EIC(jt, store_inc, store_sct);

    // Source wavefield
    int recordBornModeling = ((shot->ishot==ishotMovie) && recordMovie);

    float  *inc_sou, *sct_sou;
    if(recordBornModeling)
    {
        inc_sou   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
        sct_sou   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
    }

    FILE *fimg_sou_Hdr, *fimg_sou_Bin;
    FILE *fWavSouInc_Bin, *fWavSouInc_Hdr;
    FILE *fWavSouSct_Bin, *fWavSouSct_Hdr;
    char movieIncSouFileName[1024];
    char movieSctSouFileName[1024];
    if(recordBornModeling)
    {
        sprintf(movieIncSouFileName, "movie_inc_extBorn_sou_ishot%d", ishotMovie);
        sprintf(movieSctSouFileName, "movie_sct_extBorn_sou_ishot%d", ishotMovie);
        initFilesSmart3d(movieIncSouFileName, &fWavSouInc_Hdr, &fWavSouInc_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieSctSouFileName, &fWavSouSct_Hdr, &fWavSouSct_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
    }

    
    // Forward in time incident wavefield
    store_inc = 1;
    store_sct = 0;
    GPU_OPER_setStore_inc_sct(store_inc, store_sct);
    for(it=0; it<nt; it++)
    {
        GPU_OPER_setTimeStep(it);

        if(it%jt==0) {
            if(recordBornModeling) {
                long long int shift = nx*nz * (it/jt);
                GPU_OPER_copyPresentWavefield_secDer_GPU2CPU(inc_sou+shift);
            }
        }
        
        // Time second-derivative: part1
        GPU_OPER_wavefieldTimeSecondDerivative_part1(0);
        
        // Inject source and propagate incident wavefield
        GPU_WRAP_injectSource_Acoustic();
        GPU_WRAP_propagate_Acoustic_abs2();

        // Time second-derivative: part2
        GPU_OPER_wavefieldTimeSecondDerivative_part2(0);

        GPU_OPER_EIC_storeWavefields_GPU2GPU(it, IC_SOU_WAVE, IC_BORNSCT_IMAGE2);
    }

    // Forward in time scattered wavefield
    // store_inc = 0;
    // store_sct = 1;
    // GPU_OPER_setStore_inc_sct(store_inc, store_sct);
    for(it=0; it<nt; it++)
    {
        GPU_OPER_setTimeStep(it);
        
        GPU_WRAP_recordInjectData(2);

        if(it%jt==0) {
            if(recordBornModeling) {
                long long int shift = nx*nz * (it/jt);
                GPU_OPER_copyPresentScatteredWavefield_secDer_GPU2CPU(sct_sou+shift);
            }
        }        
        // Time second-derivative: part1
        GPU_OPER_wavefieldTimeSecondDerivative_part1(1);

        // Inject scattered-virtual source and propagate scattered
        GPU_WRAP_scatteringCondition_TAU(IC_SOU_WAVE);
        GPU_WRAP_propagate_Acoustic_Scattered_abs2();

        // Time second-derivative: part2
        GPU_OPER_wavefieldTimeSecondDerivative_part2(1);
    }

    GPU_OPER_getRecordedData_GPU2CPU(shot->seismogram);
    
    GPU_OPER_freeArrays_EIC();
    GPU_OPER_freeArrays_acoustic_ExtBornSct();
    GPU_OPER_freeArrays_acoustic();

    if(recordBornModeling)
    {
        printf("\n Saving files \n");
        for(it=0; it<nt/jt; it++)
        {
            long long int shift = it*nx*nz;
            writeSamples(inc_sou+shift, fWavSouInc_Bin, nx*nz);
            writeSamples(sct_sou+shift, fWavSouSct_Bin, nx*nz);
        }

        fclose(fWavSouInc_Bin);
        fclose(fWavSouInc_Hdr);

        fclose(fWavSouSct_Bin);
        fclose(fWavSouSct_Hdr);

        free(inc_sou);
        free(sct_sou);
    }

    
    free(sigma);
    free(sigX);
    free(sigZ);
    free(ws);
    free(wc);
    free(wstp);
    free(wav1);
    free(wav2);
    free(sct1);
    free(sct2);
    free(lap);
    free(sigmaInv);

return;
}


void migrateExtBorn_TAU_Mem_GPU(int idev, float *imageSouSide, float *imageRecSide, \
                                float *iluminSouSide, float *iluminRecSide, \
                                float *Tcig, float *vp, Shot3D *shot, \
                                MSH *msh, Eimage_MSH *eimsh, int ishotMovie, int recordMovie) {

    long long int nxb    = msh->nxb;
    long long int nyb    = msh->nyb;
    long long int nzb    = msh->nzb;
    long long int nx     = msh->nx;
    long long int ny     = msh->ny;
    long long int nz     = msh->nz;
    long long int nt     = msh->nt;
    float         dx     = msh->dx;
    float         dy     = msh->dy;
    float         dz     = msh->dz;
    float         dt     = msh->dt;
    float         ox     = msh->ox;
    float         oy     = msh->oy;
    float         oz     = msh->oz;
    float         ot     = msh->ot;
    long long int jt     = msh->jt;
    long long int jt_EIC = msh->jt_EIC;
    long long int lt0    = 0;
    long long int it;
    
    int mode_cigs = MODE_CIGS_TIME;

    size_t  err;
    float  *sigX, *sigZ, *sigma, *sigmaInv;
    sigX   = CPU_zaloc1F(nx);
    sigZ   = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);

    int store_inc = 1, store_sct = 1;

    long long int ixMin = shot->ixMin;
    long long int ixMax = shot->ixMax;
    long long int nxProp = ixMax - ixMin + 1;

    // Initializing GPU environment for current shot
    GPU_OPER_initEnvironment(procGPU, nx, nz, nt, dx, dz, dt, shot->nsou, shot->nrec, eimsh->lx0, eimsh->lz0, 0, nxb, nzb, mode_cigs);
    GPU_OPER_initTAU(eimsh->tau0, eimsh->dtau);
    GPU_OPER_alocArrays_acoustic();
    GPU_OPER_alocArrays_acoustic_ExtBornSct();
    GPU_OPER_copyTCIG_CPU2GPU(Tcig);

    GPU_OPER_zeroArrays_acoustic();
    GPU_OPER_zeroArrays_ExtBornSct();
    
    GPU_OPER_getVelP_CPU2GPU(vp);
    GPU_OPER_getSigmas_CPU2GPU(sigma, sigmaInv);
    GPU_OPER_getSourceWavelet_CPU2GPU(shot->source);
    GPU_OPER_getRecordedData_CPU2GPU(shot->seismogram);
    GPU_OPER_getSourceReceiverIndexes_CPU2GPU(shot->ixs, shot->izs, shot->ixr, shot->izr);

    GPU_OPER_alocArrays_EIC(jt, store_inc, store_sct);

    float *ilumin_sou_inc = CPU_zaloc1F(nx*nz);
    float *ilumin_sou_sct = CPU_zaloc1F(nx*nz);
    float *ilumin_rec_inc = CPU_zaloc1F(nx*nz);
    float *ilumin_rec_sct = CPU_zaloc1F(nx*nz);

    // Source wavefield
    int recordBornModeling = ((shot->ishot==ishotMovie) && recordMovie);

    float  *inc_sou, *inc_rec, *sct_sou, *sct_rec, *img_sou, *img_rec;
    if(recordBornModeling)
    {
        inc_sou   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
        inc_rec   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
        sct_sou   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
        sct_rec   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
        img_sou   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
        img_rec   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
    }

    FILE *fimg_sou_Hdr, *fimg_sou_Bin;
    FILE *fimg_rec_Hdr, *fimg_rec_Bin;
    FILE *fWavSouInc_Bin, *fWavSouInc_Hdr;
    FILE *fWavSouSct_Bin, *fWavSouSct_Hdr;
    FILE *fWavRecInc_Bin, *fWavRecInc_Hdr;
    FILE *fWavRecSct_Bin, *fWavRecSct_Hdr;
    char movieIncSouFileName[1024];
    char movieSctSouFileName[1024];
    char movieIncRecFileName[1024];
    char movieSctRecFileName[1024];
    char movieImgSouFileName[1024];
    char movieImgRecFileName[1024];
    if(recordBornModeling)
    {
        sprintf(movieIncSouFileName, "movie_inc_extBorn_sou_ishot%d", ishotMovie);
        sprintf(movieSctSouFileName, "movie_sct_extBorn_sou_ishot%d", ishotMovie);
        sprintf(movieImgSouFileName, "movie_img_extBorn_sou_ishot%d", ishotMovie);
        sprintf(movieIncRecFileName, "movie_inc_extBorn_rec_ishot%d", ishotMovie);
        sprintf(movieSctRecFileName, "movie_sct_extBorn_rec_ishot%d", ishotMovie);
        sprintf(movieImgRecFileName, "movie_img_extBorn_rec_ishot%d", ishotMovie);
        initFilesSmart3d(movieIncSouFileName, &fWavSouInc_Hdr, &fWavSouInc_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieSctSouFileName, &fWavSouSct_Hdr, &fWavSouSct_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieIncRecFileName, &fWavRecInc_Hdr, &fWavRecInc_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieSctRecFileName, &fWavRecSct_Hdr, &fWavRecSct_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieImgSouFileName, &fimg_sou_Hdr  , &fimg_sou_Bin  , nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieImgRecFileName, &fimg_rec_Hdr  , &fimg_rec_Bin  , nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
    }
    
    GPU_OPER_zeroIluminArray();
    // Forward in time incident source-wavefield
    store_inc = 1;
    store_sct = 0;
    GPU_OPER_setStore_inc_sct(store_inc, store_sct);
    for(it=0; it<nt; it++)
    {
        GPU_OPER_setTimeStep(it);

        if(it%jt==0)
        {
            GPU_WRAP_imagingCondition(IC_BORNSCT_ILUMIN_INC);
            if(recordBornModeling) {
                long long int shift = nx*nz * (it/jt);
                // GPU_OPER_copyPresentWavefield_GPU2CPU(inc_sou+shift);
                GPU_OPER_copyPresentWavefield_secDer_GPU2CPU(inc_sou+shift);
            }
        }
        
        // Time second-derivative: part1
        GPU_OPER_wavefieldTimeSecondDerivative_part1(0);
        
        // Inject source and propagate incident
        GPU_WRAP_injectSource_Acoustic();
        GPU_WRAP_propagate_Acoustic_abs2();

        // Time second-derivative: part2
        GPU_OPER_wavefieldTimeSecondDerivative_part2(0);

        GPU_OPER_EIC_storeWavefields_GPU2GPU(it, IC_SOU_WAVE, IC_BORNSCT_IMAGE2);
    }
    GPU_OPER_copyIlumin_GPU2CPU(ilumin_sou_inc);

    GPU_OPER_zeroIluminArray();
    // Forward in time scattered source-wavefield
    store_inc = 0;
    store_sct = 1;
    GPU_OPER_setStore_inc_sct(store_inc, store_sct);
    for(it=0; it<nt; it++)
    {
        GPU_OPER_setTimeStep(it);

        if(it%jt==0) 
        {
            GPU_WRAP_imagingCondition(IC_BORNSCT_ILUMIN_SCT);
            if(recordBornModeling) {
                long long int shift = nx*nz * (it/jt);
                // GPU_OPER_copyPresentScatteredWavefield_secDer_GPU2CPU(sct_sou+shift);
                GPU_OPER_copyPresentScatteredWavefield_GPU2CPU(sct_sou+shift);
            }
        }
        GPU_OPER_EIC_storeWavefields_GPU2GPU(it, IC_SOU_WAVE, IC_BORNSCT_IMAGE2);
        
        // Time second-derivative: part1
        GPU_OPER_wavefieldTimeSecondDerivative_part1(1);

        // Inject scattered-virtual source and propagate scattered
        GPU_WRAP_scatteringCondition_TAU(IC_SOU_WAVE);
        GPU_WRAP_propagate_Acoustic_Scattered_abs2();

        // Time second-derivative: part2
        GPU_OPER_wavefieldTimeSecondDerivative_part2(1);
    }
    GPU_OPER_copyIlumin_GPU2CPU(ilumin_sou_sct);

    GPU_OPER_zeroWavefieldArrays();
    GPU_OPER_zeroScatteredWavefieldArrays();

    GPU_OPER_zeroIluminArray();
    // Backward in time incident receiver-wavefield
    store_inc = 1;
    store_sct = 0;
    GPU_OPER_setStore_inc_sct(store_inc, store_sct);
    for(it=nt-1; it>=0; it--)
    {
        GPU_OPER_setTimeStep(it);

        // Imaging condition to compute gradient on receiver side (source-scattered x receiver-incident)
        GPU_WRAP_imagingCondition_smart(IC_REGULAR_IMAGE3, it);

        if(it%jt==0)
        {
            GPU_WRAP_imagingCondition(IC_BORNSCT_ILUMIN_INC);
            if(recordBornModeling)
            {
                long long int shift = nx*nz * (it/jt);
                GPU_OPER_copyPresentWavefield_GPU2CPU(inc_rec+shift);
                GPU_OPER_copyPresentScatteredWavefield_GPU2CPU(sct_rec+shift);
                GPU_OPER_copyImage_BRNSCT_GPU2CPU(img_sou+shift, img_rec+shift, 0);
            }
        }

        // Time second-derivative: part1
        GPU_OPER_wavefieldTimeSecondDerivative_part1(0);

        // Inject receivers and propagate incident wavefield
        GPU_WRAP_recordInjectData(1);
        GPU_WRAP_propagate_Acoustic_abs2();

        // Time second-derivative: part2
        GPU_OPER_wavefieldTimeSecondDerivative_part2(0);

        GPU_OPER_EIC_storeWavefields_GPU2GPU(it, IC_REC_WAVE, IC_BORNSCT_IMAGE3);
    }
    GPU_OPER_copyIlumin_GPU2CPU(ilumin_rec_inc);
    
    GPU_OPER_zeroIluminArray();
    // Backward in time scattered receiver-wavefield
    for(it=nt-1; it>=0; it--)
    {
        GPU_OPER_setTimeStep(it);

        // Imaging condition to compute gradient on source side (source-incident x receiver-scattered)
        GPU_WRAP_imagingCondition_smart(IC_REGULAR_IMAGE2, it);

        if(it%jt==0)
        {
            GPU_WRAP_imagingCondition(IC_BORNSCT_ILUMIN_SCT);
            if(recordBornModeling)
            {
                long long int shift = nx*nz * (it/jt);
                GPU_OPER_copyPresentWavefield_GPU2CPU(inc_rec+shift);
                GPU_OPER_copyPresentScatteredWavefield_GPU2CPU(sct_rec+shift);
                GPU_OPER_copyImage_BRNSCT_GPU2CPU(img_sou+shift, img_rec+shift, 0);
            }
        }

        // Inject receivers and propagate incident wavefield
        GPU_WRAP_scatteringCondition_TAU(IC_REC_WAVE);
        GPU_WRAP_propagate_Acoustic_Scattered_abs2();
    }
    GPU_OPER_copyIlumin_GPU2CPU(ilumin_rec_sct);
    GPU_OPER_zeroIluminArray();

    multiplyArrays(iluminSouSide, ilumin_sou_inc, ilumin_rec_sct, nx*nz);
    multiplyArrays(iluminRecSide, ilumin_rec_inc, ilumin_sou_sct, nx*nz);
    free(ilumin_sou_inc);
    free(ilumin_sou_sct);
    free(ilumin_rec_inc);
    free(ilumin_rec_sct);

    GPU_OPER_freeArrays_EIC();
    GPU_OPER_copyImage_BRNSCT_GPU2CPU(imageSouSide, imageRecSide, 1);
    // GPU_OPER_copyIlumin_BRNSCT_GPU2CPU(iluminSouSide, iluminRecSide, 1);
    GPU_OPER_freeArrays_acoustic_ExtBornSct();
    GPU_OPER_freeArrays_acoustic();

    if(recordBornModeling)
    {
        printf("\n Saving files \n");
        for(it=0; it<nt/jt; it++)
        {
            long long int shift = it*nx*nz;
            writeSamples(inc_sou+shift, fWavSouInc_Bin, nx*nz);
            writeSamples(sct_sou+shift, fWavSouSct_Bin, nx*nz);
            writeSamples(inc_rec+shift, fWavRecInc_Bin, nx*nz);
            writeSamples(sct_rec+shift, fWavRecSct_Bin, nx*nz);
            // writeSamples(img_sou+shift,   fimg_sou_Bin, nx*nz);
            // writeSamples(img_rec+shift,   fimg_rec_Bin, nx*nz);
        }

        fclose(fWavSouInc_Bin);
        fclose(fWavSouInc_Hdr);

        fclose(fWavSouSct_Bin);
        fclose(fWavSouSct_Hdr);

        fclose(fWavRecInc_Bin);
        fclose(fWavRecInc_Hdr);

        fclose(fWavRecSct_Hdr);
        fclose(fWavRecSct_Bin);

        fclose(fimg_sou_Bin);
        fclose(fimg_sou_Hdr);

        fclose(fimg_rec_Bin);
        fclose(fimg_rec_Hdr);

        free(inc_sou);
        free(inc_rec);
        free(sct_sou);
        free(sct_rec);
        free(img_sou);
        free(img_rec);
    }

    free(sigma);
    free(sigX);
    free(sigZ);
    free(sigmaInv);

return;
}

void migrateExtBorn_lx_lz_lt_Mem_GPU(int idev, float *imageSouSide, float *imageRecSide, \
                                     float *iluminSouSide, float *iluminRecSide, \
                                     float *eref, float *Heref, float *Veref, float *vp, Shot3D *shot, \
                                     MSH *msh, Eimage_MSH *eimsh, int ishotMovie, int recordMovie) {

    long long int nxb    = msh->nxb;
    long long int nyb    = msh->nyb;
    long long int nzb    = msh->nzb;
    long long int nx     = msh->nx;
    long long int ny     = msh->ny;
    long long int nz     = msh->nz;
    long long int nt     = msh->nt;
    float         dx     = msh->dx;
    float         dy     = msh->dy;
    float         dz     = msh->dz;
    float         dt     = msh->dt;
    float         ox     = msh->ox;
    float         oy     = msh->oy;
    float         oz     = msh->oz;
    float         ot     = msh->ot;
    long long int jt     = msh->jt;
    long long int jt_EIC = msh->jt_EIC;
    long long int lt0    = 0;
    long long int it;
    
    int mode_cigs;
    if(Heref!=NULL  &&  Veref!=NULL)    mode_cigs = MODE_CIGS_ORTH;
    else                                mode_cigs = MODE_CIGS_FULL;
    

    size_t  err;
    float  *wav1, *wav2;
    float  *sct1, *sct2;
    float  *lap, *sigX, *sigZ, *sigma, *sigmaInv, *tmp, *ws, *wc, *secDer, *wstp;
    secDer = CPU_zaloc1F(nx*nz);
    wav1   = CPU_zaloc1F(nx*nz);
    wav2   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
    ws     = CPU_zaloc1F(nx*nz);
    wstp   = CPU_zaloc1F((nx-2*nxb)*(nz-2*nzb));
    wc     = CPU_zaloc1F(nx*nz);
    sct1   = CPU_zaloc1F(nx*nz);
    sct2   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
    lap    = CPU_zaloc1F(nx*nz);
    sigX   = CPU_zaloc1F(nx);
    sigZ   = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);

    int mode, side;
    int store_inc = 1, store_sct = 1;

    mode=IC_BORNSCT_IMAGE;

    long long int ixMin = shot->ixMin;
    long long int ixMax = shot->ixMax;
    long long int nxProp = ixMax - ixMin + 1;

    // Initializing GPU environment for current shot
    GPU_OPER_initEnvironment(procGPU, nx, nz, nt, dx, dz, dt, shot->nsou, shot->nrec, eimsh->lx0, eimsh->lz0, 0, nxb, nzb, mode_cigs);
    GPU_OPER_alocArrays_acoustic();
    GPU_OPER_alocArrays_acoustic_ExtBornSct();
    if(mode_cigs == MODE_CIGS_FULL) {
        GPU_OPER_copyEImage_CPU2GPU(eref);
    }
    else if(mode_cigs == MODE_CIGS_ORTH) {
        GPU_OPER_copyOCIGS_CPU2GPU(Heref, Veref);
    }

    GPU_OPER_zeroArrays_acoustic();
    GPU_OPER_zeroArrays_ExtBornSct();
    
    GPU_OPER_getVelP_CPU2GPU(vp);
    GPU_OPER_getSigmas_CPU2GPU(sigma, sigmaInv);
    GPU_OPER_getSourceWavelet_CPU2GPU(shot->source);
    GPU_OPER_getRecordedData_CPU2GPU(shot->seismogram);
    GPU_OPER_getSourceReceiverIndexes_CPU2GPU(shot->ixs, shot->izs, shot->ixr, shot->izr);

    GPU_OPER_alocArrays_EIC(jt, store_inc, store_sct);

    // Source wavefield
    int recordBornModeling = ((shot->ishot==ishotMovie) && recordMovie);

    float  *inc_sou, *inc_rec, *sct_sou, *sct_rec, *img_sou, *img_rec;
    if(recordBornModeling)
    {
        inc_sou   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
        inc_rec   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
        sct_sou   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
        sct_rec   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
        img_sou   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
        img_rec   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
    }

    FILE *fimg_sou_Hdr, *fimg_sou_Bin;
    FILE *fimg_rec_Hdr, *fimg_rec_Bin;
    FILE *fWavSouInc_Bin, *fWavSouInc_Hdr;
    FILE *fWavSouSct_Bin, *fWavSouSct_Hdr;
    FILE *fWavRecInc_Bin, *fWavRecInc_Hdr;
    FILE *fWavRecSct_Bin, *fWavRecSct_Hdr;
    char movieIncSouFileName[1024];
    char movieSctSouFileName[1024];
    char movieIncRecFileName[1024];
    char movieSctRecFileName[1024];
    char movieImgSouFileName[1024];
    char movieImgRecFileName[1024];
    if(recordBornModeling)
    {
        sprintf(movieIncSouFileName, "movie_inc_extBorn_sou_ishot%d", ishotMovie);
        sprintf(movieSctSouFileName, "movie_sct_extBorn_sou_ishot%d", ishotMovie);
        sprintf(movieImgSouFileName, "movie_img_extBorn_sou_ishot%d", ishotMovie);
        sprintf(movieIncRecFileName, "movie_inc_extBorn_rec_ishot%d", ishotMovie);
        sprintf(movieSctRecFileName, "movie_sct_extBorn_rec_ishot%d", ishotMovie);
        sprintf(movieImgRecFileName, "movie_img_extBorn_rec_ishot%d", ishotMovie);
        initFilesSmart3d(movieIncSouFileName, &fWavSouInc_Hdr, &fWavSouInc_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieSctSouFileName, &fWavSouSct_Hdr, &fWavSouSct_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieIncRecFileName, &fWavRecInc_Hdr, &fWavRecInc_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieSctRecFileName, &fWavRecSct_Hdr, &fWavRecSct_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieImgSouFileName, &fimg_sou_Hdr  , &fimg_sou_Bin  , nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieImgRecFileName, &fimg_rec_Hdr  , &fimg_rec_Bin  , nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
    }

    // Forward
    for(it=0; it<nt; it++)
    {
        GPU_OPER_setTimeStep(it);

        // printf("\n Propagation forward in time it=%lld \n", it);

        if(it%jt==0)
        {
            long long int shift = nx*nz * (it/jt);
            GPU_WRAP_imagingCondition(IC_BORNSCT_ILUMIN);

            if(recordBornModeling)
            {
                GPU_OPER_copyPresentWavefield_secDer_GPU2CPU(inc_sou+shift);
                GPU_OPER_copyPresentScatteredWavefield_secDer_GPU2CPU(sct_sou+shift);
            }
        }
        
        side = IC_SOU_WAVE;
        GPU_OPER_EIC_storeWavefields_GPU2GPU(it, side, mode);
        
        GPU_OPER_wavefieldTimeSecondDerivative_part1(0);
        GPU_OPER_wavefieldTimeSecondDerivative_part1(1);
        GPU_WRAP_injectSource_Acoustic();
        GPU_WRAP_propagate_Acoustic_abs2();
        GPU_OPER_wavefieldTimeSecondDerivative_part2(0);
        if(mode_cigs == MODE_CIGS_FULL)
        {
            if(eimsh->lz0==0)    GPU_WRAP_scatter_Acoustic_lx_opt();
            else                 GPU_WRAP_scatter_Acoustic_lx_lz_opt();
        }
        else if(mode_cigs == MODE_CIGS_ORTH)
        {
            GPU_WRAP_scatter_Acoustic_OCIGS_opt(EIMSH_CIGS_VCIG);
        }
        GPU_WRAP_propagate_Acoustic_Scattered_abs2();
        GPU_OPER_wavefieldTimeSecondDerivative_part2(1);
    }

    GPU_OPER_zeroWavefieldArrays();
    GPU_OPER_zeroScatteredWavefieldArrays();

    // printf("\n Started propagation backward in time \n");

    // Backward
    for(it=nt-1; it>=0; it--)
    {
        GPU_OPER_setTimeStep(it);
        
        // printf("\n Propagation backward in time it=%lld \n", it);

        long long int shift = nx*nz * (it/jt);
        GPU_WRAP_imagingCondition_smart(IC_BORNSCT_IMAGE, it);

        if(it%jt==0)
        {
            if(recordBornModeling)
            {
                GPU_OPER_copyPresentWavefield_GPU2CPU(inc_rec+shift);
                GPU_OPER_copyPresentScatteredWavefield_GPU2CPU(sct_rec+shift);
                GPU_OPER_copyImage_BRNSCT_GPU2CPU(img_sou+shift, img_rec+shift, 0);
            }
        }

        GPU_OPER_wavefieldTimeSecondDerivative_part1(0);
        GPU_WRAP_recordInjectData(1);
        GPU_WRAP_propagate_Acoustic_abs2();
        GPU_OPER_wavefieldTimeSecondDerivative_part2(0);
        if(mode_cigs == MODE_CIGS_FULL)
        {
            if(eimsh->lz0==0)    GPU_WRAP_scatter_Acoustic_lx_opt();
            else                 GPU_WRAP_scatter_Acoustic_lx_lz_opt();
        }
        else if(mode_cigs == MODE_CIGS_ORTH)
        {
            GPU_WRAP_scatter_Acoustic_OCIGS_opt(EIMSH_CIGS_VCIG);
        }
        GPU_WRAP_propagate_Acoustic_Scattered_abs2();
    }

    // printf("\n Finished backprop  \n");

    GPU_OPER_freeArrays_EIC();
    GPU_OPER_copyImage_BRNSCT_GPU2CPU(imageSouSide, imageRecSide, 1);
    GPU_OPER_copyIlumin_BRNSCT_GPU2CPU(iluminSouSide, iluminRecSide, 1);
    GPU_OPER_freeArrays_acoustic_ExtBornSct();
    GPU_OPER_freeArrays_acoustic();

    if(recordBornModeling)
    {
        printf("\n Saving files \n");
        for(it=0; it<nt/jt; it++)
        {
            long long int shift = it*nx*nz;
            writeSamples(inc_sou+shift, fWavSouInc_Bin, nx*nz);
            writeSamples(sct_sou+shift, fWavSouSct_Bin, nx*nz);
            writeSamples(inc_rec+shift, fWavRecInc_Bin, nx*nz);
            writeSamples(sct_rec+shift, fWavRecSct_Bin, nx*nz);
            // writeSamples(img_sou+shift,   fimg_sou_Bin, nx*nz);
            // writeSamples(img_rec+shift,   fimg_rec_Bin, nx*nz);
        }

        fclose(fWavSouInc_Bin);
        fclose(fWavSouInc_Hdr);

        fclose(fWavSouSct_Bin);
        fclose(fWavSouSct_Hdr);

        fclose(fWavRecInc_Bin);
        fclose(fWavRecInc_Hdr);

        fclose(fWavRecSct_Hdr);
        fclose(fWavRecSct_Bin);

        fclose(fimg_sou_Bin);
        fclose(fimg_sou_Hdr);

        fclose(fimg_rec_Bin);
        fclose(fimg_rec_Hdr);

        free(inc_sou);
        free(inc_rec);
        free(sct_sou);
        free(sct_rec);
        free(img_sou);
        free(img_rec);
    }

    
    free(sigma);
    free(sigX);
    free(sigZ);
    free(ws);
    free(wc);
    free(wstp);
    free(wav1);
    free(wav2);
    free(sct1);
    free(sct2);
    free(lap);
    free(sigmaInv);

return;
}

void migrateExtBorn_lx_lz_lt_Mem_GPU_OPTIM(int idev, float *imageSouSide, float *imageRecSide, \
                                           float *imageSouSide_shot, float *imageRecSide_shot, \
                                           float *iluminSouSide, float *iluminRecSide, \
                                           float *eref, float *Heref, float *Veref, float *vp, Shot *shot, \
                                           long long int jt, long long int jt_EIC, long long int nt, \
                                           long long int nx, long long int nz, float dx, float dz, \
                                           float dt, float ox, float oz, float ot, long long int lx0, long long int lz0, long long int lt0, \
                                           int ishotMovie, int recordMovie, int imageDomain) {

    int mode_cigs;
    if(Heref!=NULL  &&  Veref!=NULL)
        mode_cigs = MODE_CIGS_ORTH;
    else 
        mode_cigs = MODE_CIGS_FULL;
    
    long long int     it;
    size_t  err;
    float  *wav1, *wav2;
    float  *sct1, *sct2;
    float  *lap, *sigX, *sigZ, *sigma, *sigmaInv, *tmp, *ws, *wc, *secDer, *wstp;
    secDer = CPU_zaloc1F(nx*nz);
    wav1   = CPU_zaloc1F(nx*nz);
    wav2   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
    ws     = CPU_zaloc1F(nx*nz);
    wstp   = CPU_zaloc1F((nx-2*nxb)*(nz-2*nzb));
    wc     = CPU_zaloc1F(nx*nz);
    sct1   = CPU_zaloc1F(nx*nz);
    sct2   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
    lap    = CPU_zaloc1F(nx*nz);
    sigX   = CPU_zaloc1F(nx);
    sigZ   = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);

    int mode, side;
    int store_inc = 1, store_sct = 1;
    if(imageDomain==0)    mode=IC_BORNSCT_IMAGE2;
    else                  mode=IC_BORNSCT_IMAGE;

    long long int ixMin = shot->ixMin;
    long long int ixMax = shot->ixMax;
    long long int nxProp = ixMax - ixMin + 1;

    // Initializing GPU environment for current shot
    GPU_OPER_zeroArrays_acoustic();
    GPU_OPER_zeroArrays_ExtBornSct();
    
    GPU_OPER_getVelP_CPU2GPU(vp);
    GPU_OPER_getSigmas_CPU2GPU(sigma, sigmaInv);
    GPU_OPER_getSourceWavelet_CPU2GPU(shot->source);
    GPU_OPER_getRecordedData_CPU2GPU(shot->seismogram);
    GPU_OPER_getSourceReceiverIndexes_CPU2GPU(shot->ixs, shot->izs, shot->ixr, shot->izr);

    GPU_OPER_alocArrays_EIC(jt, store_inc, store_sct);

    // Source wavefield
    int recordBornModeling = ((shot->ishot==ishotMovie) && recordMovie);

    float  *inc_sou, *inc_rec, *sct_sou, *sct_rec, *img_sou, *img_rec;
    if(recordBornModeling)
    {
        inc_sou   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
        inc_rec   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
        sct_sou   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
        sct_rec   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
        img_sou   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
        img_rec   = CPU_zaloc1F(nx*nz*(1+(nt/jt)));
    }

    FILE *fimg_sou_Hdr, *fimg_sou_Bin;
    FILE *fimg_rec_Hdr, *fimg_rec_Bin;
    FILE *fWavSouInc_Bin, *fWavSouInc_Hdr;
    FILE *fWavSouSct_Bin, *fWavSouSct_Hdr;
    FILE *fWavRecInc_Bin, *fWavRecInc_Hdr;
    FILE *fWavRecSct_Bin, *fWavRecSct_Hdr;
    char movieIncSouFileName[1024];
    char movieSctSouFileName[1024];
    char movieIncRecFileName[1024];
    char movieSctRecFileName[1024];
    char movieImgSouFileName[1024];
    char movieImgRecFileName[1024];
    if(recordBornModeling)
    {
        sprintf(movieIncSouFileName, "movie_inc_extBorn_sou_ishot%d", ishotMovie);
        sprintf(movieSctSouFileName, "movie_sct_extBorn_sou_ishot%d", ishotMovie);
        sprintf(movieImgSouFileName, "movie_img_extBorn_sou_ishot%d", ishotMovie);
        sprintf(movieIncRecFileName, "movie_inc_extBorn_rec_ishot%d", ishotMovie);
        sprintf(movieSctRecFileName, "movie_sct_extBorn_rec_ishot%d", ishotMovie);
        sprintf(movieImgRecFileName, "movie_img_extBorn_rec_ishot%d", ishotMovie);
        initFilesSmart3d(movieIncSouFileName, &fWavSouInc_Hdr, &fWavSouInc_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieSctSouFileName, &fWavSouSct_Hdr, &fWavSouSct_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieIncRecFileName, &fWavRecInc_Hdr, &fWavRecInc_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieSctRecFileName, &fWavRecSct_Hdr, &fWavRecSct_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieImgSouFileName, &fimg_sou_Hdr  , &fimg_sou_Bin  , nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
        initFilesSmart3d(movieImgRecFileName, &fimg_rec_Hdr  , &fimg_rec_Bin  , nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
    }

    // Forward
    for(it=0; it<nt; it++)
    {
        GPU_OPER_setTimeStep(it);

        // printf("\n Propagation forward in time it=%lld \n", it);

        if(it%jt==0)
        {
            long long int shift = nx*nz * (it/jt);
            GPU_WRAP_imagingCondition(IC_BORNSCT_ILUMIN);

            if(recordBornModeling)
            {
                GPU_OPER_copyPresentWavefield_secDer_GPU2CPU(inc_sou+shift);
                GPU_OPER_copyPresentScatteredWavefield_secDer_GPU2CPU(sct_sou+shift);
            }
        }
        
        side = IC_SOU_WAVE;
        GPU_OPER_EIC_storeWavefields_GPU2GPU(it, side, mode);
        
        GPU_OPER_wavefieldTimeSecondDerivative_part1(0);
        GPU_OPER_wavefieldTimeSecondDerivative_part1(1);
        GPU_WRAP_injectSource_Acoustic();
        GPU_WRAP_propagate_Acoustic_abs2();
        GPU_OPER_wavefieldTimeSecondDerivative_part2(0);
        if(mode_cigs == MODE_CIGS_FULL)
        {
            if(lz0==0)    GPU_WRAP_scatter_Acoustic_lx_opt();
            else          GPU_WRAP_scatter_Acoustic_lx_lz_opt();
        }
        else if(mode_cigs == MODE_CIGS_ORTH)
        {
            GPU_WRAP_scatter_Acoustic_OCIGS_opt(EIMSH_CIGS_VCIG);
        }
        GPU_WRAP_propagate_Acoustic_Scattered_abs2();
        GPU_OPER_wavefieldTimeSecondDerivative_part2(1);
    }

    GPU_OPER_zeroWavefieldArrays();
    GPU_OPER_zeroScatteredWavefieldArrays();

    // printf("\n Started propagation backward in time \n");

    // Backward
    for(it=nt-1; it>=0; it--)
    {
        GPU_OPER_setTimeStep(it);
        
        // printf("\n Propagation backward in time it=%lld \n", it);

        long long int shift = nx*nz * (it/jt);
        GPU_WRAP_imagingCondition_smart(IC_BORNSCT_IMAGE, it);

        if(it%jt==0)
        {
            if(recordBornModeling)
            {
                GPU_OPER_copyPresentWavefield_GPU2CPU(inc_rec+shift);
                GPU_OPER_copyPresentScatteredWavefield_GPU2CPU(sct_rec+shift);
                GPU_OPER_copyImage_BRNSCT_GPU2CPU(img_sou+shift, img_rec+shift, 0);
            }
        }

        GPU_OPER_wavefieldTimeSecondDerivative_part1(0);
        GPU_WRAP_recordInjectData(1);
        GPU_WRAP_propagate_Acoustic_abs2();
        GPU_OPER_wavefieldTimeSecondDerivative_part2(0);
        if(mode_cigs == MODE_CIGS_FULL)
        {
            if(lz0==0)    GPU_WRAP_scatter_Acoustic_lx_opt();
            else          GPU_WRAP_scatter_Acoustic_lx_lz_opt();
        }
        else if(mode_cigs == MODE_CIGS_ORTH)
        {
            GPU_WRAP_scatter_Acoustic_OCIGS_opt(EIMSH_CIGS_VCIG);
        }
        GPU_WRAP_propagate_Acoustic_Scattered_abs2();
    }

    // printf("\n Finished backprop  \n");

    GPU_OPER_freeArrays_EIC();

    if(recordBornModeling)
    {
        printf("\n Saving files \n");
        for(it=0; it<nt/jt; it++)
        {
            long long int shift = it*nx*nz;
            writeSamples(inc_sou+shift, fWavSouInc_Bin, nx*nz);
            writeSamples(sct_sou+shift, fWavSouSct_Bin, nx*nz);
            writeSamples(inc_rec+shift, fWavRecInc_Bin, nx*nz);
            writeSamples(sct_rec+shift, fWavRecSct_Bin, nx*nz);
            // writeSamples(img_sou+shift,   fimg_sou_Bin, nx*nz);
            // writeSamples(img_rec+shift,   fimg_rec_Bin, nx*nz);
        }

        fclose(fWavSouInc_Bin);
        fclose(fWavSouInc_Hdr);

        fclose(fWavSouSct_Bin);
        fclose(fWavSouSct_Hdr);

        fclose(fWavRecInc_Bin);
        fclose(fWavRecInc_Hdr);

        fclose(fWavRecSct_Hdr);
        fclose(fWavRecSct_Bin);

        fclose(fimg_sou_Bin);
        fclose(fimg_sou_Hdr);

        fclose(fimg_rec_Bin);
        fclose(fimg_rec_Hdr);

        free(inc_sou);
        free(inc_rec);
        free(sct_sou);
        free(sct_rec);
        free(img_sou);
        free(img_rec);
    }

    
    free(sigma);
    free(sigX);
    free(sigZ);
    free(ws);
    free(wc);
    free(wstp);
    free(wav1);
    free(wav2);
    free(sct1);
    free(sct2);
    free(lap);
    free(sigmaInv);

return;
}
void modelShotBorn_lx(float *eref, float *vp, Shot *shot, int jt, int nt, int nx, int nz, float dx, float dz, \
                   float dt, float ox, float oz, float ot, int lx0) {
    int it;
    float *wav1, *wav2, *secDer;
    float *sct1, *sct2;
    float *lap, *sigX, *sigZ, *sigma, *sigmaInv, *tmp;
    secDer = CPU_zaloc1F(nx*nz);
    wav1   = CPU_zaloc1F(nx*nz);
    wav2   = CPU_zaloc1F(nx*nz);
    sct1   = CPU_zaloc1F(nx*nz);
    sct2   = CPU_zaloc1F(nx*nz);
    lap    = CPU_zaloc1F(nx*nz);
    sigX   = CPU_zaloc1F(nx);
    sigZ   = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);

    // applyHalfDerivative2Source(shot);

    // Source wavefield
    //int recordBornModeling = (shot->ishot==(shot->nshot/2));
    int recordBornModeling = (shot->ishot==0);
    FILE *fwav_Hdr, *fwav_Bin;
    FILE *fsct_Hdr, *fsct_Bin;
    if(recordBornModeling) {
        initFilesSmart3d("movie_incident_extendedBorn" , &fwav_Hdr, &fwav_Bin, nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
        initFilesSmart3d("movie_scattered_extendedBorn", &fsct_Hdr, &fsct_Bin, nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
    }
    for(it=0; it<nt; it++)
    {
        recordData(shot, sct2, nx, nz, it);

        if(it%jt==0 && recordBornModeling) {
            writeSamples(wav2, fwav_Bin, nx*nz);
            writeSamples(sct2, fsct_Bin, nx*nz);
        }

        /*
        injectSource_Acoustic(wav1, shot->source, it, (shot->ixs)[0], (shot->izs)[0], nx, nz, dx, dz);
        scatter_Acoustic_lx(eref, sct1, 2.0f, wav2, -1.0f, wav1, nx, nz, dt, lx0);
        propagate_Acoustic_abs(wav1, wav2, lap, vp, sigma, nx, nz);
        scatter_Acoustic_lx(eref, sct1, -1.0f, wav2, 0.0f, wav1, nx, nz, dt, lx0);
        propagate_Acoustic_abs(sct1, sct2, lap, vp, sigma, nx, nz);
        */
        zeroArray(secDer,nx*nz);
        wavefieldTimeSecondDerivative(secDer, -2.0f, wav2, +1.0f, wav1, nx, nz, dt);
        injectSource_Acoustic(wav1, shot->source, it, (shot->ixs)[0], (shot->izs)[0], nx, nz, dx, dz);
        propagate_Acoustic_abs(wav1, wav2, lap, vp, sigma, nx, nz);
        wavefieldTimeSecondDerivative(secDer, +1.0f, wav2,  0.0f, wav1, nx, nz, dt);

        scatter_Acoustic_lx_opt(eref, sct1, secDer, nx, nz, dt, lx0);
        propagate_Acoustic_abs(sct1, sct2, lap, vp, sigma, nx, nz);
        
        tmp  = wav2;
        wav2 = wav1;
        wav1 = tmp;
        tmp  = sct2;
        sct2 = sct1;
        sct1 = tmp;
    }
    if(recordBornModeling) {
        fclose(fwav_Bin);
        fclose(fwav_Hdr);
        fclose(fsct_Bin);
        fclose(fsct_Hdr);
    }
    free(sigma);
    free(sigX);
    free(sigZ);
    free(wav1);
    free(wav2);
    free(sct1);
    free(sct2);
    free(lap);
    free(sigmaInv);

return;
}
void modelShotBorn_lx_lz_lt_GPU(int idev, float *eref, float *Heref, float *Veref, float *vp, Shot *shot, int jt, int nt, int nx, int nz, float dx, float dz, \
                                float dt, float ox, float oz, float ot, int lx0, int lz0, int lt0, int recordMovie, int shotToRecord, int OPTIM) {
    int mode_cigs;
    if(Heref!=NULL  &&  Veref!=NULL)
        mode_cigs = MODE_CIGS_ORTH;
    else 
        mode_cigs = MODE_CIGS_FULL;
    
    int it;
    float *wav1, *wav2, *secDer;
    float *sct1, *sct2;
    float *lap, *sigX, *sigZ, *sigma, *sigmaInv, *tmp;
    secDer = CPU_zaloc1F(nx*nz);
    wav1   = CPU_zaloc1F(nx*nz);
    wav2   = CPU_zaloc1F(nx*nz);
    sct1   = CPU_zaloc1F(nx*nz);
    sct2   = CPU_zaloc1F(nx*nz);
    lap    = CPU_zaloc1F(nx*nz);
    sigX   = CPU_zaloc1F(nx);
    sigZ   = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);

    // Initializing GPU environment for current shot
    if(!OPTIM)
    {
        GPU_OPER_initEnvironment(idev, nx, nz, nt, dx, dz, dt, shot->nsou, shot->nrec, lx0, lz0, lt0, nxb, nzb, mode_cigs);
        GPU_OPER_alocArrays_acoustic();
        GPU_OPER_alocArrays_acoustic_ExtBornSct();
        if(mode_cigs == MODE_CIGS_FULL) {
            GPU_OPER_copyEImage_CPU2GPU(eref);
        }
        else if(mode_cigs == MODE_CIGS_ORTH) {
            GPU_OPER_copyOCIGS_CPU2GPU(Heref, Veref);
        }
    }
    else
    {
        GPU_OPER_zeroArrays_acoustic();
        GPU_OPER_zeroArrays_ExtBornSct();
    }
    GPU_OPER_getVelP_CPU2GPU(vp);
    GPU_OPER_getSigmas_CPU2GPU(sigma, sigmaInv);
    GPU_OPER_getSourceWavelet_CPU2GPU(shot->source);
    GPU_OPER_getSourceReceiverIndexes_CPU2GPU(shot->ixs, shot->izs, shot->ixr, shot->izr);

    // Source wavefield
    int recordBornModeling = ((shot->ishot==shotToRecord) && recordMovie);
    FILE *fwav_Hdr, *fwav_Bin;
    FILE *fsct_Hdr, *fsct_Bin;
    if(recordBornModeling)
    {
        initFilesSmart3d("movie_incident_extendedBorn" , &fwav_Hdr, &fwav_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
        initFilesSmart3d("movie_scattered_extendedBorn", &fsct_Hdr, &fsct_Bin, nx, dx, ox, nz, dz, oz, nt/jt, dt*jt, ot);
    }
    for(it=0; it<nt; it++)
    {
        GPU_OPER_setTimeStep(it);
        GPU_WRAP_recordInjectData(2);

        if(it%jt==0 && recordBornModeling)
        {
            GPU_OPER_copyPresentWavefield_GPU2CPU(wav2);
            GPU_OPER_copyPresentScatteredWavefield_GPU2CPU(sct2);
            writeSamples(wav2, fwav_Bin, nx*nz);
            writeSamples(sct2, fsct_Bin, nx*nz);
        }
        GPU_OPER_wavefieldTimeSecondDerivative_part1(0);
        GPU_WRAP_injectSource_Acoustic();
        GPU_WRAP_propagate_Acoustic_abs2();
        GPU_OPER_wavefieldTimeSecondDerivative_part2(0);
        if(mode_cigs == MODE_CIGS_FULL)
        {
            if(lz0==0)    GPU_WRAP_scatter_Acoustic_lx_opt();
            else          GPU_WRAP_scatter_Acoustic_lx_lz_opt();
        }
        else if(mode_cigs == MODE_CIGS_ORTH)
        {
            GPU_WRAP_scatter_Acoustic_OCIGS_opt(EIMSH_CIGS_VCIG);
        }
        GPU_WRAP_propagate_Acoustic_Scattered_abs2();
    }
    if(recordBornModeling)
    {
        fclose(fwav_Bin);
        fclose(fwav_Hdr);
        fclose(fsct_Bin);
        fclose(fsct_Hdr);
    }
    free(sigma);
    free(sigX);
    free(sigZ);
    free(wav1);
    free(wav2);
    free(sct1);
    free(sct2);
    free(lap);

    GPU_OPER_getRecordedData_GPU2CPU(shot->seismogram);
    if(!OPTIM)
    {
        GPU_OPER_freeArrays_acoustic_ExtBornSct();
        GPU_OPER_freeArrays_acoustic();
    }
    free(sigmaInv);

return;
}

void modelShotBorn_ampAngle(float *ref, float *vp, Shot *shot, int jt, int nt, int nx, int nz, float dx, float dz, \
                            float dt, float ox, float oz, float ot, int lx0, int ntheta) {
    int it;
    int nlx=2*lx0+1;
    int itheta;
    float *wav1, *wav2;
    float *sct1, *sct2;
    float *lap, *sigX, *sigZ, *sigma, *sigmaInv, *tmp;
    wav1 = CPU_zaloc1F(nx*nz);
    wav2 = CPU_zaloc1F(nx*nz);
    sct1 = CPU_zaloc1F(nx*nz);
    sct2 = CPU_zaloc1F(nx*nz);
    lap  = CPU_zaloc1F(nx*nz);
    sigX = CPU_zaloc1F(nx);
    sigZ = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);

    applyHalfDerivative2Source(shot);

    // Source wavefield
    //int recordBornModeling = (shot->ishot==(shot->nshot/2));
    int recordBornModeling = (shot->ishot==0);
    FILE *fwav_Hdr, *fwav_Bin;
    FILE *fsct_Hdr, *fsct_Bin;
    if(recordBornModeling) {
        initFilesSmart3d("movie_incident_extendedBorn" , &fwav_Hdr, &fwav_Bin, nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
        initFilesSmart3d("movie_scattered_extendedBorn", &fsct_Hdr, &fsct_Bin, nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
    }
    for(it=0; it<nt; it++)
    {
        recordData(shot, sct2, nx, nz, it);

        if(it%jt==0 && recordBornModeling) {
            writeSamples(wav2, fwav_Bin, nx*nz);
            writeSamples(sct2, fsct_Bin, nx*nz);
        }

        injectSource_Acoustic(wav1, shot->source, it, (shot->ixs)[0], (shot->izs)[0], nx, nz, dx, dz);
        for(itheta=0; itheta<ntheta; itheta+=5)
            scatter_Acoustic_AmpAngle(ref+(itheta*nlx*nx*nz), sct1, 2.0f, wav2,   -1.0f,  wav1, nx, nz, dt, lx0);
        
        propagate_Acoustic_abs(wav1, wav2, lap, vp, sigma, nx, nz);
        for(itheta=0; itheta<ntheta; itheta+=5)
            scatter_Acoustic_AmpAngle(ref+(itheta*nlx*nx*nz), sct1, -1.0f, wav2, 0.0f, wav1, nx, nz, dt, lx0);
        
        propagate_Acoustic_abs(sct1, sct2, lap, vp, sigma, nx, nz);

        tmp = wav2;
        wav2=wav1;
        wav1=tmp;
        tmp = sct2;
        sct2=sct1;
        sct1=tmp;
    }
    if(recordBornModeling) {
        fclose(fwav_Bin);
        fclose(fwav_Hdr);
        fclose(fsct_Bin);
        fclose(fsct_Hdr);
    }
    free(sigma);
    free(sigX);
    free(sigZ);
    free(wav1);
    free(wav2);
    free(sct1);
    free(sct2);
    free(lap);
    free(sigmaInv);

return;
}