void modelShot_elastic(float *lambda, float *mu, float *lambda2mu, float *bc, float *bd, \
                       Shot *shot, int jt, int nt, int nx, int nz, \
                       float dx, float dz, float dt, float ox, float oz, float ot, int recordMovie, int ishotMovie) {
    int it;
    float *vx, *vz, *Sxx, *Szz, *Sxz, *lap, *sigX, *sigZ, *sigma, *sigmaInv;
    vx   = CPU_zaloc1F(nx*nz);
    vz   = CPU_zaloc1F(nx*nz);
    lap  = CPU_zaloc1F(nx*nz);
    Sxx  = CPU_zaloc1F(nx*nz);
    Szz  = CPU_zaloc1F(nx*nz);
    Sxz  = CPU_zaloc1F(nx*nz);
    sigX = CPU_zaloc1F(nx);
    sigZ = CPU_zaloc1F(nz);
    getSigmas(sigX, sigZ, nxb, OL, nx, nz, dx, dz);
    getVelSigma(&sigma, &sigmaInv, sigX, sigZ, nx, nz, dx, dt);

    applyHalfDerivative2Source(shot);

    FILE *fWavMovie_Hdr, *fWavMovie_Bin;
    if(shot->ishot==ishotMovie && recordMovie) {
        initFiles3d("./movie_modeling_elastic_vz", &fWavMovie_Hdr, &fWavMovie_Bin, nz, dz, oz, nx, dx, ox, nt/jt, dt*jt, ot);
    }
    for(it=0; it<nt; it++)
    {
        //recordData(shot, vz, nx, nz, it);
        //recordData_Elastic_Omni(shot, Sxx, Szz, nx, nz, it);
        recordData_Elastic_phi(shot, vx, vz, nx, nz, dx, dz, it);
        //recordData_Elastic_Vertical(shot, vz, nx, nz, it);
        if(it%jt==0  &&  shot->ishot==ishotMovie && recordMovie)
            writeSamples(vz, fWavMovie_Bin, nx*nz);
        
        //injectSource_elastic_Szz(Szz, shot->source, it, (shot->ixs)[0], (shot->izs)[0], nx, nz, dx, dz);
        injectSource_elastic_Omni(Sxx, Szz, shot->source, it, (shot->ixs)[0], (shot->izs)[0], nx, nz, dx, dz);
        //propagate_Elastic(vx, vz, Sxx, Szz, Sxz, lap, lambda, mu, b, nx, nz, dx, dt);
        propagate_Elastic_abs(vx, vz, Sxx, Szz, Sxz, lap, lambda, mu, lambda2mu, bc, bd, nx, nz, dx, dt, sigX, sigZ);
    }
    integrateRecTime(shot);
    if(shot->ishot==ishotMovie && recordMovie) {
        fclose(fWavMovie_Bin);
        fclose(fWavMovie_Hdr);
    }
    free(sigma);
    free(sigX);
    free(sigZ);
    free(vx);
    free(vz);
    free(Sxx);
    free(Szz);
    free(Sxz);
    free(lap);

return;
}