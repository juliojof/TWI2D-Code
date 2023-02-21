void transformVel(float *m, int nx, int nz, float dx, float dt) {
    int ix, iz, i;
    for(ix=0; ix<nx; ix++) {
        for(iz=0; iz<nz; iz++) {
            i = iz + ix*nz;
            m[i] *= dt/dx;
            m[i] *= m[i];
        }
    }
return;
}
void untransformVel(float *m, int nx, int nz, float dx, float dt) {
    int ix, iz, i;
    for(ix=0; ix<nx; ix++) {
        for(iz=0; iz<nz; iz++) {
            i = iz + ix*nz;
            m[i] = sqrtf(m[i]);
            m[i] *= dx/dt;
        }
    }
return;
}

void transformVelSlowness(float *m, long long int n) {
    int i;
    for(i=0; i<n; i++)    m[i] = 1.0/m[i];
return;
}

void model_vel2SquaredVaga(float *m, int nx, int nz) {
    int ix, iz, i;
    for(ix=0; ix<nx; ix++) {
        for(iz=0; iz<nz; iz++)
        {
            i = iz + ix*nz;
            m[i] = 1.0 / (m[i]*m[i]);
        }
    }
return;
}
void model_squaredVaga2Vel(float *m, int nx, int nz) {
    int ix, iz, i;
    for(ix=0; ix<nx; ix++) {
        for(iz=0; iz<nz; iz++)
        {
            i = iz + ix*nz;
            m[i] = 1.0 / sqrtf(m[i]);
        }
    }
return;
}

float* getVp_1Reflector(int nx, int nxb, int  nz, int nzb, float dx, float dz, float vp0, float dvp, float zref, float lref) {
    
    float *m;
    int ix, iz, i;
    float x, z;

    m = CPU_zaloc1F(nx*nz);

    for(ix=0; ix<nx; ix++) {
        x = ix*dx;
        for(iz=0; iz<nz; iz++) {
            z = (iz-nzb)*dz;
            i = iz + ix*nz;
            m[i] = vp0;
            if(zref<=z && z<=zref+lref) m[i] += dvp;
        }
    }


return m;    
}
void getVp_add1Reflector(float *m, int nx, int nxb, int  nz, int nzb, float dx, float dz, float dvp, float zref, float lref) {
    
    int ix, iz, i;
    float x, z;

    for(ix=0; ix<nx; ix++) {
        x = ix*dx;
        for(iz=0; iz<nz; iz++) {
            z = (iz-nzb)*dz;
            i = iz + ix*nz;
            if(zref<=z && z<=zref+lref) m[i] += dvp;
        }
    }
    
return;    
}

void getVp_add1DippingReflector(float *m, int nx, int nxb, int  nz, int nzb, float dx, float dz, float ox, float oz, \
                                float dvp, float xleft, float xright, float zleft, float zright, float width) {
    
    int ix, iz, i;
    float x, z;
    float ztop, zbot;

    float dz_dx = (zright-zleft)/(xright-xleft);

    for(ix=0; ix<nx; ix++)
    {
        x = ix*dx + ox;

        float deltaX = (x-xleft);
        ztop = deltaX * dz_dx + zleft;
        zbot = ztop + width;

        // if(ztop <  oz             )    ztop = oz;
        // if(ztop > (nz-2*nzb-1)*dz )    ztop = (nz-2*nzb-1)*dz;
        // if(zbot <  oz             )    zbot = oz;
        // if(zbot > (nz-2*nzb-1)*dz )    zbot = (nz-2*nzb-1)*dz;

        for(iz=0; iz<nz; iz++)
        {
            z = iz*dz + oz;
            i = iz + ix*nz;
            float deltaVp;
            
            if(z < ztop-dz)                          deltaVp = 0.0;
            else if(z<ztop  &&  fabsf(z-ztop)<dz)    deltaVp = dvp * (1.0-fabsf(z-ztop)/dz);
            else if(ztop<=z && z<=zbot)              deltaVp = dvp;
            else if(z>zbot  &&  fabsf(z-zbot)<dz)    deltaVp = dvp * (1.0-fabsf(z-zbot)/dz);
            else if(z > zbot+dz)                     deltaVp = 0.0;
            
            if(x<xleft  ||  x>xright)                deltaVp = 0.0;
            
            m[i] += deltaVp;
        }
    }

    /*
    float *mtmp = CPU_zaloc1F(nx*nz);
    memcpy(mtmp,m,nx*nz*sizeof(float));
    for(ix=0; ix<nx; ix++)
    {
        x = ix*dx + ox;

        float deltaX = (x-xleft);
        ztop = deltaX * dz_dx + zleft;
        zbot = ztop + width;

        iztop = floorf(ztop/dz);

        for(iz=0; iz<nz; iz++)
        {
            z = iz*dz + oz;
            i = iz + ix*nz;
            float deltaVp;
            if(z < ztop-dz)                          deltaVp = 0.0;
            else if(z<ztop  &&  fabsf(z-ztop)<dz)    deltaVp = dvp * (1.0-fabsf(z-ztop)/dz);
            else if(ztop<=z && z<=zbot)              deltaVp = dvp;
            else if(z>zbot  &&  fabsf(z-zbot)<dz)    deltaVp = dvp * (1.0-fabsf(z-zbot)/dz);
            if(z > zbot+dz)                          deltaVp = 0.0;
            
            m[i] += deltaVp;
        }
    }
    //*/
    
return;
}


float* getVp_1Reflector_1GaussAnomaly(int nx, int nxb, int  nz, int nzb, float dx, float dz, \
                                      float vp0, float dvp, float zref, float lref, \
                                      float dvpGauss, float gaussRadius, float x0, float z0) {
    
    float *m;
    int ix, iz, i;
    float x, z;
    float deltaX, deltaZ;
    float sigma, p, r2, ampGauss;
    p = 0.1f;
    // sigma = logf(1.0f/p)/gaussRadius*gaussRadius;
    sigma = -logf(p)/(gaussRadius*gaussRadius);

    m = CPU_zaloc1F(nx*nz);

    for(ix=0; ix<nx; ix++)
    {
        x = (ix-nxb)*dx;
        deltaX = x-x0;
        for(iz=0; iz<nz; iz++) {
            z = (iz-nzb)*dz;
            deltaZ = z-z0;
            //r2 = 0.09*deltaX*deltaX + deltaZ*deltaZ;
            r2 = deltaX*deltaX + deltaZ*deltaZ;
            
            ampGauss = dvpGauss * exp(-sigma*r2);
            if(sqrtf(r2)>gaussRadius)    ampGauss = 0.0;

            i = iz + ix*nz;
            m[i] = vp0 + ampGauss;

            if(zref<=z && z<=zref+lref) m[i] += dvp;
        }
    }


return m;    
}
void getVp_addGaussAnomaly(float *vp, int nx, int nxb, int  nz, int nzb, float dx, float dz, \
                           float dvpGauss, float gaussRadius, float x0, float z0) {    
    float *m;
    int ix, iz, i;
    float x, z;
    float deltaX, deltaZ;
    float sigma, p, r2, ampGauss;
    p = 0.1f;
    // sigma = logf(1.0f/p)/gaussRadius*gaussRadius;
    sigma = -logf(p)/(gaussRadius*gaussRadius);

    for(ix=0; ix<nx; ix++)
    {
        x = (ix-nxb)*dx;
        deltaX = x-x0;
        for(iz=0; iz<nz; iz++)
        {
            z = (iz-nzb)*dz;
            deltaZ = z-z0;
            r2 = deltaX*deltaX + deltaZ*deltaZ;

            ampGauss = dvpGauss * exp(-sigma*r2);
            if(sqrtf(r2)>gaussRadius)    ampGauss = 0.0;

            i = iz + ix*nz;
            vp[i] += ampGauss;
        }
    }

return;
}

float* getVp_VofZ(int nx, int nxb, int  nz, int nzb, float dx, float dz, \
                  float vp0, float coef, float zref, float lref) {
    
    float *m;
    int ix, iz, i;
    float x, z;
    float deltaZ;
    float dvz;

    m = CPU_zaloc1F(nx*nz);

    for(ix=0; ix<nx; ix++)
    {
        for(iz=0; iz<nz; iz++)
        {
            z = (iz-nzb)*dz;
            deltaZ = z-zref;
            
            if(z<zref)
                dvz = 0.0;
            else if(zref<=z && z<zref+lref)
                dvz = deltaZ * coef;  // coef: (m/s) / m
            if(z>=zref+lref)
                dvz = lref * coef;

            i = iz + ix*nz;
            m[i] = vp0 + dvz;
        }
    }


return m;    
}



void getLame_1Reflector(float *lambda, float *mu, float *lambda2mu, float *bc, float *bd, int nx, int nxb, int  nz, int nzb, float dx, float dz, \
                          float vp0, float dvp, float vs0, float dvs, float rho0, float drho, float zref, float lref) {
    int ix, iz, i;
    int i_d, i_x, i_z;
    float x, z;
    float vp, vs, rho;
    float *ll, *mm, *bb; 

    ll  = CPU_zaloc1F(nx*nz);
    mm  = CPU_zaloc1F(nx*nz);
    bb = CPU_zaloc1F(nx*nz);

    for(ix=0; ix<nx; ix++) {
        x = ix*dx;
        for(iz=0; iz<nz; iz++) {
            z = (iz-nzb)*dz;
            i = iz + ix*nz;
            vp  = vp0;
            vs  = vs0;
            rho = rho0;
            if(zref<=z && z<=zref+lref)
            {
                vp  += dvp;
                vs  += dvs;
                rho += drho;
            }
            ll[i] = rho*(vp*vp - 2.0f*vs*vs);
            mm[i] = rho*vs*vs;
            bb[i] = 1.0f/rho;
        }
    }
    for(ix=1; ix<nx-1; ix++) {
        for(iz=1; iz<nz-1; iz++) {
            i   =  iz +     ix   *nz;
            i_d = (iz+1) + (ix+1)*nz;
            i_z = (iz+1) +  ix   *nz;
            i_x =  iz    + (ix+1)*nz;
            lambda[i]    = 0.5*(ll[i] + ll[i_x]);
            lambda2mu[i] = lambda[i]  + 2.0*0.5*(mm[i] + mm[i_x]);
            mu[i]        = 0.5*(mm[i]+mm[i_z]);
            bc[i]        = bb[i];
            bd[i]        = 0.5*(bb[i]+bb[i_d]);
        }
    }

    free(bb);
    free(ll);
    free(mm);
    
return;
}
void getLame_1Reflector_1GaussAnomaly(float *lambda, float *mu, float *lambda2mu, float *bc, float *bd, \
                                      int nx, int nxb, int  nz, int nzb, float dx, float dz, \
                                      float vp0, float dvp, float vs0, float dvs, float rho0, float drho, float zref, float lref, \
                                      float dvpGauss, float dvsGauss, float gaussRadius, float x0, float z0) {
    
    float x, z;
    float vp, vs, rho;
    float *ll, *mm, *bb; 
    int ix, iz, i;
    int i_d, i_x, i_z;
    float deltaX, deltaZ;
    float sigma, p, r2, ampGauss;
    p = 0.1f;
    sigma = -logf(p)/(gaussRadius*gaussRadius);

    ll  = CPU_zaloc1F(nx*nz);
    mm  = CPU_zaloc1F(nx*nz);
    bb  = CPU_zaloc1F(nx*nz);

    for(ix=0; ix<nx; ix++)
    {
        x = (ix-nxb)*dx;
        deltaX = x-x0;
        for(iz=0; iz<nz; iz++) {
            z = (iz-nzb)*dz;
            deltaZ = z-z0;
            //r2 = 0.09*deltaX*deltaX + deltaZ*deltaZ;
            r2 = deltaX*deltaX + deltaZ*deltaZ;
            vp = vp0;
            vs = vs0;
            rho = rho0;
            ampGauss = dvpGauss * exp(-sigma*r2);
            if(sqrtf(r2)<=gaussRadius) {
                vp += 1.0*ampGauss;
                vs += 0.5*ampGauss;
            }
            if(zref<=z && z<=zref+lref) {
                vp  += dvp;
                vs  += dvs;
                rho += drho;
            }   

            i = iz + ix*nz;
            ll[i] = rho*(vp*vp - 2.0f*vs*vs);
            mm[i] = rho*vs*vs;
            bb[i] = 1.0f/rho;

        }
    }

    for(ix=1; ix<nx-1; ix++) {
        for(iz=1; iz<nz-1; iz++) {
            i   =  iz +     ix   *nz;
            i_d = (iz+1) + (ix+1)*nz;
            i_z = (iz+1) +  ix   *nz;
            i_x =  iz    + (ix+1)*nz;
            lambda[i]    = 0.5*(ll[i] + ll[i_x]);
            lambda2mu[i] = lambda[i]  + 2.0*0.5*(mm[i] + mm[i_x]);
            mu[i]        = 0.5*(mm[i]+mm[i_z]);
            bc[i]        = bb[i];
            bd[i]        = 0.5*(bb[i]+bb[i_d]);
        }
    }

    free(bb);
    free(ll);
    free(mm);


return;
}


float* getRicker(int nt, float dt, float fm, float t0) {
    int it;
    float t, sigma, amp;
    float *sou;
    float omega = 2*M_PI*fm;
    sou = CPU_zaloc1F(nt);
    t0 += 3.0f/fm;
    for(it=0; it<nt; it++) {
        t = it*dt;
        sigma = omega * (t-t0);
        sigma *= sigma;
        sou[it] = (1.0f-0.5*sigma) * exp(-0.25*sigma);
    }
return sou;
}

int** definePositionPSF(int nx, int nz)
{
    int **positionPSF;
    int   ix0, ix1;
    int   iz0, iz1;
    int   ix, iz, i;
    int   ix_, iz_;
    int   izIni;
    
    positionPSF = CPU_zaloc2I(nx*nz,4);

    izIni = (int) (floor(0.4*(nz-2*nzb)) + nzb);
    
    for(ix=0; ix<nx; ix++)
        for(iz=0; iz<nz; iz++)
        {
            i = iz + ix*nz;
            //*
            if     ( iz>izIni  &&  iz<nz-nzb )    iz_ = iz;
            else if( iz<izIni                )    iz_ = izIni;
            else if( iz>=nz-nzb              )    iz_ = nz-nzb-1;
            
            if     ( ix<nxb                  )    ix_ = nxb;
            else if( ix>=nx-nxb              )    ix_ = nx-nxb-1;
            else                                  ix_ = ix;
            
            ix0 = ((ix_-nxb)/jxPSF) * jxPSF + nxb;
            iz0 = ((iz_-nzb)/jzPSF) * jzPSF + nzb;
            
            //printf("\n  ix=%d  ix0=%d  iz=%d  iz0=%d", ix, ix0, iz, iz0);
            
            ix1 = ix0 + jxPSF;
            iz1 = iz0 + jzPSF;
            
            positionPSF[i][0] = ix0;
            positionPSF[i][1] = iz0;
            positionPSF[i][2] = ix1;
            positionPSF[i][3] = iz1;       //*/
            //printf("\n positionPSF[i][0]=%d  positionPSF[i][1]=%d  positionPSF[i][2]=%d  positionPSF[i][3]=%d ", \
                       positionPSF[i][0], positionPSF[i][1], positionPSF[i][2], positionPSF[i][3]); 
        }    
    

return positionPSF;
}

void createPointScatterersModel(float *m, float dv, int **positionPSF, int nx, int nz) {
    int     ix, iz, i;
    int     condition;
    float   dm;
    
    for(ix=0; ix<nx; ix++)
        for(iz=0; iz<nz; iz++) {
            i = iz + ix*nz;            
            condition = ( (ix==positionPSF[i][0])  &&  (iz==positionPSF[i][1]) );                                     
            if( condition ) dm = dv;
            else            dm = 0.0f;
            m[i] += dm;
        }
return;
}

float *extendInputModel(float *inpModel, int nx, int nz, int nxInput, int nzInput, float scalar) {

    int ix, iz;
    float *outModel = CPU_zaloc1F(nx*nz);

    for(ix=0; ix<nx; ix++)
        for(iz=0; iz<nz; iz++)
        {
            int ix0, iz0, i0, i;
            if(ix<nxb)                    ix0 = 0;
            else if((ix-nxb)<nxInput)     ix0 = (ix-nxb);
            else if((ix-nxb)>=nxInput)    ix0 = nxInput-1;

            if(iz<nzb)                    iz0 = 0;
            else if((iz-nzb)<nzInput)     iz0 = (iz-nzb);
            else if((iz-nzb)>=nzInput)    iz0 = nzInput-1;

            i  = iz  + ix  * nz;
            i0 = iz0 + ix0 * nzInput;

            outModel[i] = scalar * inpModel[i0];
        }

return outModel;
}

/*
float *extendResampInputModel_wrong(float *inpModel, int nxInput, int nzInput, float dxInput, float dzInput, float oxInput, float ozInput, \
                              int nx, int nz, float dx, float dz, float ox, float oz, int nxb, int nzb, float scalar)
{
    int ix, iz, jx, jz;
    int nx_ = nx-2*nxb;
    int nz_ = nz-2*nzb;

    int jx1 = ceilf((oxInput-ox)/dxInput);
    int jz1 = ceilf((ozInput-oz)/dzInput);
    jx1 = max(jx1, 0);
    jz1 = max(jz1, 0);


    printf("\n nxInput=%d nzInput=%d dxInput=%f dzInput=%f oxInput=%f ozInput=%f   nx=%d nz=%d ox=%f oz=%f dx=%f dz=%f   nxb=%d nzb=%d\n", \
               nxInput, nzInput, dxInput, dzInput, oxInput, ozInput, nx, nz, ox, oz, dx, dz, nxb, nzb);

    float x2 = ((nxInput-1)*dxInput) + oxInput;
    float z2 = ((nzInput-1)*dzInput) + ozInput;
    int  jx2 = floorf((x2-ox)/dx);
    int  jz2 = floorf((z2-oz)/dz);
    jx2 = min(nxInput-1, jx2);
    jz2 = min(nzInput-1, jz2);

    printf("\n  jx1=%d  jz1=%d \n", jx1, jz1);
    printf("\n  jx2=%d  jz2=%d \n", jx2, jz2);

    float *outModel = CPU_zaloc1F(nx*nz);

    for(ix=ix1; ix<=ix2; ix++)
        for(iz=iz1; iz<=iz2; iz++)
        {
            float x = ix*dx + ox;
            float z = iz*dz + oz;
            int ixa = floorf((x-oxInput)/dxInput);
            int iza = floorf((z-ozInput)/dzInput);
            int ixb = ixa+1;
            int izb = iza+1;

            float xa = ixa * dxInput + oxInput;
            float xb = ixb * dxInput + oxInput;
            float za = iza * dzInput + ozInput;
            float zb = izb * dzInput + ozInput;
            float pxa = (xb - x)/dx;
            float pxb = (x - xa)/dx;
            float pza = (xb - z)/dz;
            float pzb = (z - xa)/dz;
            float ptot_x = pxa + pxb;
            float ptot_z = pza + pzb;
            pxa /= ptot_x;
            pxb /= ptot_x;
            pza /= ptot_z;
            pzb /= ptot_z;

            float ptot = pxa*pza + pxb*pza + pxa*pzb + pxb*pzb;

            int i = (iz+nzb) + (ix+nxb) * nz;

            int i_xa_za = iza + ixa * nzInput;
            int i_xb_za = iza + ixb * nzInput;
            int i_xa_zb = izb + ixa * nzInput;
            int i_xb_zb = izb + ixb * nzInput;

            outModel[i] = scalar * ( pxa*pza*inpModel[i_xa_za] + pxb*pza*inpModel[i_xb_za] + pxa*pzb*inpModel[i_xa_zb] + pxb*pzb*inpModel[i_xb_zb] ) / ptot;
        }
    
    ix1 += nxb;
    ix2 += nxb;
    iz1 += nzb;
    iz2 += nzb;

    // Left side
    for(ix=0; ix<ix1; ix++)
        for(iz=iz1; iz<iz2; iz++)
        {
            int i  = iz + ix  * nz;
            int i1 = iz + ix1 * nz;
            outModel[i] = outModel[i1];
        }
    // Right side
    for(ix=ix2+1; ix<nx; ix++)
        for(iz=iz1; iz<iz2; iz++)
        {
            int i  = iz + ix  * nz;
            int i2 = iz + ix2 * nz;
            outModel[i] = outModel[i2];
        }
    
    for(ix=0; ix<nx; ix++)
    {
        // Upper side
        for(iz=0; iz<iz1; iz++)
        {
            int i  = iz  + ix * nz;
            int i1 = iz1 + ix * nz;
            outModel[i] = outModel[i1];
        }
        // Lower side
        for(iz=iz2+1; iz<nz; iz++)
        {
            int i  = iz  + ix * nz;
            int i2 = iz2 + ix * nz;
            outModel[i] = outModel[i2];
        }
    }

return outModel;
}
//*/

float *extendResampInputModel(float *inpModel, int nxInput, int nzInput, float dxInput, float dzInput, float oxInput, float ozInput, \
                              int nx, int nz, float dx, float dz, float ox, float oz, int nxb, int nzb, float scalar)
{
    int ix, iz;
    int nx_ = nx-2*nxb;
    int nz_ = nz-2*nzb;

    int ix1 = ceilf((oxInput-ox)/dx);
    int iz1 = ceilf((ozInput-oz)/dz);
    ix1 = max(ix1, 0);
    iz1 = max(iz1, 0);


    // printf("\n nxInput=%d nzInput=%d dxInput=%f dzInput=%f oxInput=%f ozInput=%f   nx=%d nz=%d ox=%f oz=%f dx=%f dz=%f   nxb=%d nzb=%d\n", \
               nxInput, nzInput, dxInput, dzInput, oxInput, ozInput, nx, nz, ox, oz, dx, dz, nxb, nzb);

    float x2 = ((nxInput-1)*dxInput) + oxInput;
    float z2 = ((nzInput-1)*dzInput) + ozInput;
    int  ix2 = floorf((x2-ox)/dx);
    int  iz2 = floorf((z2-oz)/dz);
    ix2 = min(nx_-1, ix2);
    iz2 = min(nz_-1, iz2);

    float *outModel = CPU_zaloc1F(nx*nz);

    for(ix=ix1; ix<=ix2; ix++)
        for(iz=iz1; iz<=iz2; iz++)
        {
            float x = ix*dx + ox;
            float z = iz*dz + oz;
            int ixa = floorf((x-oxInput)/dxInput);
            int iza = floorf((z-ozInput)/dzInput);
            int ixb = ixa+1;
            int izb = iza+1;

            float xa = ixa * dxInput + oxInput;
            float xb = ixb * dxInput + oxInput;
            float za = iza * dzInput + ozInput;
            float zb = izb * dzInput + ozInput;
            float pxa = (xb - x)/dxInput;
            float pxb = (x - xa)/dxInput;
            float pza = (zb - z)/dzInput;
            float pzb = (z - za)/dzInput;
            float ptot_x = pxa + pxb;
            float ptot_z = pza + pzb;
            pxa /= ptot_x;
            pxb /= ptot_x;
            pza /= ptot_z;
            pzb /= ptot_z;

            float ptot = pxa*pza + pxb*pza + pxa*pzb + pxb*pzb;

            int i = (iz+nzb) + (ix+nxb) * nz;

            int i_xa_za = iza + ixa * nzInput;
            int i_xb_za = iza + ixb * nzInput;
            int i_xa_zb = izb + ixa * nzInput;
            int i_xb_zb = izb + ixb * nzInput;

            outModel[i] = scalar * ( pxa*pza*inpModel[i_xa_za] + pxb*pza*inpModel[i_xb_za] + pxa*pzb*inpModel[i_xa_zb] + pxb*pzb*inpModel[i_xb_zb] ) / ptot;
        }
    
    ix1 += nxb;
    ix2 += nxb;
    iz1 += nzb;
    iz2 += nzb;

    // Left side
    for(ix=0; ix<ix1; ix++)
        for(iz=iz1; iz<=iz2; iz++)
        {
            int i  = iz + ix  * nz;
            int i1 = iz + ix1 * nz;
            outModel[i] = outModel[i1];
        }
    // Right side
    for(ix=ix2+1; ix<nx; ix++)
        for(iz=iz1; iz<=iz2; iz++)
        {
            int i  = iz + ix  * nz;
            int i2 = iz + ix2 * nz;
            outModel[i] = outModel[i2];
        }
    
    for(ix=0; ix<nx; ix++)
    {
        // Upper side
        for(iz=0; iz<iz1; iz++)
        {
            int i  = iz  + ix * nz;
            int i1 = iz1 + ix * nz;
            outModel[i] = outModel[i1];
        }
        // Lower side
        for(iz=iz2+1; iz<nz; iz++)
        {
            int i  = iz  + ix * nz;
            int i2 = iz2 + ix * nz;
            outModel[i] = outModel[i2];
        }
    }

return outModel;
}

void fillBoundariesInUpdatedModel(float *model, int nxb, int nx, int nz, float x1, float x2, float dx)
{
    int ix, iz;
    int i1, i2, i;

    int ix1 = x1/dx + nxb;
    int ix2 = x2/dx + nxb;

    for(ix=0; ix<ix1; ix++)
        for(iz=0; iz<nz; iz++)
        {
            i  = iz + ix  * nz;
            i1 = iz + ix1 * nz;
            model[i] = model[i1];
        }
    for(ix=ix2+1; ix<nx; ix++)
        for(iz=0; iz<nz; iz++)
        {
            i  = iz + ix  * nz;
            i2 = iz + ix2 * nz;
            model[i] = model[i2];
        }

return;
}

void padBoundaries(float *model, int nx, int nz, int nxb, int nzb) {

    int ix, iz;

    for(ix=0; ix<nxb; ix++)
        for(iz=nzb; iz<nz-nzb; iz++)
        {
            int i0, i;
            i  = iz + ix  * nz;
            i0 = iz + nxb * nz;
            model[i] = model[i0];
        }
    for(ix=nx-nxb; ix<nx; ix++)
        for(iz=nzb; iz<nz-nzb; iz++)
        {
            int i0, i;
            i  = iz + ix  * nz;
            i0 = iz + (nx-nxb-1) * nz;
            model[i] = model[i0];
        }
    for(ix=0; ix<nx; ix++)
        for(iz=0; iz<nzb; iz++)
        {
            int i0, i;
            i  = iz  + ix * nz;
            i0 = nzb + ix * nz;
            model[i] = model[i0];
        }
    for(ix=0; ix<nx; ix++)
        for(iz=nz-nzb; iz<nz; iz++)
        {
            int i0, i;
            i  = iz       + ix * nz;
            i0 = nz-nzb-1 + ix * nz;
            model[i] = model[i0];
        }

return;
}

void fillGradientEdges(float *gradient, int nx, int nz, float dx, float x1, float x2) {

    int ix, iz;
    int ix1, ix2;

    ix1 = ((int) ceilf(x1/dx)) + nxb;
    ix2 = ((int) ceilf(x2/dx)) + nxb;

    for(ix=0; ix<ix1; ix++) 
    {
        for(iz=0; iz<nz; iz++)
        {
            int i  = iz  + ix  * nz;
            int i1 = iz  + ix1 * nz;            
            gradient[i] = gradient[i1];
        }
    }
    for(ix=ix2+1; ix<nx; ix++) 
    {
        for(iz=0; iz<nz; iz++)
        {
            int i  = iz  + ix  * nz;
            int i2 = iz  + ix2 * nz;            
            gradient[i] = gradient[i2];
        }
    }

return;
}