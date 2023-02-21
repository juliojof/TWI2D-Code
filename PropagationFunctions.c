void eimagingCondition_lx(float *eimage, float *ws, float *wr, int nx, int nz, int lx0) {
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;
    int ix, iz, lx, nlx=2*lx0+1;
    //for(ix=lx0;ix<nx-lx0;ix++)
        // for(iz=0;iz<nz;iz++) {
    #pragma omp parallel for schedule(dynamic,2) num_threads(4) private(iz, ix, lx)
    for(ix=ixMin;ix<ixMax;ix++)
        for(iz=izMin;iz<izMax;iz++) {
            for(lx=-lx0; lx<=lx0; lx++)
            {
                int is = id2(iz,ix-lx);
                int ir = id2(iz,ix+lx);
                int i  = id3_eic_lx(lx,iz,ix);
                eimage[i] += ws[is] * wr[ir];
            }
        }
return;
}
void eimagingCondition_lx_timeDer_1stPart(float *eimage, float *ws1, float *ws2, float *wr2, int nx, int nz, int lx0) {
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;
    int ix, iz, lx, nlx=2*lx0+1;
    //for(ix=lx0;ix<nx-lx0;ix++)
        // for(iz=0;iz<nz;iz++) {
    for(ix=ixMin;ix<ixMax;ix++)
        for(iz=izMin;iz<izMax;iz++) {
            for(lx=-lx0; lx<=lx0; lx++)
            {
                int is = id2(iz,ix-lx);
                int ir = id2(iz,ix+lx);
                int i  = id3_eic_lx(lx,iz,ix);
                eimage[i] += (2.0*ws2[is]-ws1[is]) * wr2[ir];
            }
        }
return;
}
void eimagingCondition_lx_timeDer_2ndPart(float *eimage, float *ws2, float *wr1, int nx, int nz, int lx0) {
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;
    int ix, iz, lx, nlx=2*lx0+1;
    //for(ix=lx0;ix<nx-lx0;ix++)
        // for(iz=0;iz<nz;iz++) {
    for(ix=ixMin;ix<ixMax;ix++)
        for(iz=izMin;iz<izMax;iz++) {
            for(lx=-lx0; lx<=lx0; lx++)
            {
                int is = id2(iz,ix-lx);
                int ir = id2(iz,ix+lx);
                int i  = id3_eic_lx(lx,iz,ix);
                eimage[i] += ws2[is] * wr1[ir];
            }
        }
return;
}
void eimagingCondition_lx_theta(float *eimage, float *eimage_shot, float *eimage_theta, float *eimage_shot_theta, \
                                float *ws, float *wr, int nx, float dx, \
                                int nz, float dz, int lx0, int ntheta, float dtheta, float otheta) {

    int ix, iz, lx, itheta, ltheta0=ntheta/2;
    int nlx=2*lx0+1;

    int ixMin = nxb;
    int ixMax = nx-nxb;
    //int ixMin = nx/2;
    //int ixMax = ixMin+1;
    int izMin = nzb;
    int izMax = nz-nzb;

    // Extended Imaging Condition to get Lag Gather
    for(ix=ixMin;ix<ixMax;ix++)
        for(iz=izMin;iz<izMax;iz++) {
            for(lx=-lx0; lx<=lx0; lx++)
            {
                int is = id2(iz,ix-lx);
                int ir = id2(iz,ix+lx);
                int i  = id3_eic_lx(lx,iz,ix);
                eimage[i] += ws[is] * wr[ir];
                if(eimage_shot!=NULL) eimage_shot[i] += ws[is] * wr[ir];
            }
        }
    // Slant Stack to get Angle Gather
    ixMin = nx/2;
    ixMax = ixMin+1;
    int   iz_0, iz_1, i_0, i_1;
    float len, deltaZ, z, z_0, z_1;
    float pz_0, pz_1, ptot;
    for(ix=ixMin;ix<ixMax;ix++)
        for(iz=izMin;iz<izMax;iz++) 
            for(itheta=-ltheta0; itheta<=ltheta0; itheta++) {
                int   i_th  = id3_eic_th(itheta,iz,ix);
                float theta = (itheta*dtheta) * M_PI/180.0f;
                for(lx=-lx0; lx<=lx0; lx++) {
                    len    = lx*dx;
                    deltaZ = len * tanf(theta);
                    z      = dz*iz + deltaZ;
                    iz_0   = floorf(z/dz);
                    if(iz_0<nzb)          continue;
                    if(iz_0>=nz-nzb-1)    continue;
                    iz_1   = iz_0 + 1;
                    //printf("\n nz=%d  iz_0=%d  iz_1=%d  len=%f  theta=%f  tan=%f  z=%f  deltaZ=%f otheta=%f\n", \
                    nz, iz_0, iz_1, len, theta*180.0/M_PI, tanf(theta), z, deltaZ, otheta);
                    z_0    = iz_0*dz;
                    z_1    = iz_1*dz;
                    pz_0   = z_1 - z;
                    pz_1   = z - z_0;
                    ptot   = pz_0 + pz_1;
                    i_0    = id3_eic_lx(lx,iz_0,ix);
                    i_1    = id3_eic_lx(lx,iz_1,ix);
                    if(eimage_theta      !=NULL) eimage_theta[i_th]      += (pz_0*eimage_shot[i_0] + pz_1*eimage_shot[i_1])/ptot;
                    if(eimage_shot_theta !=NULL) eimage_shot_theta[i_th] += (pz_0*eimage_shot[i_0] + pz_1*eimage_shot[i_1])/ptot;
                }
            }
return;
}
void lx2theta_forward(int add, float *eimage, float *eimage_theta, \
                      int nx, float dx, int nz, float dz, \
                      int lx0, int ntheta, float dtheta, float otheta) {

    int ix, iz, lx, itheta, ltheta0=ntheta/2;
    int nlx=2*lx0+1;

    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;

    if(add == 0) {
        for(ix=ixMin;ix<ixMax;ix++)
            for(iz=izMin;iz<izMax;iz++) 
                for(itheta=-ltheta0; itheta<=ltheta0; itheta++) {
                    int   i_th  = id3_eic_th(itheta,iz,ix);
                    eimage_theta[i_th] = 0.0f;
                }
    }


    // Slant Stack to get Angle Gather
    int   iz_0, iz_1, i_0, i_1;
    float len, deltaZ, z, z_0, z_1;
    float pz_0, pz_1, ptot;
    for(ix=ixMin;ix<ixMax;ix++)
        for(iz=izMin;iz<izMax;iz++) 
            for(itheta=-ltheta0; itheta<=ltheta0; itheta++) {
                int   i_th  = id3_eic_th(itheta,iz,ix);
                float theta = (itheta*dtheta) * M_PI/180.0f;
                for(lx=-lx0; lx<=lx0; lx++) {
                    len    = lx*dx;
                    deltaZ = len * tanf(theta);
                    z      = dz*iz + deltaZ;
                    iz_0   = floorf(z/dz);
                    if(iz_0<=nzb)         continue;
                    if(iz_0>=nz-nzb-1)    continue;
                    iz_1   = iz_0 + 1;
                    z_0    = iz_0*dz;
                    z_1    = iz_1*dz;
                    pz_0   = z_1 - z;
                    pz_1   = z - z_0;
                    ptot   = pz_0 + pz_1;
                    i_0    = id3_eic_lx(lx,iz_0,ix);
                    i_1    = id3_eic_lx(lx,iz_1,ix);
                    eimage_theta[i_th] += (pz_0*eimage[i_0] + pz_1*eimage[i_1])/ptot;
                }
            }
return;
}
void lx2theta_adjoint(int add, float *eimage, float *eimage_theta, \
                      int nx, float dx, int nz, float dz, \
                      int lx0, int ntheta, float dtheta, float otheta) {

    int ix, iz, lx, itheta, ltheta0=ntheta/2;
    int nlx=2*lx0+1;

    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;

    int   iz_0, iz_1, i_0, i_1;
    float len, deltaZ, z, z_0, z_1;
    float pz_0, pz_1, ptot;

    if(add == 0) {
        for(ix=ixMin;ix<ixMax;ix++)
            for(iz=izMin;iz<izMax;iz++) 
                for(lx=-lx0; lx<=lx0; lx++) {
                    i_0 = id3_eic_lx(lx,iz,ix);
                    eimage[i_0] = 0.0f;
                }
    }

    for(ix=ixMin;ix<ixMax;ix++)
        for(iz=izMin;iz<izMax;iz++) 
            for(itheta=-ltheta0; itheta<=ltheta0; itheta++) {
                int   i_th  = id3_eic_th(itheta,iz,ix);
                float theta = (itheta*dtheta) * M_PI/180.0f;
                for(lx=-lx0; lx<=lx0; lx++) {
                    len    = lx*dx;
                    deltaZ = len * tanf(theta);
                    z      = dz*iz + deltaZ;
                    iz_0   = floorf(z/dz);
                    if(iz_0<=nzb)         continue;
                    if(iz_0>=nz-nzb-1)    continue;
                    iz_1   = iz_0 + 1;
                    z_0    = iz_0*dz;
                    z_1    = iz_1*dz;
                    pz_0   = z_1 - z;
                    pz_1   = z - z_0;
                    ptot   = pz_0 + pz_1;
                    i_0    = id3_eic_lx(lx,iz_0,ix);
                    i_1    = id3_eic_lx(lx,iz_1,ix);
                    // eimage[i_0] += pz_0*eimage_theta[i_th]/(ptot*ntheta);
                    // eimage[i_1] += pz_1*eimage_theta[i_th]/(ptot*ntheta);
                    eimage[i_0] += pz_0*eimage_theta[i_th]/(ptot*ntheta);
                    eimage[i_1] += pz_1*eimage_theta[i_th]/(ptot*ntheta);
                }
            }
return;
}
void lx2theta_adjoint_singleTheta(int add, float *eimage, float *eimage_theta, \
                      int nx, float dx, int nz, float dz, \
                      int lx0, int ntheta, float dtheta, float otheta, int itheta) {

    int ix, iz, lx, ltheta0=ntheta/2;
    int nlx=2*lx0+1;

    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;

    int   iz_0, iz_1, i_0, i_1;
    float len, deltaZ, z, z_0, z_1;
    float pz_0, pz_1, ptot;

    if(add == 0) {
        for(ix=ixMin;ix<ixMax;ix++)
            for(iz=izMin;iz<izMax;iz++) 
                for(lx=-lx0; lx<=lx0; lx++) {
                    i_0 = id3_eic_lx(lx,iz,ix);
                    eimage[i_0] = 0.0f;
                }
    }

    for(ix=ixMin;ix<ixMax;ix++)
        for(iz=izMin;iz<izMax;iz++)
        {            
            int   i_th  = id3_eic_th(itheta,iz,ix);
            float theta = (itheta*dtheta) * M_PI/180.0f;
            for(lx=-lx0; lx<=lx0; lx++) {
                len    = lx*dx;
                deltaZ = len * tanf(theta);
                z      = dz*iz + deltaZ;
                iz_0   = floorf(z/dz);
                if(iz_0<=nzb)         continue;
                if(iz_0>=nz-nzb-1)    continue;
                iz_1   = iz_0 + 1;
                z_0    = iz_0*dz;
                z_1    = iz_1*dz;
                pz_0   = z_1 - z;
                pz_1   = z - z_0;
                ptot   = pz_0 + pz_1;
                i_0    = id3_eic_lx(lx,iz_0,ix);
                i_1    = id3_eic_lx(lx,iz_1,ix);
                eimage[i_0] += pz_0*eimage_theta[i_th]/ptot;
                eimage[i_1] += pz_1*eimage_theta[i_th]/ptot;                      
            }
        }
return;
}
void scatter_Acoustic_AmpAngle(float *ref, float *sct, float coef_2, float *inc_2, float coef_1, float *inc_1, int nx, int nz, float dt, int lx0) {
    int ix, iz, lx;
    int ixMin = nxb;
    int ixMax = nx-nxb;
    // int izMin = nzb;
    // int izMax = nz-nzb;
    int izMin = nzb + 145;
    int izMax = nzb + 155;
    int nlx=2*lx0+1;
    for(ix=ixMin;ix<ixMax;ix++)
        for(iz=izMin;iz<izMax;iz++) {
            for(lx=-lx0; lx<=lx0; lx++) 
            {
                int i_inc = id2(iz,ix-lx);
                int i_sct = id2(iz,ix);
                int i     = id3_eic_lx(lx,iz,ix);
                sct[i_sct] -= ref[i] * (coef_2*inc_2[i_inc] + coef_1*inc_1[i_inc]);
            }
        }
}
void scatter_Acoustic_lx_opt(float *ref, float *sct, float *inc, int nx, int nz, float dt, int lx0) {
    int ix, iz, lx;
    int i, i_sct, i_inc;
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;
    int nlx=2*lx0+1;
#pragma omp parallel for schedule(dynamic,2) num_threads(4) private(iz, ix, i, i_inc, i_sct, lx)
    for(ix=ixMin; ix<ixMax; ix++)
        for(iz=izMin; iz<izMax; iz++) {
            for(lx=-lx0; lx<=lx0; lx++) 
            {
                i_inc   = id2(iz,ix-lx);
                i_sct   = id2(iz,ix+lx);
                i       = id3_eic_lx(lx,iz,ix);
                sct[i_sct] += ref[i] * inc[i_inc];
            }
        }
}
void scatter_Acoustic_lx(float *ref, float *sct, float coef_2, float *inc_2, float coef_1, float *inc_1, int nx, int nz, float dt, int lx0) {
    int ix, iz, lx;
    int i, i_sct, i_inc;
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;
    int nlx=2*lx0+1;
#pragma omp parallel for schedule(dynamic,2) num_threads(4) private(iz, ix, i, i_inc, i_sct, lx)
    for(ix=ixMin; ix<ixMax; ix++)
        for(iz=izMin; iz<izMax; iz++) {
            for(lx=-lx0; lx<=lx0; lx++) 
            {
                i_inc   = id2(iz,ix-lx);
                i_sct   = id2(iz,ix+lx);
                i       = id3_eic_lx(lx,iz,ix);
                sct[i_sct] -= ref[i] * (coef_2*inc_2[i_inc] + coef_1*inc_1[i_inc]);
            }
        }
}
void scatter_Acoustic(float *ref, float *sct, float coef_2, float *inc_2, float coef_1, float *inc_1, int nx, int nz, float dt) {
    int ix, iz;
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;
    for(ix=ixMin;ix<ixMax;ix++)
        for(iz=izMin;iz<izMax;iz++) {
            int i = id2(iz,ix);
            sct[i] -= ref[i] * (coef_2*inc_2[i] + coef_1*inc_1[i]);
        }
}
void wavefieldTimeSecondDerivative(float *secDer, float coef_2, float *inc_2, float coef_1, float *inc_1, int nx, int nz, float dt) {
    int ix, iz;
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;
    for(ix=ixMin;ix<ixMax;ix++)
        for(iz=izMin;iz<izMax;iz++) {
            int i = id2(iz,ix);
            secDer[i] += (coef_2*inc_2[i] + coef_1*inc_1[i])/(dt*dt);
        }
}
void imagingCondition(float *image, float *ws, float *wr, int nx, int nz) {
    int ix, iz, i;
    int ixMin = nxb;
    int ixMax = nx-nxb;
    int izMin = nzb;
    int izMax = nz-nzb;
#pragma omp parallel for schedule(dynamic,2) num_threads(4) private(iz, ix, i)
    for(ix=ixMin;ix<ixMax;ix++)
        for(iz=izMin;iz<izMax;iz++) {
            i = id2(iz,ix);
            image[i] += ws[i] * wr[i];
        }
return;
}
void injectDataTaper(Shot *shot, float *wav, int nx, int nz, float dx, float dz, int it) {
    int ir;
    int nt = shot->nt;
    int nrec = shot->nrec;
    float xedge;
    float filter;
    int nx_ = nx-2*nxb;
    int len = (int) (0.1*nx_);
    for(ir=0; ir<nrec; ir++) {
        int ixr = shot->ixr[ir];
        int izr = shot->izr[ir];
        int ixr_ = ixr-nxb;
        if(ixr_<len  &&  ixr_>=0)                  xedge = dx*(ixr_/len);
        else if((nx_-1)-ixr_<len  &&  ixr_<nx_)    xedge = dx*((nx_-1)-ixr_)/len;
        else if(ixr_<0)                            xedge = 0.0f;
        else if(ixr_>=nx_)                         xedge = 0.0f;
        else                                       xedge = 1.0f;
        filter = sin(0.5*M_PI*xedge);
        filter *= filter;
        injectSourceTaper_Acoustic(wav, (shot->seismogram)+ir*nt, filter, it, ixr, izr, nx, nz, dx, dz);
    }
return;
}
void injectData(Shot *shot, float *wav, int nx, int nz, float dx, float dz, int it) {
    int ir;
    int nt = shot->nt;
    int nrec = shot->nrec;
    for(ir=0; ir<nrec; ir++) {
        int ixr = shot->ixr[ir];
        int izr = shot->izr[ir];
        injectSource_Acoustic(wav, (shot->seismogram)+ir*nt, it, ixr, izr, nx, nz, dx, dz);
    }
return;
}
void recordData(Shot *shot, float *wav2, int nx, int nz, int it) {
    int ir;
    int nt = shot->nt;
    int nrec = shot->nrec;
    for(ir=0; ir<nrec; ir++) {
        int ixr = shot->ixr[ir];
        int izr = shot->izr[ir];
        int i = id2(izr,ixr);
        shot->seismogram[ir*nt+it] = wav2[i];
    }
return;
}
void injectSource_Acoustic(float *u, float *source, int it, int ix0, int iz0, int nx, int nz, float dx, float dz)  {
    int i = id2(iz0,ix0);
    u[i] -= source[it]/(dx*dz);
return;
}
void injectSourceTaper_Acoustic(float *u, float *source, float taper, int it, int ix0, int iz0, int nx, int nz, float dx, float dz)  {
    int i = id2(iz0,ix0);
    u[i] -= taper * source[it]/(dx*dz);
return;
}
void propagate_Acoustic(float *wav1, float *wav2, float *lap, float *vp, int nx, int nz) {
    laplacian_10thOrder_2D_Acoustic_Iso(lap, wav2, nx, nz);
    timeStep_Acoustic(wav1, wav2, lap, vp, nx, nz);
return;
}
void propagate_Acoustic_abs(float *wav1, float *wav2, float *lap, float *vp, float *sigma, int nx, int nz) {
    laplacian_10thOrder_2D_Acoustic_Iso(lap, wav2, nx, nz);
    timeStep_Acoustic_abs(wav1, wav2, lap, vp, sigma, nx, nz);
return;
}
void propagate_Acoustic_abs2(float *wav1, float *wav2, float *lap, float *vp, float *sigma, float *sigmaInv, int nx, int nz) {
    laplacian_10thOrder_2D_Acoustic_Iso(lap, wav2, nx, nz);
    timeStep_Acoustic_abs2(wav1, wav2, lap, vp, sigma, sigmaInv, nx, nz);
return;
}
void timeStep_Acoustic(float *u1, const float *u2, const float *lap, const float *vp2dt2dx2, int nx, int nz)  {
    int ix, iz;
    for(ix=0; ix<nx; ix++) {
        for(iz=0; iz<nz; iz++) {
            int i = id2(iz,ix);
            u1[i] = 2.0f*u2[i] - u1[i] + vp2dt2dx2[i] * lap[i];
        }
    }
}
void timeStep_Acoustic_abs(float *u1, const float *u2, const float *lap, const float *vp2dt2dx2, const float *sigma, int nx, int nz)  {
    int ix, iz, i;
#pragma omp parallel for schedule(dynamic,2) num_threads(4) private(iz,ix,i)
    for(ix=0; ix<nx; ix++) {
        for(iz=0; iz<nz; iz++) {
            i = id2(iz,ix);
            u1[i] = (2.0f*u2[i] - (1.0f-sigma[i])*u1[i] + vp2dt2dx2[i]*lap[i]) / (1.0f+sigma[i]);
        }
    }
}
void timeStep_Acoustic_abs2(float *u1, const float *u2, const float *lap, const float *vp2dt2dx2, \
                           const float *sigma, const float *sigmaInv, int nx, int nz)  {
    int ix, iz, i;
#pragma omp parallel for schedule(dynamic,2) num_threads(4) private(iz,ix,i)
    for(ix=0; ix<nx; ix++) {
        // #pragma clang loop vectorize_width(4) interleave_count(4)
        // #pragma unroll
        for(iz=0; iz<nz; iz++) {
            i = id2(iz,ix);
            u1[i] = (2.0f*u2[i] - (1.0f-sigma[i])*u1[i] + vp2dt2dx2[i]*lap[i]) * sigmaInv[i];
        }
    }
}
// Compute the 2D laplacian to 10th order
// Division by dx^2 is embedded in the coefficients C0 through C5.
void laplacian_10thOrder_2D_Acoustic_Iso(float *lap, const float *u2, int nx, int nz) {
    int ix, iz, i;
#pragma omp parallel for schedule(dynamic,2) num_threads(4) private(iz,ix,i)
    for(ix=OL; ix<nx-OL; ix++) {
        // #pragma clang loop vectorize(enable) interleave(enable)
        // #pragma clang loop vectorize_width(4) interleave_count(4)
        // #pragma unroll
        for(iz=OL; iz<nz-OL; iz++) {
            i = id2(iz,ix);
            lap[i] =      C5_10 * (u2[id2(iz+5,ix  )]+u2[id2(iz-5,ix  )]
                                  +u2[id2(iz  ,ix+5)]+u2[id2(iz  ,ix-5)])
                   +      C4_10 * (u2[id2(iz+4,ix  )]+u2[id2(iz-4,ix  )]
                                  +u2[id2(iz  ,ix+4)]+u2[id2(iz  ,ix-4)])
                   +      C3_10 * (u2[id2(iz+3,ix  )]+u2[id2(iz-3,ix  )]
                                  +u2[id2(iz  ,ix+3)]+u2[id2(iz  ,ix-3)])
                   +      C2_10 * (u2[id2(iz+2,ix  )]+u2[id2(iz-2,ix  )]
                                  +u2[id2(iz  ,ix+2)]+u2[id2(iz  ,ix-2)])
                   +      C1_10 * (u2[id2(iz+1,ix  )]+u2[id2(iz-1,ix  )]
                                  +u2[id2(iz  ,ix+1)]+u2[id2(iz  ,ix-1)])
                   - 4.0f*C0_10 *  u2[id2(iz  ,ix  )];
        }
    }
}

void getSigmas(float *sigmaX, float *sigmaZ, int nb, int operLen, int nx, int nz, float dx, float dz)
{

    int     ix, iz;
    float   x, z, L, sgm;
    float   reflec=0.0000001, factor=3000.0;
    

    L = (nb-operLen) * dx;
    for(ix=0; ix<nb; ix++)
    {        
	    x = (nb-ix)*dx;
        sgm = log(1.0f/reflec) * factor/L * (x/L)*(x/L);
        sigmaX[ix]      = sgm;
        sigmaX[nx-ix-1] = sgm;
	}
	
	L = (nb-operLen) * dz;
    for(iz=0; iz<nb; iz++)
    {        
	    z = (nb-iz)*dz;
        sgm = log(1.0f/reflec) * factor/L * (z/L)*(z/L);

        sigmaZ[iz]      = sgm;
        sigmaZ[nz-iz-1] = sgm;
	}

return;
}
void getVelSigma(float **sigma, float **sigmaInv, float *sigX, float *sigZ, int nx, int nz, float dx, float dt)  {
    int ix, iz;
    *sigma    = CPU_zaloc1F(nx*nz);
    *sigmaInv = CPU_zaloc1F(nx*nz);
    for(ix=0; ix<nx; ix++) {
        for(iz=0; iz<nz; iz++) {
            int i = id2(iz,ix);
            if(nzb<=iz && iz<=nz-nzb)
                (*sigma)[i]    = 0.5f * (sigX[ix] + sigZ[iz]) * dt;
            else
                (*sigma)[i] = 0.5f * (0.6*sigX[ix] + sigZ[iz]) * dt;
            
            (*sigmaInv)[i] = 1.0f / (1.0f + (*sigma)[i]);

        }
    }
}