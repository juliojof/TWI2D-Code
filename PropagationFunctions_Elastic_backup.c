void propagate_Elastic(float *vx, float *vz, float *Sxx, float *Szz, float *Sxz,
                       float *lap, float *lambda, float *mu, float *b, int nx, int nz, float dx, float dt) {
    
    // Update particle velocity
    laplacian_10thOrder_Vel_Elastic_Iso(lap, b, Sxx, Sxz, nx, nz);
    timeStep_Elastic(vx, lap, nx, nz, dt, dx);
    laplacian_10thOrder_Vel_Elastic_Iso(lap, b, Sxz, Szz, nx, nz);
    timeStep_Elastic(vz, lap, nx, nz, dt, dx);

    // Update stress tensor
    laplacian_10thOrder_Sxx_Elastic_Iso(lap, lambda, mu, vx, vz, nx, nz);
    timeStep_Elastic(Sxx, lap, nx, nz, dt, dx);
    laplacian_10thOrder_Szz_Elastic_Iso(lap, lambda, mu, vx, vz, nx, nz);
    timeStep_Elastic(Szz, lap, nx, nz, dt, dx);
    laplacian_10thOrder_Sxz_Elastic_Iso(lap,         mu, vx, vz, nx, nz);
    timeStep_Elastic(Sxz, lap, nx, nz, dt, dx);

return;
}

void laplacian_10thOrder_Vel_Elastic_Iso(float *lap, float *b, float *Sdx, float *Sdz, int nx, int nz) {
    for(int ix=OL; ix<nx-OL; ix++) {
        for(int iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            lap[i] = b[i]*(B5_10 * (Sdz[id2(iz+5,ix  )]-Sdz[id2(iz-5,ix  )]
                                   +Sdx[id2(iz  ,ix+5)]-Sdx[id2(iz  ,ix-5)])
                   +       B4_10 * (Sdz[id2(iz+4,ix  )]-Sdz[id2(iz-4,ix  )]
                                   +Sdx[id2(iz  ,ix+4)]-Sdx[id2(iz  ,ix-4)])
                   +       B3_10 * (Sdz[id2(iz+3,ix  )]-Sdz[id2(iz-3,ix  )]
                                   +Sdx[id2(iz  ,ix+3)]-Sdx[id2(iz  ,ix-3)])
                   +       B2_10 * (Sdz[id2(iz+2,ix  )]-Sdz[id2(iz-2,ix  )]
                                   +Sdx[id2(iz  ,ix+2)]-Sdx[id2(iz  ,ix-2)])
                   +       B1_10 * (Sdz[id2(iz+1,ix  )]-Sdz[id2(iz-1,ix  )]
                                   +Sdx[id2(iz  ,ix+1)]-Sdx[id2(iz  ,ix-1)]));
        }
    }
}
void laplacian_10thOrder_Sxx_Elastic_Iso(float *lap, float *lambda, float *mu, float *vx, float *vz, int nx, int nz) {
    for(int ix=OL; ix<nx-OL; ix++) {
        for(int iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            float m_vx = lambda[i] + mu[i];
            float m_vz = lambda[i];
            lap[i] = B5_10 * ( m_vz*(vz[id2(iz+5,ix  )]-vz[id2(iz-5,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+5)]-vx[id2(iz  ,ix-5)]))
                   + B4_10 * ( m_vz*(vz[id2(iz+4,ix  )]-vz[id2(iz-4,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+4)]-vx[id2(iz  ,ix-4)]))
                   + B3_10 * ( m_vz*(vz[id2(iz+3,ix  )]-vz[id2(iz-3,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+3)]-vx[id2(iz  ,ix-3)]))
                   + B2_10 * ( m_vz*(vz[id2(iz+2,ix  )]-vz[id2(iz-2,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+2)]-vx[id2(iz  ,ix-2)]))
                   + B1_10 * ( m_vz*(vz[id2(iz+1,ix  )]-vz[id2(iz-1,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+1)]-vx[id2(iz  ,ix-1)]));
        }
    }
}
void laplacian_10thOrder_Szz_Elastic_Iso(float *lap, float *lambda, float *mu, float *vx, float *vz, int nx, int nz) {
    for(int ix=OL; ix<nx-OL; ix++) {
        for(int iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            float m_vz = lambda[i] + mu[i];
            float m_vx = lambda[i];
            lap[i] = B5_10 * ( m_vz*(vz[id2(iz+5,ix  )]-vz[id2(iz-5,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+5)]-vx[id2(iz  ,ix-5)]))
                   + B4_10 * ( m_vz*(vz[id2(iz+4,ix  )]-vz[id2(iz-4,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+4)]-vx[id2(iz  ,ix-4)]))
                   + B3_10 * ( m_vz*(vz[id2(iz+3,ix  )]-vz[id2(iz-3,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+3)]-vx[id2(iz  ,ix-3)]))
                   + B2_10 * ( m_vz*(vz[id2(iz+2,ix  )]-vz[id2(iz-2,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+2)]-vx[id2(iz  ,ix-2)]))
                   + B1_10 * ( m_vz*(vz[id2(iz+1,ix  )]-vz[id2(iz-1,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+1)]-vx[id2(iz  ,ix-1)]));
        }
    }
}
void laplacian_10thOrder_Sxz_Elastic_Iso(float *lap, float *mu, float *vx, float *vz, int nx, int nz) {
    for(int ix=OL; ix<nx-OL; ix++) {
        for(int iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            float m_v = mu[i];
            lap[i] = m_v * (B5_10 * ( vx[id2(iz+5,ix  )]-vx[id2(iz-5,ix  )]
                                    + vz[id2(iz  ,ix+5)]-vz[id2(iz  ,ix-5)])
                          + B4_10 * ( vx[id2(iz+4,ix  )]-vx[id2(iz-4,ix  )]
                                    + vz[id2(iz  ,ix+4)]-vz[id2(iz  ,ix-4)])
                          + B3_10 * ( vx[id2(iz+3,ix  )]-vx[id2(iz-3,ix  )]
                                    + vz[id2(iz  ,ix+3)]-vz[id2(iz  ,ix-3)])
                          + B2_10 * ( vx[id2(iz+2,ix  )]-vx[id2(iz-2,ix  )]
                                    + vz[id2(iz  ,ix+2)]-vz[id2(iz  ,ix-2)])
                          + B1_10 * ( vx[id2(iz+1,ix  )]-vx[id2(iz-1,ix  )]
                                    + vz[id2(iz  ,ix+1)]-vz[id2(iz  ,ix-1)]));
        }
    }
return;
}
void timeStep_Elastic(float *u, float *lap, int nx, int nz, float dt, float dx) {

    float dt_dx = dt/dx;
    for(int ix=0; ix<nx; ix++) {
        for(int iz=0; iz<nz; iz++) {
            int i = id2(iz,ix);
            u[i] = u[i] + dt_dx * lap[i];
        }
    }
return;
}
void injectSource_elastic(float *sxx, float *szz, float *source, int it, int ix0, int iz0, int nx, int nz, float dx, float dz)  {
    int   i = id2(iz0,ix0);
    //sxx[i] += source[it]/(dx*dz);
    szz[i] += source[it]/(dx*dz);
return;
}

