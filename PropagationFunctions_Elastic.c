void propagate_Elastic(float *vx, float *vz, float *Sxx, float *Szz, float *Sxz,
                       float *lap, float *lambda, float *mu, float *b, \
                       int nx, int nz, float dx, float dt) {
    
    // Update particle velocity
    laplacian_40thOrder_Vx_Elastic_Iso(lap, b, Sxx, Sxz, nx, nz);
    timeStep_Elastic(vx, lap, nx, nz, dt, dx);
    laplacian_40thOrder_Vz_Elastic_Iso(lap, b, Sxz, Szz, nx, nz);
    timeStep_Elastic(vz, lap, nx, nz, dt, dx);

    // Update stress tensor
    laplacian_40thOrder_Sxx_Elastic_Iso(lap, lambda, mu, vx, vz, nx, nz);
    timeStep_Elastic(Sxx, lap, nx, nz, dt, dx);
    laplacian_40thOrder_Szz_Elastic_Iso(lap, lambda, mu, vx, vz, nx, nz);
    timeStep_Elastic(Szz, lap, nx, nz, dt, dx);
    laplacian_40thOrder_Sxz_Elastic_Iso(lap,         mu, vx, vz, nx, nz);
    timeStep_Elastic(Sxz, lap, nx, nz, dt, dx);

return;
}
void propagate_Elastic_abs(float *vx, float *vz, float *Sxx, float *Szz, float *Sxz,
                       float *lap, float *lambda, float *mu, float *lambda2mu, float *bc, \
                       float *bd, int nx, int nz, float dx, float dt, float *sigX, float *sigZ) {
    
    // Update particle velocity
    laplacian_40thOrder_Vx_Elastic_Iso(lap, bc, Sxx, Sxz, nx, nz);
    timeStep_Elastic_abs(vx, lap, nx, nz, dt, dx, sigX, sigZ);
    laplacian_40thOrder_Vz_Elastic_Iso(lap, bd, Sxz, Szz, nx, nz);
    timeStep_Elastic_abs(vz, lap, nx, nz, dt, dx, sigX, sigZ);

    // Update stress tensor
    laplacian_40thOrder_Sxx_Elastic_Iso(lap, lambda, lambda2mu, vx, vz, nx, nz);
    timeStep_Elastic(Sxx, lap, nx, nz, dt, dx);
    laplacian_40thOrder_Szz_Elastic_Iso(lap, lambda, lambda2mu, vx, vz, nx, nz);
    timeStep_Elastic(Szz, lap, nx, nz, dt, dx);
    laplacian_40thOrder_Sxz_Elastic_Iso(lap,  mu, vx, vz, nx, nz);
    timeStep_Elastic(Sxz, lap, nx, nz, dt, dx);

return;
}

void laplacian_10thOrder_Vz_Elastic_Iso(float *lap, float *b, float *Sdx, float *Sdz, int nx, int nz) {
    int ix, iz;
    for(ix=OL; ix<nx-OL; ix++) {
        for(iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            lap[i] = b[i]*(B5_10 * (Sdz[id2(iz+5,ix  )]-Sdz[id2(iz-4,ix  )]
                                   +Sdx[id2(iz  ,ix+5)]-Sdx[id2(iz  ,ix-4)])
                   +       B4_10 * (Sdz[id2(iz+4,ix  )]-Sdz[id2(iz-3,ix  )]
                                   +Sdx[id2(iz  ,ix+4)]-Sdx[id2(iz  ,ix-3)])
                   +       B3_10 * (Sdz[id2(iz+3,ix  )]-Sdz[id2(iz-2,ix  )]
                                   +Sdx[id2(iz  ,ix+3)]-Sdx[id2(iz  ,ix-2)])
                   +       B2_10 * (Sdz[id2(iz+2,ix  )]-Sdz[id2(iz-1,ix  )]
                                   +Sdx[id2(iz  ,ix+2)]-Sdx[id2(iz  ,ix-1)])
                   +       B1_10 * (Sdz[id2(iz+1,ix  )]-Sdz[id2(iz  ,ix  )]
                                   +Sdx[id2(iz  ,ix+1)]-Sdx[id2(iz  ,ix  )]));
        }
    }
}
void laplacian_20thOrder_Vz_Elastic_Iso(float *lap, float *b, float *Sdx, float *Sdz, int nx, int nz) {
    int ix, iz;
    for(ix=OL; ix<nx-OL; ix++) {
        for(iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            lap[i] = b[i]*(B10_20 * (Sdz[id2(iz+10,ix  )]-Sdz[id2(iz-9,ix  )]
                                    +Sdx[id2(iz  ,ix+10)]-Sdx[id2(iz  ,ix-9)])
                   +       B09_20 * (Sdz[id2(iz+9,ix  )]-Sdz[id2(iz-8,ix  )]
                                    +Sdx[id2(iz  ,ix+9)]-Sdx[id2(iz  ,ix-8)])
                   +       B08_20 * (Sdz[id2(iz+8,ix  )]-Sdz[id2(iz-7,ix  )]
                                    +Sdx[id2(iz  ,ix+8)]-Sdx[id2(iz  ,ix-7)])
                   +       B07_20 * (Sdz[id2(iz+7,ix  )]-Sdz[id2(iz-6,ix  )]
                                    +Sdx[id2(iz  ,ix+7)]-Sdx[id2(iz  ,ix-6)])
                   +       B06_20 * (Sdz[id2(iz+6,ix  )]-Sdz[id2(iz-5,ix  )]
                                    +Sdx[id2(iz  ,ix+6)]-Sdx[id2(iz  ,ix-5)])
                   +       B05_20 * (Sdz[id2(iz+5,ix  )]-Sdz[id2(iz-4,ix  )]
                                    +Sdx[id2(iz  ,ix+5)]-Sdx[id2(iz  ,ix-4)])
                   +       B04_20 * (Sdz[id2(iz+4,ix  )]-Sdz[id2(iz-3,ix  )]
                                    +Sdx[id2(iz  ,ix+4)]-Sdx[id2(iz  ,ix-3)])
                   +       B03_20 * (Sdz[id2(iz+3,ix  )]-Sdz[id2(iz-2,ix  )]
                                    +Sdx[id2(iz  ,ix+3)]-Sdx[id2(iz  ,ix-2)])
                   +       B02_20 * (Sdz[id2(iz+2,ix  )]-Sdz[id2(iz-1,ix  )]
                                    +Sdx[id2(iz  ,ix+2)]-Sdx[id2(iz  ,ix-1)])
                   +       B01_20 * (Sdz[id2(iz+1,ix  )]-Sdz[id2(iz  ,ix  )]
                                    +Sdx[id2(iz  ,ix+1)]-Sdx[id2(iz  ,ix  )]));
        }
    }
}
void laplacian_40thOrder_Vz_Elastic_Iso(float *lap, float *b, float *Sdx, float *Sdz, int nx, int nz) {
    int ix, iz;
    for(ix=OL; ix<nx-OL; ix++) {
        for(iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            lap[i] = b[i]*(B05_40 * (Sdz[id2(iz+5,ix  )]-Sdz[id2(iz-4,ix  )]
                                    +Sdx[id2(iz  ,ix+5)]-Sdx[id2(iz  ,ix-4)])
                   +       B04_40 * (Sdz[id2(iz+4,ix  )]-Sdz[id2(iz-3,ix  )]
                                    +Sdx[id2(iz  ,ix+4)]-Sdx[id2(iz  ,ix-3)])
                   +       B03_40 * (Sdz[id2(iz+3,ix  )]-Sdz[id2(iz-2,ix  )]
                                    +Sdx[id2(iz  ,ix+3)]-Sdx[id2(iz  ,ix-2)])
                   +       B02_40 * (Sdz[id2(iz+2,ix  )]-Sdz[id2(iz-1,ix  )]
                                    +Sdx[id2(iz  ,ix+2)]-Sdx[id2(iz  ,ix-1)])
                   +       B01_40 * (Sdz[id2(iz+1,ix  )]-Sdz[id2(iz  ,ix  )]
                                    +Sdx[id2(iz  ,ix+1)]-Sdx[id2(iz  ,ix  )]));
        }
    }
}
void laplacian_80thOrder_Vz_Elastic_Iso(float *lap, float *b, float *Sdx, float *Sdz, int nx, int nz) {
    int ix, iz;
    for(ix=OL; ix<nx-OL; ix++) {
        for(iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            lap[i] = b[i]*(B05_80 * (Sdz[id2(iz+5,ix  )]-Sdz[id2(iz-4,ix  )]
                                   +Sdx[id2(iz  ,ix+5)]-Sdx[id2(iz  ,ix-4)])
                   +       B04_80 * (Sdz[id2(iz+4,ix  )]-Sdz[id2(iz-3,ix  )]
                                   +Sdx[id2(iz  ,ix+4)]-Sdx[id2(iz  ,ix-3)])
                   +       B03_80 * (Sdz[id2(iz+3,ix  )]-Sdz[id2(iz-2,ix  )]
                                   +Sdx[id2(iz  ,ix+3)]-Sdx[id2(iz  ,ix-2)])
                   +       B02_80 * (Sdz[id2(iz+2,ix  )]-Sdz[id2(iz-1,ix  )]
                                   +Sdx[id2(iz  ,ix+2)]-Sdx[id2(iz  ,ix-1)])
                   +       B01_80 * (Sdz[id2(iz+1,ix  )]-Sdz[id2(iz  ,ix  )]
                                   +Sdx[id2(iz  ,ix+1)]-Sdx[id2(iz  ,ix  )]));
        }
    }
}
void laplacian_10thOrder_Vx_Elastic_Iso(float *lap, float *b, float *Sdx, float *Sdz, int nx, int nz) {
    int ix, iz;
    for(ix=OL; ix<nx-OL; ix++) {
        for(iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            lap[i] = b[i]*(B5_10 * (Sdz[id2(iz+4,ix  )]-Sdz[id2(iz-5,ix  )]
                                   +Sdx[id2(iz  ,ix+4)]-Sdx[id2(iz  ,ix-5)])
                   +       B4_10 * (Sdz[id2(iz+3,ix  )]-Sdz[id2(iz-4,ix  )]
                                   +Sdx[id2(iz  ,ix+3)]-Sdx[id2(iz  ,ix-4)])
                   +       B3_10 * (Sdz[id2(iz+2,ix  )]-Sdz[id2(iz-3,ix  )]
                                   +Sdx[id2(iz  ,ix+2)]-Sdx[id2(iz  ,ix-3)])
                   +       B2_10 * (Sdz[id2(iz+1,ix  )]-Sdz[id2(iz-2,ix  )]
                                   +Sdx[id2(iz  ,ix+1)]-Sdx[id2(iz  ,ix-2)])
                   +       B1_10 * (Sdz[id2(iz  ,ix  )]-Sdz[id2(iz-1,ix  )]
                                   +Sdx[id2(iz  ,ix  )]-Sdx[id2(iz  ,ix-1)]));
        }
    }
}
void laplacian_20thOrder_Vx_Elastic_Iso(float *lap, float *b, float *Sdx, float *Sdz, int nx, int nz) {
    int ix, iz;
    for(ix=OL; ix<nx-OL; ix++) {
        for(iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            lap[i] = b[i]*(B10_20 * (Sdz[id2(iz+9,ix  )]-Sdz[id2(iz-10,ix  )]
                                    +Sdx[id2(iz  ,ix+9)]-Sdx[id2(iz  ,ix-10)])
                   +       B09_20 * (Sdz[id2(iz+8,ix  )]-Sdz[id2(iz-9,ix  )]
                                    +Sdx[id2(iz  ,ix+8)]-Sdx[id2(iz  ,ix-9)])
                   +       B08_20 * (Sdz[id2(iz+7,ix  )]-Sdz[id2(iz-8,ix  )]
                                    +Sdx[id2(iz  ,ix+7)]-Sdx[id2(iz  ,ix-8)])
                   +       B07_20 * (Sdz[id2(iz+6,ix  )]-Sdz[id2(iz-7,ix  )]
                                    +Sdx[id2(iz  ,ix+6)]-Sdx[id2(iz  ,ix-7)])
                   +       B06_20 * (Sdz[id2(iz+5,ix  )]-Sdz[id2(iz-6,ix  )]
                                    +Sdx[id2(iz  ,ix+5)]-Sdx[id2(iz  ,ix-6)])
                   +       B05_20 * (Sdz[id2(iz+4,ix  )]-Sdz[id2(iz-5,ix  )]
                                    +Sdx[id2(iz  ,ix+4)]-Sdx[id2(iz  ,ix-5)])
                   +       B04_20 * (Sdz[id2(iz+3,ix  )]-Sdz[id2(iz-4,ix  )]
                                    +Sdx[id2(iz  ,ix+3)]-Sdx[id2(iz  ,ix-4)])
                   +       B03_20 * (Sdz[id2(iz+2,ix  )]-Sdz[id2(iz-3,ix  )]
                                    +Sdx[id2(iz  ,ix+2)]-Sdx[id2(iz  ,ix-3)])
                   +       B02_20 * (Sdz[id2(iz+1,ix  )]-Sdz[id2(iz-2,ix  )]
                                    +Sdx[id2(iz  ,ix+1)]-Sdx[id2(iz  ,ix-2)])
                   +       B01_20 * (Sdz[id2(iz  ,ix  )]-Sdz[id2(iz-1,ix  )]
                                    +Sdx[id2(iz  ,ix  )]-Sdx[id2(iz  ,ix-1)]));
        }
    }
}
void laplacian_40thOrder_Vx_Elastic_Iso(float *lap, float *b, float *Sdx, float *Sdz, int nx, int nz)
{
    int ix, iz;
    for(ix=OL; ix<nx-OL; ix++) {
        for(iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            lap[i] = b[i]*(B05_40 * (Sdz[id2(iz+4,ix  )]-Sdz[id2(iz-5,ix  )]
                                    +Sdx[id2(iz  ,ix+4)]-Sdx[id2(iz  ,ix-5)])
                   +       B04_40 * (Sdz[id2(iz+3,ix  )]-Sdz[id2(iz-4,ix  )]
                                    +Sdx[id2(iz  ,ix+3)]-Sdx[id2(iz  ,ix-4)])
                   +       B03_40 * (Sdz[id2(iz+2,ix  )]-Sdz[id2(iz-3,ix  )]
                                    +Sdx[id2(iz  ,ix+2)]-Sdx[id2(iz  ,ix-3)])
                   +       B02_40 * (Sdz[id2(iz+1,ix  )]-Sdz[id2(iz-2,ix  )]
                                    +Sdx[id2(iz  ,ix+1)]-Sdx[id2(iz  ,ix-2)])
                   +       B01_40 * (Sdz[id2(iz  ,ix  )]-Sdz[id2(iz-1,ix  )]
                                    +Sdx[id2(iz  ,ix  )]-Sdx[id2(iz  ,ix-1)]));
        }
    }
}
void laplacian_80thOrder_Vx_Elastic_Iso(float *lap, float *b, float *Sdx, float *Sdz, int nx, int nz) {
    int ix, iz;
    for(ix=OL; ix<nx-OL; ix++) {
        for(iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            lap[i] = b[i]*(B05_80 * (Sdz[id2(iz+4,ix  )]-Sdz[id2(iz-5,ix  )]
                                    +Sdx[id2(iz  ,ix+4)]-Sdx[id2(iz  ,ix-5)])
                   +       B04_80 * (Sdz[id2(iz+3,ix  )]-Sdz[id2(iz-4,ix  )]
                                    +Sdx[id2(iz  ,ix+3)]-Sdx[id2(iz  ,ix-4)])
                   +       B03_80 * (Sdz[id2(iz+2,ix  )]-Sdz[id2(iz-3,ix  )]
                                    +Sdx[id2(iz  ,ix+2)]-Sdx[id2(iz  ,ix-3)])
                   +       B02_80 * (Sdz[id2(iz+1,ix  )]-Sdz[id2(iz-2,ix  )]
                                    +Sdx[id2(iz  ,ix+1)]-Sdx[id2(iz  ,ix-2)])
                   +       B01_80 * (Sdz[id2(iz  ,ix  )]-Sdz[id2(iz-1,ix  )]
                                    +Sdx[id2(iz  ,ix  )]-Sdx[id2(iz  ,ix-1)]));
        }
    }
}
void laplacian_10thOrder_Sxx_Elastic_Iso(float *lap, float *lambda, float *lambda2mu, float *vx, float *vz, int nx, int nz) {
    int ix, iz;
    for(ix=OL; ix<nx-OL; ix++) {
        for(iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            float m_vx = lambda2mu[i];
            float m_vz = lambda[i];
            lap[i] = B5_10 * ( m_vz*(vz[id2(iz+4,ix  )]-vz[id2(iz-5,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+5)]-vx[id2(iz  ,ix-4)]))
                   + B4_10 * ( m_vz*(vz[id2(iz+3,ix  )]-vz[id2(iz-4,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+4)]-vx[id2(iz  ,ix-3)]))
                   + B3_10 * ( m_vz*(vz[id2(iz+2,ix  )]-vz[id2(iz-3,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+3)]-vx[id2(iz  ,ix-2)]))
                   + B2_10 * ( m_vz*(vz[id2(iz+1,ix  )]-vz[id2(iz-2,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+2)]-vx[id2(iz  ,ix-1)]))
                   + B1_10 * ( m_vz*(vz[id2(iz  ,ix  )]-vz[id2(iz-1,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+1)]-vx[id2(iz  ,ix  )]));
        }
    }
}
void laplacian_20thOrder_Sxx_Elastic_Iso(float *lap, float *lambda, float *lambda2mu, float *vx, float *vz, int nx, int nz) {
    int ix, iz;
    for(ix=OL; ix<nx-OL; ix++) {
        for(iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            float m_vx = lambda2mu[i];
            float m_vz = lambda[i];
            lap[i] = B10_20 * ( m_vz*(vz[id2(iz+9,ix  )]-vz[id2(iz-10,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+10)]-vx[id2(iz  ,ix-9)]))
                   + B09_20 * ( m_vz*(vz[id2(iz+8,ix  )]-vz[id2(iz-9,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+9)]-vx[id2(iz  ,ix-8)]))
                   + B08_20 * ( m_vz*(vz[id2(iz+7,ix  )]-vz[id2(iz-8,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+8)]-vx[id2(iz  ,ix-7)]))
                   + B07_20 * ( m_vz*(vz[id2(iz+6,ix  )]-vz[id2(iz-7,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+7)]-vx[id2(iz  ,ix-6)]))
                   + B06_20 * ( m_vz*(vz[id2(iz+5,ix  )]-vz[id2(iz-6,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+6)]-vx[id2(iz  ,ix-5)]))
                   + B05_20 * ( m_vz*(vz[id2(iz+4,ix  )]-vz[id2(iz-5,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+5)]-vx[id2(iz  ,ix-4)]))
                   + B04_20 * ( m_vz*(vz[id2(iz+3,ix  )]-vz[id2(iz-4,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+4)]-vx[id2(iz  ,ix-3)]))
                   + B03_20 * ( m_vz*(vz[id2(iz+2,ix  )]-vz[id2(iz-3,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+3)]-vx[id2(iz  ,ix-2)]))
                   + B02_20 * ( m_vz*(vz[id2(iz+1,ix  )]-vz[id2(iz-2,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+2)]-vx[id2(iz  ,ix-1)]))
                   + B01_20 * ( m_vz*(vz[id2(iz  ,ix  )]-vz[id2(iz-1,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+1)]-vx[id2(iz  ,ix  )]));
        }
    }
}
void laplacian_40thOrder_Sxx_Elastic_Iso(float *lap, float *lambda, float *lambda2mu, float *vx, float *vz, int nx, int nz) {
    int ix, iz;
    for(ix=OL; ix<nx-OL; ix++) {
        for(iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            float m_vx = lambda2mu[i];
            float m_vz = lambda[i];
            lap[i] = B05_40 * ( m_vz*(vz[id2(iz+4,ix  )]-vz[id2(iz-5,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+5)]-vx[id2(iz  ,ix-4)]))
                   + B04_40 * ( m_vz*(vz[id2(iz+3,ix  )]-vz[id2(iz-4,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+4)]-vx[id2(iz  ,ix-3)]))
                   + B03_40 * ( m_vz*(vz[id2(iz+2,ix  )]-vz[id2(iz-3,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+3)]-vx[id2(iz  ,ix-2)]))
                   + B02_40 * ( m_vz*(vz[id2(iz+1,ix  )]-vz[id2(iz-2,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+2)]-vx[id2(iz  ,ix-1)]))
                   + B01_40 * ( m_vz*(vz[id2(iz  ,ix  )]-vz[id2(iz-1,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+1)]-vx[id2(iz  ,ix  )]));
        }
    }
}
void laplacian_80thOrder_Sxx_Elastic_Iso(float *lap, float *lambda, float *lambda2mu, float *vx, float *vz, int nx, int nz) {
    int ix, iz;
    for(ix=OL; ix<nx-OL; ix++) {
        for(iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            float m_vx = lambda2mu[i];
            float m_vz = lambda[i];
            lap[i] = B05_80 * ( m_vz*(vz[id2(iz+4,ix  )]-vz[id2(iz-5,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+5)]-vx[id2(iz  ,ix-4)]))
                   + B04_80 * ( m_vz*(vz[id2(iz+3,ix  )]-vz[id2(iz-4,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+4)]-vx[id2(iz  ,ix-3)]))
                   + B03_80 * ( m_vz*(vz[id2(iz+2,ix  )]-vz[id2(iz-3,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+3)]-vx[id2(iz  ,ix-2)]))
                   + B02_80 * ( m_vz*(vz[id2(iz+1,ix  )]-vz[id2(iz-2,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+2)]-vx[id2(iz  ,ix-1)]))
                   + B01_80 * ( m_vz*(vz[id2(iz  ,ix  )]-vz[id2(iz-1,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+1)]-vx[id2(iz  ,ix  )]));
        }
    }
}
void laplacian_10thOrder_Szz_Elastic_Iso(float *lap, float *lambda, float *lambda2mu, float *vx, float *vz, int nx, int nz) {
    int ix, iz;
    for(ix=OL; ix<nx-OL; ix++) {
        for(iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            float m_vz = lambda2mu[i];
            float m_vx = lambda[i];
            lap[i] = B5_10 * ( m_vz*(vz[id2(iz+4,ix  )]-vz[id2(iz-5,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+5)]-vx[id2(iz  ,ix-4)]))
                   + B4_10 * ( m_vz*(vz[id2(iz+3,ix  )]-vz[id2(iz-4,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+4)]-vx[id2(iz  ,ix-3)]))
                   + B3_10 * ( m_vz*(vz[id2(iz+2,ix  )]-vz[id2(iz-3,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+3)]-vx[id2(iz  ,ix-2)]))
                   + B2_10 * ( m_vz*(vz[id2(iz+1,ix  )]-vz[id2(iz-2,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+2)]-vx[id2(iz  ,ix-1)]))
                   + B1_10 * ( m_vz*(vz[id2(iz  ,ix  )]-vz[id2(iz-1,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+1)]-vx[id2(iz  ,ix  )]));
        }
    }
}
void laplacian_20thOrder_Szz_Elastic_Iso(float *lap, float *lambda, float *lambda2mu, float *vx, float *vz, int nx, int nz) {
    int ix, iz;
    for(ix=OL; ix<nx-OL; ix++) {
        for(iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            float m_vz = lambda2mu[i];
            float m_vx = lambda[i];
            lap[i] = 
                   + B10_20 * ( m_vz*(vz[id2(iz+9,ix  )]-vz[id2(iz-10,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+10)]-vx[id2(iz  ,ix-9)]))
                   + B09_20 * ( m_vz*(vz[id2(iz+8,ix  )]-vz[id2(iz-9,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+9)]-vx[id2(iz  ,ix-8)]))
                   + B08_20 * ( m_vz*(vz[id2(iz+7,ix  )]-vz[id2(iz-8,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+8)]-vx[id2(iz  ,ix-7)]))
                   + B07_20 * ( m_vz*(vz[id2(iz+6,ix  )]-vz[id2(iz-7,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+7)]-vx[id2(iz  ,ix-6)]))
                   + B06_20 * ( m_vz*(vz[id2(iz+5,ix  )]-vz[id2(iz-6,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+6)]-vx[id2(iz  ,ix-5)]))
                   + B05_20 * ( m_vz*(vz[id2(iz+4,ix  )]-vz[id2(iz-5,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+5)]-vx[id2(iz  ,ix-4)]))
                   + B04_20 * ( m_vz*(vz[id2(iz+3,ix  )]-vz[id2(iz-4,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+4)]-vx[id2(iz  ,ix-3)]))
                   + B03_20 * ( m_vz*(vz[id2(iz+2,ix  )]-vz[id2(iz-3,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+3)]-vx[id2(iz  ,ix-2)]))
                   + B02_20 * ( m_vz*(vz[id2(iz+1,ix  )]-vz[id2(iz-2,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+2)]-vx[id2(iz  ,ix-1)]))
                   + B01_20 * ( m_vz*(vz[id2(iz  ,ix  )]-vz[id2(iz-1,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+1)]-vx[id2(iz  ,ix  )]));
        }
    }
}
void laplacian_40thOrder_Szz_Elastic_Iso(float *lap, float *lambda, float *lambda2mu, float *vx, float *vz, int nx, int nz) {
    int ix, iz;
    for(ix=OL; ix<nx-OL; ix++) {
        for(iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            float m_vz = lambda2mu[i];
            float m_vx = lambda[i];
            lap[i] = B05_40 * ( m_vz*(vz[id2(iz+4,ix  )]-vz[id2(iz-5,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+5)]-vx[id2(iz  ,ix-4)]))
                   + B04_40 * ( m_vz*(vz[id2(iz+3,ix  )]-vz[id2(iz-4,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+4)]-vx[id2(iz  ,ix-3)]))
                   + B03_40 * ( m_vz*(vz[id2(iz+2,ix  )]-vz[id2(iz-3,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+3)]-vx[id2(iz  ,ix-2)]))
                   + B02_40 * ( m_vz*(vz[id2(iz+1,ix  )]-vz[id2(iz-2,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+2)]-vx[id2(iz  ,ix-1)]))
                   + B01_40 * ( m_vz*(vz[id2(iz  ,ix  )]-vz[id2(iz-1,ix  )])
                              + m_vx*(vx[id2(iz  ,ix+1)]-vx[id2(iz  ,ix  )]));
        }
    }
}
void laplacian_80thOrder_Szz_Elastic_Iso(float *lap, float *lambda, float *lambda2mu, float *vx, float *vz, int nx, int nz) {
    int ix, iz;
    for(ix=OL; ix<nx-OL; ix++) {
        for(iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            float m_vz = lambda2mu[i];
            float m_vx = lambda[i];
            lap[i] = B05_80 * ( m_vz*(vz[id2(iz+4,ix  )]-vz[id2(iz-5,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+5)]-vx[id2(iz  ,ix-4)]))
                   + B04_80 * ( m_vz*(vz[id2(iz+3,ix  )]-vz[id2(iz-4,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+4)]-vx[id2(iz  ,ix-3)]))
                   + B03_80 * ( m_vz*(vz[id2(iz+2,ix  )]-vz[id2(iz-3,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+3)]-vx[id2(iz  ,ix-2)]))
                   + B02_80 * ( m_vz*(vz[id2(iz+1,ix  )]-vz[id2(iz-2,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+2)]-vx[id2(iz  ,ix-1)]))
                   + B01_80 * ( m_vz*(vz[id2(iz  ,ix  )]-vz[id2(iz-1,ix  )])
                             + m_vx*(vx[id2(iz  ,ix+1)]-vx[id2(iz  ,ix  )]));
        }
    }
}
void laplacian_10thOrder_Sxz_Elastic_Iso(float *lap, float *mu, float *vx, float *vz, int nx, int nz) {
    int ix, iz;
    for(ix=OL; ix<nx-OL; ix++) {
        for(iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            float m_v = mu[i];
            lap[i] = m_v * (B5_10 * ( vx[id2(iz+5,ix  )]-vx[id2(iz-4,ix  )]
                                    + vz[id2(iz  ,ix+4)]-vz[id2(iz  ,ix-5)])
                          + B4_10 * ( vx[id2(iz+4,ix  )]-vx[id2(iz-3,ix  )]
                                    + vz[id2(iz  ,ix+3)]-vz[id2(iz  ,ix-4)])
                          + B3_10 * ( vx[id2(iz+3,ix  )]-vx[id2(iz-2,ix  )]
                                    + vz[id2(iz  ,ix+2)]-vz[id2(iz  ,ix-3)])
                          + B2_10 * ( vx[id2(iz+2,ix  )]-vx[id2(iz-1,ix  )]
                                    + vz[id2(iz  ,ix+1)]-vz[id2(iz  ,ix-2)])
                          + B1_10 * ( vx[id2(iz+1,ix  )]-vx[id2(iz  ,ix  )]
                                    + vz[id2(iz  ,ix  )]-vz[id2(iz  ,ix-1)]));
        }
    }
return;
}
void laplacian_20thOrder_Sxz_Elastic_Iso(float *lap, float *mu, float *vx, float *vz, int nx, int nz) {
    int ix, iz;
    for(ix=OL; ix<nx-OL; ix++) {
        for(iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            float m_v = mu[i];
            lap[i] = m_v * (
                          + B10_20 * ( vx[id2(iz+10,ix  )]-vx[id2(iz-9,ix  )]
                                     + vz[id2(iz  ,ix+9)]-vz[id2(iz  ,ix-10)])
                          + B09_20 * ( vx[id2(iz+9,ix  )]-vx[id2(iz-8,ix  )]
                                     + vz[id2(iz  ,ix+8)]-vz[id2(iz  ,ix-9)])
                          + B08_20 * ( vx[id2(iz+8,ix  )]-vx[id2(iz-7,ix  )]
                                     + vz[id2(iz  ,ix+7)]-vz[id2(iz  ,ix-8)])
                          + B07_20 * ( vx[id2(iz+7,ix  )]-vx[id2(iz-6,ix  )]
                                     + vz[id2(iz  ,ix+6)]-vz[id2(iz  ,ix-7)])
                          + B06_20 * ( vx[id2(iz+6,ix  )]-vx[id2(iz-5,ix  )]
                                     + vz[id2(iz  ,ix+5)]-vz[id2(iz  ,ix-6)])
                          + B05_20 * ( vx[id2(iz+5,ix  )]-vx[id2(iz-4,ix  )]
                                     + vz[id2(iz  ,ix+4)]-vz[id2(iz  ,ix-5)])
                          + B04_20 * ( vx[id2(iz+4,ix  )]-vx[id2(iz-3,ix  )]
                                     + vz[id2(iz  ,ix+3)]-vz[id2(iz  ,ix-4)])
                          + B03_20 * ( vx[id2(iz+3,ix  )]-vx[id2(iz-2,ix  )]
                                     + vz[id2(iz  ,ix+2)]-vz[id2(iz  ,ix-3)])
                          + B02_20 * ( vx[id2(iz+2,ix  )]-vx[id2(iz-1,ix  )]
                                     + vz[id2(iz  ,ix+1)]-vz[id2(iz  ,ix-2)])
                          + B01_20 * ( vx[id2(iz+1,ix  )]-vx[id2(iz  ,ix  )]
                                     + vz[id2(iz  ,ix  )]-vz[id2(iz  ,ix-1)]));
        }
    }
return;
}
void laplacian_40thOrder_Sxz_Elastic_Iso(float *lap, float *mu, float *vx, float *vz, int nx, int nz) {
    int ix, iz;
    for(ix=OL; ix<nx-OL; ix++) {
        for(iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            float m_v = mu[i];
            lap[i] = m_v * (B05_40 * ( vx[id2(iz+5,ix  )]-vx[id2(iz-4,ix  )]
                                     + vz[id2(iz  ,ix+4)]-vz[id2(iz  ,ix-5)])
                          + B04_40 * ( vx[id2(iz+4,ix  )]-vx[id2(iz-3,ix  )]
                                     + vz[id2(iz  ,ix+3)]-vz[id2(iz  ,ix-4)])
                          + B03_40 * ( vx[id2(iz+3,ix  )]-vx[id2(iz-2,ix  )]
                                     + vz[id2(iz  ,ix+2)]-vz[id2(iz  ,ix-3)])
                          + B02_40 * ( vx[id2(iz+2,ix  )]-vx[id2(iz-1,ix  )]
                                     + vz[id2(iz  ,ix+1)]-vz[id2(iz  ,ix-2)])
                          + B01_40 * ( vx[id2(iz+1,ix  )]-vx[id2(iz  ,ix  )]
                                     + vz[id2(iz  ,ix  )]-vz[id2(iz  ,ix-1)]));
        }
    }
return;
}
void laplacian_80thOrder_Sxz_Elastic_Iso(float *lap, float *mu, float *vx, float *vz, int nx, int nz) {
    int ix, iz;
    for(ix=OL; ix<nx-OL; ix++) {
        for(iz=OL; iz<nz-OL; iz++) {
            int i = id2(iz,ix);
            float m_v = mu[i];
            lap[i] = m_v * (B05_80 * ( vx[id2(iz+5,ix  )]-vx[id2(iz-4,ix  )]
                                     + vz[id2(iz  ,ix+4)]-vz[id2(iz  ,ix-5)])
                          + B04_80 * ( vx[id2(iz+4,ix  )]-vx[id2(iz-3,ix  )]
                                     + vz[id2(iz  ,ix+3)]-vz[id2(iz  ,ix-4)])
                          + B03_80 * ( vx[id2(iz+3,ix  )]-vx[id2(iz-2,ix  )]
                                     + vz[id2(iz  ,ix+2)]-vz[id2(iz  ,ix-3)])
                          + B02_80 * ( vx[id2(iz+2,ix  )]-vx[id2(iz-1,ix  )]
                                     + vz[id2(iz  ,ix+1)]-vz[id2(iz  ,ix-2)])
                          + B01_80 * ( vx[id2(iz+1,ix  )]-vx[id2(iz  ,ix  )]
                                     + vz[id2(iz  ,ix  )]-vz[id2(iz  ,ix-1)]));
        }
    }
return;
}
void timeStep_Elastic(float *u, float *lap, int nx, int nz, float dt, float dx)
{
    float dt_dx = dt/dx;
    int ix, iz;
    for(ix=0; ix<nx; ix++) {
        for(iz=0; iz<nz; iz++) {
            int i = id2(iz,ix);
            u[i] = u[i] + dt_dx * lap[i];
        }
    }
return;
}
void timeStep_Elastic_abs(float *u, float *lap, int nx, int nz, float dt, float dx, float *sigX, float *sigZ) {

    float dt_dx = dt/dx;
    float sigma;
    int ix, iz;
    for(ix=0; ix<nx; ix++) {
        for(iz=0; iz<nz; iz++) {
            int i = id2(iz,ix);
            sigma = sigX[ix] + sigZ[iz];
            u[i] = (1.0f - dt*sigma) * u[i] + dt_dx * lap[i];
        }
    }
return;
}
void injectSource_elastic_Szz(float *szz, float *source, int it, int ix0, int iz0, int nx, int nz, float dx, float dz)  {
    int   i = id2(iz0,ix0);
    szz[i] += source[it]/(dx*dz);
return;
}
void injectSource_elastic_Omni(float *sxx, float *szz, float *source, int it, int ix0, int iz0, int nx, int nz, float dx, float dz)  {
    int   i = id2(iz0,ix0);
    //x_Szz = ix0*dx
    szz[i] += source[it]/(dx*dz);
    sxx[i] += source[it]/(dx*dz);
return;
}

void recordData_Elastic_Omni(Shot *shot, float *Sxx, float *Szz, int nx, int nz, int it) {
    int ir;
    int nt = shot->nt;
    int nrec = shot->nrec;
    for(ir=0; ir<nrec; ir++) {
        int ixr = shot->ixr[ir];
        int izr = shot->izr[ir];
        int i = id2(izr,ixr);
        shot->seismogram[ir*nt+it] = Sxx[i] + Szz[i];
    }
return;
}

void recordData_Elastic_Vertical(Shot *shot, float *vz, int nx, int nz, int it) {
    int ir;
    int nt = shot->nt;
    int nrec = shot->nrec;
    for(ir=0; ir<nrec; ir++) {
        int ixr = shot->ixr[ir];
        int izr = shot->izr[ir];
        int i = id2(izr,ixr);
        shot->seismogram[ir*nt+it] = vz[i];
    }
return;
}

void recordData_Elastic_phi(Shot *shot, float *vx, float *vz, int nx, int nz, float dx, float dz, int it) {
    float sample;
    int ir;
    int nt = shot->nt;
    int nrec = shot->nrec;
    for(ir=0; ir<nrec; ir++) {
        int ixr = shot->ixr[ir];
        int izr = shot->izr[ir];
        sample = B05_40 * ( vz[id2(izr+4,ixr  )]-vz[id2(izr-5,ixr  )]
                          + vx[id2(izr  ,ixr+5)]-vx[id2(izr  ,ixr-4)])
               + B04_40 * ( vz[id2(izr+3,ixr  )]-vz[id2(izr-4,ixr  )]
                          + vx[id2(izr  ,ixr+4)]-vx[id2(izr  ,ixr-3)])
               + B03_40 * ( vz[id2(izr+2,ixr  )]-vz[id2(izr-3,ixr  )]
                          + vx[id2(izr  ,ixr+3)]-vx[id2(izr  ,ixr-2)])
               + B02_40 * ( vz[id2(izr+1,ixr  )]-vz[id2(izr-2,ixr  )]
                          + vx[id2(izr  ,ixr+2)]-vx[id2(izr  ,ixr-1)])
               + B01_40 * ( vz[id2(izr  ,ixr  )]-vz[id2(izr-1,ixr  )]
                          + vx[id2(izr  ,ixr+1)]-vx[id2(izr  ,ixr  )]);
        shot->seismogram[ir*nt+it] = sample/dx;
    }
return;
}