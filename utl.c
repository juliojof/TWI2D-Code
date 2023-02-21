float* readInputVel(char filePathInp[], char filePathOut[], int nxInput, int nzInput, \
                  float dxInput, float dzInput, float oxInput, float ozInput, \
                  int nx, int nz, float dx, float dz, float dt, float ox, float oz, int nxb, int nzb, float velScalar) {
    
    float *vel_temp = readModel(filePathInp, nxInput, nzInput);
    if(vel_temp==NULL) {
        printf("\n WARNING: DID NOT FIND INPUT MODEL: <%s> \n", filePathInp);
        return NULL;
    }   
    float *vel = extendResampInputModel(vel_temp, nxInput, nzInput, dxInput, dzInput, oxInput, ozInput, \
                                 nx, nz, dx, dz, ox+nxb*dx, oz+nzb*dz, nxb, nzb, velScalar);
    free(vel_temp);
    padBoundaries(vel, nx, nz, nxb, nzb);
    transformVel(vel, nx, nz, dx, dt);
    MPI_Barrier(MPI_COMM_WORLD);

return vel;
}

void init_msh(MSH *msh, int nx, int nz, int nt, float dx, float dz, float dt, \
              float ox, float oz, float ot, int jt, int jt_EIC, int nxb, int nzb)
{
    msh->nxb     = nxb;
    msh->nzb     = nzb;
    msh->nx      = nx;
    msh->ny      = 1;
    msh->nz      = nz;
    msh->nt      = nt;
    msh->dx      = dx;
    msh->dy      = dx;
    msh->dz      = dz;
    msh->dt      = dt;
    msh->ox      = ox;
    msh->oy      = 0.0;
    msh->oz      = oz;
    msh->ot      = ot;
    msh->jt      = jt;
    msh->jt_EIC  = jt_EIC;
    msh->LX      = (msh->nx-1)*msh->dx;
    msh->LY      = (msh->ny-1)*msh->dy;
    msh->LZ      = (msh->nz-1)*msh->dz;
    msh->LT      = (msh->nt-1)*msh->dt;
}

void init_InvParms(InvParms *invparms, Parms *parms)
{
    invparms->n_iter = parms->n_iter;
    invparms->alpha  = parms->alpha;
}

void init_eiproc(EIProc *eiproc, Parms *parms)
{

    if(EIMSH_CIGS_TAU) {
        eiproc->tauNull    = parms->EIMSH_tauNull;
    }
    eiproc->muFocus    = parms->muFocus;
    eiproc->epsFocus   = parms->epsFocus;
    eiproc->lambdaNull = parms->lambdaNull;
    eiproc->zMin       = parms->eimageFocus_zMin;
    eiproc->zMax       = parms->eimageFocus_zMax;
    
    float muFocus, epsFocus, lambdaNull;

    eiproc->lambda1             = parms->lambda1;
    eiproc->lambda2             = parms->lambda2;
    eiproc->lambda3             = parms->lambda3;
    eiproc->lambda4             = parms->lambda4;
    eiproc->applyBandass2EImage = parms->applyBandass2EImage;

    eiproc->applyZTaper   =  parms->eimageTaper_applyZTaper;
    eiproc->z0            =  parms->eimageTaper_z0;
    eiproc->z1            =  parms->eimageTaper_z1;
    eiproc->z2            =  parms->eimageTaper_z2;
    eiproc->z3            =  parms->eimageTaper_z3;

    eiproc->applyXTaper    =  parms->eimageTaper_applyXTaper;
    eiproc->x1             =  parms->eimageTaper_x1;
    eiproc->x2             =  parms->eimageTaper_x2;
    eiproc->x3             =  parms->eimageTaper_x3;
    eiproc->x4             =  parms->eimageTaper_x4;

    eiproc->applyLXTaper   =  parms->eimageTaper_applyLXTaper;
    eiproc->z_i            =  parms->eimageTaper_z_i;
    eiproc->z_f            =  parms->eimageTaper_z_f;
    eiproc->lx1_i          =  parms->eimageTaper_lx1_i;
    eiproc->lx2_i          =  parms->eimageTaper_lx2_i;
    eiproc->lx3_i          =  parms->eimageTaper_lx3_i;
    eiproc->lx4_i          =  parms->eimageTaper_lx4_i;
    eiproc->lx1_f          =  parms->eimageTaper_lx1_f;
    eiproc->lx2_f          =  parms->eimageTaper_lx2_f;
    eiproc->lx3_f          =  parms->eimageTaper_lx3_f;
    eiproc->lx4_f          =  parms->eimageTaper_lx4_f;

    eiproc->applyTaperAngle = parms->eimageTaper_applyTaperAngle;
    eiproc->theta1          = parms->eimageTaper_theta1;
    eiproc->theta2          = parms->eimageTaper_theta2;
    eiproc->theta3          = parms->eimageTaper_theta3;
    eiproc->theta4          = parms->eimageTaper_theta4;

    eiproc->applyADFocusing     =  parms->applyADFocusing;
    eiproc->makeReciprocalADCIG =  parms->makeReciprocalADCIG;
}


void init_grdProc(GRDProc *grdProc, Parms *parms)
{
    grdProc->applyZTaper     = parms->gradProc_applyZTaper;
    grdProc->z1              = parms->gradProc_z1;
    grdProc->z2              = parms->gradProc_z2;
    grdProc->z3              = parms->gradProc_z3;
    grdProc->z4              = parms->gradProc_z4;
    grdProc->applyXTaper     = parms->gradProc_applyXTaper;
    grdProc->x1              = parms->gradProc_x1;
    grdProc->x2              = parms->gradProc_x2;
    grdProc->x3              = parms->gradProc_x3;
    grdProc->x4              = parms->gradProc_x4;
    grdProc->applySmoothing  = parms->gradProc_applySmoothing;
    grdProc->halfLengthZ     = parms->gradProc_halfLengthZ;
    grdProc->halfLengthX     = parms->gradProc_halfLengthX;
}


void init_eimsh(Eimage_MSH *eimage_msh, Parms *parms, MSH *msh)
{
    EIMSH_CIGS_ORTH = parms->EIMSH_CIGS_ORTH;
    EIMSH_CIGS_VCIG = parms->EIMSH_CIGS_VCIG;
    EIMSH_CIGS_TAU  = parms->EIMSH_CIGS_TAU;

    int    itheta;
    int    lx, lz, lt;
    
    if(EIMSH_CIGS_TAU) 
    {
        float tauMax = parms->EIMSH_tauMax;
        float dtau   = parms->EIMSH_dtau;
        int ntau  = 1 + 2 * (int) (0.5 + tauMax/dtau);
        int tau0  = ntau/2;

        eimage_msh->ntau   = ntau;
        eimage_msh->tau0   = tau0;
        eimage_msh->dtau   = dtau;
        eimage_msh->tauMax = tauMax;
    }
    
    float lambdaMax = parms->EIMSH_lambdaMax;
    float dlambda   = parms->EIMSH_dlambda;

    float dipMin    = parms->EIMSH_dipMin;
    float dipMax    = parms->EIMSH_dipMax;
    float ddip      = parms->EIMSH_ddip; 

    float azmMin    = parms->EIMSH_azmMin;
    float azmMax    = parms->EIMSH_azmMax;
    float dazm      = parms->EIMSH_dazm;

    float maxDip = max(fabsf(dipMin), fabsf(dipMax));
    float maxAzm = max(fabsf(dipMin), fabsf(dipMax));

    int    jt_lt = parms->jt_lt;
    int    lt0, nlt;
    if(jt_lt>0)    lt0   = ceilf(parms->maxSubOff_lt/(jt_lt*msh->dt));
    else           lt0   = 0;
    nlt=2*lt0+1;
    int lx0, ly0, lz0, lv0, lambda0, nlx, nly, nlz, nlv, nlambda, ndip;
    nlambda = 1 + 2 * (int) (0.5 + lambdaMax/dlambda);
    lambda0 = nlambda/2;
    ndip    = 1 +     (int) (0.5 + (dipMax-dipMin)/ddip);
    lx0   = ceilf(lambdaMax/dlambda);                           
    ly0   = ceilf(lambdaMax * sinf(M_PI*maxAzm/180.0f)/dlambda);
    nlx=2*lx0+1;
    nly=2*ly0+1;
    

    long long int lz0_tmp;
    if(parms->EIMSH_CIGS_ORTH==0)
    {
        lz0   = ceilf(lambdaMax * sinf(M_PI*maxDip/180.0f)/dlambda);
        nlz=2*lz0+1;
        lv0=0;
        nlv=2*lv0+1;
        lz0_tmp = lz0;
    }
    else
    {
        lz0 = 0.0;
        nlz=2*lz0+1;
        lv0 = ceilf(lambdaMax * sinf(M_PI*maxDip/180.0f)/dlambda);
        nlv=2*lv0+1;
        if(EIMSH_CIGS_VCIG)    lz0_tmp = lv0;
        else                   lz0_tmp = 0;
    }

    // printf("\n\n DEBUG lv0=%d  maxDip=%f  dlambda=%f  lambdaMax=%f\n\n", lv0, maxDip, dlambda, lambdaMax);

    int    nx_LambdaDip = (msh->nx-2*msh->nxb) * cosf(M_PI*maxDip/180.0f) + (msh->nz-2*msh->nzb) * sinf(M_PI*maxDip/180.0f);
    int    nz_LambdaDip = (msh->nx-2*msh->nxb) * sinf(M_PI*maxDip/180.0f) + (msh->nz-2*msh->nzb) * cosf(M_PI*maxDip/180.0f);
    float  ox_LambdaDip = ((msh->nx-2*msh->nxb)*msh->dx)/2.0;
    float  oz_LambdaDip = ((msh->nz-2*msh->nzb)*msh->dz)/2.0;

    eimage_msh->azmMin    = parms->EIMSH_azmMin;
    eimage_msh->azmMax    = parms->EIMSH_azmMax;
    eimage_msh->dazm      = parms->EIMSH_dazm;

    eimage_msh->dipMin    = parms->EIMSH_dipMin;
    eimage_msh->dipMax    = parms->EIMSH_dipMax;
    eimage_msh->ddip      = parms->EIMSH_ddip;

    eimage_msh->dlambda   = parms->EIMSH_dlambda;
    eimage_msh->lambdaMax = parms->EIMSH_lambdaMax;

    eimage_msh->nlambda   = nlambda;
    eimage_msh->ndip      = ndip;

    eimage_msh->nx_LambdaDip = nx_LambdaDip;
    eimage_msh->nz_LambdaDip = nz_LambdaDip;
    eimage_msh->ox_LambdaDip = ox_LambdaDip;
    eimage_msh->oz_LambdaDip = oz_LambdaDip;

    eimage_msh->dlx = parms->EIMSH_dlambda;
    eimage_msh->dlz = parms->EIMSH_dlambda;
    eimage_msh->dlv = parms->EIMSH_dlambda;
    eimage_msh->dly = parms->EIMSH_dlambda;
    eimage_msh->dlt = parms->EIMSH_dlambda;
    eimage_msh->lx0 = lx0;
    eimage_msh->lz0 = lz0;
    eimage_msh->lv0 = lv0;
    eimage_msh->ly0 = ly0;
    eimage_msh->lt0 = lt0;
    eimage_msh->nlx = nlx;
    eimage_msh->nlz = nlz;
    eimage_msh->nlv = nlv;
    eimage_msh->nly = nly;
    eimage_msh->nlt = nlt;

    float  otheta  = -parms->EIMSH_thetaMax;
    float  dtheta  =  parms->EIMSH_dtheta;
    int    ltheta0 =  parms->EIMSH_thetaMax/dtheta;
    int    ntheta  =  2*ltheta0+1;

    eimage_msh->ltheta0 = ltheta0;
    eimage_msh->ntheta  = ntheta;
    eimage_msh->dtheta  = dtheta;
    eimage_msh->otheta  = otheta;
}
