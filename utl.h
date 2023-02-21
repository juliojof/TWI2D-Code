float* readInputVel(char filePathInp[], char filePathOut[], int nxInput_vel, int nzInput_vel, \
                  float dxInput_vel, float dzInput_vel, float oxInput_vel, float ozInput_vel, \
                  int nx, int nz, float dx, float dz, float dt, float ox, float oz, int nxb, int nzb, float velScalar);

void init_msh(MSH *msh, int nx, int nz, int nt, float dx, float dz, float dt, \
              float ox, float oz, float ot, int jt, int jt_EIC, int nxb, int nzb);

void init_InvParms(InvParms *invparms, Parms *parms);

void init_eiproc(EIProc *eiproc, Parms *parms);

void init_grdProc(GRDProc *grdProc, Parms *parms);

void init_eimsh(Eimage_MSH *eimage_msh, Parms *parms, MSH *msh);