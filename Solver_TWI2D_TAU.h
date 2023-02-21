void Solver_TWI2D_TAU(Dataset3D *dataset, float *vel, InvParms *invparms, Parms *parms, MSH *msh, \
                      Eimage_MSH *eimsh, EIProc *eiproc, GRDProc *grdProc, int iproc, int nproc);

double ComputeResiduals_TAU(Dataset3D *dataset, float *vel, float *Tcig_dif, float *image, float *ilumin, \
                            MSH *msh, Eimage_MSH *eimsh, EIProc *eiproc, Parms *parms, int performEBM, int iproc, int nproc);

double ComputeResiduals_TAU_LS(Dataset3D *dataset, float *vel, float *Tcig_dif, float *image, float *ilumin, \
                               MSH *msh, Eimage_MSH *eimsh, EIProc *eiproc, Parms *parms, int performEBM, int iproc, int nproc, int suffix);

void CIG_Preprocessing_TAU(float *Tcig, EIProc *eiproc, MSH *msh, Eimage_MSH *eimsh, int iproc, int procGPU);

void ComputeGradient_TAU(Dataset3D *dataset, float *gradient, float *ilumin, float *vel, float *Tcig_dif, \
                         MSH *msh, Eimage_MSH *eimsh, EIProc *eiproc, Parms *parms, GRDProc *grdProc, \
                         int iter, int iproc, int nproc);

float LineSearch_TAU(Dataset3D *dataset, float *vel, float *gradient, float *Tcig_dif, float *image, float *ilumin, \
                     MSH *msh, Eimage_MSH *eimsh, EIProc *eiproc, Parms *parms, float alpha_max, double objFunc0, int iproc, int nproc);

void Apply_FocusingOperator_TAU(float *Tcig_ref, float *Tcig_foc, MSH *msh, \
                                Eimage_MSH *eimsh, EIProc *eiproc, int iproc, int nproc, int procGPU);

