void Solver_TWI2D(Dataset3D *dataset, float *vel, InvParms *invparms, Parms *parms, MSH *msh, \
                  Eimage_MSH *eimsh, EIProc *eiproc, GRDProc *grdProc, int iproc, int nproc);
                  
void updateVel(float *vel, float *gradient, MSH *msh, float alpha_optm);

float LineSearch_ORTH(Dataset3D *dataset, float *vel, float *gradient, float *Hocig_dif, float *Vocig_dif, float *image, \
                      MSH *msh, Eimage_MSH *eimsh, EIProc *eiproc, Parms *parms, float alpha_max, double objFunc0, int iproc, int nproc);

void ComputeGradient_ORTH(Dataset3D *dataset, float *gradient, float *vel, float *Hocig_dif, float *Vocig_dif, \
                          MSH *msh, Eimage_MSH *eimsh, EIProc *eiproc, Parms *parms, GRDProc *grdProc, int iter, int iproc, int nproc);

double ComputeResiduals_ORTH(Dataset3D *dataset, float *vel, float *Hocig_dif, float *Vocig_dif, float *image, \
                             MSH *msh, Eimage_MSH *eimsh, EIProc *eiproc, Parms *parms, int iproc, int nproc);

void CIG_Preprocessing_ORTH(float *Hocig, float *Vocig, EIProc *eiproc, MSH *msh, Eimage_MSH *eimsh, int iproc, int procGPU);

void Apply_FocusingOperator_ORTH(float *Hocig_ref, float *Hocig_foc, float *Vocig_ref, float *Vocig_foc, \
                                 MSH *msh, Eimage_MSH *eimsh, EIProc *eiproc, int iproc, int nproc, int procGPU);

