

void Solver_TL_TWI2D(Dataset3D *dataset_base, Dataset3D *dataset_moni, float *vel, float *Dvel, InvParms *invparms, Parms *parms, 
                     MSH *msh, Eimage_MSH *eimsh, EIProc *eiproc, GRDProc *grdProc, int acquisitionGeom, float offMax, float perc, \
                     float agc_time_window, int iproc, int nproc);


float LineSearch_TL_ORTH(Dataset3D *dataset_moni, float *vel_base, float *Dvel, float *gradient, float *image, \
                         float *Hocig_dif, float *Hocig_base, float *Hocig_moni, \
                         float *Vocig_dif, float *Vocig_base, float *Vocig_moni, \
                         MSH *msh, Eimage_MSH *eimsh, EIProc *eiproc, Parms *parms, float alpha_max, double objFunc0, \
                         int acquisitionGeom, float offMax, float perc, float agc_time_window, int iproc, int nproc);


void ComputeGradient_TL_ORTH(Dataset3D *dataset, float *gradient, float *vel, float *Hocig_dif, float *Vocig_dif, \
                             MSH *msh, Eimage_MSH *eimsh, EIProc *eiproc, Parms *parms, GRDProc *grdProc, \
                             int acquisitionGeom, float offMax, float perc, float agc_time_window, int iter, int iproc, int nproc);


void ComputeFocusedCIG_TL_ORTH(Dataset3D *dataset, float *vel, float *Hocig, float *Vocig, float *image, \
                               MSH *msh, Eimage_MSH *eimsh, EIProc *eiproc, Parms *parms, \
                               int acquisitionGeom, float offMax, float perc, float agc_time_window, int iproc, int nproc, int isBaseline);

double ComputeResiduals_TL(float *Hocig_dif, float *Hocig_base, float *Hocig_moni, \
                           float *Vocig_dif, float *Vocig_base, float *Vocig_moni, \
                           MSH *msh, Eimage_MSH *eimsh, int iproc, int nproc);

void Apply_FocusingOperator_TL_ORTH(float *Hocig, float *Vocig, MSH *msh, Eimage_MSH *eimsh, EIProc *eiproc, int iproc, int nproc, int procGPU);