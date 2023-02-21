#include "shot.h"

void modelShot_GPU(float *vp, Shot *shot, int jt, int nt, int nx, int nz, float dx, float dz, float dt, \
                   float ox, float oz, float ot, int recordMovie, int ishotMovie);