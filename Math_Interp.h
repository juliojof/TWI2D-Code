#ifndef MATH_INTERP_H
#define MATH_INTERP_H

float* interSinc1D_2P(float alpha);
float* interSinc1D_4P(float alpha);
float* interSinc1D_8P(float alpha);
float* interSinc2D_8P(float alpha_x, float alpha_y);
void   interSinc3D_8P(float *coef_x, float *coef_y, float *coef_z, float alpha_x, float alpha_y, float alpha_z);
int interpTraceLinear(float *inp, float *out, int nin, float din, int nou, float dou);
float* interp1DVolume3D(float *inp, int *n1, int *n2, int *n3, int interpDim, float ratio, float shift);

#endif