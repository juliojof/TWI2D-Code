

// #include "Math_Includes.h"



// Hicks interpolation
float* interSinc1D_2P(float alpha) {
    int r = 2;
    int diam = 2*r;
    float b = 1.24f;
    float x, sum = 0.0f;
    float *coef = CPU_zaloc1F(diam);
    int i;
    int r_ = r-1;
    for(i=0; i<diam; i++)
    {
        x = 1.0*(i-r_) - alpha;
        coef[i] = sinc(x) * kaiser(x,r,b);
        sum += coef[i];
    }
    for(i=0; i<diam; i++) {
        coef[i] /= sum;
    }
return coef;
}
float* interSinc1D_4P(float alpha) {
    int r = 2;
    int diam = 2*r;
    float b = 2.94f;
    float x, sum = 0.0f;
    float *coef = CPU_zaloc1F(diam);
    int i;
    int r_ = r-1;
    for(i=0; i<diam; i++)
    {
        x = 1.0*(i-r_) - alpha;
        coef[i] = sinc(x) * kaiser(x,r,b);
        sum += coef[i];
    }
    for(i=0; i<diam; i++) {
        coef[i] /= sum;
    }
return coef;
}
float* interSinc1D_8P(float alpha) {
    int r = 4;
    int diam = 2*r;
    float b = 6.31f;
    float x, sum = 0.0f;
    float *coef = CPU_zaloc1F(diam);
    int i;
    int r_ = r-1;
    for(i=0; i<diam; i++)
    {
        x = 1.0*(i-r_) - alpha;
        coef[i] = sinc(x) * kaiser(x,r,b);
        sum += coef[i];
    }
    for(i=0; i<diam; i++) {
        coef[i] /= sum;
    }
return coef;
}
float* interSinc2D_8P(float alpha_x, float alpha_y) {
    int r = 4;
    int diam = 2*r;
    float b = 6.31f;
    float x, y, sum = 0.0f;
    float *coef = CPU_zaloc1F(diam*diam);
    float coef1[8];
    float coef2[8];
    int i1, i2, i;
    int r_ = r-1;

    for(i=0; i<diam; i++) {
        x = 1.0f*(i-r_) + alpha_x;
        coef1[i] = sinc(x) * kaiser(x,r,b);
    }
    for(i=0; i<diam; i++) {
        y = 1.0f*(i-r_) + alpha_y;
        coef2[i] = sinc(y) * kaiser(y,r,b);
    }
    for(i2=0; i2<diam; i2++) {
        for(i1=0; i1<diam; i1++) {
            i = i1 + i2*diam;
            coef[i] = coef1[i1]*coef2[i2];
            sum += coef[i];
        }
    }
    for(i=0; i<diam*diam; i++)    coef[i] /= sum;

return coef;
}

void interSinc3D_8P(float *coef_x, float *coef_y, float *coef_z, float alpha_x, float alpha_y, float alpha_z)
{
    int r = 4;
    int diam = 2*r;
    float b = 6.31f;
    float x, y, z, sum = 0.0f;
    int i1, i2, i3, i;
    int r_ = r-1;

    for(i=0; i<diam; i++)  {
        coef_x[i] = 0.0;
        coef_y[i] = 0.0;
        coef_z[i] = 0.0;
    }

    for(i=0; i<diam; i++) {
        x = 1.0f*(i-r_) - alpha_x;
        coef_x[i] = sinc(x) * kaiser(x,r,b);
    }
    for(i=0; i<diam; i++) {
        y = 1.0f*(i-r_) - alpha_y;
        coef_y[i] = sinc(y) * kaiser(y,r,b);
    }
    for(i=0; i<diam; i++) {
        z = 1.0f*(i-r_) - alpha_z;
        coef_z[i] = sinc(z) * kaiser(z,r,b);
    }
    for(i3=0; i3<diam; i3++)
        for(i2=0; i2<diam; i2++)
            for(i1=0; i1<diam; i1++)
            {
                sum += coef_x[i1]*coef_y[i2]*coef_z[i3];
            }
    for(i=0; i<diam; i++)    coef_x[i] /= sum;

return;
}


int interpTraceLinear(float *inp, float *out, int nin, float din, int nou, float dou) {
    int iin, iou;
    if(nou*dou>nin*din) {
        printf("\n\n Interpolated trace is longer than input trace. This function does not extrapolate. \n\n");
        return -1;
    }
    for(iou=0; iou<nou; iou++) {
        float z = iou*dou;
        iin = floor(z/din);
        float zin = iin*din;
        float p1 = (z-zin)/din;
        float p0 = 1.0f - p1;
        if(0<=iin && iin<nou-1)
            out[iou] = p0*inp[iin] + p1*inp[iin+1];
    }
return 0;
}

float* interp1DVolume3D(float *inp, int *n1, int *n2, int *n3, int interpDim, float ratio, float shift)
{
    // ratio is newSamplingRatio/OriginalSamplingRatio
    // 0<=shift<1 is the offset within one cell corresponding to ((oxOut-oxInp)/dxOut) - floorf((oxOut-oxInp)/dxOut);

    int i1, i2, i3, i;
    int i1_, i2_, i3_;

    int n1_ = *n1;
    int n2_ = *n2;
    int n3_ = *n3;
    if(interpDim==1)    n1_ = (int) ( ((*n1)-1)/ratio + 1 );
    if(interpDim==2)    n2_ = (int) ( ((*n2)-1)/ratio + 1 );
    if(interpDim==3)    n3_ = (int) ( ((*n3)-1)/ratio + 1 );

    int szNew = n1_ * n2_ * n3_;
    float *out = CPU_zaloc1F(szNew);

    int r, r_, diam;
    r = 1;
    r_ = r-1;
    diam = 2*r;

    // Loop over output indexes
    for(i3_=0; i3_<n3_; i3_++)
        for(i2_=0; i2_<n2_; i2_++)
            for(i1_=0; i1_<n1_; i1_++)
            {
                int ia;
                float p, alpha;
                if(interpDim==1)    p = ((float) i1_) *ratio + shift;
                if(interpDim==2)    p = ((float) i2_) *ratio + shift;
                if(interpDim==3)    p = ((float) i3_) *ratio + shift;
                ia    = (int) floorf(p);
                alpha = p - (1.0*ia);

                // float *coef = interSinc1D_8P(alpha);
                float *coef = interSinc1D_2P(alpha);

                int iout = i1_ + i2_ * n1_ + i3_ * n1_*n2_;
                for(i=0; i<diam; i++)
                {
                    int i1 = i1_;
                    int i2 = i2_;
                    int i3 = i3_;
                    if(interpDim==1)    i1 = ia - r_ + i;
                    if(interpDim==2)    i2 = ia - r_ + i;
                    if(interpDim==3)    i3 = ia - r_ + i;

                    if(0<=i1 && i1<(*n1)  && \
                       0<=i2 && i2<(*n2)  && \
                       0<=i3 && i3<(*n3) )
                       {
                           int iinp = i1 + i2 * (*n1) + i3 * (*n1) * (*n2);
                           out[iout] += inp[iinp] * coef[i];
                       }
                }

                free(coef);
            }

    *n1 = n1_;
    *n2 = n2_;
    *n3 = n3_;

return out;
}