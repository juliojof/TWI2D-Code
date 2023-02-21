#include <stdio.h>
#include <math.h>

int main() {
    double z = 700.0;
    double x = 2000.0;
    printf("\n\n atan2(%g,%g)=%g \n\n", x, z, (180.0/M_PI) * atan2(x,z));

    z = 3000.0;
    x = 8000.0;
    printf("\n\n atan2(%g,%g)=%g \n\n", x, z, (180.0/M_PI) * atan2(x,z));

return 0;
}