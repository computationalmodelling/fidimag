#include "clib.h"
#include "math.h"

void sd_step (double *spin, double *spin_last, double *field, double *mxH,
              double *mxmxH, double *mxmxH, double *tau,
              int n) {

    // A step for the SD method
    
    double mxH_sq;
    double factor_plus;
    double factor_minus;
    double new_spin[3];


    // Update the field -------------------------------------------------------
    for (int i = 0; i < n; i++) {
        spin_idx = 3 * i; 
        
        mxH_sq = 0;
        for (int j = 0; j < 3; j++) {
            mxH_sq += mxH[spin_idx + j] * mxH[spin_idx + j]
        }

        factor_plus = 4 + tau[i] * tau[i] * mxH_sq;
        factor_minus = 4 - tau[i] * tau[i] * mxH_sq;

        for (int j = 0; j < 3; j++) {
            new_spin[j] = factor_minus * spin[spin_idx + j] 
                          - 4 * tau[i] * mxmxH[spin_idx + j];

            spin[spin_idx + j] = new_spin[j] / factor_plus;
        }



}
