#include "clib.h"
#include "fidimag_random.h"

inline double dot(double a[3], double b[3]){
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

inline double dmi_energy_site(double a[3], double b[3], int i){
    double res=0;
    switch (i) {
        case 0:
           res=-(a[1] * b[2] - a[2] * b[1]);break;//-x
        case 1:
           res= (a[1] * b[2] - a[2] * b[1]); break;//+x
        case 2:
           res = -(a[2] * b[0] - a[0] * b[2]); break;//-y
        case 3:
           res = (a[2] * b[0] - a[0] * b[2]); break;//+y
        case 4:
           res = -(a[0] * b[1] - a[1] * b[0]); break;//-z
        case 5:
           res = (a[0] * b[1] - a[1] * b[0]); break;//-z
        default:
        break;
    }
    return res;

}

inline double cubic_energy_site(double *m, double Kc){

    double mx2 = m[0]*m[0];
    double my2 = m[1]*m[1];
    double mz2 = m[2]*m[2];
    
    return -Kc*(mx2*mx2+my2*my2+mz2*mz2);
    
}

/*
 * n is the total spin number
 */
double compute_energy_difference(double *spin, double *new_spin, int *ngbs, double J, double D, double *h, double Kc, int i, int n){
    
    int id_nn = 6 * i;
    double energy1 = -dot(&spin[3*i], &h[3*i]); //zeeman energy
    double energy2 = -dot(&new_spin[3*i], &h[3*i]);
    
    energy1 += cubic_energy_site(&spin[3*i], Kc); //cubic anisotropy energy
    energy2 += cubic_energy_site(&new_spin[3*i], Kc);
    
    for (int j = 0; j < 6; j++) {
        int k = ngbs[id_nn + j];
        if (k >= 0) {
            energy1 -= J*dot(&spin[3*i], &spin[3*k]);//exchange energy
            energy1 += D*dmi_energy_site(&spin[3*i], &spin[3*k], j); //DMI energy
            
            energy2 -= J*dot(&new_spin[3*i], &spin[3*k]);
            energy2 += D*dmi_energy_site(&new_spin[3*i], &spin[3*k], j);
            
        }
    }

    return energy2-energy1;

}


void run_step_mc(mt19937_state *state, double *spin, double *new_spin, int *ngbs, double J, double D, double *h, double Kc, int n, double T){

    double delta_E, r;
    int update=0;
    
    uniform_random_sphere(state, new_spin, n);
    
    for(int i=0;i<n;i++){
        int j=3*i;
        delta_E = compute_energy_difference(&spin[0], &new_spin[0], &ngbs[0], J, D, &h[0], Kc, i, n);
        
        if (delta_E<0) {update=1;}
        else{
            r = random_double_half_open(state);
            if (r<exp(-delta_E/T)) update=1;
        }
    
        if(update){
            spin[j]=new_spin[j];
            spin[j+1]=new_spin[j+1];
            spin[j+2]=new_spin[j+2];
        }
        
        update = 0;
    
    }

}
