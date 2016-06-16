#include "clib.h"

void uniform_random_sphere(double *phi, double *ct, int n){
	for (int i=0;i<n;i++){
		phi[i] = single_random()*2*M_PI;
		ct[i] = 2*single_random()-1;
	}
}

/*
 * n is the total spin number
 */
void random_spin_uniform(double *spin, int n){
    for (int i=0;i<n;i++){
        int j=3*i;
        double phi= single_random()*2*M_PI;
        double ct = 2*single_random()-1;
        double st = sqrt(1-ct*ct);
        spin[j] = st*cos(phi); //mx = sin(theta)*cos(phi)
        spin[j+1] = st*sin(phi); //my = sin(theta)*sin(phi)
        spin[j+2] = ct; //mz = cos(theta)
    }

}

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

/*
 * n is the total spin number
 */
double compute_energy_difference(double *spin, double *new_spin, int *ngbs, double J, double D, double *h, int i, int n){
    
    int id_nn = 6 * i;
    double energy1 = -dot(&spin[3*i], &h[3*i]);
    double energy2 = -dot(&new_spin[3*i], &h[3*i]);
    
    for (int j = 0; j < 6; j++) {
        int k = ngbs[id_nn + j];
        if (k >= 0) {
            energy1 -= J*dot(&spin[3*i], &spin[3*k]);
            energy1 += D*dmi_energy_site(&spin[3*i], &spin[3*k], j);
            
            energy2 -= J*dot(&new_spin[3*i], &spin[3*k]);
            energy2 += D*dmi_energy_site(&new_spin[3*i], &spin[3*k], j);
            
        }
    }

    return energy2-energy1;

}


void run_step_mc(double *spin, double *new_spin, int *ngbs, double J, double D, double *h, int n, double T){
    
    random_spin_uniform(&new_spin[0], n);
    
    double ed, r;
    int update=0;
    for(int i=0;i<n;i++){
        int j=3*i;
        ed = compute_energy_difference(&spin[0], &new_spin[0], &ngbs[0], J, D, &h[0], i, n);
        
        if (ed<0) {update=1;}
        else{
            r= single_random();
            if (r<exp(-ed/T)) update=1;
        }
    
        if(update){
            spin[j]=new_spin[j];
            spin[j+1]=new_spin[j+1];
            spin[j+2]=new_spin[j+2];
        }
        
        update = 0;
    
    }

}
