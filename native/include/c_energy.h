#include<cmath>

class Energy {
public:
    bool set_up;
    Energy();
    int nx, ny, nz, n;
    double dx, dy, dz;
    double unit_length;
    double *spin;
    double *Ms;
    double *Ms_inv;
    double *field;
    double *energy;
    double *coordinates;
    int *ngbs;
    virtual void compute_field(double t) = 0;
    double compute_energy();
    void setup(int nx, int ny, int nz, double dx, double dy, double dz, double unit_length, double *spin, double *Ms, double *Ms_inv);
};

class Exchange : public Energy {
public:
    Exchange(double *A) : A(A) {
        set_up = false;
    }
    double *A;
    void compute_field(double t);
};
