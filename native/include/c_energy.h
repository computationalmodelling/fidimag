#include<cmath>

class Energy {
public:
    Energy();
    bool set_up;
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
    double compute_energy();
    void setup(int nx, int ny, int nz, double dx, double dy, double dz, double unit_length, double *spin, double *Ms, double *Ms_inv);
    virtual void compute_field(double t) = 0;
};

class ExchangeEnergy : public Energy {
public:
    ExchangeEnergy() {};
    void init(double *A) {
        set_up = false;
        this->A = A;
    }
    double *A;
    void compute_field(double t);
};
