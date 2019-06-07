#include<cmath>
#include<iostream>

class Energy {
public:
    Energy() {std::cout << "In A; at " << this << "\n";};
    ~Energy() {std::cout << "Killing A\n";};
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
    void setup(int nx, int ny, int nz, double dx, double dy, double dz,
               double unit_length, double *spin, double *Ms, double *Ms_inv,
               double *coordinates, double *ngbs, 
               double *energy, double *field
               );
    virtual void compute_field(double t) {};
};

class ExchangeEnergy : public Energy {
public:
    ExchangeEnergy() {std::cout << "In B; at " << this << "\n";};
    void init(double *A) {
        this->set_up = false;
        this->A = A;
    }
    double *A;
    void compute_field(double t);
};
