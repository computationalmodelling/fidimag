#include<cmath>
#include<iostream>

class Energy {
public:
    Energy() {};
    virtual ~Energy() {std::cout << "Killing A\n";};
    int interaction_id;
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
               double *coordinates, int *ngbs, 
               double *energy, double *field
               );
    virtual void compute_field(double t) {};
};

class ExchangeEnergy : public Energy {
public:
    ExchangeEnergy() {
        std::cout << "Instatiating ExchangeEnergy class; at " << this << "\n";
        this->interaction_id = 1;
    };
    ~ExchangeEnergy() {std::cout << "Killing B\n";};
    void init(double *A) {
        this->set_up = false;
        this->A = A;
    }
    double *A;
    void compute_field(double t);
};
