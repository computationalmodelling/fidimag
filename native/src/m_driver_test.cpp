#include "m_driver_test.h"
#include "m_driver.h"
#include "c_vectormath.h"
#include<cmath>

struct TestParams {
    double foo;
};

void test_f(double * x, std::vector<double>& dxdt, unsigned int n,
            double t, TestParams * params) {

    // Function: x' = x ((2 / (e^t + 1)) - 1)
    for (int i = 0; i < n; ++i) {
        dxdt[i] = x[i] * ((2 / (std::exp(t) + 1)) - 1); 
    }

}

void test_integrator_RK4 (void) {

    int N = 10;
    double dt = 1e-3;
    TestParams * testp;
    testp->foo = 1.;
    // TestParams * testp = &(TestParams) {.foo = 1.}; -> ?
    double x[N];
    std::vector<double> dxdt(N);
    

    IntegratorRK4<TestParams> integratorrk4(N);
    integratorrk4._setup(N, dt);

    // For now set initial m manually:
    for (int i = 0; i < N; ++i) {
        x[i] = 1.0;
    }

    integratorrk4.integration_step(test_f, x, dxdt, testp);

}
