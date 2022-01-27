#include "m_driver_test.h"
#include "c_vectormath.h"
#include<cmath>

void test_f(double * x, std::vector<double>& dxdt, unsigned int n,
            double t, TestParams * params) {

    // Function: x' = x ((2 / (e^t + 1)) - 1)
    // Solution:  x(t) = 12 * e^t / (e^t + 1)^2
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
    
    // This particular template is declared in m_driver.cpp for now
    IntegratorRK4<TestParams> integratorrk4(N);
    integratorrk4.setup(N, dt);

    // For now set initial m manually:
    for (int i = 0; i < N; ++i) {
        x[i] = 1.0;
    }

    int int_steps = 20;
    for (int i = 0; i < N; ++i) {
        integratorrk4.integration_step(test_f, x, dxdt, testp);
    }
}
