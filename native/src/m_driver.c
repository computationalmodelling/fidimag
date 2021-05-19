#include "m_driver.h"

void Integrator_RK4::_setup(int N) {
    this->rk_steps(N * 4, 0.0);
}
