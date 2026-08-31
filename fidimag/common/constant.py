from math import pi

# The pre-2019 definition, which was exact until the SI redefinition turned
# mu_0 into a measured quantity. OOMMF 2.x uses the CODATA value instead,
# 12.5663706127e-7 (pkg/nb/evoc.h), which is smaller by 1.32033e-10 in
# relative terms. Every field carrying a factor of 1 / mu_0 inherits that
# difference exactly, so it is the floor on how closely the two codes can
# agree; see the OOMMF section of the installation docs. It is far below any
# modelling error, and changing it here would move every result in the last
# few digits
mu_0 = 4 * pi * 1e-7
mu_B = 9.27400949e-24
k_B = 1.3806505e-23
c_e = 1.602176565e-19
eV = 1.602176565e-19
meV = 1.602176565e-22
m_e = 9.10938291e-31
g_e = 2.0023193043737
h_bar = 1.05457172647e-34
h = h_bar * 2. * pi
gamma = g_e * mu_B / h_bar
mu_s_1 = g_e * mu_B * 1.0  # for S=1, 1.856952823077189e-23
h_bar_gamma = h_bar * gamma
