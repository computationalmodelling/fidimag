import numpy as np

class Constant(object):
    mu_0 = 4*np.pi*1e-7
    mu_B = 9.27400949e-24
    k_B = 1.3806505e-23
    c_e = 1.602176565e-19
    m_e = 9.10938291e-31
    g_e = 2.0023193043737
    h_bar = 1.05457172647e-34
    h = h_bar*2.*np.pi
    gamma = g_e*mu_B/h_bar
    mu_s_1 = g_e*mu_B*1.0 #for S=1
    