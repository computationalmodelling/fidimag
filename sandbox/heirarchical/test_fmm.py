import numpy as np
import fidimag
from fidimag.atomistic.energy import Energy
from fidimag.atomistic import Demag, DemagFull
import fidimag.extensions.fmm as fmm
import fidimag.extensions.bh as bh
import time

order = 9
ncrit = 48
theta=0.0

Nx=20
Ny=16
Nz=15

N = Nx*Ny*Nz

class DemagFMM(Energy): 
    def __init__(self, order, ncrit, theta, name="DemagFMM"):
        self.name = name
        assert order > 0, "Order must be 1 or higher"
        self.order = order
        assert ncrit >= 2, "ncrit must be greater than 1."
        self.ncrit = ncrit
        assert theta >= 0.0, "theta must be >= 0.0"
        self.theta = theta

    def setup(self, mesh, spin, mu_s, mu_s_inv):
        super(DemagFMM, self).setup(mesh, spin, mu_s, mu_s_inv)
        self.n = mesh.n
        print(mesh.coordinates)
        self.m_temp = sim.spin.copy()
        self.m_temp[0::3] *= self.mu_s
        self.m_temp[1::3] *= self.mu_s
        self.m_temp[2::3] *= self.mu_s
        self.fmm = fmm.FMM(self.n, self.ncrit, self.order, mesh.coordinates * mesh.unit_length, self.m_temp)

    def compute_field(self, t=0, spin=None):
        self.m_temp[:] = spin if spin is not None else self.spin
        self.m_temp[0::3] *= self.mu_s
        self.m_temp[1::3] *= self.mu_s
        self.m_temp[2::3] *= self.mu_s

        self.field[:] = 0.0
        #self.fmm.set(self.m_temp)
        self.fmm.compute_field(self.theta, self.field)
        self.field *= 1e-7
        return self.field


# class DemagBH(Energy):
#     def __init__(self, order, ncrit, theta, name="DemagBH"):
#         self.name = name
#         assert order > 0, "Order must be 1 or higher"
#         self.order = order
#         assert ncrit >= 2, "ncrit must be greater than 1."
#         self.ncrit = ncrit
#         assert theta >= 0.0, "theta must be >= 0.0"
#         self.theta = theta

#     def setup(self, mesh, spin, mu_s, mu_s_inv):
#         super(Demag, self).setup(mesh, spin, mu_s, mu_s_inv)
#         self.n = mesh.n
#         print(mesh.coordinates)
#         self.fmm = fmm.FMM(self.n, self.ncrit, self.order, mesh.coordinates * mesh.unit_length, self.spin)
#         self.m_temp = np.zeros_like(sim.spin)

#     def compute_field(self, t=0, spin=None):
#         m = spin.copy() if spin is not None else self.spin.copy()
#         m[0::3] *= self.mu_s
#         m[1::3] *= self.mu_s
#         m[2::3] *= self.mu_s

#         self.field[:] = 0.0
#         #self.fmm.set(m)
#         self.fmm.compute_field(self.theta, self.field)
#         self.field *= 1e-7
#         return self.field


    
mesh = fidimag.common.CuboidMesh(nx=Nx, ny=Ny, nz=Nz, dx=1, dy=1, dz=1, unit_length=1e-9)
N = mesh.n
print("Nparticles = {}".format(N))

demagfmm = DemagFMM(order, ncrit, theta=theta)
demagfft = Demag()
demagfull = DemagFull()

sim = fidimag.atomistic.Sim(mesh)
sim.set_m((0, 0, 1), normalise=True)
print(f"Bohr Magnetom = {fidimag.common.constant.mu_B}")
sim.set_mu_s(fidimag.common.constant.mu_B)

sim.add(demagfft)
sim.add(demagfmm)
sim.add(demagfull)



start = time.time()
demagfft.compute_field()
end = time.time()
t = end - start
print('fft t = {}'.format(t))
print('{}\n\n'.format(demagfft.field))

start = time.time()

demagfmm.fmm.numcells()

print(f"N * bohr = {fidimag.common.constant.mu_B * N}")
demagfmm.compute_field()
err = (demagfft.field - demagfmm.field) / demagfft.field

errx = (err[0::3]*err[0::3])**0.5
erry = (err[1::3]*err[1::3])**0.5
errz = (err[2::3]*err[2::3])**0.5

print(f'errx = {np.mean(errx)}')
print(f'erry = {np.mean(erry)}')
print(f'errz = {np.mean(erry)}')


# demagfmm.fmm.checkM()

start = time.time()
demagfull.compute_field()
end = time.time()
t = end - start
print('direct t = {}, {}\n'.format(t, demagfull.field))

Exrel_err = 0
Eyrel_err = 0
Ezrel_err = 0

direct = np.zeros(3*N)
demagfmm.fmm.compute_field_exact(direct)
direct *= 1e-7

for i in range(N):
    exerr = (direct[3*i + 0] - demagfmm.field[3*i + 0])/direct[3*i + 0]
    eyerr = (direct[3*i + 1] - demagfmm.field[3*i + 1])/direct[3*i + 1]
    ezerr = (direct[3*i + 2] - demagfmm.field[3*i + 2])/direct[3*i + 2]

    Exrel_err += np.sqrt(exerr * exerr);
    Eyrel_err += np.sqrt(eyerr * eyerr);
    Ezrel_err += np.sqrt(ezerr * ezerr);
    

Exrel_err /= N
Eyrel_err /= N
Ezrel_err /= N

print(f'err = {Exrel_err}, {Eyrel_err}, {Ezrel_err}')


