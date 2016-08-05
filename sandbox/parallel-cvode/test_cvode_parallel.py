import fidimag
import time
d = 100
t = 10
dx = dy = 2.5

dz = 5
mesh = fidimag.common.CuboidMesh(nx=int(d/dx), ny=int(d/dy), nz=int(t/dz), dx=dx, dy=dy, dz=dz, unit_length=1e-9)

def Ms_function(Ms):
    def wrapped_function(pos):
        x, y, z = pos[0], pos[1], pos[2]
    
        r = ((x-d/2.)**2 + (y-d/2.)**2)**0.5  # distance from the centre
    
        if r <= d/2:
            # Mesh point is inside the disk.
            return Ms
        else:
            # Mesh point is outside the disk.
            return 0
    return wrapped_function

def init_m(pos):
    x,y,z = pos
    x0, y0 = d/2., d/2.
    r = ((x-x0)**2 + (y-y0)**2)**0.5
    
    if r<10:
        return (0,0, 1)
    elif r<30:
        return (0,0, -1)
    elif r<60:
        return (0, 0, 1)
    else:
        return (0, 0, -1)

# FeGe material paremeters.
Ms = 3.84e5  # saturation magnetisation (A/m)
A = 8.78e-12  # exchange energy constant (J/m)
D = 1.58e-3  # Dzyaloshinkii-Moriya energy constant (J/m**2)
alpha = 1  # Gilbert damping
gamma = 2.211e5  # gyromagnetic ration (m/As)

# Create simulation object.
sim = fidimag.micro.Sim(mesh, integrator='sundials_openmp')
sim.Ms = Ms_function(Ms)
sim.alpha = alpha
sim.gamma = gamma

# Add energies.
sim.add(fidimag.micro.UniformExchange(A=A))
sim.add(fidimag.micro.DMI(D=D))
#sim.add(Demag())

# Initialise the system.
sim.set_m(init_m)
a = time.time()
sim.run_until(5e-9)
b = time.time()
print('elpased time: %g'%(b - a))
