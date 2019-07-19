Version 3.0 Alpha
------------------
* Changes to the helper functions init_scalar and init_vector (fidimag/common/helper.py)
  which allow users to pass additional parameters. These are then used within the sim
  classes to allow 
  
  For example:

  ```
  
  mesh = CuboidMesh(nx=10, ny=1, nz=1, unit_length=1e-9)
  sim = Sim(mesh, Ms)

  def init_domain_wall(pos, domain_centre)
      x, y, z = pos

      if x < domain_centre:
          return (0, 0, 1)

      else:
          return (0, 0, -1)
   
  # Place a domain wall at 5nm
  sim.set_m(init_domain_wall, 5)
  # Place a domain wall at 3nm
  sim.set_m(init_domain_wall, 3)

  ```

* Setting currents is now more general, and is standardised across the simulation types.
  This allows us to use more general functions for setting the current.
  Previously, the current function was set as:
  ```
  sim(mesh, driver='llg_stt')
  sim.driver.jx = 1e14 # A / m^2
  sim.driver.update_j_fun = lambda t: np.sin(t)
  ```
  with the actual current used being multiplicative:

  $ jx * sin(t) $

  For the current-perpendicular to the plane STT ('llg_stt_cpp') driver
   we would now change this to 

  ```
  sim.driver(mesh, driver='llg_stt_cpp')
  sim.driver.j_function = lambda t: 1e14 * np.sin(t)
  ```
  and for the standard STT driver:

  ```
  sim.drive(mesh, driver='llg_stt')
  sim.driver.jz_function = lambda t: 1e14 * np.sin(t)
  # Can also set:
  # sim.driver.jx_function = ...
  # sim.driver.jy_function = ...

* Similarly, the TimeZeeman interaction is also no longer multiplicative;
  you can have an arbitrary function of the form:
 
  def time_function(pos, t):
      x, y, z = pos
      # some calculation of Bx(pos, t), By(pos, t), Bz(pos, t)
      return (Bx, By, Bz)
  zee = TimeZeeman(np.array([0, 0, 0]), time_function)
  sim.add(zee)

* You can now remove energy classes from the simulation.

  This can be useful in cases where you have an interaction
  but no longer need to calculate it because the simulation has 
  reached a certain point, e.g. an applied field has been turned off.
 
  In the data table, the entries corresponding to a given interaction 
  will be zero everywhere once the interaction is removed.
 
  
  For example:

  ```
  sim.add(Zeeman((0, 0, B), name='Zeeman'))
  
  sim.run_until(1e-9)
  sim.remove('Zeeman')
  ```
