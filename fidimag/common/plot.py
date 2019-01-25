import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
try:
    from collections.abc import Iterable
except:
    from collections import Iterable
import numpy as np

def plot_m(sim, component='x', filename=None, figsize=(10, 5),
           z=0.0, cbar=True, ncbarticks=5,
           cmap='RdBu', bgcolor='w', scale_by_Ms=False, axis_units='nm',
           axes_pad=0.25):
    """
    Plotting function for the magnetisation.
    Inputs
    ------
    sim: simulation object of type Simulation or NormalModeSimulation

    component: str
        'x', 'y', 'z', 'all' or 'angle'

    filename: str (None)
        File to save the figure to - does not save if not specified.

    figsize: 2-tuple
        Matplotlib figure size specifiers in inches.
        If not 'all', the x-size is 1/3 of what is specified,
        such that consistency in size is maintained between
        'all' and other plots.

    extent: None, float, or 4-tuple
        If None, the bounds are calculated from the
        finite difference mesh.

        If a single number, the calculated extent is
        multiplied by that - e.g. if extent=1.1, if the
        mesh bounds are [-50, 50, -50, 50] then the new
        bounds are [-55, 55, -55, 55].

        If specified, directly, must be in
        the format:

        [-ve x bounds of plot, +ve x-bounds of plot,
         -ve y bounds of plot, +ve y-bounds of plot]

    z: float
        Height at which to sample the field.

    cbar: Boolean
        Plot a colorbar, or not...

    ncbarticks:
        Number of values listed on colorbar axis.
        Ignored if no colorbar.

    cmap:
        Matplotlib colormap name. For magnetisation,
        divergent colormaps, like RdBu, tend to work
        well.

        For spin angles, a cyclic map like 'hsv' works
        better.
        See the full list here:

        https://matplotlib.org/examples/color/colormaps_reference.html

    bcolor: str
        Color specifier for background. Areas outside of
        the mesh are set to this color. 'w' for white,
        'k' for black, or use a hexadecimal color code,
        or a tuple of RGB values.

    axis_units: str
        Units for the x and y axis labels, in case
        nm does not make sense.
        
    axes_pad: float
        Adjust the padding between the plots - necessary
        if the 
    """
    mesh = sim.mesh
    
    if not (mesh.z0 <= z <= mesh.z0 + mesh.Lz):
        raise ValueError("Z value outside of mesh!")
        
    valid_z = np.array(sorted(list(set(mesh.coordinates[:, 2]))))
    layer = np.argmin(np.abs(valid_z - z))
    n_layer = mesh.nz
    
    m = sim.spin.copy()
    # Number of spins in a layer
    n_layer = int(mesh.nx) * int(mesh.ny)
    m.shape = (-1, 3)
    mx = m[:, 0][layer*n_layer:(layer+1)*n_layer]
    my = m[:, 1][layer*n_layer:(layer+1)*n_layer]
    mz = m[:, 2][layer*n_layer:(layer+1)*n_layer]
    Ms = sim.Ms.copy()[layer*n_layer:(layer+1)*n_layer]

    mx[Ms == 0.0] = np.nan
    my[Ms == 0.0] = np.nan
    mz[Ms == 0.0] = np.nan

    mx = mx.reshape(mesh.ny, mesh.nx)
    my = my.reshape(mesh.ny, mesh.nx)
    mz = mz.reshape(mesh.ny, mesh.nx)
    Ms = Ms.reshape(mesh.ny, mesh.nx)
    

    
    extent = [mesh.x0, mesh.x0+mesh.Lx,
              mesh.y0, mesh.y0+mesh.Ly]

    components = ['x', 'y', 'z']

    if scale_by_Ms is True:
        mx *= Ms
        my *= Ms
        mz *= Ms

    m = [mx, my, mz]

    if component in ['x', 'y', 'z', 'angle']:
        fig = plt.figure(figsize=(figsize[0]/3, figsize[1]))
        # Have to use ImageGrid in order to get the Colorbar
        # to scale in size with the subplots!
        grid = ImageGrid(fig, 111,
                         nrows_ncols=(1, 1),
                         axes_pad=0.15,
                         share_all=True,
                         cbar_location="right",
                         cbar_mode="single",
                         cbar_size="7%",
                         cbar_pad=0.15,
                         )
        ax = grid[0]
        ax.set_xlabel('$x$ ({})'.format(axis_units))
        ax.set_ylabel('$y$ ({})'.format(axis_units))

        # Note: By default, imshow plots like a matrix rather
        # than a Cartesian axis, and so below we have to set
        # origin = 'lower' everywhere.

        if component is 'angle':
            theta = np.arctan2(my, mx)
            with np.testing.suppress_warnings() as sup:
                sup.filter(RuntimeWarning)
                theta[theta < 0] += 2*np.pi

            cmap_edited = plt.get_cmap(cmap)
            # Edit the colormap and set bad values to bgcolor
            cmap_edited.set_bad(color=bgcolor, alpha=1.0)
            # Here we set the bounds between 0 and 2*pi for the angle,
            # though there's no reason why it couldn't be -pi and pi
            # really.
            im = ax.imshow(theta, origin='lower',
                           extent=extent, vmin=0,
                           vmax=2*np.pi, cmap=cmap_edited)
            ax.set_title('$xy$ angle')

        else:
            cmap_edited = plt.get_cmap(cmap)
            cmap_edited.set_bad(color=bgcolor, alpha=1.0)
            if scale_by_Ms is True:
                vmin = -np.max(Ms)
                vmax = np.max(Ms)
            else:
                vmin = -1.0
                vmax = 1.0

            im = ax.imshow(m[components.index(component)], origin='lower',
                           extent=extent, vmin=vmin, vmax=vmax,
                           cmap=cmap_edited)
            if scale_by_Ms is True:
                ax.set_title('$M_{}$'.format(component))
            else:
                ax.set_title('$m_{}$'.format(component))

    elif component is 'all':
        fig = plt.figure(figsize=figsize)
        grid = ImageGrid(fig, 111,
                         nrows_ncols=(1, 3),
                         axes_pad=axes_pad,
                         share_all=True,
                         cbar_location="right",
                         cbar_mode="single",
                         cbar_size="7%",
                         cbar_pad=0.15,
                         )

        if scale_by_Ms is True:
            vmin = -np.max(Ms)
            vmax = np.max(Ms)
        else:
            vmin = -1.0
            vmax = 1.0

        for ax, comp, label in zip(grid, [mx, my, mz], ['x', 'y', 'z']):
            im = ax.imshow(comp, origin='lower', extent=extent,
                           vmin=vmin, vmax=vmax, cmap=cmap)
            if scale_by_Ms is True:
                ax.set_title('$M_{}$'.format(label))
            else:
                ax.set_title('$m_{}$'.format(label))
                
            ax.set_xlabel('$x$ ({})'.format(axis_units))
            ax.set_ylabel('$y$ ({})'.format(axis_units))

    else:
        raise ValueError("Component is not valid")

    if cbar is True:
        if component is 'angle':
            # Some special handling to print \pi
            # rather than the numbers!
            cbar = ax.cax.colorbar(im,
                                   ticks=np.linspace(0, 2*np.pi, ncbarticks))

            cbarlabels = ['${:.1f} \pi$'.format(x/(np.pi))
                          if x != 0.0 else '0.0'
                          for x in np.linspace(0, 2*np.pi, ncbarticks)]
            cbar.ax.set_yticklabels(cbarlabels)

        else:
            cbar = ax.cax.colorbar(im,
                                   ticks=np.linspace(vmin, vmax, ncbarticks))

        if scale_by_Ms:
            cbar.ax.set_ylabel('A / m', rotation=0)
        ax.cax.toggle_label(True)

    if filename:
        fig.savefig(filename, dpi=1)
