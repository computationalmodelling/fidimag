from ipywidgets import *
import textwrap


def GUI():
    try:
        a = _widget()
        return a.GUI
    except:
        print("This must be called within a Jupyter Notebook.\n"
              "Requires IPython >= 4.2 and ipywidgets >= 5")
    

class _widget:

    def __init__(self):
        # Trying to keep MVC logic:
        # Functions provide logic to handle button presses, changing values
        # Only on click of "Generate code" style
        # These functions represent Control
        # Todo: Dictionary as data model?

        self.dict = {}

        def on_dimension_change(b):
            if self.property_dim.value == 1:
                self.property_ny.disabled = True
                self.property_nz.disabled = True
                self.property_dy.disabled = True
                self.property_dz.disabled = True
                self.property_periodicity.options = {'x': 1}
            elif self.property_dim.value == 2:
                self.property_ny.disabled = False
                self.property_nz.disabled = True
                self.property_dy.disabled = False
                self.property_dz.disabled = True
                self.property_periodicity.options = {'x': 1, 'y': 2}
            elif self.property_dim.value == 3:
                self.property_ny.disabled = False
                self.property_nz.disabled = False
                self.property_dy.disabled = False
                self.property_dz.disabled = False
                self.property_periodicity.options = {'x': 1, 'y': 2, 'z': 3}

        def on_change_anis(c):
            if not self.property_anis_enabled.value:
                self.property_K.disabled = True
                self.property_anis_axis.disabled = True
            else:
                self.property_K.disabled = False
                self.property_anis_axis.disabled = False

        def on_change_exch(c):
            if not self.property_exch_enabled.value:
                self.property_exchconst.disabled = True
            else:
                self.property_exchconst.disabled = False

        def on_change_DMI(c):
            if not self.property_DMI_enabled.value:
                self.property_DMI_const.disabled = True
            else:
                self.property_DMI_const.disabled = False

        def on_select_init_mag(c):
            # Not implementing yet but left here so don't have to rewrite
            if dim.value == 1:
                dim.options = ['Constant', 'Domain Wall']
            if dim.value == 2 or dim.value == 3:
                dim.options = ['Constant', 'Domain Wall', 'Vortex', 'Skyrmion']

        def on_change_simtype(c):
            if self.property_simtype.value == 'Relax':
                self.property_rununtil.disabled = True
            else:
                self.property_rununtil.disabled = False

        def normpress(c):
            r = np.sqrt(self.property_mx.value**2 + self.property_my.value**2 +
                        self.property_mz.value**2)
            self.property_mx.value /= r
            self.property_my.value /= r
            self.property_mz.value /= r

        def on_click_getcode(c):
            self.code_output.value = "import fidimag\n"

        def update_dictionary():
            self.dict = {}
            internal_objects = dir(self)
            for item in internal_objects:
                if "property" in item:
                    var = getattr(self, item)
                    self.dict[item[9:]] = var.value

        def get_code(c):
            update_dictionary()
            code = 'import fidimag'
            code += assemble_mesh_code()
            code += assemble_interactions_code()
            code += assemble_properties_code()
            self.code_output.value = code

        def assemble_mesh_code():
            update_dictionary()
            dim = self.dict['dim']

            periodic = self.dict['periodicity']
            px, py, pz = False, False, False
            if 1 in periodic or periodic == 1:
                px = True
            if 2 in periodic:
                py = True
            if 3 in periodic:
                pz = True

            if dim == 1:
                nx = self.dict['nx']
                ny = 1
                nz = 1
            elif dim == 2:
                nx = self.dict['nx']
                ny = self.dict['ny']
                nz = 1
            elif dim == 3:
                nx = self.dict['nx']
                ny = self.dict['ny']
                nz = self.dict['nz']

            code = textwrap.dedent("""\

                   mesh = fidimag.common.CuboidMesh(dx={}, dy={}, dz={},
                                                    nx={}, ny={}, nz={},
                                                    periodicity=({},{},{})
                                                    )

                   sim = fidimag.micro.Sim(mesh)

                   """)

            return code.format(self.dict['dx'],
                               self.dict['dy'],
                               self.dict['dz'],
                               nx, ny, nz,
                               px, py, pz
                               )

        def assemble_interactions_code():
            axis = self.dict['anis_axis']
            ax_x, ax_y, ax_z = 0, 0, 0
            if 1 in axis or axis == 1:
                ax_x = 1
            if 2 in axis or axis == 2:
                ax_y = 1
            if 3 in axis or axis == 3:
                ax_z = 1
            code = 'sim.Ms = {}'.format(self.dict['Ms'])
            if self.dict['exch_enabled']:
                code += textwrap.dedent("""\

                    exchange = fidimag.micro.UniformExchange({})
                    sim.add(exchange)

                    """).format(self.dict['exchconst'])

            if self.dict['demag_enabled']:
                code += textwrap.dedent("""\

                    demag = fidimag.micro.Demag()
                    sim.add(demag)

                    """)
            if self.dict['anis_enabled']:
                code += textwrap.dedent("""\

                    anis = fidimag.micro.UniaxialAnisotropy(Ku={},
                                                            axis=({},{},{}))
                    sim.add(anis)

                    """).format(self.dict['K'],
                                ax_x, ax_y, ax_z)
            return code

        def assemble_properties_code():
            code = 'sim.do_precession = {}'.format(
                self.dict['precession_enabled'])
            if self.dict['simtype'] == 'Relax':
                code += textwrap.dedent("""\

                    sim.relax(dt={}, stopping_dmdt={}, max_steps={},
                              save_m_steps={},
                              save_vtk_steps={})

                    """).format(self.dict['dt'],
                                self.dict['stoppingdmdt'],
                                self.dict['maxsteps'],
                                self.dict['saveevery'],
                                self.dict['saveevery'])
            elif self.dict['simtype'] == 'Run until':
                code += textwrap.dedent("""\

                    sim.run_until({})
                    """).format(self.dict['rununtil'])
            return code

        # View:

        # Page 1 : Mesh Shape

        self.label_dim = Label("Dimensions")
        self.property_dim = IntSlider(value=1, min=1, max=3, step=1)
        self.label_lengths = Label("Mesh Lengths")
        self.property_nx = IntText("1", description="nx")
        self.property_ny = IntText("1", description="ny", disabled=True)
        self.property_nz = IntText("1", description="nz", disabled=True)
        self.property_scale = FloatText(description="scale (m)", value=1e-9)
        self.label_discretisation = Label("Discretisation")
        self.property_dx = FloatText("1", description="dx")
        self.property_dy = FloatText("1", description="dy", disabled=True)
        self.property_dz = FloatText("1", description="dz", disabled=True)

        self.property_dim.observe(on_dimension_change)
        self.page1 = widgets.Box((self.label_dim, self.property_dim,
                                  self.label_lengths, self.property_nx,
                                  self.property_ny, self.property_nz,
                                  self.label_discretisation, self.property_dx,
                                  self.property_dy, self.property_dz,
                                  ))

        # Page 2 : Interactions

        self.Mslab = Label("Saturation Magnetization (A/m)")
        self.property_Ms = FloatText(value=8e5, readout_format='.5e')
        self.Msbox = Box([self.Mslab, self.property_Ms])
        self.exchlab = Label("Exchange (J/m)")
        self.property_exch_enabled = Checkbox(
            value=True, description="Enabled")
        self.property_exchconst = FloatText(value=13e-12, readout_format='.5e')
        self.property_exch_enabled.observe(on_change_exch)
        self.exchbox = Box(
            [HBox([self.exchlab, self.property_exch_enabled]),
             self.property_exchconst])
        self.demaglab = Label("Demagnetisation")
        self.property_demag_enabled = Checkbox(
            value=False, description="Enabled")
        self.demagbox = HBox([self.demaglab, self.property_demag_enabled])
        self.anislab = Label("Uniaxial Anisotropy")
        self.property_anis_enabled = Checkbox(
            value=False, description="Enabled")
        self.property_K = FloatText(value=1.3e-11, disabled=True)
        self.property_anis_axis = SelectMultiple(description="Axis",
                                                 options={
                                                     'x': 1, 'y': 2, 'z': 3},
                                                 disabled=True)
        self.property_anis_enabled.observe(on_change_anis)
        self.anisbox = Box([HBox([self.anislab, self.property_anis_enabled]),
                            self.property_K, self.property_anis_axis])
        self.DMIlab = Label("DMI")
        self.property_DMI_enabled = Checkbox(
            value=False, description="Enabled")
        self.property_DMI_const = FloatText(value=1e-3, disabled=True)
        self.property_DMI_enabled.observe(on_change_DMI)
        self.DMIbox = Box([HBox([self.DMIlab, self.property_DMI_enabled]),
                           self.property_DMI_const])
        self.page2 = Box([self.Msbox, self.exchbox, self.demagbox,
                          self.anisbox, self.DMIbox])

        # Page 3 - Initial Magnetisation
        # Put this in later after  are working
        # self.init_mag = Select(options=['Constant', 'Domain Wall']) #
        # 'Vortex', 'Skyrmion', for 2D and 3D.

        # Page 3 : Simulation Parameters
        self.label_precession = Label("Allow Precession")
        self.property_precession_enabled = Checkbox(value=True)
        self.label_dt = Label("dt (s)")
        self.property_dt = FloatText(value=1e-9)

        # self.init_mag.observe
        self.m0 = Box([FloatText(description="m0_x", min=-1, max=1, value=1),
                       FloatText(description="m0_y", min=-1, max=1, value=0),
                       FloatText(description="m0_z", min=-1, max=1, value=0)])

        self.property_stoppingdmdt = FloatText(
            description="Stopping dm/dt", value=0.01)
        self.property_saveevery = IntText(
            description="Save after every .. steps", value=100)
        self.property_maxsteps = IntText(description="Maximum number of steps")
        self.property_periodicity = SelectMultiple(description="Periodicity",
                                                   options={"x": 1})
        self.property_mx = FloatText(value=0.0)
        self.property_my = FloatText(value=0.0)
        self.property_mz = FloatText(value=0.0)
        self.norm = Button(description="Normalise m0")
        self.norm.on_click(normpress)
        self.label_simtype = Label("Simulation Type")
        self.property_simtype = RadioButtons(options=['Relax', 'Run until'])
        self.property_simtype.observe(on_change_simtype)
        self.property_rununtil = FloatText(value=10e-9, disabled=True)
        self.get_code_button = Button(description="Get Code!")
        self.page3 = Box([HBox([self.label_precession,
                                self.property_precession_enabled]),
                          HBox([self.label_dt, self.property_dt]),
                          self.property_stoppingdmdt,
                          self.property_saveevery, self.property_maxsteps,
                          self.property_periodicity,
                          self.label_simtype, self.property_simtype,
                          self.property_rununtil])

        self.code_output = Textarea(value="")
        self.code_output.width = '600px'
        self.code_output.height = '400px'

        self.page4 = Box([self.get_code_button, self.code_output])

        self.get_code_button.on_click(get_code)

        # Header, assemble final layout

        self.maintitle = Label("FIDIMAG")
        self.maintitle.layout.padding = '10px'
        self.tabs = Accordion([self.page1, self.page2, self.page3, self.page4])
        self.tabs.set_title(0, 'Mesh')
        self.tabs.set_title(1, 'Interactions')
        self.tabs.set_title(2, 'Simulation Parameters')
        self.tabs.set_title(3, 'Generate Code')
        self.GUI = Box((self.maintitle, self.tabs))

    def get_mesh_code(self):
        pass

    def print_dir(self):
        print(dir(self))
