import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fidimag.common.fileio import DataReader

# labels on plot
mathmode = lambda txt: "$" + txt + "$"
rm = lambda txt: "\mathrm{" + txt + "}"
label = lambda name, mi: mathmode(rm(rm(name) + "\; m_" + rm(mi)))
fidi = mathmode(rm("fidimag"))
m = lambda i: label("fidimag", i)

for datafile, descr in (
        ("sim_sundials_J0", "w/o Jacobian"),
        ("sim_sundials_diag", "Jacobian approx."),
        ("sim_sundials_J1", "with Jacobian")):
    data = DataReader(datafile + ".txt")
    total_wall_time = (data["real_time"][-1] - data["real_time"][0]) / 60.0

    fig, axes = plt.subplots(2, figsize=(10, 8), sharex=True)

    # magnetisation dynamics
    axes[0].plot(data["time"] * 1e9, data["m_x"], "bx", label=m("x"))
    axes[0].plot(data["time"] * 1e9, data["m_y"], "gx", label=m("y"))
    axes[0].plot(data["time"] * 1e9, data["m_z"], "rx", label=m("z"))
    axes[0].set_ylim((-0.3, 1.1))
    axes[0].set_ylabel("unit magnetisation (1)")
    axes[0].legend()
    axes[0].set_title("{}, wall time {:.2} minutes".format(descr, total_wall_time))

    # number of RHS evaluations
    axes[1].plot(data["time"] * 1e9, data["rhs_evals"], label=fidi)
    axes[1].legend(loc=0)
    axes[1].set_ylabel("# of RHS evaluations")
    axes[1].set_xlabel("time (ns)")
    axes[1].set_xlim((0, 0.3))

    fig.tight_layout()
    fig.savefig(datafile + ".png")
    plt.close(fig)
