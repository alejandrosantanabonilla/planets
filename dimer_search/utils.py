import ase
from ase.io import read, write
from ase.optimize.minimahopping import MinimaHopping
from tblite.ase import TBLite
from ase.optimize.minimahopping import MHPlot

def relax(atoms, tblite_params=None, mh_params=None, totalsteps=20, output_filename_prefix="relaxed"):
    """Relaxes the molecule using TBLite and MinimaHopping with options."""

    if tblite_params is None:
        tblite_params = {}

    if mh_params is None:
        mh_params = {}

    # Default TBLite Parameters
    default_tblite_params = {
        "method": "GFN2-xTB",
        "charge": None,
        "multiplicity": None,
        "accuracy": 1.0,
        "electronic_temperature": 300.0,
        "max_iterations": 250,
        "initial_guess": "sad",
        "mixer_damping": 0.4,
        "electric_field": None,
        "spin_polarization": None,
        "cache_api": True,
        "verbosity": 1,
    }
    tblite_params = {**default_tblite_params, **tblite_params}

    # Default MinimaHopping Parameters
    default_mh_params = {
        "T0": 1000.0,  # K, initial MD ‘temperature’
        "beta1": 1.1,  # temperature adjustment parameter
        "beta2": 1.1,  # temperature adjustment parameter
        "beta3": 1.0 / 1.1,  # temperature adjustment parameter
        "Ediff0": 0.5,  # eV, initial energy acceptance threshold
        "alpha1": 0.98,  # energy threshold adjustment parameter
        "alpha2": 1.0 / 0.98,  # energy threshold adjustment parameter
        "mdmin": 2,  # criteria to stop MD simulation (no. of minima)
        "logfile": "hop.log",  # text log
        "minima_threshold": 0.5,  # A, threshold for identical configs
        "timestep": 1.0,  # fs, timestep for MD simulations
        "minima_traj": "minima.traj",  # storage file for minima list
        "fmax": 0.05,  # eV/A, max force for optimizations
    }
    mh_params = {**default_mh_params, **mh_params}

    print(f"Relaxing molecule with TBLite parameters: {tblite_params}, MinimaHopping parameters: {mh_params}, totalsteps={totalsteps}")

    calculator = TBLite(**tblite_params)
    atoms.set_calculator(calculator)

    hop = MinimaHopping(atoms=atoms, **mh_params)
    hop(totalsteps=totalsteps)

    mhplot = MHPlot()
    mhplot.save_figure(f"{output_filename_prefix}_summary.png")

    traj_filename = str("minima.traj")
    traj = read(traj_filename)
    write(f"{output_filename_prefix}_minima.xyz", traj[-1], format="xyz")

    #return traj[-1]
