import numpy as np

import ase.units as units

def default_options():
    """
    Return degault bader options
    
    Returns
    -------
    opts : bader.options_mod.options_obj
        populated with default options
    """
    opts = options_mod.options_obj()
    # default values
    opts.out_opt = opts.out_chgcar4
    opts.in_opt = opts.in_auto
    # print options
    opts.vac_flag = False
    opts.weight_flag = False
    opts.vacval = 1E-3
    opts.print_all_atom = False
    opts.print_all_bader = False
    opts.print_sel_atom = False
    opts.print_sel_bader = False
    opts.print_sum_atom = False
    opts.print_sum_bader = False
    opts.print_bader_index = False
    opts.print_atom_index = False
    # end of print options
    opts.bader_opt = opts.bader_neargrid
    opts.quit_opt = opts.quit_known
    opts.refine_edge_itrs = -1
    opts.bader_flag = True
    opts.voronoi_flag = False
    opts.dipole_flag = False
    opts.ldos_flag = False
    opts.verbose_flag = False
    opts.badertol = 1.0e-4
    opts.stepsize = 0.0
    opts.ref_flag = False
    return opts


def atoms_to_ions(atoms, rho):
    """
    Convert ASE Atoms object to bader ions_obj object
    
    Given an ase.atoms.Atoms instance `atoms` and and
    associated charge density array `rho` (e.g. computed
    using `get_pseudo_density()` method of a calculator,
    return bader.ions_mod.ions_obj and bader.charge_mod.charge_obj
    objects suitable for performing a Bader analysis.

    Parameters
    ----------
    atoms : ase.atoms.Atoms

    rho : array_like

    Returns
    -------
    ions : bader.ions_mod.ions_obj

    chg  : bader.charge_mod.charge_mod.charge_obj
    """
    
    nions = len(atoms)
    types = atoms.get_chemical_symbols()
    niontypes = len(set(types))

    chg = charge_mod.charge_obj(rho.shape)
    chg.nrho = np.product(chg.npts)
    chg.i_npts = 1.0/chg.npts
    
    ions = ions_mod.ions_obj(nions, niontypes)
    ions.num_ion = [types.count(typ) for typ in set(types)]
    ions.atomic_num = atoms.get_atomic_numbers()
    
    ions.lattice = atoms.get_cell() / units.Bohr
    ions.dir2car = ions.lattice.T
    ions.car2dir = np.linalg.inv(ions.lattice)

    # convert from ASE units as output to get_pseudo_density() to atomic units
    chg.rho = rho * units.Bohr**3 * np.linalg.det(ions.lattice)

    for i in range(3):
        chg.lat2car[:, i] = ions.dir2car[:, i]/chg.npts[i]
    chg.car2lat = np.linalg.inv(chg.lat2car)

    # origin of lattice is at chg(1,1,1)
    chg.org_lat[:] = [1., 1., 1.]
    chg.org_car[:] = [0., 0., 0.]

    # ion positions in grid points
    ions.r_car = atoms.get_positions() / units.Bohr
    ions.r_dir = atoms.get_scaled_positions()
    for i in range(nions):
        ions.r_lat[i, :] = np.dot(chg.car2lat, ions.r_car[i, :])
        ions.r_lat[i, :] += chg.org_lat
        charge_mod.pbc_r_lat(np.asfortranarray(ions.r_lat[i, :]), chg.npts)

    # distance between neighboring points
    dlat = np.zeros(3)
    for d1 in [-1, 0, 1]:
        dlat[0] = d1
        for d2 in [-1, 0, 1]:
            dlat[1] = d2
            for d3 in [-1, 0, 1]:
                dlat[2] = d3
                dcar = np.dot(chg.lat2car, dlat)
                chg.lat_dist[d1+1,d2+1,d3+1] = np.sqrt((dcar*dcar).sum())
                if (d1 == 0) and (d2 == 0) and (d3 == 0):
                    chg.lat_i_dist[d1+1,d2+1,d3+1] = 0.0
                else:
                    chg.lat_i_dist[d1+1,d2+1,d3+1] = 1/chg.lat_dist[d1+1,d2+1,d3+1]
                        
    return ions, chg


def bader(atoms, rho, full_output=False, **kwargs):
    """
    Perform a Bader analysis on charge density `rho` and structure `atoms`.

    Additional keyword arguments can be used to specify options to Bader code
    (elements of bader.options_mod.options_obj).

    Parameters
    ----------
    atoms : ase.atoms.Atoms
       Atoms object used to define atom positions and unit cell

    rho : array_like
       Charge density array on a real-space grid

    full_output : bool
        If True, return `ions`, `chg` and `opts` in addition to `bdr`
       
    Returns
    -------
    bdr : bader.bader_mod.bader_obj

    ions : bader.ions_mod.ions_obj

    chg : bader.charge_mod.charge_obj

    opts : bader.option_mod.Option_Obj
    """

    ions, chg = atoms_to_ions(atoms, rho)
    opts = default_options()
    for (key, value) in kwargs.items():
        if not hasattr(opts, key):
            raise KeyError('Unknown Bader option {0}'.format(key))
        setattr(opts, key, value)

    bdr = bader_mod.bader_obj()        
    bader_mod.bader_calc(bdr,ions, chg, opts)
    bader_mod.bader_mindist(bdr, ions, chg)

    if full_output:
        return bdr, ions, chg, opts
    else:
        return bdr

    
