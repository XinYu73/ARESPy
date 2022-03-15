from bader import *

opts = default_options()
opts.chargefile = 'CHGCAR'

ions = ions_mod.ions_obj(0, 0)
chgval = charge_mod.charge_obj([0, 0, 0])
bdr = bader_mod.bader_obj()

io_mod.read_charge(ions,chgval,opts)
bader_mod.bader_calc(bdr,ions,chgval,opts)
bader_mod.bader_mindist(bdr,ions,chgval)
bader_mod.bader_output(bdr,ions,chgval)


