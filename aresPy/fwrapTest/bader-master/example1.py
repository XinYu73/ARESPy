from bader import *

opts = default_options()
opts.chargefile = 'CHGCAR'

ions = ions_mod.Ions_Obj(0, 0)
chgval = charge_mod.Charge_Obj([0, 0, 0])
bdr = bader_mod.Bader_Obj()

io_mod.read_charge(ions,chgval,opts)
bader_mod.bader_calc(bdr,ions,chgval,opts)
bader_mod.bader_mindist(bdr,ions,chgval)
bader_mod.bader_output(bdr,ions,chgval)


