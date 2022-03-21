import numpy as np
import Example

d1 = Example.Dual_Num_Auto_Diff.DUAL_NUM()
d2 = Example.Dual_Num_Auto_Diff.DUAL_NUM()
d3 = Example.Mcyldnad.cyldnad(d1, d2)

cylin = Example.Mcyldnad.cyldnad
print(cylin.__doc__)

# Specify radius (r)
d1.x_ad_ = 3
# Specify that we want dv/dr, where v is cylinder volume and r is cylinder radius
d1.xp_ad_ = np.array((1., 0.))
print("d1:", d1)

# Specify height (h)
d2.x_ad_ = 5
# Specify that we want dv/dh, where h is cylinder height
d2.xp_ad_ = np.array((0, 1.))
print("d2:", d2)

d3 = Example.Mcyldnad.cyldnad(d1, d2)
# Print computed v, dv/dr, dv/dh (thanks to dual numbers)
print("result:", d3)