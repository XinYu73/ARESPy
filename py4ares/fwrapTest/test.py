import f2pytest
import add
import numpy as np

a = np.array([1, 2, 3], np.float32)
f2pytest.fmodule.fast_reverse(a, 2)
print(a)
print(add.zadd(a[0], a[1]))
print(add.zadd.__doc__)
print(add.zadd(a, a))