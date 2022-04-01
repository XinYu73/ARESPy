import numpy as np
from time import *
import sys
import my

n1 = int(3)
n2 = int(3)
n3 = int(3)

A = np.random.rand(n1,n2)
B = np.random.rand(n2,n3)
C = np.random.rand(1,3)
D = np.zeros(n1)
#-------
# Case 1
#-------
print(my.__doc__)
beg = time()
print(A)
my.my(A,C)
print(C)
end = time()
