import numpy as np
from time import *
import sys
import matrixMult

n1 = int(4)
n2 = int(2)
n3 = int(5)

A = np.random.rand(n1,n2)
B = np.random.rand(n2,n3)

#-------
# Case 1
#-------
beg = time()
AB = matrixMult.matrixmult(A,B,1)
print(AB)
end = time()

print ('Loop1: time for','AB'+str(np.shape(AB)),'=','A'+str(np.shape(A)),'B'+str(
np.shape(B)),'is', end - beg,'s')

#-------
# Case 2
#-------
beg = time()
AB = matrixMult.matrixmult(A,B,2)
end = time()

print ('Loop2: time for','AB'+str(np.shape(AB)),'=','A'+str(np.shape(A)),'B'+str(np.shape(B)),'is', end - beg,'s')

#-------
# Case 3
#-------
beg = time()
AB = matrixMult.matrixmult(A,B,3)
end = time()

print ('matmul function: time for','AB'+str(np.shape(AB)),'=','A'+str(np.shape(A)),'B'+str(np.shape(B)),'is', end - beg,'s')
