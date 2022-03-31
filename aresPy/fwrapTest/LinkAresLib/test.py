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
D=my.my(A,C)
print(D)
end = time()
""" 
print ('Loop1: time for','AB'+str(np.shape(AB)),'=','A'+str(np.shape(A)),'B'+str(
np.shape(B)),'is', end - beg,'s')

#-------
# Case 2
#-------
beg = time()
AB = my.matrixmult(A,B,C,2)
end = time()

print ('Loop2: time for','AB'+str(np.shape(AB)),'=','A'+str(np.shape(A)),'B'+str(np.shape(B)),'is', end - beg,'s')

#-------
# Case 3
#-------
beg = time()
AB = my.matrixmult(A,B,C,3)
end = time()

print ('matmul function: time for','AB'+str(np.shape(AB)),'=','A'+str(np.shape(A)),'B'+str(np.shape(B)),'is', end - beg,'s')
 """