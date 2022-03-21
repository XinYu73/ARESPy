from numpy import f2py
with open("fib2.f90") as sourcefile:
    sourcecode = sourcefile.read()
f2py.compile(sourcecode, modulename='fib4')
import fib4
print(fib4.fib.__doc__)