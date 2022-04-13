# Question

## 关于

```fortran
c file: dot.f 
FUNCTION dot(n, x, y) 
INTEGER n, i 
DOUBLE PRECISION dot, x(n), y(n) 
dot = 0d0 
DO 10 i=1,n 
    dot = dot + x(i) * y(i) 
10 CONTINUE 
END
```

```python
sh> python 
>>> import foo 
>>> result = foo.dot([1,2], [3, 4]) 
>>> print result 
11.0
```
