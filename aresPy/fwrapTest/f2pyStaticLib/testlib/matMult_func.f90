      subroutine matMult_func(A, B, C, n1, n2, n3)

      INTEGER n1, n2, n3
      DOUBLE PRECISION A(n1,n2)
      DOUBLE PRECISION B(n2,n3)
      DOUBLE PRECISION C(n1,n3)

      C = matmul(A,B)

      RETURN

      end subroutine matMult_func
!end