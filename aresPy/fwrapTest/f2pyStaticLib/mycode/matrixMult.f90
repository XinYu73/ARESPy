      subroutine matrixMult(A, B, C, n1, n2, n3, iCase)

      INTEGER, intent(in) :: n1, n2, n3
      INTEGER, intent(in) :: iCase
      DOUBLE PRECISION, intent(in) :: A(n1,n2)
      DOUBLE PRECISION, intent(in) :: B(n2,n3)
      DOUBLE PRECISION, intent(out) :: C(n1,n3)

      SELECT CASE (iCase)
         CASE (1)
            CALL matMult_loop1(A, B, C, n1, n2, n3)
         CASE (2)
            CALL matMult_loop2(A, B, C, n1, n2, n3)
         CASE (3)
            CALL matMult_func (A, B, C, n1, n2, n3)
      END SELECT

      RETURN

      end subroutine matrixMult
