      subroutine matMult_loop2(A, B, C, n1, n2, n3)

      INTEGER n1, n2, n3
      DOUBLE PRECISION A(n1,n2)
      DOUBLE PRECISION B(n2,n3)
      DOUBLE PRECISION C(n1,n3)

      INTEGER i, j, k

      do j = 1, n3
         do i = 1, n1
            C(i, j) = 0
         enddo
         do k = 1, n2
            do i = 1, n1
               C(i, j) = C(i, j) + A(i, k)*B(k, j)
            enddo
         enddo
      enddo

      RETURN

      end subroutine matMult_loop2
