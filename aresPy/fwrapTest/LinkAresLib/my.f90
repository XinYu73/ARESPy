subroutine matrixMult(n, a, d)
    INTEGER, intent(in) :: n
    DOUBLE PRECISION, intent(in) :: a(n,n)
    DOUBLE PRECISION, intent(inout) :: d(n)
    integer :: i,info,lda,ldb ,nrhs
    integer ,dimension(n) :: ipiv
    CALL dgesv(n,nrhs,a,ldb,ipiv,d,ldb,info)
    write(*,*) d
end subroutine matrixMult