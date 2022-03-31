subroutine matrixMult(n, A, D)
     external dgesv
    INTEGER, intent(in) :: n
    DOUBLE PRECISION, intent(in) :: a(n,n)
    DOUBLE PRECISION, intent(inout) :: D(n)
    integer :: i,info,lda,ldb ,nrhs
    integer ,dimension(n) :: ipiv
    CALL dgesv(n,nrhs,a,ldb,ipiv,D,ldb,info)
    RETURN
end subroutine matrixMult