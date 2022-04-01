subroutine my(a, b)
implicit none
integer,parameter :: n = 3
real(kind=8), dimension(n,n),intent(in) :: a
real(kind=8),dimension(n), intent(inout) :: b
integer :: info,lda,ldb,nrhs
integer ,dimension(n) :: ipiv
nrhs=1
lda=n
ldb=n
ipiv=(/0.0,0.0,0.0/)
write(*,*) b
call dgesv(n,nrhs,a,ldb,ipiv,b,ldb,info)
write(*,*) b
RETURN
end subroutine my