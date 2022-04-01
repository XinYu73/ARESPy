program test
implicit none
integer,parameter :: n = 3
real(kind=8), dimension(n) :: x,b
real(kind=8), dimension(n,n) :: a,xx,yy
integer :: i,info,lda,ldb ,nrhs
integer ,dimension(n) :: ipiv

a = reshape((/1.0,45.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0/),(/n,n/))
xx = reshape((/1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0/),(/n,n/))
yy = reshape((/1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0/),(/n,n/))
b=(/2.0,4.1,5.5/)
x=b
nrhs=1
lda=n
ldb=n
call dgesv(n,nrhs,a,ldb,ipiv,x,ldb,info)
write(*,*) x
write(*,*) a
write(*,*) b
write(*,*) "RRRRRRRRRR"
a = reshape((/1.0,45.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0/),(/n,n/))
call my(a,b)
write(*,*) "RRRRRRRRRR"
end program test

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