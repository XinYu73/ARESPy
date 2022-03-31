program test
implicit none
integer,parameter :: n= 3
real(kind=8), dimension(n) :: x,b
real(kind=8), dimension(n,n) :: a,xx,yy
integer :: i,info,lda,ldb ,nrhs
integer ,dimension(n) :: ipiv

a = reshape((/1.0,45.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0/),(/n,n/))
xx = reshape((/1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0/),(/n,n/))
yy = reshape((/1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0/),(/n,n/))
b=(/2.0,2.0,3.5/)
write(*,*) 'j(1) = '
x=b
nrhs=1
lda=n
ldb=n
call dgesv(n,nrhs,a,ldb,ipiv,x,ldb,info)
write(*,*) x
end program test