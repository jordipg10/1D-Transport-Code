! LU decomposition
subroutine LU(A,L,U)
    implicit none
    real(kind=8), intent(in) :: A(:,:) ! square matrix
    real(kind=8), intent(out) :: L(:,:), U(:,:)
    integer(kind=4) :: n,i,j
    real(kind=8) :: factor
    n=size(A,1)
    U=A
    L=0d0
    do j=1,n
        L(j,j)=1d0
        do i=j+1,n
            factor=U(i,j)/U(j,j)
            L(i,j)=factor
            U(i,1:n)=U(i,1:n)-factor*U(j,1:n)
        end do
    end do
    !do i=2,n
    !    do j=1,i-1
    !        if (U(i,j)/=0d0) error stop "U is not upper triangular"
    !    end do
    !end do
end subroutine LU