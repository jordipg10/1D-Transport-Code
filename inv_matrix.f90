! Inverse of square matrix using LU decomposition
subroutine inv_matrix(A,inv)
    use matrices_m
    use vectors_m
    use metodos_sist_lin_m
    implicit none
    real(kind=8), intent(in) :: A(:,:)
    real(kind=8), intent(out) :: inv(:,:)
    
    integer(kind=4) :: n,j,i
    real(kind=8), parameter :: epsilon=1d-6
    real(kind=8), allocatable :: id(:,:), prod_A_invA(:,:), id_col(:), inv_col(:)
    logical :: nzdiag
    
    if (size(A,1)/=size(A,2)) error stop "Matrix must be square (inv_matrix)"
    n=size(A,1)
    if (n==2) then
        inv(1,1)=A(2,2)
        inv(1,2)=-A(1,2)
        inv(2,1)=-A(2,1)
        inv(2,2)=A(1,1)
        inv=inv/det(A)
    else
        allocate(inv_col(n))
        id=id_matrix(n)
        nzdiag=.true.
        do i=1,n
            if (abs(A(i,i))<epsilon) then
                nzdiag=.false.
                exit
            end if
        end do
        if (nzdiag==.true.) then
            do j=1,n
                id_col=id(1:n,j)
                call LU_lin_syst(A,id_col,inv_col) ! LU decomposition
                inv(1:n,j)=inv_col
            end do
        else
            do j=1,n
                id_col=id(1:n,j)
                call Gauss_Jordan(A,id_col,inv_col) ! LU decomposition
                inv(1:n,j)=inv_col
            end do
        end if
        prod_A_invA=matmul(A,inv)
        do i=1,n
            if (inf_norm_vec(prod_A_invA(i,:)-id(i,:))>=epsilon) then
                error stop "Error in inverse matrix"
            end if
        end do
    end if
end subroutine