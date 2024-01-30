subroutine eigenvectors_tridiag_toeplitz_matrix(A)
    use eigenvectors_m
    implicit none
    class(tridiag_Toeplitz_matrix_c) :: A
    
    integer(kind=4) :: i,j
    real(kind=8) :: L2_norm_vj
    real(kind=8), parameter :: pi=4d0*atan(1d0)
    
    if (A%sub*A%super<=0d0) error stop "a*c must be positive"
    allocate(A%eigenvectors(A%dim,A%dim))
    do j=1,A%dim
        do i=1,A%dim
            A%eigenvectors(i,j)=(A%sub/A%super)**((i-1d0)/2d0)*sin((A%dim-j+1)*pi*i/(A%dim+1d0))
        end do
        L2_norm_vj=p_norm_vec(A%eigenvectors(:,j),2)
        A%eigenvectors(:,j)=A%eigenvectors(:,j)/L2_norm_vj
    end do
end subroutine