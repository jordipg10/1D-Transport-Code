subroutine LU_lin_syst(A,b,x) ! Ax=b
! Solves linear system using LU decomposition
    use vectors_m
    use matrices_m, only : LU
    use metodos_sist_lin_m, only : forward_substitution,backward_substitution
    implicit none
    real(kind=8), intent(in) :: A(:,:) ! square matrix (A=LU)
    real(kind=8), intent(in) :: b(:) ! vector
    real(kind=8), intent(out) :: x(:) ! solution of linear system
    
    real(kind=8), allocatable :: L(:,:), U(:,:), y(:)
    integer(kind=4) :: n
    real(kind=8), parameter :: tol=1d-9
    
    n=size(b)
    allocate(L(n,n),U(n,n),y(n))
    call LU(A,L,U)
    call forward_substitution(L,b,y)
    call backward_substitution(U,y,x)
    if (inf_norm_vec(matmul(A,x)-b)>=tol) then
        error stop "Wrong solution in LU_lin_syst"
    end if
end subroutine