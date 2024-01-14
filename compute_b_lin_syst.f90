subroutine compute_b_lin_syst(this,theta,conc_old,b,k)
! A*c^(k+1)=b
! A=Id-theta*E
! B=(Id+(1-theta)*E) (tridiagonal, negative semi-definite)
! b=B*c^k+f
    use BCs_subroutines_m
    implicit none
    
    class(PDE_1D_transient_c), intent(in) :: this
    real(kind=8), intent(in) :: theta
    real(kind=8), intent(in) :: conc_old(:)
    real(kind=8), intent(out) :: b(:)
    integer(kind=4), intent(in), optional :: k
    
    integer(kind=4) :: i,n
    type(tridiag_matrix_c) :: B_mat
    
    n=size(conc_old)
    if (n/=this%spatial_discr%Num_targets) error stop "Dimension error in subroutine 'compute_b_lin_syst'"
    
    call B_mat%allocate_matrix(n)
    call this%compute_B_mat(theta,B_mat,k)
    
    b(1)=B_mat%diag(1)*conc_old(1)+B_mat%super(1)*conc_old(2)
    do i=2,n-1
        b(i)=B_mat%sub(i-1)*conc_old(i-1)+B_mat%diag(i)*conc_old(i)+B_mat%super(i)*conc_old(i+1)
    end do
    b(n)=B_mat%sub(n-1)*conc_old(n-1)+B_mat%diag(n)*conc_old(n)
    b=b+this%compute_f(k)
end subroutine
    