subroutine compute_B_mat(this,theta,B_mat,k)
! B=(Id+(1-theta)*E) (tridiagonal, negative semi-definite)
    use PDE_transient_m
    implicit none
    
    class(PDE_1D_transient_c) :: this
    real(kind=8), intent(in) :: theta
    type(tridiag_matrix_c), intent(out) :: B_mat
    integer(kind=4), intent(in), optional :: k
    
    type(tridiag_matrix_c) :: E_mat
    integer(kind=4) :: n
    real(kind=8) :: B_norm_inf,B_norm_1
    
    n=this%spatial_discr%Num_targets
    
    call E_mat%allocate_matrix(n)
    call this%compute_E_mat(E_mat,k)
    B_mat%sub=(1d0-theta)*E_mat%sub
    B_mat%diag=1d0+(1d0-theta)*E_mat%diag
    B_mat%super=(1d0-theta)*E_mat%super
    

end subroutine 
