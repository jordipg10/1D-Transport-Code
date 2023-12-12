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
    !print *, this%props%get_props()
    
    !E_mat%sub=this%trans_mat%sub
    !E_mat%diag=this%trans_mat%diag
    !E_mat%super=this%trans_mat%super
    !select type (time_discr=>this%time_discr)
    !type is (time_discr_homog_c)
    !    E_mat%sub=E_mat%sub*time_discr%Delta_t
    !    E_mat%diag=E_mat%diag*time_discr%Delta_t
    !    E_mat%super=E_mat%super*time_discr%Delta_t
    !type is (time_discr_heterog_c)
    !    E_mat%sub=E_mat%sub*time_discr%Delta_t(k)
    !    E_mat%diag=E_mat%diag*time_discr%Delta_t(k)
    !    E_mat%super=E_mat%super*time_discr%Delta_t(k)
    !end select
    !E_mat%diag=E_mat%diag/this%F_mat%diag
    call E_mat%allocate_matrix(n)
    call this%compute_E_mat(E_mat,k)
    B_mat%sub=(1d0-theta)*E_mat%sub
    B_mat%diag=1d0+(1d0-theta)*E_mat%diag
    B_mat%super=(1d0-theta)*E_mat%super
    
    !B_norm_inf=B_mat%compute_norm_inf()
    !B_norm_1=B_mat%compute_norm_1()
    !print *, this%props%get_props()
end subroutine 