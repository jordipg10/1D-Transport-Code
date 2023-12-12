! Computes mixing ratios matrix
subroutine compute_mixing_ratios(this,theta,k)
    use BCs_subroutines_m
    implicit none
    
    class(diffusion_1D_transient_c) :: this
    real(kind=8), intent(in) :: theta
    integer(kind=4), intent(in), optional :: k
    
    integer(kind=4) :: i,n
    real(kind=8) :: lambda
    real(kind=8), allocatable :: r(:)
    
    type(tridiag_matrix_c) :: B_mat
    
    n=this%spatial_discr%Num_targets
    
    call this%mixing_ratios%allocate_matrix(n)

    call this%compute_B_mat(theta,B_mat,k)
    this%mixing_ratios%sub=B_mat%sub
    this%mixing_ratios%diag=B_mat%diag
    this%mixing_ratios%super=B_mat%super
    
    select type (time=>this%time_discr)
    type is (time_discr_homog_c)
        this%mixing_ratios%vector=this%mixing_ratios%vector*time%Delta_t/this%F_mat%diag
    type is (time_discr_heterog_c)
        this%mixing_ratios%vector=this%mixing_ratios%vector*time%Delta_t(k)/this%F_mat%diag
    end select
   
end subroutine 