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
    
! We compute transport arrays
    call this%allocate_trans_mat()
    call this%compute_trans_mat_PDE()
    call this%compute_source_term_PDE()
    call this%F_mat%allocate_matrix(n-this%spatial_discr%targets_flag)
    call this%compute_F_mat_PDE()
! We impose BCs
    if (this%BCs%BCs_label(1)==1 .and. this%BCs%BCs_label(2)==1) then
        call Dirichlet_BCs_PDE(this)
    else if (this%BCs%BCs_label(1)==2 .and. this%BCs%BCs_label(2)==2) then
        call Neumann_homog_BCs(this)
    else if (this%BCs%BCs_label(1)==3 .and. this%BCs%BCs_label(2)==2) then
        call Robin_Neumann_homog_BCs(this)
    else
        error stop "Boundary conditions not implemented yet"
    end if

    call this%compute_B_mat(theta,B_mat,k)
    this%mixing_ratios%sub=B_mat%sub
    this%mixing_ratios%diag=B_mat%diag
    this%mixing_ratios%super=B_mat%super
    
    select type (this)
    type is (transport_1D_transient_c)
        select type (time=>this%time_discr)
        type is (time_discr_homog_c)
            this%mixing_ratios%vector=this%tpt_props_heterog%source_term_flag*this%tpt_props_heterog%source_term*time%Delta_t/this%F_mat%diag
        type is (time_discr_heterog_c)
            this%mixing_ratios%vector=this%tpt_props_heterog%source_term_flag*this%tpt_props_heterog%source_term*time%Delta_t(k)/this%F_mat%diag
        end select
   end select
end subroutine 