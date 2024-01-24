! Computes mixing ratios matrix with uniform time stepping
subroutine compute_mixing_ratios_Delta_t_homog(this,theta)
    use BCs_subroutines_m
    implicit none
    
    class(PDE_1D_transient_c) :: this
    real(kind=8), intent(in) :: theta
    
    integer(kind=4) :: i,n
    real(kind=8) :: lambda
    real(kind=8), allocatable :: r(:)
    
    type(tridiag_matrix_c) :: E_mat
    
    n=this%spatial_discr%Num_targets
    
    !call this%mixing_ratios%allocate_matrix(n)
    call this%allocate_arrays_PDE_1D()
! We compute PDE arrays
    !call this%allocate_trans_mat()
    call this%compute_trans_mat_PDE()
    call this%compute_source_term_PDE()
    !call this%allocate_F_mat()
    call this%compute_F_mat_PDE()
    !call this%allocate_B_mat()
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
! We compute arrays for linear system
    call this%compute_E_mat(E_mat)
    call this%compute_B_mat(theta,E_mat)
    call this%compute_A_mat(theta,E_mat)
    call this%compute_f_vec()
    
    !print *, this%B_mat%diag
    !print *, this%A_mat%diag
end subroutine 