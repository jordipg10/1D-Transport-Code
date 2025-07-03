! Computes mixing ratios matrix with uniform time stepping
subroutine compute_mixing_ratios_Delta_t_homog(this,theta,A_mat_lumped)
    use BCs_subroutines_m
    implicit none
    
    class(PDE_1D_transient_c) :: this
    real(kind=8), intent(in) :: theta
    type(diag_matrix_c), intent(out), optional :: A_mat_lumped
    
    integer(kind=4) :: i,n
    real(kind=8) :: lambda
    real(kind=8), allocatable :: r(:)
    
    type(tridiag_matrix_c) :: E_mat
    
    
    if (theta<0d0 .or. theta>1d0) error stop "Time weighting factor must be in [0,1]"
    
    n=this%spatial_discr%Num_targets
    
    !call this%mixing_ratios%allocate_matrix(n)
    !call this%allocate_arrays_PDE_1D()
! We compute PDE arrays
    !call this%allocate_trans_mat()
    call this%compute_trans_mat_PDE()
    call this%compute_source_term_PDE()
    !call this%allocate_F_mat()
    call this%compute_F_mat_PDE()
    !call this%allocate_B_mat()
! We impostoich_mat_react_zone BCs
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
    
    if (present(A_mat_lumped)) then
        call this%compute_lumped_A_mat(A_mat_lumped)
    end if
    !print *, this%A_mat%diag
    !print *, this%B_mat%diag
    !print *, A_mat_lumped%diag
end subroutine 