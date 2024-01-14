! Calls subroutines that perform 1D PDE computations
subroutine solve_PDE_1D(this,Time_out,output)
    use BCs_subroutines_m
    implicit none
! Variables
    class(PDE_1D_c) :: this ! 1D PDE object
    real(kind=8), intent(in) :: Time_out(:)
    real(kind=8), intent(out) :: output(:,:)
    
    integer(kind=4) :: k,n,j,i
    real(kind=8) :: sum,B_norm_inf
    real(kind=8), allocatable :: Delta_r(:)
    real(kind=8), parameter :: tol=1d-12
    type(tridiag_matrix_c) :: B_mat,A_mat
    
    n=this%spatial_discr%Num_targets
    
! We compute arrays
    call this%allocate_trans_mat()
    call this%compute_trans_mat_PDE()
    call this%compute_source_term_PDE()
    select type (this)
    class is (diffusion_1D_c)
        if (this%sol_method==1) then
            call this%solve_PDE_1D_stat()
        end if
    class is (PDE_1D_transient_c)
        call this%allocate_F_mat()
        call this%compute_F_mat_PDE()
        if (this%sol_method==1) then
            select type (time_discr=>this%time_discr)
            type is (time_discr_homog_c)
                if (time_discr%int_method==1) then
                    call this%solve_PDE_EE_Delta_t_homog(Time_out,output)
                else if (time_discr%int_method==2 .or. time_discr%int_method==3) then
                    call this%solve_PDE_EI_Delta_t_homog(1d0,Time_out,output)
                else if (time_discr%int_method==4) then
                    call this%solve_PDE_EI_Delta_t_homog(5d-1,Time_out,output)
                else if (time_discr%int_method==5) then
                    call this%solve_PDE_RKF45(time_discr%Delta_t,tol)
                end if
            type is (time_discr_heterog_c)
                if (time_discr%int_method==1) then
                    call this%solve_PDE_EE_Delta_t_heterog(Time_out,output)
                end if
            end select
        else if (this%sol_method==2) then
            select type (this)
            type is (diffusion_1D_transient_c)
                call this%compute_A_mat_ODE(A_mat)
                call A_mat%compute_eigenvalues()
                call A_mat%compute_eigenvectors()
                if (minval(A_mat%eigenvalues)<=0d0) then
                    error stop "Eigenvalues are not positive"
                end if
                allocate(A_mat%eigenvectors(n,n))
                call A_mat%check_eigenvectors_tridiag_sym_matrix(tol)
            end select
        end if
    end select
end subroutine