subroutine solve_PDE_1D_stat(this)
    ! Solves 1D steady-state PDE
    
    ! this: 1D PDE object
        
    use transport_m
    use metodos_sist_lin_m
    implicit none
    
    ! Variables
    class(PDE_1D_c) :: this

    integer(kind=4) :: n,i,icol,k,out_freq,conc_r_flag,source_term_flag,Num_output
    real(kind=8) :: Time,MBE
    real(kind=8), parameter :: tol_Thomas=1d-9
    real(kind=8), allocatable :: conc_old(:),conc_new(:)
    type(tridiag_matrix_c) :: E_mat,B_mat

    procedure(mass_balance_error_ADE_stat_Dirichlet_discharge), pointer :: p_MBE=>null()
    
    n=this%spatial_discr%Num_targets
    
    select type (this)
    type is (transport_1D_c)
            allocate(this%conc(n))
            call Thomas(this%trans_mat,-this%source_term_PDE,tol_Thomas,this%conc) ! we solve linear system Tc=-g
            if (maxval(this%tpt_props_heterog%source_term)<0d0 .and. this%BCs%evap==.false.) then ! discharge
                p_MBE=>mass_balance_error_ADE_stat_Dirichlet_discharge
            else
                error stop "Mass balance error not implemented yet"
            end if
            MBE=abs(p_MBE(this))
    end select
end subroutine 