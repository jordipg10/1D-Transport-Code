module transport_transient_m
    use diffusion_transient_m
    use transport_properties_heterog_m
    use transport_stab_params_m
    implicit none
    save
    type, public, extends(diffusion_1D_transient_c) :: transport_1D_transient_c ! 1D transient transport subclass
        type(tpt_props_heterog_c) :: tpt_props_heterog      ! properties
        type(stab_params_tpt_c) :: stab_params_tpt          ! stability parameters
    contains
        procedure, public :: set_conc_r_flag=>set_conc_r_flag_tpt
        procedure, public :: compute_F_mat_PDE=>compute_F_mat_tpt
        procedure, public :: compute_trans_mat_PDE=>compute_trans_mat_tpt_transient
        procedure, public :: set_stab_params_tpt
        procedure, public :: initialise_PDE=>initialise_transport_1D_transient
        procedure, public :: write_PDE_1D=>write_transport_1D_transient
        procedure, public :: set_tpt_props_heterog_obj
        procedure, public :: mass_balance_error_ADE_trans_PMF_evap
        procedure, public :: mass_balance_error_ADE_trans_Dirichlet_evap
        procedure, public :: mass_balance_error_ADE_trans_PMF_recharge
        procedure, public :: mass_balance_error_ADE_trans_Dirichlet_recharge
        procedure, public :: mass_balance_error_ADE_trans_PMF_discharge
        procedure, public :: mass_balance_error_ADE_trans_Dirichlet_discharge
        procedure, public :: check_Delta_t
    end type
    
    interface
        subroutine compute_F_mat_tpt(this)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c) :: this
        end subroutine
        
        subroutine compute_trans_mat_tpt_transient(this)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c) :: this
        end subroutine
        
        subroutine initialise_transport_1D_transient(this)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c) :: this
        end subroutine
        
        subroutine write_transport_1D_transient(this,Time_out,output)
            import transport_1D_transient_c
            import props_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this ! transport object
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(in) :: output(:,:)
        end subroutine
        
        function mass_balance_error_ADE_trans_Dirichlet_recharge(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
        
        function mass_balance_error_ADE_trans_Dirichlet_discharge(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
        
        function mass_balance_error_ADE_trans_Dirichlet_evap(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
        
        function mass_balance_error_ADE_trans_PMF_recharge(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
        
        function mass_balance_error_ADE_trans_PMF_discharge(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
        
        function mass_balance_error_ADE_trans_PMF_evap(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
    end interface
    
    contains
      
        subroutine set_stab_params_tpt(this,stab_params_tpt)
            implicit none
            class(transport_1D_transient_c) :: this
            class(stab_params_tpt_c), intent(in) :: stab_params_tpt
            this%stab_params_tpt=stab_params_tpt
        end subroutine
        
        subroutine set_tpt_props_heterog_obj(this,tpt_props_heterog)
            implicit none
            class(transport_1D_transient_c) :: this
            class(tpt_props_heterog_c), intent(in) :: tpt_props_heterog
            this%tpt_props_heterog=tpt_props_heterog
        end subroutine
        
        subroutine set_conc_r_flag_tpt(this)
            implicit none
            class(transport_1D_transient_c) :: this
            
            integer(kind=4) :: i
            
            allocate(this%conc_r_flag(this%spatial_discr%Num_targets))
            this%conc_r_flag=0
            do i=1,this%spatial_discr%Num_targets
                if (this%tpt_props_heterog%source_term(i)>0) then
                    this%conc_r_flag(i)=1
                end if
            end do
        end subroutine
        
        subroutine check_Delta_t(this)
            implicit none
            class(transport_1D_transient_c) :: this
            select type (time=>this%time_discr)
            type is (time_discr_homog_c)
                if (time%Delta_t>this%stab_params_tpt%Delta_t_crit) then
                    !print *, "Critical time step: ", this%stab_params_tpt%Delta_t_crit
                    !error stop "You must reduce time step to have stability"
                end if
            end select
        end subroutine
end module