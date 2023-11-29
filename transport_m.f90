! 1D steady-state transport equation:
! 0 = -q*𝜕c/𝜕x + D*(𝜕^2)c/𝜕(c^2) + r*(c_r-c)
module transport_m
    use diffusion_m
    use transport_properties_heterog_m
    !use transport_properties_homog_m
    !use BCs_m
    implicit none
    save
    type, public, extends(diffusion_1D_c) :: transport_1D_c
        type(tpt_props_heterog_c) :: tpt_props_heterog
        integer(kind=4), allocatable :: conc_r_flag(:)      ! 1 if r>0
                                                            ! 0 otherwise
    contains
    ! Set
        !procedure, public :: set_spatial_discr
        !procedure, public :: set_BCs
        !procedure, public :: set_source_terms => set_source_terms_transport_1D
        !procedure, public :: set_conc_ext
        !procedure, public :: set_parameters => set_parameters_pointer
        procedure, public :: set_tpt_props_heterog_obj
        procedure, public :: set_conc_r_flag
        !procedure, public :: set_source_term_flag
        !procedure, public :: set_conc_init => set_conc_init_transport_1D
        !procedure, public :: set_time_discr
    ! Get
        !procedure, public :: get_spatial_discr
        !procedure, public :: get_props
    ! Computations
        procedure, public :: compute_trans_mat_PDE=>compute_trans_mat_tpt
        !procedure, public :: main_transport_1D
        !procedure, public :: solve_transport_1D
        !procedure, public :: solve_transport_1D_EE_Delta_t_homog
        !procedure, public :: solve_transport_1D_EE_Delta_t_heterog
        !!procedure, public :: transport_CFD_EE_1D
        !procedure, public :: transport_iteration_EE_1D
        !!procedure, public :: transport_EFD_EE_1D
        !!procedure, public :: transport_upwind_EE_1D
        !procedure, public :: solve_transport_1D_nonEE_Delta_t_homog
        !procedure, public :: solve_transport_1D_nonEE_Delta_t_heterog
        !procedure, public :: compute_trans_mat
        !procedure, public :: compute_mixing_ratios
        !procedure, public :: compute_A_lin_syst
        !procedure, public :: compute_b_mat
        !procedure, public :: compute_b_lin_syst
        !procedure, public :: compute_g
        !procedure, public :: compute_F_mat
        !procedure, public :: solve_transport_1D
        procedure, public :: mass_balance_error_ADE_stat_Dirichlet_discharge
        !procedure, public :: mass_balance_error_Dirichlet_recharge
    ! Abstract
        procedure, public :: initialise_PDE=>initialise_transport_1D
        procedure, public :: write_PDE_1D=>write_transport_1D
    end type
!*****************************************************************************************************************************
    !type, public, extends(transport_1D) :: transport_1D_cst_parameters
    !    type(cst_parameters) :: cst_parameters ! stationary 1D transport with constant parameters
    !end type
    !
    !type, public, extends(transport_1D) :: transport_1D_var_parameters
    !    type(var_parameters) :: var_parameters ! stationary 1D transport with variable parameters
    !end type
!*****************************************************************************************************************************
    !type, public, extends(transport_1D_c) :: transport_1D_transient_c ! transient subclass
    !    real(kind=8), allocatable :: conc_init(:) ! Initial concentration
    !    
    !    class(time_discr_c), pointer :: time_discr ! Time discretisation polymorphic variable
    !    !class(stab_parameters_c), pointer :: stab_parameters
    !contains
    !    procedure, public :: set_conc_init => set_conc_init_transport_1D_transient
    !    !procedure, public :: set_time_discr
    !end type
    
    !type, public, extends(transport_1D_transient) :: transport_1D_transient_cst_parameters
    !    type(cst_parameters_transient) :: cst_parameters_transient ! transient 1D transport with constant parameters
    !end type
    !
    !type, public, extends(transport_1D_transient) :: transport_1D_transient_var_parameters
    !    type(var_parameters_transient) :: var_parameters_transient ! transient 1D transport with variable parameters
    !end type
!*****************************************************************************************************************************
    !type, public, extends(transport_1D_transient_cst_parameters) :: transport_1D_spatial_discr
    !    !private
    !    integer(kind=4) :: scheme   ! Spatial discretisation scheme
    !                                ! 1: CFD
    !                                ! 2: Upwind
    !end type
    
    !type, public, extends(transport_1D_spatial_discr) :: transport_1D_time_int
    !    !private
    !    integer(kind=4) :: int_method   ! Time integration method
    !                                    ! 1: Euler explicit
    !                                    ! 2: Euler implicit
    !                                    ! 3: Crank-Nicolson
    !                                    ! 4: 𝜃-method
    !end type
!*****************************************************************************************************************************
    interface
        subroutine compute_trans_mat_tpt(this)
            import transport_1D_c
            implicit none
            class(transport_1D_c) :: this
            !real(kind=8), intent(out), optional :: T_sub(:),T_diag(:),T_super(:)
            !integer(kind=4), intent(in), optional :: k
        end subroutine
        
        !subroutine main_transport_1D(this,anal_sol)
        !! Calls subroutines that perform transport computations and writes results
        !    import transport_1D_c 
        !    implicit none
        !    class(transport_1D_c), intent(in) :: this
        !    real(kind=8), external, optional :: anal_sol
        !end subroutine
        !
        !subroutine solve_transport_1D(this)
        !    import transport_1D_c 
        !    implicit none
        !    class(transport_1D_c) :: this ! transport object
        !    !integer(kind=4), intent(in) :: Num_output ! size(Time_out)
        !    !real(kind=8), intent(in) :: Time_out(:)
        !    !real(kind=8), intent(out) :: conc_out(:,:)
        !    !real(kind=8), intent(out) :: theta
        !    !real(kind=8), external, optional :: anal_sol
        !end subroutine
        !
        !
        !
        !subroutine solve_transport_1D_EE_Delta_t_homog(this,Num_output,Time_out,conc_out,anal_sol)
        !    import transport_1D_c
        !    implicit none
        !    class(transport_1D_c) :: this
        !    integer(kind=4), intent(in) :: Num_output
        !    real(kind=8), intent(in) :: Time_out(:)
        !    real(kind=8), intent(out) :: conc_out(:,:)
        !    real(kind=8), external, optional :: anal_sol
        !end subroutine
        !
        !subroutine solve_transport_1D_EE_Delta_t_heterog(this,Num_output,Time_out,conc_out,anal_sol)
        !    import transport_1D_c
        !    implicit none
        !    class(transport_1D_c) :: this
        !    integer(kind=4), intent(in) :: Num_output
        !    real(kind=8), intent(in) :: Time_out(:)
        !    real(kind=8), intent(out) :: conc_out(:,:)
        !    real(kind=8), external, optional :: anal_sol
        !end subroutine 
        !            
        !subroutine transport_CFD_EE_1D(this,conc_old,conc_new,k)
        !    import transport_1D_c
        !    implicit none
        !    class(transport_1D_c), intent(in) :: this
        !    real(kind=8), intent(in) :: conc_old(:)
        !    !integer(kind=4), intent(in) :: conc_star_flag
        !    !integer(kind=4), intent(in) :: source_term_flag
        !    real(kind=8), intent(inout) :: conc_new(:)
        !    integer(kind=4), intent(in), optional :: k
        !end subroutine transport_CFD_EE_1D
        !
        !subroutine transport_iteration_EE_1D(this,b_mat_sub,b_mat_diag,b_mat_super,conc_old,conc_new,k,anal_sol)
        !    import transport_1D_c
        !    implicit none
        !    class(transport_1D_c), intent(in) :: this
        !    real(kind=8), intent(in) :: b_mat_sub(:),b_mat_diag(:),b_mat_super(:)
        !    real(kind=8), intent(in) :: conc_old(:)
        !    real(kind=8), intent(inout) :: conc_new(:)
        !    integer(kind=4), intent(in), optional :: k
        !    real(kind=8), external, optional :: anal_sol
        !end subroutine
        !
        !!subroutine transport_CFD_EE_1D_var_param_flux_nodes(this,conc_old,conc_new,k)
        !!    import transport_1D_c
        !!    implicit none
        !!    class(transport_1D_c), intent(in) :: this
        !!    real(kind=8), intent(in) :: conc_old(:)
        !!    real(kind=8), intent(inout) :: conc_new(:)
        !!    integer(kind=4), intent(in) :: k
        !!end subroutine transport_CFD_EE_1D_var_param_flux_nodes
        !
        !subroutine transport_EFD_EE_1D(this,conc_old,conc_new,k)
        !    import transport_1D_c
        !    implicit none
        !    class(transport_1D_c), intent(in) :: this
        !    real(kind=8), intent(in) :: conc_old(:)
        !    real(kind=8), intent(inout) :: conc_new(:)
        !    integer(kind=4), intent(in), optional :: k
        !end subroutine 
        !
        !subroutine transport_upwind_EE_1D(this,conc_old,conc_new,k)
        !    import transport_1D_c
        !    implicit none
        !    class(transport_1D_c), intent(in) :: this
        !    real(kind=8), intent(in) :: conc_old(:)
        !    real(kind=8), intent(inout) :: conc_new(:)
        !    integer(kind=4), intent(in), optional :: k
        !end subroutine 
        !
        !subroutine solve_transport_1D_nonEE_Delta_t_homog(this,theta,Num_output,Time_out,conc_out)
        !    import transport_1D_c
        !    implicit none
        !    class(transport_1D_c) :: this
        !    real(kind=8), intent(in) :: theta
        !    integer(kind=4), intent(in) :: Num_output
        !    real(kind=8), intent(in) :: Time_out(:)
        !    real(kind=8), intent(out) :: conc_out(:,:)
        !end subroutine
        !
        !subroutine solve_transport_1D_nonEE_Delta_t_heterog(this,theta,Num_output,Time_out,conc_out)
        !    import transport_1D_c
        !    implicit none
        !    class(transport_1D_c) :: this
        !    real(kind=8), intent(in) :: theta
        !    integer(kind=4), intent(in) :: Num_output
        !    real(kind=8), intent(in) :: Time_out(:)
        !    real(kind=8), intent(out) :: conc_out(:,:)
        !end subroutine
        !
        !subroutine compute_trans_mat(this,E_sub,E_diag,E_super,k)
        !    import transport_1D_c
        !    implicit none
        !    class(transport_1D_c) :: this
        !    real(kind=8), intent(out), optional :: E_sub(:),E_diag(:),E_super(:)
        !    integer(kind=4), intent(in), optional :: k
        !end subroutine
        !
        !subroutine compute_mixing_ratios(this,theta,k)
        !    import transport_1D_c
        !    implicit none
        !    class(transport_1D_c) :: this
        !    real(kind=8), intent(in), optional :: theta
        !    integer(kind=4), intent(in), optional :: k
        !end subroutine
        !
        !subroutine compute_A_lin_syst(this,theta,E_sub,E_diag,E_super,A_sub,A_diag,A_super)
        !    import transport_1D_c
        !    implicit none
        !    class(transport_1D_c), intent(in) :: this
        !    real(kind=8), intent(in) :: theta
        !    real(kind=8), intent(in) :: E_sub(:),E_diag(:),E_super(:)
        !    real(kind=8), intent(out) :: A_sub(:),A_diag(:),A_super(:)
        !end subroutine
        !
        !subroutine compute_b_mat(this,theta,E_sub,E_diag,E_super,b_mat_sub,b_mat_diag,b_mat_super)
        !    import transport_1D_c
        !    implicit none
        !    class(transport_1D_c), intent(in) :: this
        !    real(kind=8), intent(in) :: theta
        !    real(kind=8), intent(in) :: E_sub(:),E_diag(:),E_super(:)
        !    real(kind=8), intent(out) :: b_mat_sub(:),b_mat_diag(:),b_mat_super(:)
        !end subroutine 
        !
        !subroutine compute_b_lin_syst(this,b_mat_sub,b_mat_diag,b_mat_super,conc_old,b,k)
        !    import transport_1D_c
        !    implicit none
        !    class(transport_1D_c), intent(in) :: this
        !    real(kind=8), intent(in) :: b_mat_sub(:),b_mat_diag(:),b_mat_super(:),conc_old(:)
        !    real(kind=8), intent(out) :: b(:)
        !    integer(kind=4), intent(in), optional :: k
        !end subroutine
        !
        !function compute_g(this,k,conc_old,anal_sol) result(g)
        !    import transport_1D_c
        !    implicit none
        !    class(transport_1D_c), intent(in) :: this
        !    integer(kind=4), intent(in), optional :: k
        !    real(kind=8), intent(in), optional :: conc_old(:)
        !    real(kind=8), external, optional :: anal_sol
        !    real(kind=8), allocatable :: g(:)
        !end function
        !
        !function compute_F_mat(this) result(F_mat)
        !    import transport_1D_c
        !    implicit none
        !    class(transport_1D_c), intent(in) :: this
        !    real(kind=8), allocatable :: F_mat(:)
        !end function
        
        subroutine initialise_transport_1D(this)
            import transport_1D_c
            implicit none
            class(transport_1D_c) :: this
        end subroutine
        
        subroutine write_transport_1D(this,Time_out,output)
            import transport_1D_c
            implicit none
            class(transport_1D_c), intent(in) :: this
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(in) :: output(:,:)
        end subroutine
        
        function mass_balance_error_ADE_stat_Dirichlet_discharge(this) result(mass_bal_err)
            import transport_1D_c
            implicit none
            class(transport_1D_c), intent(in) :: this
            real(kind=8) :: mass_bal_err
        end function
    
    end interface
!*****************************************************************************************************************************
    contains
        !subroutine set_spatial_discr(this,spatial_discr_obj)
        !    implicit none
        !    class(transport_1D_c) :: this
        !    class(spatial_discr_c), intent(in), target :: spatial_discr_obj
        !    this%spatial_discr=>spatial_discr_obj
        !end subroutine 
        
        !subroutine set_BCs(this,BCs_obj)
        !    implicit none
        !    class(transport_1D_c) :: this
        !    class(BCs_t), intent(in) :: BCs_obj
        !    this%BCs=BCs_obj
        !end subroutine
        
        !subroutine set_conc_ext(this,conc_ext)
        !    class(transport_1D_c) :: this
        !    real(kind=8), intent(in) :: conc_ext(:)
        !    if (size(conc_ext)/=this%spatial_discr%Num_targets) error stop "Dimension error in external concentration"
        !    this%conc_ext=conc_ext
        !end subroutine 
        
        !subroutine set_parameters_pointer(this,parameters)
        !    implicit none
        !    class(transport_1D_c) :: this
        !    class(parameters_c), intent(in), target  :: parameters
        !    this%parameters=>parameters
        !end subroutine
        
        subroutine set_tpt_props_heterog_obj(this,tpt_props_heterog)
            implicit none
            class(transport_1D_c) :: this
            class(tpt_props_heterog_c), intent(in)  :: tpt_props_heterog
            this%tpt_props_heterog=tpt_props_heterog
        end subroutine
        
        subroutine set_conc_r_flag(this)
            implicit none
            class(transport_1D_c) :: this
            
            integer(kind=4) :: i
            
            allocate(this%conc_r_flag(this%spatial_discr%Num_targets))
            this%conc_r_flag=0
            do i=1,this%spatial_discr%Num_targets
                if (this%tpt_props_heterog%source_term(i)>0) then
                    this%conc_r_flag(i)=1
                end if
            end do
        end subroutine 
        
        !subroutine set_source_term_flag(this)
        !    implicit none
        !    class(transport_1D_c) :: this
        !    !integer(kind=4), intent(out) :: source_term_flag
        !    this%props%source_term_flag=1
        !    if (this%props%source_term<0 .and. this%BCs%evap==.false.) then
        !        this%props%source_term_flag=0
        !    end if
        !end subroutine 
        
        !subroutine set_conc_init_transport_1D(this,conc_init)
        !    implicit none
        !    class(transport_1D_c) :: this
        !    real(kind=8), intent(in) :: conc_init(:)
        !    
        !    select type (this)
        !    type is (transport_1D_transient_c)
        !        call this%set_conc_init(conc_init)
        !    end select
        !end subroutine 
        
        !subroutine set_conc_init_transport_1D_transient(this,conc_init)
        !    implicit none
        !    class(transport_1D_transient_c) :: this
        !    real(kind=8), intent(in) :: conc_init(:)
        !    if (this%spatial_discr%Num_targets_defined==.true.) then
        !        if (size(conc_init)/=this%spatial_discr%Num_targets) error stop "Dimension error in initial concentration"
        !    else
        !        this%spatial_discr%Num_targets=size(conc_init)
        !        this%spatial_discr%Num_targets_defined=.true.
        !    end if
        !    this%conc_init=conc_init
        !end subroutine
        
        
        
       
        
      
        
        !subroutine set_time_discr(this,time_discr_obj)
        !    implicit none
        !    class(transport_1D_c) :: this
        !    class(time_discr_c), intent(in), target :: time_discr_obj
        !    select type (this)
        !    type is (transport_1D_transient_c)
        !        this%time_discr=>time_discr_obj
        !    end select
        !end subroutine 
        
        !function get_spatial_discr(this) result(spatial_discr)
        !    implicit none
        !    class(transport_1D_c), intent(in) :: this
        !    class(spatial_discr_c), pointer :: spatial_discr
        !    spatial_discr=>this%spatial_discr
        !end function
        !
        !function get_props(this) result(props)
        !    implicit none
        !    class(transport_1D_c), intent(in) :: this
        !    class(props_c), pointer :: props
        !    props=>this%props
        !end function
end module 