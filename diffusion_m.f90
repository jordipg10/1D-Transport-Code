! Diffusion equation:
module diffusion_m
    use PDE_m
    !use PDE_bis_m
    !use time_discr_m
    !use matrices_m
    use diff_props_heterog_m
    implicit none
    save
    type, public, extends(PDE_1D_c) :: diffusion_1D_c ! 1D diffusion equation class
        real(kind=8), allocatable :: conc(:) ! concentrations (c)
        real(kind=8), allocatable :: conc_ext(:) ! (c_e)
        type(diff_props_heterog_c) :: diff_props_heterog
    contains
    ! Set
        
        !procedure, public :: set_spatial_discr
        !procedure, public :: set_BCs
        !procedure, public :: set_source_terms => set_source_terms_transport_1D
        procedure, public :: set_conc_ext
        !procedure, public :: set_parameters => set_parameters_pointer
        !procedure, public :: set_props
        !procedure, public :: set_source_term_flag
        !procedure, public :: set_conc_init => set_conc_init_transport_1D
        !procedure, public :: set_time_discr
        procedure, public :: set_diff_props_heterog
    ! Get
        !procedure, public :: get_spatial_discr
        !procedure, public :: get_props
    ! Computations
        !procedure, public :: compute_trans_mat_diff_bis
        !procedure, public :: compute_trans_mat_tpt
        !procedure, public :: compute_trans_mat
        procedure, public :: compute_trans_mat_PDE=>compute_trans_mat_diff
        !procedure, public :: compute_source_term_PDE=>compute_g
        !procedure, public :: compute_E_mat_diff
        !procedure, public :: compute_E_mat_tpt
        !procedure, public :: compute_E_mat=>compute_E_mat_diff
        !procedure, public :: iteration_EE_1D
        procedure, public :: initialise_PDE=>initialise_diffusion_1D
        !procedure, public :: solve_PDE_1D=>solve_diffusion
        !procedure, public :: main_PDE=>main_diffusion
        !procedure, public :: solve_transport_1D_EE_Delta_t_heterog
        !procedure, public :: transport_CFD_EE_1D
        !procedure, public :: transport_EFD_EE_1D
        !procedure, public :: transport_upwind_EE_1D
        !procedure, public :: solve_transport_1D_nonEE_Delta_t_homog
        !procedure, public :: solve_transport_1D_nonEE_Delta_t_heterog
        procedure, public :: write_PDE_1D=>write_diffusion_1D
        !procedure, public :: compute_A_lin_syst
        !procedure, public :: compute_b_lin_syst
        !procedure, public :: compute_g
        !procedure, public :: compute_A_mat_ODE
        !procedure, public :: compute_b_ODE        
        !procedure, public :: solve_diffusion_1D
        !procedure, public :: prod_total_sym_mat
    end type
!*****************************************************************************************************************************
    interface
        
        
        !subroutine solve_diffusion_1D(this)
        !    import diffusion_1D_c
        !    class(diffusion_1D_c) :: this
        !    !real(kind=8), intent(in) :: Time_out(:)
        !    !real(kind=8), intent(out) :: output(:,:)
        !    !integer(kind=4), intent(in), optional :: opcion
        !    !real(kind=8), external, optional :: anal_sol
        !end subroutine
        
        subroutine initialise_diffusion_1D(this)
            import diffusion_1D_c
            class(diffusion_1D_c) :: this
        end subroutine
        
        subroutine compute_trans_mat_diff(this)
            import diffusion_1D_c
            implicit none
            class(diffusion_1D_c) :: this
        end subroutine
        
      
        
        !subroutine compute_E_mat_diff(this,E_sub,E_diag,E_super,k)
        !    import diffusion_c
        !    implicit none
        !    class(diffusion_c) :: this
        !    real(kind=8), intent(out) :: E_sub(:),E_diag(:),E_super(:)
        !    integer(kind=4), intent(in), optional :: k
        !end subroutine
        
       subroutine write_diffusion_1D(this,Time_out,output)
            import diffusion_1D_c
            import props_c
            class(diffusion_1D_c), intent(in) :: this
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(in) :: output(:,:)
            !class(props_c), intent(in), optional :: props
            !integer(kind=4), intent(in), optional :: opcion
        end subroutine
      
      
        
        
        
    end interface
!*****************************************************************************************************************************
    contains
        
        
        !subroutine set_spatial_discr(this,spatial_discr_obj)
        !    implicit none
        !    class(diffusion_c) :: this
        !    class(spatial_discr_c), intent(in), target :: spatial_discr_obj
        !    this%spatial_discr=>spatial_discr_obj
        !end subroutine 
        !
        !subroutine set_BCs(this,BCs_obj)
        !    implicit none
        !    class(diffusion_c) :: this
        !    class(BCs_t), intent(in) :: BCs_obj
        !    this%BCs=BCs_obj
        !end subroutine
        
        subroutine set_conc_ext(this,conc_ext)
            class(diffusion_1D_c) :: this
            real(kind=8), intent(in) :: conc_ext(:)
            if (size(conc_ext)/=this%spatial_discr%Num_targets) error stop "Dimension error in external concentration"
            this%conc_ext=conc_ext
        end subroutine 
        
        !subroutine set_parameters_pointer(this,parameters)
        !    implicit none
        !    class(transport_1D_c) :: this
        !    class(parameters_c), intent(in), target  :: parameters
        !    this%parameters=>parameters
        !end subroutine
        
        !subroutine set_props(this,props_obj)
        !    implicit none
        !    class(diffusion_c) :: this
        !    class(props_c), intent(in), target  :: props_obj
        !    this%props=>props_obj
        !end subroutine
        
      
        
        !subroutine set_source_term_flag(this)
        !    implicit none
        !    class(diffusion_c) :: this
        !    integer(kind=4) :: i
        !    allocate(this%props%source_term_flag(this%spatial_discr%Num_targets))
        !    this%props%source_term_flag=1
        !    do i=1,this%spatial_discr%Num_targets
        !        if (this%props%source_term(i)<0 .and. this%BCs%evap==.false.) then
        !            this%props%source_term_flag(i)=0
        !        end if
        !    end do
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
        
        subroutine set_diff_props_heterog(this,diff_props_heterog)
            implicit none
            class(diffusion_1D_c) :: this
            class(diff_props_heterog_c), intent(in) :: diff_props_heterog
            this%diff_props_heterog=diff_props_heterog
        end subroutine
        
       
        
      
        
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
        !    class(diffusion_c), intent(in) :: this
        !    class(spatial_discr_c), pointer :: spatial_discr
        !    spatial_discr=>this%spatial_discr
        !end function
        !
        !function get_props(this) result(props)
        !    implicit none
        !    class(diffusion_c), intent(in) :: this
        !    class(props_c), pointer :: props
        !    props=>this%props
        !end function
end module 