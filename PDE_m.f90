! PDE module
! T*c+g=0
module PDE_m
    !use spatial_discr_m
    use properties_m
    !use spatial_discr_1D_m
    !use spatial_discr_rad_m
    !use transport_properties_homog_m
    !use transport_properties_heterog_m
    !use properties_viena_m
    use matrices_m
    implicit none
    save
    type, public, abstract :: PDE_1D_c ! superclass
        class(spatial_discr_c), pointer :: spatial_discr ! spatial discretisation (polymorphic variable)
        type(BCs_t) :: BCs ! Boundary conditions
        type(tridiag_matrix_c) :: trans_mat ! Transition matrix (T)
        real(kind=8), allocatable :: source_term_PDE(:) ! g
        logical :: dimensionless
        integer(kind=4) :: sol_method   ! 1: Numerical
                                        ! 2: Eigendecomposition
    contains
    ! Set
        procedure, public :: set_spatial_discr
        procedure, public :: set_BCs
        !procedure, public :: set_conc_ext
        !procedure, public :: set_props
        !procedure, public :: set_props_homog
        !procedure, public :: set_props_heterog
        !procedure, public :: set_conc_star_flag
        !procedure, public :: set_source_term_flag
        !procedure, public :: set_conc_init
        !procedure, public :: set_time_discr
        procedure, public :: set_sol_method
    ! Allocate
        procedure, public :: allocate_trans_mat
        procedure, public :: allocate_source_term_PDE
    ! Update
        procedure, public :: update_trans_mat
    ! Get
        !procedure, public :: get_spatial_discr
        !procedure, public :: get_props
        !procedure, public :: get_trans_mat
    ! Computations
        procedure(compute_trans_mat_PDE), public, deferred :: compute_trans_mat_PDE
        procedure(initialise_PDE), public, deferred :: initialise_PDE
        !procedure(solve_PDE_1D), public, deferred :: solve_PDE_1D
        !procedure(main_PDE), public, deferred :: main_PDE
        procedure(write_PDE_1D), public, deferred :: write_PDE_1D
        !procedure, public :: compute_trans_mat_diff
        !procedure, public :: compute_trans_mat_tpt
        !procedure, public :: compute_trans_mat
        !procedure, public :: compute_E_mat_diff
        !procedure, public :: compute_E_mat_tpt
        !procedure, public :: compute_E_mat
        !procedure, public :: compute_Q_mat
        !procedure, public :: compute_mixing_ratios
        !procedure, public :: iteration_EE_1D
        
        !procedure(solve_PDE), public, deferred :: solve_PDE
        !procedure, public :: solve_transport_1D_EE_Delta_t_homog
        !procedure, public :: solve_transport_1D_EE_Delta_t_heterog
        !procedure, public :: transport_CFD_EE_1D
        !procedure, public :: transport_EFD_EE_1D
        !procedure, public :: transport_upwind_EE_1D
        !procedure, public :: solve_transport_1D_nonEE_Delta_t_homog
        !procedure, public :: solve_transport_1D_nonEE_Delta_t_heterog
        !procedure, public :: write_transport_1D
        !procedure, public :: compute_A_lin_syst
        !procedure, public :: compute_b_mat
        !procedure, public :: compute_b_lin_syst
        procedure, public :: compute_source_term_PDE
        procedure, public :: solve_PDE_1D
        procedure, public :: solve_write_PDE_1D
        procedure, public :: main_PDE
        !procedure, public :: compute_F_mat_tpt
        !procedure, public :: compute_F_mat
        !procedure, public :: compute_A_mat_ODE
        !procedure, public :: compute_b_ODE        
        !procedure, public :: prod_total_conc
        !procedure, public :: prod_total_sym_mat
        procedure, public :: solve_PDE_1D_stat
    end type
!****************************************************************************************************************************************************
    abstract interface
        subroutine compute_trans_mat_PDE(this)
            import PDE_1D_c
            class(PDE_1D_c) :: this
        end subroutine
        
        subroutine initialise_PDE(this)
            import PDE_1D_c
            class(PDE_1D_c) :: this
        end subroutine
        
        !subroutine solve_PDE_1D(this,Time_out,output,anal_sol)
        !    import PDE_1D_c
        !    class(PDE_1D_c) :: this
        !    real(kind=8), intent(in) :: Time_out(:)
        !    real(kind=8), intent(out) :: output(:,:)
        !    real(kind=8), external, optional :: anal_sol
        !end subroutine
        
        !subroutine main_PDE(this)
        !    import PDE_1D_c
        !    class(PDE_1D_c) :: this
        !end subroutine
        
        subroutine write_PDE_1D(this,Time_out,output)
            import PDE_1D_c
            import props_c
            class(PDE_1D_c), intent(in) :: this
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(in) :: output(:,:)
            !class(props_c), intent(in), optional :: props
            !integer(kind=4), intent(in), optional :: opcion
            !real(kind=8), intent(in), optional :: theta
            !real(kind=8), external, optional :: anal_sol
        end subroutine
    end interface
!****************************************************************************************************************************************************
    interface
        subroutine compute_source_term_PDE(this,k,anal_sol)
            import PDE_1D_c
            class(PDE_1D_c) :: this
            integer(kind=4), intent(in), optional :: k
            real(kind=8), external, optional :: anal_sol
        end subroutine
        
        subroutine solve_PDE_1D(this,Time_out,output)
            import PDE_1D_c
            class(PDE_1D_c) :: this
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(out) :: output(:,:)
            !integer(kind=4), intent(in), optional :: opcion
            !real(kind=8), external, optional :: anal_sol
        end subroutine
        
        subroutine solve_PDE_1D_stat(this)
            import PDE_1D_c
            class(PDE_1D_c) :: this
        end subroutine
        
        subroutine solve_write_PDE_1D(this,Time_out)
            import PDE_1D_c
            class(PDE_1D_c) :: this
            real(kind=8), intent(in) :: Time_out(:)
            !integer(kind=4), intent(in), optional :: opcion
            !real(kind=8), external, optional :: anal_sol
        end subroutine
        
        subroutine main_PDE(this)
            import PDE_1D_c
            class(PDE_1D_c) :: this
            !integer(kind=4), intent(in), optional :: opcion
        end subroutine
    end interface
!****************************************************************************************************************************************************
    contains
        subroutine set_spatial_discr(this,spatial_discr_obj)
            implicit none
            class(PDE_1D_c) :: this
            class(spatial_discr_c), intent(in), target :: spatial_discr_obj
            this%spatial_discr=>spatial_discr_obj
        end subroutine 
        
        subroutine set_BCs(this,BCs_obj)
            implicit none
            class(PDE_1D_c) :: this
            class(BCs_t), intent(in) :: BCs_obj
            this%BCs=BCs_obj
        end subroutine
        
        !subroutine set_props(this,props_obj)
        !    implicit none
        !    class(PDE_1D_c) :: this
        !    class(props_c), intent(in), target :: props_obj
        !    if (size(props_obj%source_term)/=this%spatial_discr%Num_targets) error stop "Dimension error in source term"
        !    this%props=>props_obj
        !end subroutine
        
        !subroutine set_props(this,props_obj)
        !    implicit none
        !    class(PDE_1D_c) :: this
        !    class(tpt_props_homog_c), intent(in) :: props_obj
        !    if (size(props_obj%source_term)/=this%spatial_discr%Num_targets) error stop "Dimension error in source term"
        !    this%props=props_obj
        !end subroutine
        
        !subroutine set_props(this,props_obj)
        !    implicit none
        !    class(PDE_1D_c) :: this
        !    class(tpt_props_heterog_c), intent(in) :: props_obj
        !    if (size(props_obj%source_term)/=this%spatial_discr%Num_targets) error stop "Dimension error in source term"
        !    this%props=props_obj
        !end subroutine
        
        !subroutine set_props(this,props_obj)
        !    implicit none
        !    class(PDE_1D_c) :: this
        !    class(diff_props_homog_c), intent(in) :: props_obj
        !    if (size(props_obj%source_term)/=this%spatial_discr%Num_targets) error stop "Dimension error in source term"
        !    this%props=props_obj
        !end subroutine
        
        subroutine allocate_trans_mat(this)
            implicit none
            class(PDE_1D_c) :: this
            !class(tridiag_sym_matrix_c), intent(in), target  :: trans_mat
            call this%trans_mat%allocate_matrix(this%spatial_discr%Num_targets)
            !this%trans_mat=>trans_mat
        end subroutine
        
        !subroutine update_trans_mat(this,trans_mat)
        !    implicit none
        !    class(PDE_1D_c) :: this
        !    class(tridiag_sym_matrix_c), intent(in), target  :: trans_mat
        !    this%trans_mat=>trans_mat
        !    !call this%trans_mat%allocate_matrix(this%spatial_discr%Num_targets)
        !end subroutine
        
        subroutine update_trans_mat(this,trans_mat)
            implicit none
            class(PDE_1D_c) :: this
            class(tridiag_matrix_c), intent(in)  :: trans_mat
            this%trans_mat=trans_mat
        end subroutine
        
        subroutine allocate_source_term_PDE(this,source_term_PDE)
            implicit none
            class(PDE_1D_c) :: this
            real(kind=8), intent(in)  :: source_term_PDE(:)
            allocate(this%source_term_PDE(this%spatial_discr%Num_targets))
        end subroutine
        
        !function get_spatial_discr(this) result(spatial_discr)
        !    implicit none
        !    class(PDE_1D_c), intent(in) :: this
        !    class(spatial_discr_c), pointer :: spatial_discr
        !    spatial_discr=>this%spatial_discr
        !end function
        
        !function get_props(this) result(props)
        !    implicit none
        !    class(PDE_1D_c), intent(in) :: this
        !    class(props_c), pointer :: props
        !    props=>this%props
        !end function
        
        !subroutine get_trans_mat(this,trans_mat)
        !    implicit none
        !    class(PDE_1D_c), intent(in) :: this
        !    class(diag_matrix_c), intent(out) :: trans_mat
        !    trans_mat%diag=this%trans_mat%diag
        !end subroutine
        
        subroutine set_sol_method(this,method)
            implicit none
            class(PDE_1D_c) :: this
            integer(kind=4), intent(in) :: method
            if (method<0 .or. method>2) error stop "Solution method not implemented"
            this%sol_method=method
        end subroutine
end module 