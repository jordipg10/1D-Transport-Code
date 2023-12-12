module diffusion_transient_m
    use PDE_transient_m
    !use time_discr_m
    use diff_stab_params_m
    implicit none
    save
    type, public, extends(PDE_1D_transient_c) :: diffusion_1D_transient_c ! transient subclass
        real(kind=8), allocatable :: conc(:) ! concentration (c)
        real(kind=8), allocatable :: conc_ext(:) ! (c_e)
        integer(kind=4), allocatable :: conc_star_flag(:)   ! 1 if r>0
                                                            ! 0 otherwise
        real(kind=8), allocatable :: conc_init(:) ! initial concentration (c_0)
        type(tridiag_matrix_vec_c) :: mixing_ratios
        type(diff_props_heterog_c) :: diff_props_heterog        ! properties
        type(stab_params_diff_c) :: stab_params_diff            ! stability parameters
    contains
        procedure, public :: set_conc_init
        procedure, public :: set_conc_ext
        procedure, public :: set_conc_star_flag=>set_conc_star_flag_diff
        !procedure, public :: set_source_term_flag
        procedure, public :: compute_trans_mat_PDE=>compute_trans_mat_diff
        !procedure, public :: compute_source_term_PDE=>compute_g
        procedure, public :: compute_F_mat_PDE=>compute_F_mat_diff
        !procedure, public :: compute_E_mat=>compute_E_mat_diff
        !procedure, public :: compute_A_mat_ODE=>compute_A_mat_ODE_diff
        procedure, public :: initialise_PDE=>initialise_diffusion_transient
        procedure, public :: compute_mixing_ratios
        !procedure, public :: compute_b_ODE=>compute_b_ODE_diff
        procedure, public :: set_stab_params_diff
        !procedure, public :: initialise_props
        !procedure, public :: update_props
        procedure, public :: update_conc_ext
        procedure, public :: prod_total_conc
        !procedure, public :: solve_PDE_1D=>solve_diffusion_transient
        !procedure, public :: solve_diffusion_EE_Delta_t_homog
        !procedure, public :: solve_diffusion_EE_Delta_t_heterog
        !procedure, public :: main_PDE=>main_diffusion_transient
        procedure, public :: write_PDE_1D=>write_diffusion_transient
        !procedure, public :: solve_write_diffusion_transient
        !procedure, public :: compute_c_tilde
        procedure, public :: set_diff_props_heterog
    end type
!****************************************************************************************************************************************************
    interface
        subroutine compute_F_mat_diff(this)
            import diffusion_1D_transient_c
            implicit none
            class(diffusion_1D_transient_c) :: this
        end subroutine
        
        subroutine compute_trans_mat_diff(this)
            import diffusion_1D_transient_c
            implicit none
            class(diffusion_1D_transient_c) :: this
            !integer(kind=4), intent(in), optional :: k
        end subroutine
        
        subroutine initialise_diffusion_transient(this)
            import diffusion_1D_transient_c
            implicit none
            class(diffusion_1D_transient_c) :: this
        end subroutine
        
        !subroutine compute_g(this,k,anal_sol)
        !    import diffusion_transient_c
        !    implicit none
        !    class(diffusion_transient_c) :: this
        !    integer(kind=4), intent(in), optional :: k
        !    real(kind=8), external, optional :: anal_sol
        !end subroutine
        
        !subroutine compute_E_mat_diff(this,E_mat,k)
        !    import diffusion_transient_c
        !    implicit none
        !    class(diffusion_transient_c) :: this
        !    class(matrix_c), pointer, intent(out) :: E_mat
        !    integer(kind=4), intent(in), optional :: k
        !end subroutine
        
        !subroutine compute_A_mat_ODE_diff(this,A_mat)
        !    import diffusion_transient_c
        !    import matrix_c
        !    implicit none
        !    class(diffusion_transient_c), intent(in) :: this
        !    class(matrix_c), pointer, intent(out) :: A_mat
        !end subroutine
        
        !function compute_Q_mat(this,k) result(Q_mat)
        !    import diffusion_transient_c
        !    implicit none
        !    class(diffusion_transient_c) :: this
        !    integer(kind=4), intent(in), optional :: k
        !    real(kind=8), allocatable :: Q_mat(:)
        !end function
        
        subroutine compute_mixing_ratios(this,theta,props,k)
            import diffusion_1D_transient_c
            import props_c
            implicit none
            class(diffusion_1D_transient_c) :: this
            !class(matrix_c), intent(in) :: B_mat
            real(kind=8), intent(in) :: theta
            class(props_c), intent(in), optional :: props
            integer(kind=4), intent(in), optional :: k
        end subroutine
        
        !subroutine prod_total_conc(this,k)
        !    import diffusion_transient_c
        !    implicit none
        !    class(diffusion_transient_c) :: this
        !    integer(kind=4), intent(in), optional :: k
        !end subroutine
        
        subroutine prod_total_conc(this,A_mat,time)
            import diffusion_1D_transient_c
            import tridiag_matrix_c
            implicit none
            class(diffusion_1D_transient_c) :: this
            real(kind=8), intent(in), optional :: time
            class(tridiag_matrix_c), intent(in) :: A_mat
            !real(kind=8), intent(in), optional :: Time_out(:)
            !real(kind=8), intent(out), optional :: output(:,:)
        end subroutine
        
        !subroutine solve_diffusion_transient(this,Time_out,output,anal_sol)
        !    import diffusion_transient_c
        !    class(diffusion_transient_c) :: this
        !    real(kind=8), intent(in) :: Time_out(:)
        !    real(kind=8), intent(out) :: output(:,:)
        !    real(kind=8), external, optional :: anal_sol
        !end subroutine
        !
        !subroutine solve_diffusion_EE_Delta_t_homog(this,Time_out,conc_out,anal_sol)
        !    import diffusion_transient_c
        !    class(diffusion_transient_c) :: this
        !    real(kind=8), intent(in) :: Time_out(:)
        !    real(kind=8), intent(out) :: conc_out(:,:)
        !    real(kind=8), external, optional :: anal_sol
        !end subroutine
        !
        !subroutine solve_diffusion_EE_Delta_t_heterog(this,Time_out,conc_out,anal_sol)
        !    import diffusion_transient_c
        !    class(diffusion_transient_c) :: this
        !    real(kind=8), intent(in) :: Time_out(:)
        !    real(kind=8), intent(out) :: conc_out(:,:)
        !    real(kind=8), external, optional :: anal_sol
        !end subroutine
        
        !subroutine main_diffusion_transient(this)
        !    import diffusion_transient_c
        !    class(diffusion_transient_c) :: this
        !end subroutine
        
        subroutine write_diffusion_transient(this,Time_out,output)
            import diffusion_1D_transient_c
            import props_c
            implicit none
            class(diffusion_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(in) :: output(:,:)
            !class(props_c), intent(in), optional :: props
            !integer(kind=4), intent(in), optional :: opcion
        end subroutine
        
        subroutine solve_write_diffusion_transient(this,Time_out)
            import diffusion_1D_transient_c
            implicit none
            class(diffusion_1D_transient_c) :: this
            real(kind=8), intent(in) :: Time_out(:)
        end subroutine
        
        function compute_c_tilde(this,j,conc,conc_star,mixing_ratios) result(c_tilde)
            import diffusion_1D_transient_c
            import tridiag_matrix_vec_c
            implicit none
            class(diffusion_1D_transient_c), intent(in) :: this
            integer(kind=4), intent(in) :: j
            real(kind=8), intent(in) :: conc(:,:)
            real(kind=8), intent(in) :: conc_star(:)
            class(tridiag_matrix_vec_c), intent(in) :: mixing_ratios
            real(kind=8), allocatable :: c_tilde(:)
        end function
    end interface
!****************************************************************************************************************************************************
    contains
        subroutine set_conc_init(this,conc_init)
            implicit none
            class(diffusion_1D_transient_c) :: this
            real(kind=8), intent(in) :: conc_init(:)
            if (this%spatial_discr%Num_targets_defined==.true.) then
                if (size(conc_init)/=this%spatial_discr%Num_targets) error stop "Dimension error in initial concentration"
            else
                this%spatial_discr%Num_targets=size(conc_init)
                this%spatial_discr%Num_targets_defined=.true.
            end if
            this%conc_init=conc_init
        end subroutine
    
        !subroutine set_time_discr(this,time_discr_obj)
        !    implicit none
        !    class(diffusion_transient_c) :: this
        !    class(time_discr_c), intent(in), target :: time_discr_obj
        !    this%time_discr=>time_discr_obj
        !end subroutine
        !
        subroutine set_stab_params_diff(this,stab_params_diff)
            implicit none
            class(diffusion_1D_transient_c) :: this
            type(stab_params_diff_c), intent(in) :: stab_params_diff
            this%stab_params_diff=stab_params_diff
        end subroutine
        !
        !subroutine initialise_props(this,props_conc_init)
        !    implicit none
        !    class(diffusion_transient_c) :: this
        !    class(diff_props_c), intent(in), target :: props_conc_init
        !    this%props=>props_conc_init
        !end subroutine
        !
        !subroutine update_props(this,props_new)
        !    implicit none
        !    class(diffusion_transient_c) :: this
        !    class(diff_props_c), intent(in), target :: props_new
        !    this%props=>props_new
        !end subroutine
        
        subroutine set_conc_ext(this,conc_ext)
            class(diffusion_1D_transient_c) :: this
            real(kind=8), intent(in) :: conc_ext(:)
            if (size(conc_ext)/=this%spatial_discr%Num_targets) error stop "Dimension error in external concentration"
            this%conc_ext=conc_ext
        end subroutine 
        
        subroutine update_conc_ext(this,conc_ext_new)
            implicit none
            class(diffusion_1D_transient_c) :: this
            real(kind=8), intent(in) :: conc_ext_new
            this%conc_ext=conc_ext_new
        end subroutine
        
        subroutine set_conc_star_flag_diff(this)
            implicit none
            class(diffusion_1D_transient_c) :: this
            integer(kind=4) :: i
            allocate(this%conc_star_flag(this%spatial_discr%Num_targets))
            this%conc_star_flag=0
            do i=1,this%spatial_discr%Num_targets
                if (this%diff_props_heterog%source_term(i)>0) then
                    this%conc_star_flag(i)=1
                end if
            end do
        end subroutine
        
        !subroutine set_source_term_flag(this)
        !    implicit none
        !    class(diffusion_transient_c) :: this
        !    integer(kind=4) :: i
        !    allocate(this%props%source_term_flag(this%spatial_discr%Num_targets))
        !    this%props%source_term_flag=1
        !    do i=1,this%spatial_discr%Num_targets
        !        if (this%props%source_term(i)<0 .and. this%BCs%evap==.false.) then
        !            this%props%source_term_flag(i)=0
        !        end if
        !    end do
        !end subroutine
        
        subroutine set_diff_props_heterog(this,diff_props_heterog)
            implicit none
            class(diffusion_1D_transient_c) :: this
            class(diff_props_heterog_c), intent(in) :: diff_props_heterog
            this%diff_props_heterog=diff_props_heterog
        end subroutine
end module