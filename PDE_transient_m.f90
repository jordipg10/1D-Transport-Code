! Transient PDE module
! F*dc/dt=T*c+g
module PDE_transient_m
    use PDE_m
    use time_discr_m
    use stability_parameters_m
    use char_params_m
    implicit none
    save
    type, public, abstract, extends(PDE_1D_c) :: PDE_1D_transient_c ! subclass
        class(time_discr_c), pointer :: time_discr ! time discretisation (polymorphic variable)
        class(char_params_c), pointer :: char_params ! characteristic parameters (polymorphic variable)
        type(diag_matrix_c) :: F_mat ! F
    contains
    ! Set
        procedure, public :: set_time_discr
        !procedure, public :: set_stab_params
        !procedure, public :: compute_char_params
        procedure, public :: set_char_params
        !procedure, public :: initialise_props
        !procedure, public :: update_props
        !procedure, public :: set_spatial_discr
        !procedure, public :: set_BCs
        !procedure, public :: set_source_terms => set_source_terms_transport_1D
        !procedure, public :: set_conc_ext
        !procedure, public :: set_parameters => set_parameters_pointer
        !procedure, public :: set_props
        !procedure, public :: set_conc_star_flag
        !procedure, public :: set_source_term_flag
        !procedure, public :: set_conc_init => set_conc_init_transport_1D
        !procedure, public :: set_time_discr
        
    ! Get
        !procedure, public :: get_spatial_discr
        !procedure, public :: get_props
        !procedure, public :: get_time_discr
    ! Computations
        procedure(compute_F_mat_PDE), public, deferred :: compute_F_mat_PDE
        !procedure(compute_A_mat_ODE), public, deferred :: compute_A_mat_ODE
        !procedure(initialise_PDE_1D_transient), public, deferred :: initialise_PDE_1D_transient
        !procedure(solve_PDE_1D_transient), public, deferred :: solve_PDE_1D_transient
        !procedure, public :: main_PDE_1D_transient
        !procedure, public :: compute_trans_mat_diff
        !procedure, public :: compute_trans_mat_tpt
        !procedure, public :: compute_trans_mat
        !procedure, public :: compute_E_mat_diff
        !procedure, public :: compute_E_mat_tpt
        procedure, public :: compute_E_mat
        !procedure, public :: compute_mixing_ratios
        !procedure, public :: iteration_EE_1D
        !procedure, public :: main_transport_1D
        !procedure, public :: solve_transport_1D_EE_Delta_t_homog
        !procedure, public :: solve_transport_1D_EE_Delta_t_heterog
        !procedure, public :: transport_CFD_EE_1D
        !procedure, public :: transport_EFD_EE_1D
        !procedure, public :: transport_upwind_EE_1D
        !procedure, public :: solve_transport_1D_nonEE_Delta_t_homog
        !procedure, public :: solve_transport_1D_nonEE_Delta_t_heterog
        !procedure, public :: write_transport_1D
        procedure, public :: compute_B_mat
        procedure, public :: compute_A_mat_lin_syst
        procedure, public :: compute_f
        procedure, public :: compute_b_lin_syst
        procedure, public :: compute_A_mat_ODE
        procedure, public :: compute_b_ODE     
        !procedure, public :: solve_PDE_1D_transient
        !procedure, public :: prod_total_sym_mat
        procedure, public :: solve_PDE_EE_Delta_t_homog
        procedure, public :: solve_PDE_EE_Delta_t_heterog
        procedure, public :: solve_PDE_EI_Delta_t_homog
        procedure, public :: solve_PDE_RKF45
        procedure, public :: compute_k_RKF45
        !procedure, public :: compute_A_conc_mob
        !procedure, public :: compute_b_conc_mob
        !procedure, public :: compute_conc_imm_MRMT
        !procedure, public :: solve_PDE_EI_Delta_t_homog_MRMT_bis
    end type
!*****************************************************************************************************************************
    abstract interface
        
        subroutine compute_F_mat_PDE(this)
            import PDE_1D_transient_c
            implicit none
            class(PDE_1D_transient_c) :: this
        end subroutine
        
        !subroutine initialise_PDE_1D_transient(this)
        !    import PDE_1D_transient_c
        !    implicit none
        !    class(PDE_1D_transient_c) :: this
        !end subroutine
        
        !subroutine compute_A_mat_ODE(this,A_mat)
        !    import PDE_1D_transient_c
        !    import matrix_c
        !    implicit none
        !    class(PDE_1D_transient_c), intent(in) :: this
        !    class(matrix_c), pointer, intent(out) :: A_mat
        !end subroutine
        
        !subroutine solve_PDE_1D_transient(this,Time_out,output,anal_sol)
        !    import PDE_1D_transient_c
        !    class(PDE_1D_transient_c), intent(in) :: this
        !    real(kind=8), intent(in) :: Time_out(:)
        !    real(kind=8), intent(out) :: output(:,:)
        !    real(kind=8), external, optional :: anal_sol
        !end subroutine
    end interface
!*****************************************************************************************************************************
    interface
        !subroutine compute_char_params(this,k)
        !    import PDE_1D_transient_c
        !    class(PDE_1D_transient_c) :: this
        !    integer(kind=4), intent(in), optional :: k
        !end subroutine
        
        subroutine compute_E_mat(this,E_mat,k)
            import PDE_1D_transient_c
            import tridiag_matrix_c
            class(PDE_1D_transient_c) :: this
            type(tridiag_matrix_c), intent(out) :: E_mat
            integer(kind=4), intent(in), optional :: k
        end subroutine
        
        subroutine compute_B_mat(this,theta,B_mat,k)
            import PDE_1D_transient_c
            import tridiag_matrix_c
            implicit none
            class(PDE_1D_transient_c) :: this
            real(kind=8), intent(in) :: theta
            type(tridiag_matrix_c), intent(out) :: B_mat
            integer(kind=4), intent(in), optional :: k
            !class(tridiag_matrix_c), intent(in), optional :: E_mat
        end subroutine
        
        !subroutine compute_B_mat_bis(this,theta,B_mat,k)
        !    import PDE_1D_transient_c
        !    import tridiag_matrix_c
        !    implicit none
        !    class(PDE_1D_transient_c) :: this
        !    real(kind=8), intent(in) :: theta
        !    class(tridiag_matrix_c), intent(out) :: B_mat
        !    integer(kind=4), intent(in), optional :: k
        !end subroutine
        
        subroutine compute_A_mat_lin_syst(this,theta,A_mat,k)
            import PDE_1D_transient_c
            import tridiag_matrix_c
            implicit none
            class(PDE_1D_transient_c), intent(in) :: this
            !class(tridiag_matrix_c), intent(in) :: E_mat
            real(kind=8), intent(in) :: theta
            type(tridiag_matrix_c), intent(out) :: A_mat
            integer(kind=4), intent(in), optional :: k
        end subroutine
        
        function compute_f(this,k) result(f)
            import PDE_1D_transient_c
            implicit none
            class(PDE_1D_transient_c) :: this
            integer(kind=4), intent(in), optional :: k
            real(kind=8), allocatable :: f(:)
        end function
        
        subroutine compute_b_lin_syst(this,theta,conc_old,b,k)
            import PDE_1D_transient_c
            import tridiag_matrix_c
            implicit none
            class(PDE_1D_transient_c), intent(in) :: this
            !class(tridiag_matrix_c), intent(in) :: B_mat
            real(kind=8), intent(in) :: theta
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(inout) :: b(:)
            integer(kind=4), intent(in), optional :: k
            !real(kind=8), allocatable :: b(:)
        !end function
        end subroutine
        
        subroutine compute_A_mat_ODE(this,A_mat)
            import PDE_1D_transient_c
            import tridiag_matrix_c
            implicit none
            class(PDE_1D_transient_c), intent(in) :: this
            type(tridiag_matrix_c), intent(out) :: A_mat
        end subroutine
        
        
        function compute_b_ODE(this) result(b)
            import PDE_1D_transient_c
            implicit none
            class(PDE_1D_transient_c), intent(in) :: this
            real(kind=8), allocatable :: b(:)
        end function
        
        subroutine solve_PDE_EE_Delta_t_homog(this,Time_out,output)
            import PDE_1D_transient_c
            class(PDE_1D_transient_c) :: this
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(out) :: output(:,:)
        end subroutine
        
        subroutine solve_PDE_EE_Delta_t_heterog(this,Time_out,output)
            import PDE_1D_transient_c
            class(PDE_1D_transient_c) :: this
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(out) :: output(:,:)
        end subroutine
        
        subroutine solve_PDE_EI_Delta_t_homog(this,theta,Time_out,output)
            import PDE_1D_transient_c
            class(PDE_1D_transient_c) :: this
            real(kind=8), intent(in) :: theta
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(out) :: output(:,:)
        end subroutine
        
        !subroutine solve_PDE_EI_Delta_t_homog_MRMT_bis(this,theta,Time_out,output)
        !    import PDE_1D_transient_c
        !    class(PDE_1D_transient_c) :: this
        !    real(kind=8), intent(in) :: theta
        !    real(kind=8), intent(in) :: Time_out(:)
        !    real(kind=8), intent(out) :: output(:,:)
        !end subroutine
        !
        !subroutine compute_A_conc_mob(this,theta,alpha,prob,Delta_t,phi_mob,phi_imm,A)
        !    import PDE_1D_transient_c
        !    import tridiag_matrix_c
        !    implicit none
        !    class(PDE_1D_transient_c), intent(in) :: this
        !    real(kind=8), intent(in) :: theta ! time weighting factor
        !    real(kind=8), intent(in) :: alpha(:) ! exchange rates
        !    real(kind=8), intent(in) :: prob(:) ! probabilities
        !    real(kind=8), intent(in) :: Delta_t ! time step
        !    real(kind=8), intent(in) :: phi_mob ! mobile porosity
        !    real(kind=8), intent(in) :: phi_imm ! immmobile porosity
        !    class(tridiag_matrix_c), intent(out) :: A ! A*c_mob^(k+1)=b
        !end subroutine
        !
        !subroutine compute_b_conc_mob(this,theta,alpha,prob,Delta_t,phi_mob,phi_imm,conc_mob_old,conc_imm_old,b)
        !    import PDE_1D_transient_c
        !    implicit none
        !    class(PDE_1D_transient_c), intent(in) :: this
        !    real(kind=8), intent(in) :: theta ! time weighting factor
        !    real(kind=8), intent(in) :: alpha(:) ! exchange rates
        !    real(kind=8), intent(in) :: prob(:) ! probabilities
        !    real(kind=8), intent(in) :: Delta_t ! time step
        !    real(kind=8), intent(in) :: phi_mob ! mobile porosity
        !    real(kind=8), intent(in) :: phi_imm ! immmobile porosity
        !    real(kind=8), intent(in) :: conc_mob_old(:) ! c_mob^k
        !    real(kind=8), intent(in) :: conc_imm_old(:) ! c_imm^k
        !    real(kind=8), intent(out) :: b(:) ! A*c_mob^(k+1)=b
        !end subroutine
        !
        !subroutine compute_conc_imm_MRMT(this,theta,conc_imm_old,conc_mob_old,conc_mob_new,alpha,Delta_t,conc_imm_new)
        !    import PDE_1D_transient_c
        !    implicit none
        !    class(PDE_1D_transient_c), intent(in) :: this
        !    real(kind=8), intent(in) :: theta ! time weighting factor
        !    real(kind=8), intent(in) :: conc_imm_old(:) ! c_imm^k
        !    real(kind=8), intent(in) :: conc_mob_old(:) ! c_m^k
        !    real(kind=8), intent(in) :: conc_mob_new(:) ! c_m^(k+1)
        !    real(kind=8), intent(in) :: alpha(:) ! exchange rates
        !    real(kind=8), intent(in) :: Delta_t ! time step
        !    real(kind=8), intent(out) :: conc_imm_new(:) ! c_imm^(k+1)
        !end subroutine
        
        subroutine solve_PDE_RKF45(this,Delta_t_init,tolerance)
            import PDE_1D_transient_c
            class(PDE_1D_transient_c) :: this
            real(kind=8), intent(in) :: Delta_t_init
            real(kind=8), intent(in) :: tolerance
            !real(kind=8), intent(in) :: Time_out(:)
            !real(kind=8), intent(out) :: output(:,:)
        end subroutine
        
        subroutine update_time_step_RKF45(Delta_t_old,tolerance,conc_RK4,conc_RK5,Delta_t_new)
            implicit none
            real(kind=8), intent(in) :: Delta_t_old
            real(kind=8), intent(in) :: tolerance
            real(kind=8), intent(in) :: conc_RK4(:)
            real(kind=8), intent(in) :: conc_RK5(:)
            real(kind=8), intent(out) :: Delta_t_new
        end subroutine
        
        function compute_k_RKF45(this,Delta_t,conc_RK4) result(k)
            import PDE_1D_transient_c
            implicit none
            class(PDE_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: conc_RK4(:)
            real(kind=8), allocatable :: k(:,:)
        end function
    end interface
!*****************************************************************************************************************************
    contains
        subroutine set_time_discr(this,time_discr_obj)
            implicit none
            class(PDE_1D_transient_c) :: this
            class(time_discr_c), intent(in), target :: time_discr_obj
            this%time_discr=>time_discr_obj
        end subroutine
        
        !subroutine set_stab_params(this,stab_params_obj)
        !    implicit none
        !    class(PDE_1D_transient_c) :: this
        !    class(stab_params_c), intent(in), target :: stab_params_obj
        !    this%stab_params=>stab_params_obj
        !end subroutine
        
        subroutine set_char_params(this,char_params_obj)
            implicit none
            class(PDE_1D_transient_c) :: this
            class(char_params_c), intent(in), target :: char_params_obj
            this%char_params=>char_params_obj
        end subroutine
        
        !subroutine initialise_props(this,props_conc_init)
        !    implicit none
        !    class(PDE_1D_transient_c) :: this
        !    class(props_c), intent(in), target :: props_conc_init
        !    this%props=>props_conc_init
        !end subroutine
        
        !subroutine update_props(this,props_new)
        !    implicit none
        !    class(PDE_1D_transient_c) :: this
        !    class(props_c), intent(in), target :: props_new
        !    this%props=>props_new
        !end subroutine
        
       ! subroutine get_time_discr(this,time_discr)
       !     implicit none
       !     class(PDE_1D_transient_c), intent(in) :: this
       !     class(time_discr_c), intent(out) :: time_discr
       !     time_discr%Final_time=this%time_discr%Final_time
       !     time_discr%Num_time=this%time_discr%Num_time
       !     select type (time=>this%time_discr)
       !     type is (time_discr_homog_c)
       !        ! time_discr%Delta_t=time%Delta_t
       !     end select
       !     !time_discr=time
       !     !time_discr=this%time_discr
       !end subroutine
        
        
end module