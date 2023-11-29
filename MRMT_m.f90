! Multirate mass transfer module
module MRMT_m
    use PDE_model_m
    use transport_transient_m
    use mob_zone_m
    use imm_zone_m
    implicit none
    save
    type, public, extends(PDE_model_c) :: MRMT_c
        
        type(mob_zone_c) :: mob_zone ! mobile zone
        integer(kind=4) :: n_imm ! number of immobile zones
        class(imm_zone_c), allocatable :: imm_zones(:) ! immobile zones
        !real(kind=8) :: mob_por ! mobile porosity
        !real(kind=8) :: imm_por ! immobile porosity
        !real(kind=8), allocatable :: exchange_rates(:) ! (alpha)
        !real(kind=8), allocatable :: res_times(:) ! residence times (tau)
        !real(kind=8), allocatable :: probs(:) ! probabilities
        !real(kind=8), allocatable :: flux(:) ! flux of immobile zones
    contains
        !procedure, public :: set_PDE
        !procedure, public :: set_mob_zone
        procedure, public :: set_n_imm
        !procedure, public :: set_imm_zones
        procedure, public :: allocate_imm_zones
        procedure, public :: compute_A_mat_conc_mob
        procedure, public :: compute_b_conc_mob
        procedure, public :: compute_conc_imm_MRMT
        procedure, public :: solve_PDE_EI_Delta_t_homog_MRMT
        procedure, public :: check_imm_zones
    end type
    
    interface
        subroutine solve_PDE_EI_Delta_t_homog_MRMT(this,theta,Time_out,output)
            import MRMT_c
            class(MRMT_c) :: this
            real(kind=8), intent(in) :: theta
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(out) :: output(:,:)
        end subroutine
        
        subroutine compute_A_mat_conc_mob(this,theta,Delta_t,A_mat)
            import MRMT_c
            import tridiag_matrix_c
            implicit none
            class(MRMT_c), intent(in) :: this
            real(kind=8), intent(in) :: theta ! time weighting factor
            !real(kind=8), intent(in) :: alpha(:) ! exchange rates
            !real(kind=8), intent(in) :: prob(:) ! probabilities
            real(kind=8), intent(in) :: Delta_t ! time step
            !real(kind=8), intent(in) :: phi_mob ! mobile porosity
            !real(kind=8), intent(in) :: phi_imm ! immmobile porosity
            class(tridiag_matrix_c), intent(out) :: A_mat ! A*c_mob^(k+1)=b
        end subroutine
        
        subroutine compute_b_conc_mob(this,theta,Delta_t,conc_mob_old,conc_imm_old,b)
            import MRMT_c
            implicit none
            class(MRMT_c), intent(in) :: this
            real(kind=8), intent(in) :: theta ! time weighting factor
            !real(kind=8), intent(in) :: alpha(:) ! exchange rates
            !real(kind=8), intent(in) :: prob(:) ! probabilities
            real(kind=8), intent(in) :: Delta_t ! time step
            !real(kind=8), intent(in) :: phi_mob ! mobile porosity
            !real(kind=8), intent(in) :: phi_imm ! immmobile porosity
            real(kind=8), intent(in) :: conc_mob_old(:) ! c_mob^k
            real(kind=8), intent(in) :: conc_imm_old(:) ! c_imm^k
            real(kind=8), intent(out) :: b(:) ! A*c_mob^(k+1)=b
        end subroutine
        
        subroutine compute_conc_imm_MRMT(this,theta,conc_imm_old,conc_mob_old,conc_mob_new,Delta_t,conc_imm_new)
            import MRMT_c
            implicit none
            class(MRMT_c), intent(in) :: this
            real(kind=8), intent(in) :: theta ! time weighting factor
            real(kind=8), intent(in) :: conc_imm_old(:) ! c_imm^k
            real(kind=8), intent(in) :: conc_mob_old(:) ! c_m^k
            real(kind=8), intent(in) :: conc_mob_new(:) ! c_m^(k+1)
            !real(kind=8), intent(in) :: alpha(:) ! exchange rates
            real(kind=8), intent(in) :: Delta_t ! time step
            real(kind=8), intent(out) :: conc_imm_new(:) ! c_imm^(k+1)
        end subroutine

    end interface
    
    contains
        !subroutine set_PDE(this,PDE)
        !    implicit none
        !    class(MRMT_c) :: this
        !    class(diffusion_1D_transient_c), intent(in), target :: PDE
        !    this%PDE=>PDE
        !end subroutine
        
        subroutine set_n_imm(this,n_imm)
            implicit none
            class(MRMT_c) :: this
            integer(kind=4), intent(in) :: n_imm
            this%n_imm=n_imm
        end subroutine
        
        subroutine allocate_imm_zones(this)
            implicit none
            class(MRMT_c) :: this
            allocate(this%imm_zones(this%n_imm))
        end subroutine
        
        subroutine check_imm_zones(this)
            implicit none
            class(MRMT_c) :: this
            
            integer(kind=4) :: i
            real(kind=8) :: prob_tot=0d0
            real(kind=8), parameter :: epsilon=1d-12
            
            do i=1,this%n_imm
                if (this%imm_zones(i)%prob<0d0 .or. this%imm_zones(i)%prob>1d0) error stop "Error in probabilities"
                prob_tot=prob_tot+this%imm_zones(i)%prob
            end do
            if (abs(prob_tot-1d0)>=epsilon) error stop "Probabilities must sum 1"
        end subroutine
        
       
end module