﻿! 1D steady-state transport equation:
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
        procedure, public :: set_tpt_props_heterog_obj
        procedure, public :: set_conc_r_flag
    ! Get

    ! Computations
        procedure, public :: compute_trans_mat_PDE=>compute_trans_mat_tpt
        procedure, public :: mass_balance_error_ADE_stat_Dirichlet_discharge
    ! Abstract
        procedure, public :: initialise_PDE=>initialise_transport_1D
        procedure, public :: write_PDE_1D=>write_transport_1D
    end type
!*****************************************************************************************************************************
    interface
        subroutine compute_trans_mat_tpt(this)
            import transport_1D_c
            implicit none
            class(transport_1D_c) :: this
            !real(kind=8), intent(out), optional :: T_sub(:),T_diag(:),T_super(:)
            !integer(kind=4), intent(in), optional :: k
        end subroutine
        
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
end module 
