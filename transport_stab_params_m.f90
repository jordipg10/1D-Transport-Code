module transport_stab_params_m
    use diff_stab_params_m
    use transport_properties_heterog_m
    use spatial_discr_1D_m
    
    implicit none
    save
    type, public, extends(stab_params_diff_c) :: stab_params_tpt_c ! transport stability parameters subclass
        real(kind=8) :: alpha ! advection stability parameter (alpha=q*Delta_t/(2*phi*Delta_x))
        real(kind=8) :: Peclet ! Pe=|q|*Delta_x/D
    contains
        procedure, public :: compute_stab_params=>compute_stab_params_tpt
    end type
    
    contains
        subroutine compute_stab_params_tpt(this,props_obj,mesh_size,time_step)
            implicit none
            class(stab_params_tpt_c) :: this
            class(props_c), intent(in) :: props_obj
            real(kind=8), intent(in) :: mesh_size
            real(kind=8), intent(in) :: time_step
            
            real(kind=8) :: Delta_t_crit1,Delta_t_crit2
            real(kind=8), parameter :: epsilon=1d-9
            
            call compute_stab_params_diff(this,props_obj,mesh_size,time_step)
            
            select type (props_obj)
            type is (tpt_props_heterog_c)
                this%alpha=inf_norm_vec(props_obj%flux)*time_step/(2d0*minval(props_obj%porosity)*mesh_size)
                this%Peclet=inf_norm_vec(props_obj%flux)*mesh_size/minval(props_obj%dispersion)
            end select
            if (this%alpha>=5d-1) print *, "Unstable advection", this%alpha
            if (this%Peclet>2d0) print *, "Peclet too large", this%Peclet
        end subroutine
end module