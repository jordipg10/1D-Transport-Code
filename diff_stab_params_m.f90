module diff_stab_params_m
    use stability_parameters_m
    use diff_props_heterog_m
    use spatial_discr_rad_m
    use spatial_discr_1D_m
    use vectors_m
    implicit none
    save
    type, public, extends(stab_params_c) :: stab_params_diff_c ! diffusion equation stability parameters subclass
        real(kind=8) :: beta ! dispersion stability parameter (beta=D*Delta_t/phi*Delta_x^2)
    contains
        procedure, public :: compute_stab_params=>compute_stab_params_diff
    end type
    
    contains
        subroutine compute_stab_params_diff(this,props_obj,mesh_size,time_step)
            implicit none
            class(stab_params_diff_c) :: this
            !real(kind=8), intent(in) :: phi
            !real(kind=8), intent(in) :: D
            !real(kind=8), intent(in) :: Delta_r_D
            !real(kind=8), intent(in) :: Delta_t_D
            class(props_c), intent(in) :: props_obj
            real(kind=8), intent(in) :: mesh_size
            real(kind=8), intent(in) :: time_step
            
            !real(kind=8), allocatable :: mesh_size(:)
            !mesh_size=spatial_discr_obj%get_mesh_size(1)
            
            select type (props=>props_obj)
            class is (diff_props_heterog_c)
                this%beta=maxval(props%dispersion)*time_step/(minval(props%porosity)*mesh_size**2)
            end select
            if (this%beta>=5d-1) print *, "Unstable diffusion", this%beta
        end subroutine
end module