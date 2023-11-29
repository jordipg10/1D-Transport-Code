module transport_stab_params_m
    use diff_stab_params_m
    use transport_properties_heterog_m
    use spatial_discr_1D_m
    
    implicit none
    save
    type, public, extends(stab_params_diff_c) :: stab_params_tpt_c ! homogeneous transport stability parameters subclass
        real(kind=8) :: alpha ! advection stability parameter (alpha=q*Delta_t/(2*phi*Delta_x))
        real(kind=8) :: Peclet ! Pe=|q|*Delta_x/D
    contains
        procedure, public :: compute_stab_params=>compute_stab_params_tpt
        procedure, public :: get_stab_params => get_stab_params_tpt_homog
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
                !this%Delta_x_crit=2*props_obj%dispersion/props_obj%flux
                !select type (spatial_discr_obj)
                !type is (mesh_1D_Euler_homog_c)
                !    if (abs(spatial_discr_obj%Delta_x-this%Delta_x_crit)<epsilon) then
                !        this%Delta_t_crit=2*props_obj%porosity*props_obj%retardo*(this%Delta_x_crit**2)/(2*props_obj%dispersion+props_obj%source_term(1)*this%Delta_x_crit**2)
                !    else if (spatial_discr_obj%Delta_x<this%Delta_x_crit) then
                !        this%Delta_t_crit=2*props_obj%porosity*props_obj%retardo*(2*props_obj%dispersion-sqrt(4*props_obj%dispersion**2-(props_obj%flux**2)*(spatial_discr_obj%Delta_x**2)))/(props_obj%flux**2)
                !    else
                !        Delta_t_crit1=props_obj%porosity*props_obj%retardo*spatial_discr_obj%Delta_x/props_obj%flux
                !        Delta_t_crit2=2*props_obj%porosity*props_obj%retardo/(props_obj%source_term(1)+2*props_obj%dispersion/(spatial_discr_obj%Delta_x**2)+props_obj%flux/spatial_discr_obj%Delta_x)
                !        this%Delta_t_crit=minval([Delta_t_crit1,Delta_t_crit2])
                !    end if
                !    select type (time_discr_obj)
                !    type is (time_discr_homog_c)
                !        this%alpha=props_obj%velocity*time_discr_obj%Delta_t/(2*props_obj%retardo*spatial_discr_obj%Delta_x)
                !        this%Peclet=abs(props_obj%velocity)*spatial_discr_obj%Delta_x/(2*props_obj%dispersion)
                !    end select
                !end select
            end select
            if (this%alpha>=5d-1) print *, "Unstable advection", this%alpha
            if (this%Peclet>2d0) print *, "Peclet too large", this%Peclet
        end subroutine
        
        function get_stab_params_tpt_homog(this) result(stab_params)
            implicit none
            class(stab_params_tpt_c) :: this
            real(kind=8), allocatable :: stab_params(:)
            allocate(stab_params(5))
            stab_params(1)=this%Delta_x_crit
            stab_params(2)=this%Delta_t_crit
            stab_params(3)=this%beta
            stab_params(4)=this%alpha
            stab_params(5)=this%Peclet            
        end function
        
        !function get_courant(this) result(courant)
        !    implicit none
        !    class(stab_params_homog_c) :: this
        !    real(kind=8) :: courant
        !    courant=this%courant
        !end function
        
        !function get_beta(this) result(beta)
        !    implicit none
        !    class(stab_params_homog_c) :: this
        !    real(kind=8) :: beta
        !    beta=this%beta
        !end function
end module