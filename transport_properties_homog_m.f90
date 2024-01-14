module transport_properties_homog_m
    use diff_props_homog_m
    implicit none
    save
    type, public, extends(diff_props_homog_c) :: tpt_props_homog_c ! homogeneous transport properties subclass
        ! physical properties
        real(kind=8) :: flux ! q
        real(kind=8) :: velocity ! v=q/phi
    contains
        procedure, public :: set_tpt_props_homog_c
        procedure, public :: read_props=>read_tpt_props_homog
        procedure, public :: get_props=>get_props_tpt_homog
        !procedure, public :: get_dispersion=>get_dispersion_tpt_homog
        procedure, public :: compute_retardo=>compute_retardo_homog
    end type
    
    contains
        subroutine set_tpt_props_homog_c(this,porosity,dispersion,flux)
            implicit none
            class(tpt_props_homog_c) :: this
            !real(kind=8), intent(in) :: source_term,retardo
            real(kind=8), intent(in) :: porosity,dispersion
            real(kind=8), intent(in) :: flux
            !this%source_term=source_term
            !this%source_term_flag=1 ! by default
            !this%retardo=retardo
            this%porosity=porosity
            this%flux=flux
            this%velocity=this%flux/this%porosity
            this%dispersion=dispersion
        end subroutine
        
        subroutine read_tpt_props_homog(this,filename,spatial_discr)
            implicit none
            class(tpt_props_homog_c) :: this
            character(len=*), intent(in) :: filename
            class(spatial_discr_c), intent(in) :: spatial_discr
            
            open(unit=1,file=filename,status='old',action='read')
            !read(1,"(/,F10.2)") this%source_term
            !this%source_term_flag=1 ! by default
            !read(1,*) this%retardo
            read(1,"(/,F10.2)") this%porosity
            read(1,*) this%dispersion
            read(1,*) this%flux
            !this%velocity=this%flux/this%porosity
            close(1)
        end subroutine
        
        function get_props_tpt_homog(this,i) result(props)
            implicit none
            class(tpt_props_homog_c) :: this
            integer(kind=4), intent(in), optional :: i
            real(kind=8), allocatable :: props(:,:)
            allocate(props(1,5))
            props(1,:)=[this%porosity,this%retardo,this%dispersion,this%flux,this%velocity]
        end function
        
        !function get_dispersion_tpt_homog(this) result(dispersion)
        !    implicit none
        !    class(tpt_props_homog_c) :: this
        !    real(kind=8), allocatable :: dispersion(:)
        !    allocate(dispersion(1))
        !    dispersion=this%dispersion
        !end function
        
         subroutine compute_retardo_homog(this,filename)
            implicit none
            class(tpt_props_homog_c) :: this
            character(len=*), intent(in) :: filename
            
            real(kind=8) :: rho_d
            real(kind=8) :: K_d
            
            open(unit=1,file=filename,status='old',action='read')
            read(1,"(/,F10.5)") rho_d
            read(1,*) K_d
            close(1)
            this%retardo=1d0+rho_d*K_d/this%porosity
            print *, this%retardo
         end subroutine
         

end module