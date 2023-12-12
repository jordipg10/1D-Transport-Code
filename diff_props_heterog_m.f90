module diff_props_heterog_m
    use properties_m
    implicit none
    save
    type, public, extends(props_c) :: diff_props_heterog_c ! heterogeneous diffusion properties subclass
        ! physical properties
        real(kind=8), allocatable :: porosity(:)   ! phi
        real(kind=8), allocatable :: dispersion(:) ! D
    contains
        procedure, public :: set_props_diff_heterog
        procedure, public :: read_props=>read_props_diff_heterog
        !procedure, public :: compute_retardo=>compute_retardo_heterog
        !procedure, public :: get_dispersion
    end type
    
    contains
        subroutine set_props_diff_heterog(this,porosity,dispersion)
            implicit none
            class(diff_props_heterog_c) :: this
            !real(kind=8), intent(in) :: source_term,retardo
            real(kind=8), intent(in) :: porosity(:),dispersion(:)
            !this%source_term=source_term
            !this%retardo=retardo
            if (size(porosity)/=size(dispersion)) error stop "Dimensions of porosity and dispersion must be the same"
            this%porosity=porosity
            this%dispersion=dispersion
        end subroutine
        
        subroutine read_props_diff_heterog(this,filename,spatial_discr)
            implicit none
            class(diff_props_heterog_c) :: this
            character(len=*), intent(in) :: filename
            class(spatial_discr_c), intent(in), optional :: spatial_discr
            
            real(kind=8), parameter :: epsilon=1d-12
            real(kind=8) :: phi,D,r
            integer(kind=4) :: flag
            
            open(unit=1,file=filename,status='old',action='read')
            read(1,"(/,F10.2)") flag
            if (flag==.true.) then
                backspace(1)
                read(1,*) flag, r
                allocate(this%source_term(spatial_discr%Num_targets-spatial_discr%targets_flag))
                this%source_term=r
                this%source_term_order=0
            else
                read(1,*) this%source_term
                if (size(this%source_term)/=spatial_discr%Num_targets-spatial_discr%targets_flag) error stop "Dimension error in source term"
            end if
            read(1,*) flag
            if (flag==.true.) then
                backspace(1)
                read(1,*) flag, phi
                allocate(this%porosity(spatial_discr%Num_targets-spatial_discr%targets_flag))
                this%porosity=phi
            end if
            read(1,*) flag
            if (flag==.true.) then
                backspace(1)
                read(1,*) flag, D
                allocate(this%dispersion(spatial_discr%Num_targets-spatial_discr%targets_flag))
                this%dispersion=D
            end if
            close(1)
        end subroutine
        
      
        
        
end module