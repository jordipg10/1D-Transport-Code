! Diffusion equation:
module diffusion_m
    use PDE_m
    use diff_props_heterog_m
    implicit none
    save
    type, public, extends(PDE_1D_c) :: diffusion_1D_c ! 1D diffusion equation class
        real(kind=8), allocatable :: conc(:) ! concentrations (c)
        real(kind=8), allocatable :: conc_ext(:) ! (c_e)
        type(diff_props_heterog_c) :: diff_props_heterog
    contains
    ! Set
        procedure, public :: set_conc_ext
        procedure, public :: set_diff_props_heterog
    ! Get

    ! Computations
        procedure, public :: compute_trans_mat_PDE=>compute_trans_mat_diff
        procedure, public :: initialise_PDE=>initialise_diffusion_1D
        procedure, public :: write_PDE_1D=>write_diffusion_1D
    end type
!*****************************************************************************************************************************
    interface
        
        subroutine initialise_diffusion_1D(this)
            import diffusion_1D_c
            class(diffusion_1D_c) :: this
        end subroutine
        
        subroutine compute_trans_mat_diff(this)
            import diffusion_1D_c
            implicit none
            class(diffusion_1D_c) :: this
        end subroutine
        
        
       subroutine write_diffusion_1D(this,Time_out,output)
            import diffusion_1D_c
            import props_c
            class(diffusion_1D_c), intent(in) :: this
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(in) :: output(:,:)
            !class(props_c), intent(in), optional :: props
            !integer(kind=4), intent(in), optional :: opcion
        end subroutine
      
      
        
        
        
    end interface
!*****************************************************************************************************************************
    contains
        
        
      
        
        subroutine set_conc_ext(this,conc_ext)
            class(diffusion_1D_c) :: this
            real(kind=8), intent(in) :: conc_ext(:)
            if (size(conc_ext)/=this%spatial_discr%Num_targets) error stop "Dimension error in external concentration"
            this%conc_ext=conc_ext
        end subroutine 
        
       
        
        subroutine set_diff_props_heterog(this,diff_props_heterog)
            implicit none
            class(diffusion_1D_c) :: this
            class(diff_props_heterog_c), intent(in) :: diff_props_heterog
            this%diff_props_heterog=diff_props_heterog
        end subroutine
        
       
        
      
        
      
end module 
