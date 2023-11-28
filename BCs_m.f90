module BCs_m
    implicit none
    save
    type, public :: BCs_t ! Boundary conditions
        integer(kind=4) :: BCs_label    ! 1: Dirichlet constant
                                        ! 2: Neumann homogeneous
                                        ! 3: Dirichlet analytical
                                        ! 4: Cauchy
        integer(kind=4) :: option   ! 1: 1 and 0
                                    ! 2: symmetry in T
        logical :: evap ! evaporation
        !real(kind=8), external :: BCs_fct
        real(kind=8), allocatable :: Dirichlet_BCs(:)
        integer(kind=4) :: Dirichlet_BC_location ! 0: center
                                                 ! 1: perimeter
    contains
        procedure, public :: set_BCs_label
        procedure, public :: set_option
        procedure, public :: set_evap
        procedure, public :: set_Dirichlet_BCs
        procedure, public :: set_Dirichlet_BC_location
        procedure, public :: read_BCs
        procedure, public :: read_Dirichlet_BCs
    end type
    
    contains
        subroutine set_BCs_label(this,BCs)
            implicit none
            class(BCs_t) :: this
            integer(kind=4), intent(in) :: BCs
            if (BCs>3) error stop "BCs not implemented yet"
            this%BCs_label=BCs
        end subroutine
        
        subroutine set_option(this,option)
            implicit none
            class(BCs_t) :: this
            integer(kind=4), intent(in) :: option
            if (option>2) error stop "BCs option not implemented yet"
            this%option=option
        end subroutine 
        
        subroutine set_evap(this,evap)
            implicit none
            class(BCs_t) :: this
            logical, intent(in) :: evap
            this%evap=evap
        end subroutine 
        
        subroutine set_Dirichlet_BCs(this,Dirichlet_BCs)
            implicit none
            class(BCs_t) :: this
            real(kind=8), intent(in) :: Dirichlet_BCs(:)
            this%Dirichlet_BCs=Dirichlet_BCs
        end subroutine
        
        subroutine set_Dirichlet_BC_location(this,Dirichlet_BC_location)
            implicit none
            class(BCs_t) :: this
            integer(kind=4), intent(in) :: Dirichlet_BC_location
            this%Dirichlet_BC_location=Dirichlet_BC_location
        end subroutine
        
        subroutine read_BCs(this,filename)
            implicit none
            class(BCs_t) :: this
            character(len=*), intent(in) :: filename
            open(unit=1,file=filename,status='old',action='read')
            read(1,"(/,I10)") this%BCs_label
            read(1,"(L10)") this%evap
            read(1,*) this%option
            close(1)
        end subroutine
        
        subroutine read_Dirichlet_BCs(this,filename)
            implicit none
            class(BCs_t) :: this
            character(len=*), intent(in) :: filename
            allocate(this%Dirichlet_BCs(2))
            open(unit=1,file=filename,status='old',action='read')
            read(1,"(/,F10.2)") this%Dirichlet_BCs(1)
            read(1,*) this%Dirichlet_BCs(2)
            close(1)
        end subroutine
end module