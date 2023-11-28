module spatial_discr_1D_m
    use spatial_discr_m
    implicit none
    save
    type, public, extends(spatial_discr_c) :: mesh_1D_Euler_homog_c
        real(kind=8) :: Delta_x ! Mesh size
    contains
        procedure, public :: set_Delta_x_homog
        procedure, public :: read_mesh=>read_mesh_homog
        procedure, public :: get_mesh_size=>get_Delta_x_homog
        procedure, public :: compute_measure=>compute_measure_homog
        procedure, public :: compute_Delta_x
        procedure, private :: compute_Num_targets
        procedure, public :: get_dim=>get_dim_1D_homog
        procedure, public :: refine_mesh=>refine_mesh_homog
    end type
    
    type, public, extends(spatial_discr_c) :: mesh_1D_Euler_heterog_c
        real(kind=8), allocatable :: Delta_x(:) ! Target sizes
    contains
        procedure, public :: set_Delta_x_heterog
        procedure, public :: read_mesh=>read_mesh_heterog
        procedure, public :: get_mesh_size=>get_Delta_x_heterog
        procedure, public :: compute_measure=>compute_measure_heterog
        procedure, public :: get_dim=>get_dim_1D_heterog
        procedure, public :: refine_mesh=>refine_mesh_heterog
    end type
    
    interface
        subroutine refine_mesh_homog(this,conc,conc_ext,rel_tol)
            import mesh_1D_Euler_homog_c
            implicit none
            class(mesh_1D_Euler_homog_c) :: this
            real(kind=8), intent(inout), allocatable :: conc(:,:) ! Num_columns=Num_targets
            real(kind=8), intent(inout), allocatable :: conc_ext(:,:) ! Num_columns=Num_targets
            real(kind=8), intent(in) :: rel_tol ! relative tolerance
            !integer(kind=4), intent(out) :: n_new
        end subroutine
        
        subroutine refine_mesh_heterog(this,conc,conc_ext,rel_tol)
            import mesh_1D_Euler_heterog_c
            implicit none
            class(mesh_1D_Euler_heterog_c) :: this
            real(kind=8), intent(inout), allocatable :: conc(:,:) ! Num_columns=Num_targets
            real(kind=8), intent(inout), allocatable :: conc_ext(:,:) ! Num_columns=Num_targets
            real(kind=8), intent(in) :: rel_tol ! relative tolerance
            !integer(kind=4), intent(out) :: n_new
        end subroutine
    end interface
    
    contains
        subroutine set_Delta_x_homog(this,Delta_x)
            implicit none
            class(mesh_1D_Euler_homog_c) :: this
            real(kind=8), intent(in) :: Delta_x
            this%Delta_x=Delta_x
        end subroutine
        
        subroutine set_Delta_x_heterog(this,Delta_x)
            implicit none
            class(mesh_1D_Euler_heterog_c) :: this
            real(kind=8), intent(in) :: Delta_x(:)
            this%Delta_x=Delta_x
        end subroutine
        
        subroutine read_mesh_homog(this,filename)
            implicit none
            class(mesh_1D_Euler_homog_c) :: this
            character(len=*), intent(in) :: filename
            open(unit=1,file=filename,status='old',action='read')
            read(1,*) this%scheme
            read(1,*) this%targets_flag
            read(1,*) this%measure
            read(1,*) this%Num_targets
            read(1,*) this%adapt_ref
            close(1)
            this%Num_targets_defined=.true.
            call this%compute_Delta_x()
        end subroutine
        
        subroutine read_mesh_heterog(this,filename)
            implicit none
            class(mesh_1D_Euler_heterog_c) :: this
            character(len=*), intent(in) :: filename
            open(unit=1,file=filename,status='old',action='read')
            read(1,"(/,I10)") this%scheme
            read(1,*) this%targets_flag
            read(1,*) this%Num_targets
            allocate(this%Delta_x(this%Num_targets))
            read(1,*) this%Delta_x
            !read(1,*) this%adapt_ref
            close(1)
            call this%compute_measure()
        end subroutine
        
        !subroutine read_mesh_1D(this,filename)
        !    implicit none
        !    class(spatial_discr_c) :: this
        !    character(len=*), intent(in) :: filename
        !    
        !    real(kind=8), allocatable :: Delta_x(:)
        !    
        !    open(unit=1,file=filename,status='old',action='read')
        !    read(1,*) this%Num_targets
        !    select type (this)
        !    type is (mesh_1D_Euler_heterog_c)
        !        allocate(Delta_x(this%Num_targets)) ! size(Delta_x_vec)=Num_targets
        !        read(1,*) Delta_x
        !        allocate(this%Delta_x(this%Num_targets)) 
        !        this%Delta_x=Delta_x
        !    end select
        !    close(1)
        !    this%Num_targets_defined=.true.
        !end subroutine 
        !
        !subroutine set_mesh_1D(this,Delta_x,Num_targets)
        !    implicit none
        !    class(spatial_discr_c) :: this
        !    real(kind=8), intent(in) :: Delta_x
        !    integer(kind=4), intent(in), optional :: Num_targets
        !    select type (this)
        !    type is (mesh_1D_Euler_homog_c)
        !        this%Delta_x=Delta_x
        !        if (this%Num_targets_defined==.false. .and. present(Num_targets)) then
        !            this%Num_targets=Num_targets
        !            this%Num_targets_defined=.true.
        !        else if (this%Num_targets_defined==.true. .and. present(Num_targets)) then
        !            error stop "Num_targets already defined"
        !        else if (this%Num_targets==.false. .and. (.not. present(Num_targets))) then
        !            error stop "Num_targets missing"
        !        end if
        !    end select
        !end subroutine
        !
        !subroutine set_mesh_1D_Euler_homog(this,Delta_x,Num_targets)
        !    implicit none
        !    class(mesh_1D_Euler_homog_c) :: this
        !    real(kind=8), intent(in) :: Delta_x
        !    integer(kind=4), intent(in), optional :: Num_targets
        !    this%Delta_x=Delta_x
        !    if (this%Num_targets_defined==.false. .and. present(Num_targets)) then
        !        this%Num_targets=Num_targets
        !        this%Num_targets_defined=.true.
        !    else if (this%Num_targets_defined==.true. .and. present(Num_targets)) then
        !        error stop "Num_targets already defined"
        !    else if (this%Num_targets==.false. .and. (.not. present(Num_targets))) then
        !        error stop "Num_targets missing"
        !    end if
        !end subroutine
        !
        !subroutine set_mesh_1D_Euler_heterog(this,Delta_x,Num_targets)
        !    implicit none
        !    class(mesh_1D_Euler_heterog_c) :: this
        !    real(kind=8), intent(in) :: Delta_x(:)
        !    integer(kind=4), intent(in), optional :: Num_targets
        !    this%Delta_x=Delta_x
        !    if (this%Num_targets_defined==.true. .and. size(Delta_x)/=this%Num_targets) then
        !        error stop "Dimension error in heterogeneous mesh"
        !    else if (this%Num_targets_defined==.false.) then
        !        this%Num_targets=size(Delta_x)
        !        this%Num_targets_defined=.true.
        !    end if
        !end subroutine 
        !
        !!subroutine set_var_Delta_x(this,file_Delta_x)
        !!    implicit none
        !!    class(mesh_transport_1D) :: this
        !!    character(len=*), intent(in) :: file_Delta_x
        !!    
        !!    integer(kind=4) :: i,j,Num_elements
        !!    real(kind=8), allocatable :: Delta_x_vec(:)
        !!    !integer(kind=4), allocatable :: n_vec(:)
        !!    ! We assume file_Delta_x contains element sizes
        !!    select type (this)
        !!    type is (heterog_mesh_transport_1D)
        !!        open(unit=2,file=file_Delta_x,status='old',action='read')
        !!        read(2,*) Num_elements
        !!        allocate(Delta_x_vec(Num_elements)) ! size(Delta_x_vec)=Num_elements
        !!        read(2,*) Delta_x_vec
        !!        allocate(this%Delta_x(Num_elements)) 
        !!        this%Delta_x=Delta_x_vec
        !!        !print *, this%Delta_x
        !!        close(2)
        !!    class default
        !!        error stop "Wrong subclass"
        !!    end select
        !!end subroutine set_var_Delta_x
        !
        !subroutine set_measure(this,measure)
        !    implicit none
        !    class(spatial_discr_c) :: this
        !    real(kind=8), intent(in) :: measure
        !    this%measure=measure
        !    !select type (this)
        !    !type is (homog_mesh_transport_1D)
        !    !    this%measure=this%Num_elements*this%Delta_x
        !    !type is (heterog_mesh_transport_1D)
        !    !    this%measure=sum(this%Delta_x)
        !    !end select
        !end subroutine set_measure
        !
        !function get_Num_targets(this)
        !    implicit none
        !    class(spatial_discr_c) :: this
        !    integer(kind=4) :: get_Num_targets
        !    get_Num_targets=this%Num_targets
        !end function 
        !
        !function get_Delta_x(this) result(Delta_x)
        !    implicit none
        !    class(spatial_discr_c) :: this
        !    real(kind=8), allocatable :: Delta_x(:)
        !    select type (this)
        !    type is (mesh_1D_Euler_homog_c)
        !        allocate(Delta_x(1))
        !        Delta_x=get_Delta_x_homog(this)
        !    type is (mesh_1D_Euler_heterog_c)
        !        Delta_x=get_Delta_x_heterog(this)
        !    end select
        !end function
        !
        function get_Delta_x_homog(this,i) result(Delta_x)
            implicit none
            class(mesh_1D_Euler_homog_c) :: this
            integer(kind=4), intent(in), optional :: i
            real(kind=8) :: Delta_x
            Delta_x=this%Delta_x
        end function
        
        function get_Delta_x_heterog(this,i) result(Delta_x)
            implicit none
            class(mesh_1D_Euler_heterog_c) :: this
            integer(kind=4), intent(in), optional :: i
            real(kind=8) :: Delta_x
            if (present(i)) then
                Delta_x=this%Delta_x(i)
            else
                error stop "Heterogeneous mesh"
            end if
        end function
        
        subroutine compute_measure_homog(this)
            implicit none
            class(mesh_1D_Euler_homog_c) :: this
            this%measure=(this%Num_targets-this%targets_flag)*this%Delta_x
        end subroutine
        
        subroutine compute_measure_heterog(this)
            implicit none
            class(mesh_1D_Euler_heterog_c) :: this
            this%measure=sum(this%Delta_x)
        end subroutine
        
        subroutine compute_Delta_x(this)
            implicit none
            class(mesh_1D_Euler_homog_c) :: this
            this%Delta_x=this%measure/(this%Num_targets-this%targets_flag)
        end subroutine
        
        subroutine compute_Num_targets(this)
            implicit none
            class(mesh_1D_Euler_homog_c) :: this
            this%Num_targets=this%measure/this%Delta_x + this%targets_flag
        end subroutine
        
        function get_dim_1D_homog(this) result(dim)
            implicit none
            class(mesh_1D_Euler_homog_c) :: this
            integer(kind=4) :: dim
            dim=1
        end function
        
        function get_dim_1D_heterog(this) result(dim)
            implicit none
            class(mesh_1D_Euler_heterog_c) :: this
            integer(kind=4) :: dim
            dim=1
        end function
        
        
end module
