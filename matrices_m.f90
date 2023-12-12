module matrices_m
    implicit none
    save
    
    type, public :: matrix_c
    contains
        procedure, public :: allocate_matrix        
        procedure, public :: prod_mat_vec
        procedure, public :: prod_mat_mat
        procedure, public :: get_diag
        procedure, public :: get_sub
        procedure, public :: get_super
        procedure, public :: compute_norm_inf
        procedure, public :: compute_norm_1   
    end type
    
    type, public, extends(matrix_c) :: sq_matrix_c
        real(kind=8), allocatable :: eigenvalues(:)
        real(kind=8), allocatable :: eigenvectors(:,:)
    contains
        procedure, public :: compute_eigenvalues
        procedure, public :: compute_eigenvectors
    end type
!****************************************************************************************************************************************************
    type, public, extends(sq_matrix_c) :: tridiag_sym_Toeplitz_matrix_c
        real(kind=8) :: sub
        real(kind=8) :: diag
    end type
    
    type, public, extends(tridiag_sym_Toeplitz_matrix_c) :: tridiag_Toeplitz_matrix_c
        real(kind=8) :: super
    contains
        procedure, public :: set_tridiag_Toeplitz_matrix
    end type
!****************************************************************************************************************************************************
    type, public, extends(sq_matrix_c) :: diag_matrix_c
        real(kind=8), allocatable :: diag(:)
    contains
        procedure, public :: set_diag_matrix
    end type
    
    type, public, extends(diag_matrix_c) :: tridiag_sym_matrix_c
        real(kind=8), allocatable :: sub(:)
    contains
        procedure, public :: check_eigenvectors_tridiag_sym_matrix
    end type

    type, public, extends(tridiag_sym_matrix_c) :: tridiag_matrix_c
        real(kind=8), allocatable :: super(:)
    contains
        procedure, public :: set_tridiag_matrix
        procedure, public :: compute_transpose_tridiag_matrix
    end type

    type, public, extends(tridiag_matrix_c) :: tridiag_matrix_vec_c
        real(kind=8), allocatable :: vector(:)
    end type
!****************************************************************************************************************************************************
    interface
        subroutine compute_eigenvalues(this)
            import sq_matrix_c
            implicit none
            class(sq_matrix_c) :: this
        end subroutine
        
        subroutine compute_eigenvectors(this)
            import sq_matrix_c
            implicit none
            class(sq_matrix_c) :: this
        end subroutine
        
        function id_matrix(n)
        ! Identity nxn matrix
            implicit none
            integer(kind=4), intent(in) :: n
            integer(kind=4) :: i
            real(kind=8) :: id_matrix(n,n)
        end function
        
        function norm_mat_inf(A) result(norm)
            implicit none
            real(kind=8), intent(in) :: A(:,:)
            real(kind=8) :: norm
        end function
        
        function norm_mat_1(A) result(norm)
            implicit none
            real(kind=8), intent(in) :: A(:,:)
            real(kind=8) :: norm
        end function
        
        function det(A)
        ! Determinant of square matrix using LU decomposition
            implicit none
            real(kind=8), intent(in) :: A(:,:)
            real(kind=8) :: det
        end function det
        
        subroutine inv_matrix(A,inv)
        ! Inverse of square matrix using LU decomposition
            implicit none
            real(kind=8), intent(in) :: A(:,:)
            real(kind=8), intent(out) :: inv(:,:)
        end subroutine inv_matrix
        
        function prod_mat_vec(this,b) result(x) ! Ab=x
            import matrix_c
            implicit none
            class(matrix_c), intent(in) :: this
            real(kind=8), intent(in) :: b(:)
            real(kind=8), allocatable :: x(:)
        end function
        
        function prod_mat_mat(this,B_mat) result(C_mat) ! AB=C
            import matrix_c
            implicit none
            class(matrix_c), intent(in) :: this
            real(kind=8), intent(in) :: B_mat(:,:)
            real(kind=8), allocatable :: C_mat(:,:)
        end function 
        
        subroutine LU(A,L,U)
            implicit none
            real(kind=8), intent(in) :: A(:,:)
            real(kind=8), intent(out) :: L(:,:), U(:,:)
        end subroutine LU
        
        subroutine Cholesky(A,L)
            implicit none
            real(kind=8), intent(in) :: A(:,:)
            real(kind=8), intent(out) :: L(:,:)
        end subroutine Cholesky
        
        subroutine QR_Householder(A,Q,R)
            implicit none
            real(kind=8), intent(in) :: A(:,:) ! square matrix
            real(kind=8), intent(out) :: Q(:,:), R(:,:)
        end subroutine QR_Householder
        
        function Householder(x)
            implicit none
            real(kind=8), intent(in) :: x(:)
            real(kind=8), allocatable :: Householder(:,:)
        end function Householder
        
        subroutine Gram_Schmidt_mat(v,u)
            implicit none
            real(kind=8), intent(in) :: v(:,:)
            real(kind=8), intent(out) :: u(:,:)
        end subroutine
        
        subroutine non_zeros(A,n,ind)
            implicit none
            real(kind=8), intent(in) :: A(:,:)
            integer(kind=4), intent(out) :: n
            integer(kind=4), intent(out) :: ind(:,:)
        end subroutine
        
        function prod_AT_A(A) result(AT_A)
            implicit none
            real(kind=8), intent(in) :: A(:,:)
            real(kind=8), allocatable :: AT_A(:,:)
        end function
        
        subroutine check_eigenvectors_tridiag_sym_matrix(this,tolerance)
            import tridiag_sym_matrix_c
            implicit none
            class(tridiag_sym_matrix_c), intent(in) :: this
            real(kind=8), intent(in) :: tolerance
        end subroutine
        
        function prod_diag_tridiag_mat(A,B) result(C)
            import diag_matrix_c
            import tridiag_matrix_c
            class(diag_matrix_c), intent(in) :: A
            class(tridiag_matrix_c), intent(in) :: B
            type(tridiag_matrix_c) :: C
        end function
        
        function prod_tridiag_diag_mat(A,B) result(C)
            import diag_matrix_c
            import tridiag_matrix_c
            class(tridiag_matrix_c), intent(in) :: A
            class(diag_matrix_c), intent(in) :: B
            type(tridiag_matrix_c) :: C
        end function
    end interface
!****************************************************************************************************************************************************
    contains
        subroutine allocate_matrix(this,n)
            implicit none
            class(matrix_c) :: this
            integer(kind=4), intent(in) :: n
            select type (this)
            class is (diag_matrix_c)
                allocate(this%diag(n))
                select type (this)
                class is (tridiag_sym_matrix_c)
                    allocate(this%sub(n-1))
                    select type (this)
                    type is (tridiag_matrix_c)
                        allocate(this%super(n-1))
                        select type (this)
                        type is (tridiag_matrix_vec_c)
                            allocate(this%vector(n))
                        end select
                    end select
                end select
            end select
        end subroutine
        
        subroutine set_diag_matrix(this,diag)
            implicit none
            class(diag_matrix_c) :: this
            real(kind=8), intent(in) :: diag(:)
            this%diag=diag
        end subroutine
        
        subroutine set_tridiag_Toeplitz_matrix(this,sub,diag,super)
            implicit none
            class(tridiag_Toeplitz_matrix_c) :: this
            real(kind=8), intent(in) :: sub,diag,super
            this%sub=sub
            this%diag=diag
            this%super=super
        end subroutine
        
        subroutine set_tridiag_matrix(this,sub,diag,super)
            implicit none
            class(tridiag_matrix_c) :: this
            real(kind=8), intent(in) :: sub(:),diag(:),super(:)
            if (size(sub)/=size(super) .or. size(diag)/=(size(sub)+1)) error stop "Dimension error in set_tridiag_matrix" 
            this%sub=sub
            this%diag=diag
            this%super=super
        end subroutine
        
        function get_diag(this) result(diag)
            implicit none
            class(matrix_c), intent(in) :: this
            real(kind=8), allocatable :: diag(:)
            select type (this)
            class is (diag_matrix_c)
                diag=this%diag
            end select
        end function
        
        function get_sub(this) result(sub)
            implicit none
            class(matrix_c), intent(in) :: this
            real(kind=8), allocatable :: sub(:)
            select type (this)
            class is (tridiag_sym_matrix_c)
                sub=this%sub
            end select
        end function
        
        function get_super(this) result(super)
            implicit none
            class(matrix_c), intent(in) :: this
            real(kind=8), allocatable :: super(:)
            select type (this)
            class is (tridiag_matrix_c)
                super=this%super
            end select
        end function
        
        function compute_norm_inf(this) result(norm)
            implicit none
            class(matrix_c), intent(in) :: this
            real(kind=8) :: norm
            
            integer(kind=4) :: i,n
            real(kind=8) :: norm_i
            
            select type (this)
            type is (diag_matrix_c)
                n=size(this%diag)
                norm=abs(this%diag(1))
                do i=2,n
                    norm_i=abs(this%diag(1))
                    if (norm_i>norm) then
                        norm=norm_i
                    end if
                end do
            type is (tridiag_matrix_c)
                n=size(this%diag)
                norm=abs(this%diag(1))+abs(this%super(1))
                do i=2,n-1
                    norm_i=abs(this%sub(i-1))+abs(this%diag(i))+abs(this%super(i))
                    if (norm_i>norm) then
                        norm=norm_i
                    end if
                end do
                norm_i=abs(this%sub(n-1))+abs(this%diag(n))
                if (norm_i>norm) then
                    norm=norm_i
                end if
            end select
        end function
        
        function compute_norm_1(this) result(norm)
            implicit none
            class(matrix_c), intent(in) :: this
            real(kind=8) :: norm
            
            type(tridiag_matrix_c) :: transpose
            
            select type (this)
            type is (tridiag_matrix_c)
                call this%compute_transpose_tridiag_matrix(transpose)
                norm=transpose%compute_norm_inf()
            end select
        end function
        
        subroutine compute_transpose_tridiag_matrix(this,transpose)
            implicit none
            class(tridiag_matrix_c), intent(in) :: this
            class(tridiag_matrix_c), intent(out) :: transpose
            transpose%sub=this%super
            transpose%diag=this%diag
            transpose%super=this%sub
        end subroutine

end module