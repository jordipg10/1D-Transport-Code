! Computes mixing ratios matrix
subroutine compute_mixing_ratios(this,theta,k)
    use BCs_subroutines_m
    !use properties_viena_m
    implicit none
    
    class(diffusion_1D_transient_c) :: this
    real(kind=8), intent(in) :: theta
    !class(tridiag_matrix_c), intent(in) :: B_mat
    integer(kind=4), intent(in), optional :: k
    
    integer(kind=4) :: i,n
    real(kind=8) :: lambda
    real(kind=8), allocatable :: r(:)
    
    !type(tridiag_matrix_vec_c) :: Lambda_mat
    type(tridiag_matrix_c) :: B_mat
    
    n=this%spatial_discr%Num_targets
    
    call this%mixing_ratios%allocate_matrix(n)

    !allocate(this%mixing_ratios(n,2*n),E_sub(n-1),E_diag(n),E_super(n-1))
    
    !select type (this)
    !type is (diffusion_c)
    !    error stop "Wrong subclass"
    !type is (transport_1D_c)
    !    error stop "Wrong subclass"
    !end select
    !select type (this)
    !class is (diffusion_transient_c)
    !if (present(theta)) then
    !    theta_star=theta
    !else
    !    theta_star=0d0
    !end if
    call this%compute_B_mat(theta,B_mat,k)
    this%mixing_ratios%sub=B_mat%sub
    this%mixing_ratios%diag=B_mat%diag
    this%mixing_ratios%super=B_mat%super
    !select type (E_mat)
    !type is (tridiag_matrix_c)
    !    !Lambda_mat%sub=(1d0-theta)*E_mat%get_sub()
    !    !Lambda_mat%diag=(1d0-theta)*E_mat%get_diag()
    !    !Lambda_mat%super=(1d0-theta)*E_mat%get_super()
    !    Lambda_mat%super=(1d0-theta)*E_mat%super
    !end select
    !allocate(this%mixing_ratios%vector(n))
    !r=this%props%get_source_term()
    
    !if (present(props)) then
    !    !print *, props%source_term
    !    this%mixing_ratios%vector=this%conc_star_flag*props%source_term
    !end if
    select type (time=>this%time_discr)
    type is (time_discr_homog_c)
        this%mixing_ratios%vector=this%mixing_ratios%vector*time%Delta_t/this%F_mat%diag
    type is (time_discr_heterog_c)
        this%mixing_ratios%vector=this%mixing_ratios%vector*time%Delta_t(k)/this%F_mat%diag
    end select
    !print *, this%mixing_ratios%vector
    
    !select type (props=>this%props)
    !type is (props_viena_c)
    !    select type (time=>this%time_discr)
    !    type is (time_discr_homog_c)
    !        ! Euler fully implicit
    !        lambda=(props%caudal/props%volume)/((props%porosity/time%Delta_t)+props%caudal/props%volume)
    !        this%mixing_ratios%sub=0d0
    !        this%mixing_ratios%diag=1d0-lambda
    !        this%mixing_ratios%super=0d0
    !        this%mixing_ratios%vector=lambda
    !    end select
    !end select
end subroutine 