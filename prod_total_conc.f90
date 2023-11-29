! Computes concentrations using eigendecomposition
subroutine prod_total_conc(this,A_mat,time)
    use diffusion_transient_m
    implicit none
    class(diffusion_1D_transient_c) :: this
    class(tridiag_matrix_c), intent(in) :: A_mat
    !real(kind=8), intent(in) :: Time_out(:)
    !real(kind=8), intent(out) :: output(:,:)
    real(kind=8), intent(in), optional :: time
    !real(kind=8), intent(in) :: y(:)
    !integer(kind=4), intent(in), optional :: k
    
    real(kind=8), allocatable :: y0(:),b(:),c0(:),g(:)
    real(kind=8) :: sumj,sumk1,sumk2,sumk,t,conc_lim_n
    integer(kind=4) :: k,i,j,n
    !type(tridiag_matrix_c) :: A_mat
    n=this%spatial_discr%Num_targets
    
    !if (present(k)) then
    !    select type (time=>this%time_discr)
    !    type is (time_discr_homog_c)
    !        t=k*time%Delta_t
    !    type is (time_discr_heterog_c)
    !        t=sum(time%Delta_t(1:k))
    !    end select
    !else
    !    t=this%time_discr%Final_time
    !    call this%compute_A_mat_ODE(A_mat)
    !    y0=sqrt(this%F_mat%diag)*this%conc_init
    !    b=this%source_term_PDE/sqrt(this%F_mat%diag)
    !    this%conc=prod_total_sym_mat_bis(A_mat,y0,b,t)/sqrt(this%F_mat%diag)
    !end if
    
    !call this%compute_A_mat_ODE(A_mat)
    
    y0=sqrt(this%F_mat%diag)*this%conc_init
    b=this%source_term_PDE/sqrt(this%F_mat%diag)
    
    if (.not. allocated(this%conc)) then
        allocate(this%conc(n))
    end if
    !call A_mat%compute_eigenvalues_eigenvectors()
    !print *, A_mat%eigenvalues
    !allocate(A_mat%eigenvectors(this%spatial_discr%Num_targets,this%spatial_discr%Num_targets))
    !call eigenvectors_tridiag_sym_matrix(A_mat%diag,A_mat%sub,A_mat%eigenvalues,A_mat%eigenvectors)
    !call A_mat%check_eigenvectors_tridiag_sym_matrix(tol)
    if (present(time)) then
        do i=1,n
            sumj=0d0
            do j=1,n
                sumk1=0d0
                do k=1,n
                    sumk1=sumk1+A_mat%eigenvectors(k,j)*sqrt(this%F_mat%diag(k))*this%conc_init(k)
                end do
                sumk2=0d0
                do k=1,n
                    sumk2=sumk2+A_mat%eigenvectors(k,j)*this%source_term_PDE(k)/sqrt(this%F_mat%diag(k))
                end do
                sumj=sumj+A_mat%eigenvectors(i,j)*(exp(-A_mat%eigenvalues(j)*time)*sumk1+((1d0-exp(-A_mat%eigenvalues(j)*time))/A_mat%eigenvalues(j))*sumk2)
            end do
            this%conc(i)=sumj/sqrt(this%F_mat%diag(i))
        end do
    else
        do i=1,n
            sumj=0d0
            do j=1,n
                sumk=0d0
                do k=1,n
                    sumk=sumk+A_mat%eigenvectors(k,j)*this%source_term_PDE(k)/sqrt(this%F_mat%diag(k))
                end do
                sumj=sumj+A_mat%eigenvectors(i,j)*sumk/A_mat%eigenvalues(j)
            end do
            this%conc(i)=sumj/sqrt(this%F_mat%diag(i))
        end do
    end if
end subroutine