subroutine compute_E_mat(this,E_mat,k)
! E=Delta_t(k)*inv(F)*T (tridiagonal, negative semi-definite)
! rows sum = 0 if r=0
    use PDE_transient_m
    implicit none
    
    class(PDE_1D_transient_c), intent(in) :: this
    type(tridiag_matrix_c), intent(out) :: E_mat
    integer(kind=4), intent(in), optional :: k
    
    integer(kind=4) :: j,n
    !type(tridiag_sym_matrix_c), target :: E_mat_sym
    !type(tridiag_matrix_c), target :: E_mat_non_sym
    
    n=this%spatial_discr%Num_targets
    
    !select type (trans_mat=>this%trans_mat)
    !type is (tridiag_sym_matrix_c)
        !print *, this%trans_mat%diag
    !print *, this%props%get_props()
    if (this%spatial_discr%adapt_ref==1) then
        call this%compute_trans_mat_PDE()
        call this%compute_F_mat_PDE()
    end if
        E_mat%sub=this%trans_mat%sub
        E_mat%diag=this%trans_mat%diag
        E_mat%super=this%trans_mat%super
        select type (time_discr=>this%time_discr)
        type is (time_discr_homog_c)
            E_mat%sub=E_mat%sub*time_discr%Delta_t
            E_mat%diag=E_mat%diag*time_discr%Delta_t
            E_mat%super=E_mat%super*time_discr%Delta_t
        type is (time_discr_heterog_c)
            E_mat%sub=E_mat%sub*time_discr%Delta_t(k)
            E_mat%diag=E_mat%diag*time_discr%Delta_t(k)
            E_mat%super=E_mat%super*time_discr%Delta_t(k)
        end select
        !print *, this%F_mat%diag
        do j=1,n-1
            E_mat%super(j)=E_mat%super(j)/this%F_mat%diag(j)
            E_mat%sub(j)=E_mat%sub(j)/this%F_mat%diag(j+1)
            E_mat%diag(j)=E_mat%diag(j)/this%F_mat%diag(j)
        end do
        E_mat%diag(n)=E_mat%diag(n)/this%F_mat%diag(n)
        !E_mat=>E_mat_sym
        !print *, E_mat%diag
        !E_mat%diag=E_mat_sym%diag
        !E_mat%sub=E_mat_sym%sub
    !type is (tridiag_matrix_c)
    !    E_mat_non_sym%sub=trans_mat%sub
    !    E_mat_non_sym%diag=trans_mat%diag
    !    E_mat_non_sym%super=trans_mat%super
    !    select type (time_discr=>this%time_discr)
    !    type is (time_discr_homog_c)
    !        E_mat_non_sym%sub=E_mat_non_sym%sub*time_discr%Delta_t
    !        E_mat_non_sym%diag=E_mat_non_sym%diag*time_discr%Delta_t
    !        E_mat_non_sym%super=E_mat_non_sym%super*time_discr%Delta_t
    !    type is (time_discr_heterog_c)
    !        E_mat_non_sym%sub=E_mat_non_sym%sub*time_discr%Delta_t(k)
    !        E_mat_non_sym%diag=E_mat_non_sym%diag*time_discr%Delta_t(k)
    !        E_mat_non_sym%super=E_mat_non_sym%super*time_discr%Delta_t(k)
    !    end select
    !    E_mat_non_sym%diag=E_mat_non_sym%diag/this%F_mat%diag
    !    E_mat=>E_mat_non_sym
    !end select
    !print *, E_mat%get_diag()
    !E_mat_sym%sub=this%trans_mat%sub
    !E_mat_sym%diag=this%trans_mat%diag
    !select type (time_discr=>this%time_discr)
    !type is (time_discr_homog_c)
    !    E_mat_sym%sub=E_mat_sym%sub*time_discr%Delta_t
    !    E_mat_sym%diag=E_mat_sym%diag*time_discr%Delta_t
    !type is (time_discr_heterog_c)
    !    E_mat_sym%sub=E_mat_sym%sub*time_discr%Delta_t(k)
    !    E_mat_sym%diag=E_mat_sym%diag*time_discr%Delta_t(k)
    !end select
    !E_mat=>E_mat_sym
    !select type (trans_mat=>this%trans_mat)
    !type is (tridiag_matrix_c)
    !    E_mat_non_sym%super=trans_mat%super
    !    select type (time_discr=>this%time_discr)
    !    type is (time_discr_homog_c)
    !        E_mat_non_sym%super=E_mat_non_sym%super*time_discr%Delta_t
    !    type is (time_discr_heterog_c)
    !        E_mat_non_sym%super=E_mat_non_sym%super*time_discr%Delta_t(k)
    !    end select
    !    
    !end select
    
end subroutine 