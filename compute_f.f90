function compute_f(this,k) result(f)
! f=Delta_t(k)*inv(F)*g
    use transport_transient_m
    implicit none
    class(PDE_1D_transient_c) :: this
    integer(kind=4), intent(in), optional :: k
    real(kind=8), allocatable :: f(:)
    
    !select type (this)
    !class is (diffusion_1D_transient_c)
        f=this%source_term_PDE
        select type (time=>this%time_discr)
        type is (time_discr_homog_c)
            f=f*time%Delta_t/this%F_mat%diag
        type is (time_discr_heterog_c)
            f=f*time%Delta_t(k)/this%F_mat%diag
        end select
    !end select
    !select type (this)
    !type is (transport_1D_transient_c)
    !    select type (time=>this%time_discr)
    !    type is (time_discr_homog_c)
    !        Q_mat=Q_mat*time%Delta_t/this%compute_F_mat_tpt()
    !    type is (time_discr_heterog_c)
    !        Q_mat=Q_mat*time%Delta_t(k)/this%compute_F_mat_tpt()
    !    end select
    !type is (diffusion_transient_c)
        
    !end select
end function