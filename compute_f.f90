function compute_f(this,k) result(f)
! f=Delta_t(k)*inv(F)*g
    use transport_transient_m
    implicit none
    class(PDE_1D_transient_c) :: this
    integer(kind=4), intent(in), optional :: k
    real(kind=8), allocatable :: f(:)
    
    f=this%source_term_PDE
    select type (time=>this%time_discr)
    type is (time_discr_homog_c)
        f=f*time%Delta_t/this%F_mat%diag
    type is (time_discr_heterog_c)
        f=f*time%Delta_t(k)/this%F_mat%diag
    end select
end function