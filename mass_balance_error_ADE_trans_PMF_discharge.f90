function mass_balance_error_ADE_trans_PMF_discharge(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
    use transport_transient_m
    implicit none
    class(transport_1D_transient_c), intent(in) :: this
    real(kind=8), intent(in) :: conc_old(:)
    real(kind=8), intent(in) :: conc_new(:)
    !real(kind=8), intent(in) :: q_inf
    !real(kind=8), intent(in) :: q_out
    real(kind=8), intent(in) :: Delta_t
    real(kind=8), intent(in) :: Delta_x
    !integer(kind=4), intent(in) :: k
    real(kind=8) :: mass_bal_err
    
    real(kind=8) :: sum_time,sum_r,inf_term,out_term,dummy_inf,D,q_1
    integer(kind=4) :: i,n
    
    n=size(conc_old)
    
    sum_time=0d0
    
    do i=1,n
        sum_time=sum_time+this%tpt_props_heterog%porosity(i)*(conc_new(i)-conc_old(i))
    end do
    sum_time=sum_time*Delta_x
    
    sum_r=0d0
    do i=1,n
        sum_r=sum_r+this%tpt_props_heterog%source_term(i)*conc_old(i)
    end do
    sum_r=sum_r*Delta_x*Delta_t
    
    !inf_term=conc_old(1)*(this%BCs%flux_inf-this%tpt_props_heterog%flux(1)-2d0*this%tpt_props_heterog%dispersion(1)/Delta_x)+this%BCs%conc_inf*(this%tpt_props_heterog%flux(1)+2d0*this%tpt_props_heterog%dispersion(1)/Delta_x)
    !inf_term = 2d0*this%tpt_props_heterog%dispersion(1)*(this%BCs%conc_inf-conc_old(1))/Delta_x + this%tpt_props_heterog%flux(1)*(this%BCs%conc_inf-conc_old(1)) + this%BCs%flux_inf*conc_old(1)
    !out_term = 2d0*this%tpt_props_heterog%dispersion(n)*(this%BCs%conc_out-conc_old(n))/Delta_x - this%tpt_props_heterog%flux(n)*(this%BCs%conc_out-conc_old(n)) - this%BCs%flux_out*conc_old(n)
    !out_term=conc_old(n)*(this%tpt_props_heterog%flux(n)-this%BCs%flux_out-2d0*this%tpt_props_heterog%dispersion(n)/Delta_x)-this%BCs%conc_out*(this%tpt_props_heterog%flux(n)-2d0*this%tpt_props_heterog%dispersion(n)/Delta_x)
    
    inf_term = this%BCs%flux_inf*this%BCs%conc_inf
    out_term = - this%BCs%flux_out*conc_old(n)
    
    inf_term=inf_term*Delta_t
    out_term=out_term*Delta_t
    
    mass_bal_err=sum_time-sum_r-inf_term-out_term
end function