function mass_balance_error_ADE_CFD_PMF_quad_flow_recharge(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
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
    
    real(kind=8) :: sum_time,sum_r,sum_inner,inf_term,out_term,dummy_inf,D,q_1,mass_flux_inf,mass_flux_out,d2q_dx2
    integer(kind=4) :: i,n
    
    n=size(conc_old)
    
    d2q_dx2=2d0 ! q(x)=x^2+1
    
    q_1=this%tpt_props_heterog%flux(1)
    D=this%tpt_props_heterog%dispersion(1) ! constant disperision
    
    dummy_inf=conc_old(1)*(2d0*D-this%BCs%flux_inf*Delta_x)/(this%BCs%flux_inf*Delta_x+2d0*D) + this%BCs%conc_inf*2d0*this%BCs%flux_inf*Delta_x/(this%BCs%flux_inf*Delta_x+2d0*D)
    mass_flux_inf=D*(dummy_inf-conc_old(1))/Delta_x + this%BCs%flux_inf*(dummy_inf+conc_old(1))/2d0
    mass_flux_out=this%BCs%flux_out*conc_old(n)
    
    sum_time=0d0
    do i=1,n
        sum_time=sum_time+this%tpt_props_heterog%porosity(i)*(conc_new(i)-conc_old(i))
    end do
    sum_time=sum_time*Delta_x
    
    sum_r=0d0
    do i=1,n
        sum_r=sum_r+this%tpt_props_heterog%source_term(i)*this%conc_ext(i)
    end do
    sum_r=sum_r*Delta_x*Delta_t
    
    sum_inner=0d0
    do i=2,n-1
        sum_inner = sum_inner + conc_old(i)*(this%tpt_props_heterog%source_term(i+1)-2d0*this%tpt_props_heterog%source_term(i)+this%tpt_props_heterog%source_term(i-1))
    end do
    sum_inner=sum_inner*Delta_x*Delta_t/4d0
    
    inf_term=mass_flux_inf + Delta_x*(conc_old(1)*(this%tpt_props_heterog%source_term(2)-2d0*this%tpt_props_heterog%source_term(1))+dummy_inf*this%tpt_props_heterog%source_term(1))/4d0 - (Delta_x**2)*d2q_dx2*(dummy_inf+conc_old(1))/1.6d1
    out_term=-mass_flux_out + Delta_x*conc_old(n)*(this%tpt_props_heterog%source_term(n-1)-this%tpt_props_heterog%source_term(n))/4d0 + (Delta_x**2)*d2q_dx2*conc_old(n)/8d0
        
    inf_term=inf_term*Delta_t
    out_term=out_term*Delta_t
    
    mass_bal_err=sum_time-sum_r-sum_inner-inf_term-out_term
end function