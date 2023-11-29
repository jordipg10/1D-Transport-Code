function mass_balance_error_ADE_stat_Dirichlet_discharge(this) result(mass_bal_err)
    use transport_m
    implicit none
    class(transport_1D_c), intent(in) :: this
    !real(kind=8), intent(in) :: conc(:)
    !real(kind=8), intent(in) :: conc_new(:)
    !real(kind=8), intent(in) :: q_inf
    !real(kind=8), intent(in) :: q_out
    !real(kind=8), intent(in) :: Delta_t
    !real(kind=8), intent(in) :: Delta_x
    !integer(kind=4), intent(in) :: k
    real(kind=8) :: mass_bal_err
    
    real(kind=8) :: sum_r,inf_term,out_term
    integer(kind=4) :: i,n
    
    n=this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
    
    !sum_time=0d0
    !
    !do i=1,n
    !    sum_time=sum_time+this%tpt_props_heterog%porosity(i)*(conc_new(i)-conc_old(i))
    !end do
    !sum_time=sum_time*Delta_x!/Delta_t
    sum_r=0d0
    select type (mesh=>this%spatial_discr)
    type is (mesh_1D_Euler_homog_c) 
        do i=1,n
            sum_r=sum_r+this%tpt_props_heterog%source_term(i)*this%conc(i)
        end do
        sum_r=sum_r*mesh%Delta_x
    
        !inf_term=this%BCs%conc_inf*(this%BCs%flux_inf+2d0*this%tpt_props_heterog%dispersion(1)/Delta_x) - conc_old(1)*2d0*this%tpt_props_heterog%dispersion(1)/Delta_x
        !out_term=-this%BCs%conc_out*(this%BCs%flux_out-2d0*this%tpt_props_heterog%dispersion(n)/Delta_x) - conc_old(n)*2d0*this%tpt_props_heterog%dispersion(n)/Delta_x
    
        inf_term = 2d0*this%tpt_props_heterog%dispersion(1)*(this%BCs%conc_inf-this%conc(1))/mesh%Delta_x + this%BCs%flux_inf*this%BCs%conc_inf
        out_term = 2d0*this%tpt_props_heterog%dispersion(n)*(this%BCs%conc_out-this%conc(n))/mesh%Delta_x - this%BCs%flux_out*this%BCs%conc_out
    
        !inf_term=inf_term*Delta_t
        !out_term=out_term*Delta_t
    end select
    mass_bal_err=-sum_r-inf_term-out_term
end function