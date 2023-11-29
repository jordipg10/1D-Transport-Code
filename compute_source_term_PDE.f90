subroutine compute_source_term_PDE(this,anal_sol)
! g=rc*
! g=BCs_fct if BCs_label=3
    !use BCs_subroutines_m
    !use transport_properties_homog_m
    !use transport_properties_heterog_m
    use transport_m
    use transport_transient_m
    implicit none
    class(PDE_1D_c) :: this
    !real(kind=8), intent(in), optional :: conc_old(:)
    real(kind=8), external, optional :: anal_sol
    
    !real(kind=8), allocatable :: g(:), porosity(:)
    real(kind=8) :: anal_BCs(2)

    select type (this)
    type is (transport_1D_c)
        this%source_term_PDE=this%conc_r_flag*this%tpt_props_heterog%source_term*this%conc_ext
    type is (transport_1D_transient_c)
        this%source_term_PDE=this%conc_r_flag*this%tpt_props_heterog%source_term*this%conc_ext
    type is (diffusion_1D_c)
        this%source_term_PDE=this%diff_props_heterog%source_term*this%conc_ext
    type is (diffusion_1D_transient_c)
        this%source_term_PDE=this%diff_props_heterog%source_term*this%conc_ext
    end select
    !if (this%BCs%option==1) then
        
    !if (this%BCs%BCs_label==3) then
    !    call Dirichlet_anal_BCs_1D(this,anal_sol,anal_BCs,k)
    !    g(1)=anal_BCs(1)-conc_old(1)
    !    g(this%spatial_discr%Num_targets)=anal_BCs(2)-conc_old(this%spatial_discr%Num_targets)
    !end if
end subroutine