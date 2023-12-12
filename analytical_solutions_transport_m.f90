module analytical_solutions_transport_m
    use transport_transient_m
    use transport_m
    use transport_properties_heterog_m
    use special_fcts_m
    implicit none
    save
    contains
        
    ! Solucion analitica para transporte 1D estacionario con: D=cste, q(x)=-x,  BCs: c(0)=1, c(L)=0
        function anal_sol_tpt_1D_stat_flujo_lin(this,x) result(conc)
            implicit none
            class(transport_1D_c), intent(in) :: this
            real(kind=8), intent(in) :: x
            real(kind=8) :: conc
            
            real(kind=8) :: D
            real(kind=8), parameter :: pi=4d0*atan(1d0)
            
            D=this%tpt_props_heterog%dispersion(1) ! asumimos D cste

            if (this%dimensionless==.true.) then
                conc=-erf(x/sqrt(2d0))/erf(1d0/sqrt(2d0)) + 1
            else
                conc=-erf(x/sqrt(2d0*D))/erf(this%spatial_discr%measure/sqrt(2d0*D)) + 1d0
            end if
        end function
        
    ! Derivada solucion analitica para transporte 1D estacionario con: D=cste, q(x)=-x,  BCs: c(0)=1, c(L)=0
        function der_anal_sol_tpt_1D_stat_flujo_lin(this,x) result(der_conc)
            implicit none
            class(transport_1D_c), intent(in) :: this
            real(kind=8), intent(in) :: x
            real(kind=8) :: der_conc
            
            real(kind=8) :: D,L 
            real(kind=8), parameter :: pi=4d0*atan(1d0)
            
            L=this%spatial_discr%measure ! 1D
            D=this%tpt_props_heterog%dispersion(1) ! asumimos D cste
            der_conc=(-sqrt(2d0)/(erf(this%spatial_discr%measure/sqrt(2d0*D))*sqrt(pi*D)))*exp(-(x**2)/(2d0*D)) ! derivative analytical solution
        end function

        
    ! Derivada solucion analitica para transporte 1D estacionario con: D=cste, q(x)=-x^2,  BCs: c(0)=1, c(L)=0
        function der_anal_sol_tpt_1D_stat_flujo_cuad(this,x) result(der_conc)
            implicit none
            class(transport_1D_c), intent(in) :: this
            real(kind=8), intent(in) :: x
            real(kind=8) :: der_conc
            
            real(kind=8) :: D,L
            real(kind=8), parameter :: pi=4d0*atan(1d0),incompl_gamma_term=0.327336564991358
            
            D=this%tpt_props_heterog%dispersion(1) ! asumimos D cste
            L=this%spatial_discr%measure ! 1D

            der_conc=-exp(-(x**3)/(3d0*D))*(3d0**(2d0/3d0))/(gamma(1d0/3d0)-incompl_gamma_term)
        end function
end module