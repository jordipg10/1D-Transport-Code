module analytical_solutions_transport_m
    use transport_transient_m
    use transport_m
    !use transport_properties_homog_m
    use transport_properties_heterog_m
    use special_fcts_m
    implicit none
    save
    contains    
        
        
    ! Solucion analitica para transporte 1D estacionario con: D=cste, q(x)=-x,  BCs:    c(0)=1, c(L)=0
        function anal_sol_tpt_1D_stat_flujo_lin(this,x) result(conc)
            implicit none
            class(transport_1D_c), intent(in) :: this
            real(kind=8), intent(in) :: x
            !integer(kind=4), intent(in) :: opcion
            real(kind=8) :: conc
            
            real(kind=8) :: D
            real(kind=8), parameter :: pi=4d0*atan(1d0)
            
            D=this%tpt_props_heterog%dispersion(1) ! asumimos D cste
            !c1=1d0
            !c2=c1
            !c=1d0
            !c=-c3
            !a=c+(1d0/(exp(c1/D)*sqrt(D*pi/2d0)*erf(1d0/sqrt(2d0*D))))
            !K=1d0
            if (this%dimensionless==.true.) then
                !if (opcion==1) then
                    conc=-erf(x/sqrt(2d0))/erf(1d0/sqrt(2d0)) + 1
                !else if (opcion==2) then
                !    conc=erf(x/sqrt(2d0))/erf(1d0/sqrt(2d0))
                !else
                !    error stop "Opcion no implementada todavia"
                !end if
            else
                !if (opcion==1) then
                    conc=-erf(x/sqrt(2d0*D))/erf(this%spatial_discr%measure/sqrt(2d0*D)) + 1d0
                !else if (opcion==2) then
                !    conc=erf(x/sqrt(2d0*D))/erf(this%spatial_discr%measure/sqrt(2d0*D))
                !else
                !    error stop "Opcion no implementada todavia"
                !end if
            end if
        end function
        
    ! Derivada solucion analitica para transporte 1D estacionario con: D=cste, q(x)=-x,  BCs:    c(0)=1, c(L)=0
        function der_anal_sol_tpt_1D_stat_flujo_lin(this,x) result(der_conc)
            implicit none
            class(transport_1D_c), intent(in) :: this
            real(kind=8), intent(in) :: x
            !integer(kind=4), intent(in) :: opcion
            real(kind=8) :: der_conc
            
            real(kind=8) :: D,L 
            real(kind=8), parameter :: pi=4d0*atan(1d0)
            
            L=this%spatial_discr%measure ! 1D
            D=this%tpt_props_heterog%dispersion(1) ! asumimos D cste
            !c1=1d0
            !c2=c1
            !c=1d0
            !c=-c3
            !a=c+(1d0/(exp(c1/D)*sqrt(D*pi/2d0)*erf(1d0/sqrt(2d0*D))))
            !K=1d0
            if (this%dimensionless==.true.) then
                !if (opcion==1) then
                    !conc=-erf(x/sqrt(2d0))/erf(1d0/sqrt(2d0)) + 1
                !else if (opcion==2) then
                !    conc=erf(x/sqrt(2d0))/erf(1d0/sqrt(2d0))
                !else
                !    error stop "Opcion no implementada todavia"
                !end if
            else
                !if (opcion==1) then
                    der_conc=(-sqrt(2d0)/(erf(this%spatial_discr%measure/sqrt(2d0*D))*sqrt(pi*D)))*exp(-(x**2)/(2d0*D)) ! derivative analytical solution
                !else if (opcion==2) then
                !    conc=erf(x/sqrt(2d0*D))/erf(this%spatial_discr%measure/sqrt(2d0*D))
                !else
                !    error stop "Opcion no implementada todavia"
                !end if
            end if
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
            !c1=1d0
            !c2=c1
            !c=1d0
            !c=-c3
            !a=c+(1d0/(exp(c1/D)*sqrt(D*pi/2d0)*erf(1d0/sqrt(2d0*D))))
            !K=1d0
            if (this%dimensionless==.true.) then
                
            else
                der_conc=-exp(-(x**3)/(3d0*D))*(3d0**(2d0/3d0))/(gamma(1d0/3d0)-incompl_gamma_term)
            end if
        end function
        
    ! Solucion analitica para transporte 1D transitorio con: D=cste, q(x)=-x, BCs=exp(-t_D/2)
        function anal_sol_tpt_1D_trans(this,x,t) result(conc)
            implicit none
            class(transport_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: x
            real(kind=8), intent(in) :: t
            real(kind=8) :: conc
            
            real(kind=8) :: K1,K2
            real(kind=8), parameter :: pi=4d0*atan(1d0)
            complex(kind=8) :: i=(0d0,1d0),i2=(-1d0,0d0)
            
            !ix=i*x
            
            
            !ix(1)=0d0
            !ix(2)=x
            
            !D=this%tpt_props_heterog%dispersion(1) ! asumimos D cste
            !c1=1d0
            !c2=c1
            !c=1d0
            !c=-c3
            !a=c+(1d0/(exp(c1/D)*sqrt(D*pi/2d0)*erf(1d0/sqrt(2d0*D))))
            !K=1d0
            
            K2=exp(0.25)*(1d0-parab_cyl_fct(-5d-1,(1d0,0d0))/parab_cyl_fct(-5d-1,(0d0,0d0)))/(parab_cyl_fct(-15d-1,i)-parab_cyl_fct(-5d-1,(1d0,0d0))*parab_cyl_fct(-15d-1,(0d0,0d0))/parab_cyl_fct(-5d-1,(0d0,0d0)))
            K1=(1d0-K2*parab_cyl_fct(-15d-1,(0d0,0d0)))/parab_cyl_fct(-5d-1,(0d0,0d0))
            if (this%dimensionless==.true.) then
                conc=(K1*parab_cyl_fct(-5d-1,-i2*x)+K2*parab_cyl_fct(-15d-1,i*x))!*exp(-0.25*x**2)*exp(-t/2)
            end if
        end function
end module