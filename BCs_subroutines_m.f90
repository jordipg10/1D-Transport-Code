module BCs_subroutines_m
    !use analytical_solutions_transport_m
    use spatial_discr_rad_m
    use transport_m
    use transport_transient_m
    implicit none
    save
    contains
         !subroutine Dirichlet_cst_BCs_1D(this,sub,diag,super,nu,k)
         !! Imposes constant Dirichlet boundary conditions
         !   implicit none
         !   class(PDE_1D_c), intent(in) :: this
         !   !class(spatial_discr_c), intent(in) :: mesh
         !   real(kind=8), intent(inout) :: sub(:),diag(:),super(:)
         !   real(kind=8), intent(in), optional :: nu
         !   integer(kind=4), intent(in), optional :: k
         !   
         !   if (size(sub)/=this%spatial_discr%Num_targets-1) then
         !       error stop "Dimension error in subdiagonal array"
         !   else if (size(super)/=size(sub)) then
         !       error stop "Dimension error in superdiagonal array"
         !   else if (size(diag)/=size(super)+1) then
         !       error stop "Dimension error in diagonal array"
         !   end if
         !   sub(this%spatial_discr%Num_targets-1)=0d0
         !   super(1)=0d0
         !   diag(1)=1d0
         !   diag(this%spatial_discr%Num_targets)=1d0
         !end subroutine
         
         subroutine Dirichlet_BCs_PDE(this)
         ! Imposes Dirichlet boundary conditions in transition matrix & source term
            implicit none
            class(PDE_1D_c) :: this
            !real(kind=8), intent(in), optional :: t ! time
            !class(PDE_bis_c) :: this
            !real(kind=8), intent(inout) :: g(:) ! source term
            !real(kind=8), intent(in) :: BCs(:)
            
            integer(kind=4) :: n,opcion,n_flux
            !type(tridiag_sym_matrix_c), target :: trans_mat_sym_BCs
            real(kind=8) :: r_n,r_n_12
            
            n=this%spatial_discr%Num_targets
            
            select type (this)
            type is (transport_1D_c)
                select type (mesh=>this%spatial_discr)
                type is (mesh_1D_Euler_homog_c)
                    if (mesh%targets_flag==1) then
                        this%trans_mat%diag(1)=-1d0
                        this%trans_mat%super(1)=0d0
                        this%trans_mat%diag(n)=-1d0
                        this%trans_mat%sub(n-1)=0d0
                        this%source_term_PDE(1)=this%BCs%conc_inf
                        this%source_term_PDE(1)=this%BCs%conc_out
                    else if (mesh%targets_flag==0 .and. mesh%scheme==1) then ! CFDS
                        this%trans_mat%diag(1)=this%trans_mat%diag(1)-this%tpt_props_heterog%flux(1)/(2d0*mesh%Delta_x) - 3d0*this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
                        this%trans_mat%super(1)=-this%tpt_props_heterog%flux(1)/(2d0*mesh%Delta_x) + this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
                        this%trans_mat%diag(n)=this%trans_mat%diag(n)+this%tpt_props_heterog%flux(n)/(2d0*mesh%Delta_x) - 3d0*this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                        this%trans_mat%sub(n-1)=this%tpt_props_heterog%flux(n)/(2d0*mesh%Delta_x) + this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                        this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*(this%tpt_props_heterog%flux(1)+2d0*this%tpt_props_heterog%dispersion(1)/mesh%Delta_x)/mesh%Delta_x
                        this%source_term_PDE(n)=this%source_term_PDE(n)-this%BCs%conc_out*(this%tpt_props_heterog%flux(n)-2d0*this%tpt_props_heterog%dispersion(n)/mesh%Delta_x)/mesh%Delta_x
                    else if (mesh%targets_flag==0 .and. mesh%scheme==2) then ! IFDS
                        this%trans_mat%diag(1)=this%trans_mat%diag(1)-this%tpt_props_heterog%flux(1)/mesh%Delta_x + this%tpt_props_heterog%flux(2)/(2d0*mesh%Delta_x) - 3d0*this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
                        this%trans_mat%super(1)=-this%tpt_props_heterog%flux(2)/(2d0*mesh%Delta_x) + this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
                        this%trans_mat%diag(n)=this%trans_mat%diag(n)+this%tpt_props_heterog%flux(n+1)/mesh%Delta_x - this%tpt_props_heterog%flux(n)/(2d0*mesh%Delta_x) - 3d0*this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                        this%trans_mat%sub(n-1)=this%tpt_props_heterog%flux(n)/(2d0*mesh%Delta_x) + this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                        this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*(this%tpt_props_heterog%flux(1)+2d0*this%tpt_props_heterog%dispersion(1)/mesh%Delta_x)/mesh%Delta_x
                        this%source_term_PDE(n)=this%source_term_PDE(n)-this%BCs%conc_out*(this%tpt_props_heterog%flux(n+1)-2d0*this%tpt_props_heterog%dispersion(n)/mesh%Delta_x)/mesh%Delta_x
                    else
                        error stop "BCs not implemented yet"
                    end if
                end select
            type is (transport_1D_transient_c)
                select type (mesh=>this%spatial_discr) 
                type is (mesh_1D_Euler_homog_c)
                    if (mesh%targets_flag==1) then
                        this%trans_mat%diag(1)=-1d0
                        this%trans_mat%super(1)=0d0
                        this%trans_mat%diag(n)=-1d0
                        this%trans_mat%sub(n-1)=0d0
                        this%source_term_PDE(1)=this%conc_init(1)
                        this%source_term_PDE(n)=this%conc_init(n)
                    else if (mesh%targets_flag==0 .and. mesh%scheme==1) then ! CFDS
                        this%trans_mat%diag(1)=this%trans_mat%diag(1)-this%tpt_props_heterog%flux(1)/(2d0*mesh%Delta_x) - 3d0*this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
                        this%trans_mat%super(1)=-this%tpt_props_heterog%flux(1)/(2d0*mesh%Delta_x) + this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
                        this%trans_mat%diag(n)=this%trans_mat%diag(n)+this%tpt_props_heterog%flux(n)/(2d0*mesh%Delta_x) - 3d0*this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                        this%trans_mat%sub(n-1)=this%tpt_props_heterog%flux(n)/(2d0*mesh%Delta_x) + this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                        this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*(this%tpt_props_heterog%flux(1)+2d0*this%tpt_props_heterog%dispersion(1)/mesh%Delta_x)/mesh%Delta_x
                        this%source_term_PDE(n)=this%source_term_PDE(n)-this%BCs%conc_out*(this%tpt_props_heterog%flux(n)-2d0*this%tpt_props_heterog%dispersion(n)/mesh%Delta_x)/mesh%Delta_x
                    else if (mesh%targets_flag==0 .and. mesh%scheme==2) then ! IFDS
                        this%trans_mat%diag(1)=this%trans_mat%diag(1)-this%tpt_props_heterog%flux(1)/mesh%Delta_x + this%tpt_props_heterog%flux(2)/(2d0*mesh%Delta_x) - 3d0*this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
                        this%trans_mat%super(1)=-this%tpt_props_heterog%flux(2)/(2d0*mesh%Delta_x) + this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
                        this%trans_mat%diag(n)=this%trans_mat%diag(n)+this%tpt_props_heterog%flux(n+1)/mesh%Delta_x - this%tpt_props_heterog%flux(n)/(2d0*mesh%Delta_x) - 3d0*this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                        this%trans_mat%sub(n-1)=this%tpt_props_heterog%flux(n)/(2d0*mesh%Delta_x) + this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                        this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*(this%tpt_props_heterog%flux(1)+2d0*this%tpt_props_heterog%dispersion(1)/mesh%Delta_x)/mesh%Delta_x
                        this%source_term_PDE(n)=this%source_term_PDE(n)-this%BCs%conc_out*(this%tpt_props_heterog%flux(n+1)-2d0*this%tpt_props_heterog%dispersion(n)/mesh%Delta_x)/mesh%Delta_x
                    else
                        error stop "BCs not implemented yet"
                    end if
                end select
            type is (diffusion_1D_transient_c)
                select type (mesh=>this%spatial_discr)
                type is (spatial_discr_rad_c)
                    if (this%dimensionless==.true.) then
                        this%trans_mat%diag(n)=this%trans_mat%diag(n)-2d0/mesh%Delta_r(n) ! dimensionless
                        this%source_term_PDE(n)=2d0/mesh%Delta_r(n) ! dimensionless
                    end if
                end select
            end select
         end subroutine
         
         !subroutine Dirichlet_anal_BCs_fund_sol_tpt_1D(this,anal_BCs,k)
         !! Computes analytical Dirichlet boundary conditions using fundamental solution transport equation 1D
         !   implicit none
         !   class(PDE_1D_c), intent(in) :: this
         !   real(kind=8), intent(out) :: anal_BCs(2)
         !   integer(kind=4), intent(in), optional :: k
         !   select type (this)
         !   type is (transport_1D_transient_c)
         !       anal_BCs(1)=fund_sol_tpt_eqn_1D(this,1,k)
         !       anal_BCs(2)=fund_sol_tpt_eqn_1D(this,this%spatial_discr%Num_targets,k)
         !   end select
         !end subroutine 
         
         !subroutine Neumann_homog_BCs_EE_1D(this,sub,diag,super)
         !! Imposes Neumann homogeneous boundary conditions for the Euler explicit method
         !   implicit none
         !   class(transport_1D_c), intent(in) :: this
         !   real(kind=8), intent(inout) :: sub(:),diag(:),super(:)
         !   if (size(sub)/=this%spatial_discr%Num_targets-1) then
         !       error stop "Dimension error in subdiagonal array"
         !   else if (size(super)/=size(sub)) then
         !       error stop "Dimension error in superdiagonal array"
         !   else if (size(diag)/=size(super)+1) then
         !       error stop "Dimension error in diagonal array"
         !   end if
         !   select type (this)
         !   type is (transport_1D_transient_c)
         !       select type (parameters=>this%parameters)
         !       type is (parameters_homog_transient_c)
         !           super(1)=2*parameters%beta
         !           sub(this%spatial_discr%Num_targets-1)=super(1)
         !           select type (time_discr=>this%time_discr)
         !           type is (time_discr_homog_c)
         !               diag(1)=1d0-2*parameters%beta-parameters%source_term_flag*time_discr%Delta_t*parameters%source_term
         !           !type is (time_discr_heterog_c)
         !           !    diag(1)=1d0-2*parameters%beta-parameters%source_term_flag*time_discr%Delta_t(k)*parameters%source_term
         !           end select
         !           diag(this%spatial_discr%Num_targets)=diag(1)
         !       end select
         !   end select
         !end subroutine 
        
        !subroutine Neumann_homog_BCs_EE_1D(this,conc_old,conc_new,k)
        !! Computes Neumann homogeneous boundary conditions for the Euler Explicit method
        !    implicit none
        !    class(transport_1D_transient_c), intent(in) :: this
        !    real(kind=8), intent(in) :: conc_old(:)
        !    real(kind=8), intent(inout) :: conc_new(:)
        !    integer(kind=4), intent(in), optional :: k
        !    
        !    !if (this%BCs==1) error stop "This subroutine is for Neumann boundary conditions"
        !    
        !
        !        select type (time_discr=>this%time_discr)
        !        type is (time_discr_homog_c)
        !            select type (parameters=>this%parameters)
        !            type is (parameters_homog_transient_c)
        !                !print *, parameters%source_term_flag*time_discr%Delta_t*parameters%source_term
        !                conc_new(1)=conc_old(2)*parameters%beta*2+conc_old(1)*(1-2*parameters%beta-parameters%source_term_flag*time_discr%Delta_t*parameters%source_term)+this%conc_star_flag*time_discr%Delta_t*parameters%source_term*this%conc_ext(1)
        !                !print *, conc_new(1)
        !                conc_new(this%spatial_discr%Num_targets)=conc_old(this%spatial_discr%Num_targets-1)*parameters%beta*2+conc_old(this%spatial_discr%Num_targets)*(1-2*parameters%beta-parameters%source_term_flag*time_discr%Delta_t*parameters%source_term)+this%conc_star_flag*time_discr%Delta_t*parameters%source_term*this%conc_ext(this%spatial_discr%Num_targets)
        !            !type is (parameters_homog_transient_discr_heterog_c)
        !            !    conc_new(1)=conc_old(2)*parameters%beta(1,k)*2+conc_old(1)*(1-2*parameters%beta(1,k)-parameters%source_term_flag*time_discr%Delta_t*parameters%source_term)+this%conc_star_flag*time_discr%Delta_t*parameters%source_term*this%conc_ext(1)
        !            !    conc_new(this%spatial_discr%Num_targets)=conc_old(this%spatial_discr%Num_targets-1)*parameters%beta(this%spatial_discr%Num_targets,k)*2+conc_old(this%spatial_discr%Num_targets)*(1-2*parameters%beta(this%spatial_discr%Num_targets,k)-parameters%source_term_flag*time_discr%Delta_t*parameters%source_term)+this%conc_star_flag*time_discr%Delta_t*parameters%source_term*this%conc_ext(this%spatial_discr%Num_targets)
        !            !type is (parameters_heterog_transient_discr_homog_c)
        !            !    conc_new(1)=conc_old(2)*parameters%beta(1)*2+conc_old(1)*(1-2*parameters%beta(1)-parameters%source_term_flag*time_discr%Delta_t*parameters%source_term)+this%conc_star_flag*time_discr%Delta_t*parameters%source_term*this%conc_ext(1)
        !            !    conc_new(this%spatial_discr%Num_targets)=conc_old(this%spatial_discr%Num_targets-1)*parameters%beta(this%spatial_discr%Num_targets)*2+conc_old(this%spatial_discr%Num_targets)*(1-2*parameters%beta(this%spatial_discr%Num_targets)-parameters%source_term_flag*time_discr%Delta_t*parameters%source_term)+this%conc_star_flag*time_discr%Delta_t*parameters%source_term*this%conc_ext(this%spatial_discr%Num_targets)
        !            !type is (parameters_heterog_transient_c)
        !            !    conc_new(1)=conc_old(2)*parameters%beta(1,k)*2+conc_old(1)*(1-2*parameters%beta(1,k)-parameters%source_term_flag*time_discr%Delta_t*parameters%source_term)+this%conc_star_flag*time_discr%Delta_t*parameters%source_term*this%conc_ext(1)
        !            !    conc_new(this%spatial_discr%Num_targets)=conc_old(this%spatial_discr%Num_targets-1)*parameters%beta(this%spatial_discr%Num_targets,k)*2+conc_old(this%spatial_discr%Num_targets)*(1-2*parameters%beta(this%spatial_discr%Num_targets,k)-parameters%source_term_flag*time_discr%Delta_t*parameters%source_term)+this%conc_star_flag*time_discr%Delta_t*parameters%source_term*this%conc_ext(this%spatial_discr%Num_targets)
        !            end select
        !        type is (time_discr_heterog_c)
        !            !select type (parameters=>this%parameters)
        !            !type is (parameters_homog_transient_discr_heterog_c)
        !            !    conc_new(1)=conc_old(2)*parameters%beta(1,k)*2+conc_old(1)*(1-2*parameters%beta(1,k)-parameters%source_term_flag*time_discr%Delta_t(k)*parameters%source_term)+this%conc_star_flag*time_discr%Delta_t(k)*parameters%source_term*this%conc_ext(1)
        !            !    conc_new(this%spatial_discr%Num_targets)=conc_old(this%spatial_discr%Num_targets-1)*parameters%beta(this%spatial_discr%Num_targets,k)*2+conc_old(this%spatial_discr%Num_targets)*(1-2*parameters%beta(this%spatial_discr%Num_targets,k)-parameters%source_term_flag*time_discr%Delta_t(k)*parameters%source_term)+this%conc_star_flag*time_discr%Delta_t(k)*parameters%source_term*this%conc_ext(this%spatial_discr%Num_targets)
        !            !type is (parameters_heterog_transient_c)
        !            !    conc_new(1)=conc_old(2)*parameters%beta(1,k)*2+conc_old(1)*(1-2*parameters%beta(1,k)-parameters%source_term_flag*time_discr%Delta_t(k)*parameters%source_term)+this%conc_star_flag*time_discr%Delta_t(k)*parameters%source_term*this%conc_ext(1)
        !            !    conc_new(this%spatial_discr%Num_targets)=conc_old(this%spatial_discr%Num_targets-1)*parameters%beta(this%spatial_discr%Num_targets,k)*2+conc_old(this%spatial_discr%Num_targets)*(1-2*parameters%beta(this%spatial_discr%Num_targets,k)-parameters%source_term_flag*time_discr%Delta_t(k)*parameters%source_term)+this%conc_star_flag*time_discr%Delta_t(k)*parameters%source_term*this%conc_ext(this%spatial_discr%Num_targets)
        !            !end select
        !        end select
        !end subroutine Neumann_homog_BCs_EE_1D
        
        !subroutine Neumann_EE_BCs_1D_var_param(this,conc_old,conc_new,k)
        !! Computes Neumann homogeneous boundary conditions for the Euler Explicit method with constant parameters
        !    implicit none
        !    class(transport_1D_c), intent(in) :: this
        !    real(kind=8), intent(in) :: conc_old(:)
        !    real(kind=8), intent(inout) :: conc_new(:)
        !    integer(kind=4), intent(in), optional :: k
        !    
        !    !if (this%BCs==1) error stop "This subroutine is for Neumann boundary conditions"
        !    
        !    select type (this)
        !    type is (transport_1D_transient)
        !        select type (parameters=>this%parameters)
        !        type is (var_parameters_transient)
        !            !if (this%BCs==2 .and. this%int_method/=1) then
        !            !    error stop "This subroutine is for the Euler Explicit method"
        !            !else if (this%BCs==2 .and. this%int_method==1) then
        !                ! Neumann homogeneous
        !                conc_new(1)=conc_old(2)*parameters%beta(1,k)*2+conc_old(1)*(1-2*parameters%beta(1,k)-parameters%source_term_flag*parameters%Delta_t(k)*parameters%source_term)+this%conc_star_flag*parameters%Delta_t(k)*parameters%source_term*this%conc_ext(1)
        !                conc_new(this%spatial_discr%Num_targets)=conc_old(this%spatial_discr%Num_targets-1)*parameters%beta(this%spatial_discr%Num_targets,k)*2+conc_old(this%spatial_discr%Num_targets)*(1-2*parameters%beta(this%spatial_discr%Num_targets,k)-parameters%source_term_flag*parameters%Delta_t(k)*parameters%source_term)+this%conc_star_flag*parameters%Delta_t(k)*parameters%source_term*this%conc_ext(this%spatial_discr%Num_targets)
        !            !else
        !            !    error stop "Boundary conditions not implemented yet"
        !            !end if
        !        end select
        !    end select
        !end subroutine Neumann_EE_BCs_1D_var_param
        
         subroutine Neumann_homog_BCs(this)
         ! Imposes Neumann homogeneous boundary conditions in transition matrix
            implicit none
            class(PDE_1D_c) :: this
            !class(PDE_bis_c) :: this
            !real(kind=8), intent(inout) :: sub(:),diag(:),super(:)
            !real(kind=8), intent(in), optional :: nu
            !integer(kind=4), intent(in), optional :: k
            
            integer(kind=4) :: n
            
            n=this%spatial_discr%Num_targets

            select type (this)
            type is (transport_1D_transient_c)
                select type (mesh=>this%spatial_discr)
                type is (mesh_1D_Euler_homog_c)
                    if (mesh%scheme==1 .and. mesh%targets_flag==0) then ! CFD
                        this%trans_mat%super(1)=-this%tpt_props_heterog%flux(1)/(2d0*mesh%Delta_x)+this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
                        this%trans_mat%sub(n-1)=this%tpt_props_heterog%flux(n)/(2d0*mesh%Delta_x)+this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                        this%trans_mat%diag(1)=this%trans_mat%diag(1)-this%trans_mat%super(1)
                        this%trans_mat%diag(n)=this%trans_mat%diag(n)-this%trans_mat%sub(n-1)
                    else if (mesh%scheme==2 .and. mesh%targets_flag==0) then ! IFD
                        this%trans_mat%super(1)=-this%tpt_props_heterog%flux(2)/(2d0*mesh%Delta_x)+this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
                        this%trans_mat%sub(n-1)=this%tpt_props_heterog%flux(n)/(2d0*mesh%Delta_x)+this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                        this%trans_mat%diag(1)=this%trans_mat%diag(1)-this%trans_mat%super(1)
                        this%trans_mat%diag(n)=this%trans_mat%diag(n)-this%trans_mat%sub(n-1)
                    end if
                type is (mesh_1D_Euler_heterog_c)
                    !this%trans_mat%super(1)=2d0*this%this%diff_props_homog%dispersion/(mesh%Delta_x(1)**2)
                    !this%trans_mat%sub(n-1)=2d0*this%this%diff_props_homog%dispersion/(mesh%Delta_x(n)**2)
                    !this%trans_mat%diag(1)=-this%trans_mat%super(1)
                    !this%trans_mat%diag(n)=-this%trans_mat%sub(n-1)
                end select
            end select
         end subroutine 
         
         subroutine Robin_Neumann_homog_BCs(this)
         ! Imposes Robin BC inflow & Neumann homogeneous BC outflow in transition matrix & sink/source term
            implicit none
            class(PDE_1D_c) :: this
            !class(PDE_bis_c) :: this
            !real(kind=8), intent(inout) :: sub(:),diag(:),super(:)
            !real(kind=8), intent(in), optional :: nu
            !integer(kind=4), intent(in), optional :: k
            
            integer(kind=4) :: n
            real(kind=8) :: q_1,q_32,q_n,q_inf,q_out,D_1,D_n,c_inf
            
            n=this%spatial_discr%Num_targets
            
            select type (this)
            type is (transport_1D_transient_c)
                q_1=this%tpt_props_heterog%flux(1)
                q_32=this%tpt_props_heterog%flux(2)
                q_n=this%tpt_props_heterog%flux(n)
                q_inf=this%BCs%flux_inf
                q_out=this%BCs%flux_out
                D_1=this%tpt_props_heterog%dispersion(1)
                D_n=this%tpt_props_heterog%dispersion(n)
                c_inf=this%BCs%conc_inf
                select type (mesh=>this%spatial_discr)
                type is (mesh_1D_Euler_homog_c)
                    if (mesh%scheme==1 .and. mesh%targets_flag==0) then ! CFDS
                        this%trans_mat%super(1)=-q_1/(2d0*mesh%Delta_x)+D_1/(mesh%Delta_x**2)
                        this%trans_mat%sub(n-1)=q_n/(2d0*mesh%Delta_x)+D_n/(mesh%Delta_x**2)
                        this%trans_mat%diag(1)=this%trans_mat%diag(1)+(2d0*D_1-q_inf*mesh%Delta_x)*(q_1*mesh%Delta_x+2*D_1)/(2*q_inf*mesh%Delta_x**3+4*D_1*mesh%Delta_x**2) - 2*D_1/(mesh%Delta_x**2)
                        this%trans_mat%diag(n)=this%trans_mat%diag(n)-this%trans_mat%sub(n-1)
                        this%source_term_PDE(1)=this%source_term_PDE(1)+q_inf*c_inf*(q_1*mesh%Delta_x+2*D_1)/(q_inf*mesh%Delta_x**2+2*D_1*mesh%Delta_x)
                    else if (mesh%scheme==2 .and. mesh%targets_flag==0) then ! IFDS
                        this%trans_mat%super(1)=-q_32/(2d0*mesh%Delta_x)+D_1/(mesh%Delta_x**2)
                        this%trans_mat%sub(n-1)=this%tpt_props_heterog%flux(n)/(2d0*mesh%Delta_x)+this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                        this%trans_mat%diag(1)=this%trans_mat%diag(1)-(q_inf**2+D_1*(3*q_inf*mesh%Delta_x+2*D_1)/(mesh%Delta_x**2))/(q_inf*mesh%Delta_x+2*D_1) + q_32/(2*mesh%Delta_x)
                        this%trans_mat%diag(n)=this%trans_mat%diag(n)-this%trans_mat%sub(n-1)
                        this%source_term_PDE(1)=this%source_term_PDE(1)+q_inf*c_inf/mesh%Delta_x
                    end if
                type is (mesh_1D_Euler_heterog_c)
                    !this%trans_mat%super(1)=2d0*this%this%diff_props_homog%dispersion/(mesh%Delta_x(1)**2)
                    !this%trans_mat%sub(n-1)=2d0*this%this%diff_props_homog%dispersion/(mesh%Delta_x(n)**2)
                    !this%trans_mat%diag(1)=-this%trans_mat%super(1)
                    !this%trans_mat%diag(n)=-this%trans_mat%sub(n-1)
                end select
            end select
         end subroutine 
       
        
end module