! Writes data and results of 1D diffusion model
subroutine write_diffusion_1D(this,Time_out,output)
    use diffusion_m
    use spatial_discr_1D_m
    !use time_discr_m
    implicit none
    ! Variables
    class(diffusion_1D_c), intent(in) :: this ! transport object
    real(kind=8), intent(in) :: Time_out(:)
    real(kind=8), intent(in) :: output(:,:)
    !real(kind=8), intent(in) :: Time_out(:)
    !real(kind=8), intent(in) :: output(:,:)
    !class(props_c), intent(in), optional :: props
    !real(kind=8), intent(in), optional :: theta

    !integer(kind=4) :: Num_output ! size(Time_out)
    real(kind=8), allocatable :: props_mat(:,:),Delta_t(:),stab_params(:)
    integer(kind=4) :: i,j,k
    
    !Num_output=size(Time_out)
    open(unit=1,file='diffusion_1D.out',status='unknown')!,position='append')!,iostat=1,iomsg='error')
    write(1,"(2x,'Equation:',5x,'0 = T*c + g',/)")
    select type (mesh=>this%spatial_discr)
    type is (mesh_1D_Euler_homog_c)
        write(1,"(2x,'Length of domain:',F15.5/)") mesh%measure
        write(1,"(2x,'Number of cells:',I5/)") this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
        write(1,"(2x,'Mesh size:',F15.5/)") mesh%Delta_x
    end select
    !write(1,"(2x,'Critical mesh size:',F15.5/)") this%stab_params%Delta_x_crit
    if (this%BCs%BCs_label(1)==1 .and. this%BCs%BCs_label(2)==1) then
        !write(1,"(2x,'Boundary conditions:',10x,'Dirichlet',4x,F15.5/)") Dirichlet_BC
        write(1,"(2x,'Boundary conditions:',10x,'Dirichlet constant',/)")
    else if (this%BCs%BCs_label(1)==2 .and. this%BCs%BCs_label(2)==2) then
        write(1,"(2x,'Boundary conditions:',10x,'Neumann homogeneous',/)")
    end if
    if (this%spatial_discr%scheme==1) then
        write(1,"(2x,'Scheme:',10x,'CFD',/)")
    else if (this%spatial_discr%scheme==2) then
        write(1,"(2x,'Scheme:',10x,'IFD',/)")
    else if (this%spatial_discr%scheme==3) then
        write(1,"(2x,'Scheme:',10x,'Upwind',/)")
    end if
    !write(1,"(2x,'Time step:',ES15.5/)") this%time_discr%get_Delta_t()
    !write(1,"(2x,'Critical time step:',ES15.5/)") this%stab_params%Delta_t_crit
    !write(1,"(2x,'Final time:',ES15.5/)") this%time_discr%Final_time
    !write(1,"(2x,'Number of time steps:',I10/)") this%time_discr%Num_time
    !if (this%time_discr%int_method==1) then
    !    write(1,"(2x,'Integration method:',10x,'Euler explicit',/)")
    !else if (this%time_discr%int_method==2) then
    !    write(1,"(2x,'Integration method:',10x,'Euler fully implicit',/)")
    !else if (this%time_discr%int_method==3) then
    !    write(1,"(2x,'Integration method:',10x,'Crank-Nicolson',/)")
    !else if (this%time_discr%int_method==4) then
    !    write(1,"(2x,'Integration method:',10x,'RKF45',/)")
    !end if
    !select type (props)
    !type is (tpt_props_homog_c)
        write(1,"(2x,'Properties:'/)")
        !props_mat=props%get_props()
        !write(1,"(10x,'Porosity:',ES15.5,10x,'Retardo',ES15.5,10x,'Dispersion:',ES15.5,10x,'Flux:',ES15.5/)") props(1,1), props(1,2), props(1,3), props(1,4)
        write(1,"(10x,'Dispersion:'/)")
        !do i=1,size(props_mat,1)
        !write(1,"(10x,*(ES15.5)/)") this%tpt_props_heterog%dispersion, this%tpt_props_heterog%flux
        !end do
        !if (this%time_discr%int_method==1) then
        !    stab_params=this%stab_params%get_stab_params()
        !    write(1,"(2x,'Stability parameters:'/)")
        !    write(1,"(10x,'Critical time step:',ES15.5/)") stab_params(2)
        !    !write(1,"(10x,'Critical mesh size:',F15.5,10x,'Critical time step:',ES15.5/)") stab_params(1), stab_params(2)
        !    write(1,"(10x,'Beta:',ES15.5,10x,'Alpha:',ES15.5/)") stab_params(3), stab_params(4)
        !    write(1,"(10x,'Peclet:',ES15.5/)") stab_params(5)
        !end if
    !type is (tpt_props_heterog_c)
    !    !write(1,*) props%flux
    !end select
    !select type (spatial_discr=>this%spatial_discr)
    !type is (mesh_1D_Euler_homog_c)
    !    write(1,"(2x,'Delta x',F20.2)") spatial_discr%Delta_x
    !!type is (mesh_1D_Euler_heterog_c)
    !!    write(1,"(2x,'Delta x',F10.2)") spatial_discr%Delta_x
    !end select
    !write(1,"(2x,'F:'/)") 
    !do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
    !    write(1,"(2x,F15.5)") this%F_mat%diag(i)
    !end do
    write(1,"(/,2x,'Transition matrix T (with BCs):'/)") 
    write(1,"(17x,2F15.5)") this%trans_mat%diag(1), this%trans_mat%super(1)    
    do i=2,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag-1
        write(1,"(2x,3F15.5)") this%trans_mat%sub(i-1), this%trans_mat%diag(i), this%trans_mat%super(i)
    end do
    write(1,"(2x,2F15.5/)") this%trans_mat%sub(this%spatial_discr%Num_targets-1), this%trans_mat%diag(this%spatial_discr%Num_targets)
    !write(1,"(2x,'r:'/)")
    !do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
    !    write(1,"(2x,F15.5)") this%props%source_term(i)
    !end do
    !write(1,"(/,2x,'External concentration:'/)")
    !do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
    !    write(1,"(2x,F15.5)") this%conc_ext(i)
    !end do
    write(1,"(/,2x,'Source term g:'/)")
    do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
        write(1,"(2x,F15.5)") this%source_term_PDE(i)
    end do
    !call write_tridiag_sym_matrix(this,filename)
    !write(1,"(2x,'Matrix A:'/)")
    !    do i=1,Num_targets-1
    !        write(1,"(2x,2F15.5/)") A_mat%diag(i), A_mat%sub(i)
    !    end do
    !    write(1,"(2x,F15.5/)") A_mat%diag(Num_targets)
    ! write(1,"(2x,'Eigenvalues:',/)")
    !    write(1,"(2x,5F15.5/)") lambda
    !    allocate(K(Num_targets,2))
    !    call Gerschgorin_tridiag_sym_mat(A_diag,A_super,K,nK)
    !    write(1,"(2x,'Gerschgorin disks:',/)")
    !    do i=1,nK
    !        write(1,"(2x,2F15.5/)") (K(i,j), j=1,2)
    !    end do
    !     write(1,"(2x,'Eigenvectors (by columns):',/)")
    !    do i=1,Num_targets
    !        write(1,"(2x,5F15.5/)") (eigenvectors(i,j), j=1,Num_targets)
    !    end do
       
    !    write(1,"(2x,'Initial concentration (using total product):',/)")
    !    !write(1,"(2x,'y(t=0):',/)")
    !    if (opcion_BCs==1) then
    !        y=prod_total(lambda,eigenvectors,c0,b,my_homog_time_discr,0)
    !    else
    !        y=prod_total_sym_mat(lambda,eigenvectors,c0,b,my_homog_time_discr,0)
    !    end if
    !    call my_diff_transient%prod_total_conc(y)
    !    write(1,"(2x,5F15.5/)") my_diff_transient%conc
    !    !write(1,"(2x,F15.5,9F15.5/)") y0(1),prod_total_y(lambda_star,eigenvectors_star,y0(2:Num_targets),b(2:Num_targets),my_homog_time_discr,0)
    !    write(1,"(2x,'Final concentration (using total product):',/)")
    !    !write(1,"(2x,'y(t=T):',/)")
    !    deallocate(my_diff_transient%conc)
    !    if (opcion_BCs==1) then
    !        y=prod_total(lambda,eigenvectors,c0,b,my_homog_time_discr)
    !    else
    !        y=prod_total_sym_mat(lambda,eigenvectors,c0,b,my_homog_time_discr)
    !    end if
    !    call my_diff_transient%prod_total_conc(y)
    !    write(1,"(2x,5F15.5/)") my_diff_transient%conc
    !if (size(output,1)==this%spatial_discr%Num_targets-this%spatial_discr%targets_flag) then
        !if (this%time_discr%int_method<4) then
            write(1,"(/,2x,'Cell'/)")
            do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
                write(1,"(2x,I5,ES15.5)") i, this%conc(i)
            end do
    !    else
    !        write(1,"(/,2x,'Cell',2ES20.5/)") Time_out(1), this%time_discr%Final_time
    !        do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
    !            write(1,"(2x,I5,2ES20.10)") i, this%conc_init(i), this%conc(i)
    !        end do
    !    end if
    !else
    !    write(1,"(/,2x,'Mobile zone:'/)")
    !    write(1,"(10x,'Cell',3ES20.5/)") (Time_out(k), k=1,Num_output)
    !    do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
    !        write(1,"(10x,I4,3F20.5)") i,(output(i,k), k=1,Num_output)
    !    end do
    !    write(1,"(/,2x,'Immobile zones:'/)")
    !    !write(1,"(/,10x,'Immobile zones',/)")
    !    do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
    !        write(1,"(10x,I4,3F20.5)") i,(output(this%spatial_discr%Num_targets-this%spatial_discr%targets_flag+i,k), k=1,Num_output)
    !    end do
    !end if
    !select type (spatial_discr=>this%spatial_discr)
    !type is (mesh_1D_Euler_homog_c)
    !    write(1,"(8x,5F22.10/)") (sum(spatial_discr%Delta_x*conc_out(:,k)), k=1,Num_output)
    !type is (mesh_1D_Euler_heterog_c)
    !    write(1,"(8x,5F22.10/)") (sum(spatial_discr%Delta_x(:)*conc_out(:,k)), k=1,Num_output)
    !end select
    !write(1,"(2x,'Fundamental solution transport: (rows->targets, columns->time steps)'/)")
    !do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
    !    write(1,"(2x,5F15.5)") (fund_sol_tpt_eqn_1D(this,i,j-1), j=1,5)!this%time_discr%Num_time) ! Num_time=10
    !end do
    rewind(1)
    close(1)
end subroutine