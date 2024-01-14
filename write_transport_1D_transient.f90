! Writes data and results of 1D transient transport equation
subroutine write_transport_1D_transient(this,Time_out,output)
    use analytical_solutions_transport_m
    implicit none
    ! Variables
    class(transport_1D_transient_c), intent(in) :: this ! 1D transient transport object
    real(kind=8), intent(in) :: Time_out(:)
    real(kind=8), intent(in) :: output(:,:)

    integer(kind=4) :: Num_output
    real(kind=8), allocatable :: Delta_t(:),stab_params(:)
    integer(kind=4) :: i,j,k,n,n_flux
    character(len=256) :: file_out
    
    n=this%spatial_discr%Num_targets
    Num_output=size(Time_out)
    if (this%spatial_discr%scheme==1) then
        if (this%dimensionless==.true.) then
            write(file_out,"('transport_1D_trans_CFDS_adim.out')")
        else
            write(file_out,"('transport_1D_trans_CFDS.out')")
        end if
    else if (this%spatial_discr%scheme==2) then
        if (this%dimensionless==.true.) then
            write(file_out,"('transport_1D_trans_IFDS_adim.out')")
        else
            write(file_out,"('transport_1D_trans_IFDS.out')")
        end if
    end if
    open(unit=1,file=file_out,status='unknown')
    write(1,"(2x,'Equation:',5x,'F*dc/dt = T*c + g',/)")
    select type (mesh=>this%spatial_discr)
    type is (mesh_1D_Euler_homog_c)
        write(1,"(2x,'Length of domain:',F15.5/)") mesh%measure
        write(1,"(2x,'Number of cells:',I5/)") this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
        write(1,"(2x,'Mesh size:',F15.5/)") mesh%Delta_x
    end select
    if (this%BCs%BCs_label(1)==1 .and. this%BCs%BCs_label(2)==1) then
        write(1,"(2x,'Boundary conditions:',10x,'Dirichlet',10x,2F15.5,/)") this%BCs%conc_inf, this%BCs%conc_out
    else if (this%BCs%BCs_label(1)==2 .and. this%BCs%BCs_label(2)==2) then
        write(1,"(2x,'Boundary conditions:',10x,'Neumann homogeneous',/)")
    else if (this%BCs%BCs_label(1)==3 .and. this%BCs%BCs_label(2)==2) then
        write(1,"(2x,'Boundary conditions:',10x,'Robin inflow, Neumann homogeneous outflow',/)")
    end if
    if (this%spatial_discr%scheme==1) then
        write(1,"(2x,'Scheme:',10x,'CFD',/)")
    else if (this%spatial_discr%scheme==2) then
        write(1,"(2x,'Scheme:',10x,'IFD',/)")
    else if (this%spatial_discr%scheme==3) then
        write(1,"(2x,'Scheme:',10x,'Upwind',/)")
    end if
    write(1,"(2x,'Time step:',ES15.5/)") this%time_discr%get_Delta_t()
    write(1,"(2x,'Final time:',ES15.5/)") this%time_discr%Final_time
    write(1,"(2x,'Number of time steps:',I10/)") this%time_discr%Num_time
    if (this%time_discr%int_method==1) then
        write(1,"(2x,'Integration method:',10x,'Lagr explicit',/)")
    else if (this%time_discr%int_method==2) then
        write(1,"(2x,'Integration method:',10x,'Lagr fully implicit',/)")
    else if (this%time_discr%int_method==3) then
        write(1,"(2x,'Integration method:',10x,'Crank-Nicolson',/)")
    else if (this%time_discr%int_method==4) then
        write(1,"(2x,'Integration method:',10x,'RKF45',/)")
    end if

    if (this%dimensionless==.false.) then
        write(1,"(2x,'Properties:'/)")
        write(1,"(10x,'Dispersion:',/)")
        do i=1,n
            write(1,"(20x,ES15.5)") this%tpt_props_heterog%dispersion(i)
        end do
        write(1,"(/,10x,'Flux:'/)")
        n_flux=size(this%tpt_props_heterog%flux)
        do i=1,n_flux
            write(1,"(20x,ES15.5)") this%tpt_props_heterog%flux(i)
        end do
    end if
    if (this%tpt_props_heterog%homog_flag==.true.) then
        write(1,"(2x,'Stability parameters:'/)")
        write(1,"(10x,'Critical time step:',/)")
        write(1,"(20x,ES15.5)") this%stab_params_tpt%Delta_t_crit
        write(1,"(/,10x,'Peclet:'/)")
        write(1,"(20x,ES15.5)") this%stab_params_tpt%Peclet
    end if
    write(1,"(2x,'F:'/)") 
    do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
        write(1,"(2x,F15.5)") this%F_mat%diag(i)
    end do
    write(1,"(/,2x,'Transition matrix T (with BCs):'/)") 
    write(1,"(17x,2F15.5)") this%trans_mat%diag(1), this%trans_mat%super(1)    
    do i=2,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag-1
        write(1,"(2x,3F15.5)") this%trans_mat%sub(i-1), this%trans_mat%diag(i), this%trans_mat%super(i)
    end do
    write(1,"(2x,2F15.5/)") this%trans_mat%sub(this%spatial_discr%Num_targets-1), this%trans_mat%diag(this%spatial_discr%Num_targets)
    
    write(1,"(/,2x,'Source term g:'/)")
    do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
        write(1,"(2x,F15.5)") this%source_term_PDE(i)
    end do
   
    if (size(output,1)==this%spatial_discr%Num_targets-this%spatial_discr%targets_flag) then
        if (this%time_discr%int_method<4) then
            write(1,"(/,2x,'Cell',*(ES20.5)/)") (Time_out(k), k=1,Num_output)
            do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
                write(1,"(2x,I4,*(F20.5))") i,(output(i,k), k=1,Num_output)
            end do
        else
            write(1,"(/,2x,'Cell',2ES20.5/)") Time_out(1), this%time_discr%Final_time
            do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
                write(1,"(2x,I5,2ES20.10)") i, this%conc_init(i), this%conc(i)
            end do
        end if
    else
        write(1,"(/,2x,'Mobile zone:'/)")
        write(1,"(10x,'Cell',3ES20.5/)") (Time_out(k), k=1,Num_output)
        do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
            write(1,"(10x,I4,3F20.5)") i,(output(i,k), k=1,Num_output)
        end do
        write(1,"(/,2x,'Immobile zones:'/)")
        do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
            write(1,"(10x,I4,3F20.5)") i,(output(this%spatial_discr%Num_targets-this%spatial_discr%targets_flag+i,k), k=1,Num_output)
        end do
    end if
    rewind(1)
    close(1)
end subroutine