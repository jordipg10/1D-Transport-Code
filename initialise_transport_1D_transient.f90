subroutine initialise_transport_1D_transient(this)
    use BCs_subroutines_m
    use transport_stab_params_m
    implicit none

    class(transport_1D_transient_c) :: this
    
    type(tpt_props_heterog_c) :: my_props_tpt
    class(spatial_discr_c), pointer :: my_mesh=>null()
    type(mesh_1D_Euler_homog_c), target :: my_homog_mesh
    type(mesh_1D_Euler_heterog_c), target :: my_heterog_mesh
    class(time_discr_c), pointer :: my_time_discr=>null()
    type(time_discr_homog_c), target :: my_homog_time_discr
    type(time_discr_heterog_c), target :: my_heterog_time_discr
    type(BCs_t) :: my_BCs
    type(stab_params_tpt_c) :: my_stab_params_tpt
    type(char_params_tpt_c) :: my_char_params_tpt
    
    real(kind=8) :: q0,Delta_t,Delta_x,theta,measure,Final_time,x_1,x_2,x_3,x
    real(kind=8), allocatable :: c0(:),c_e(:),source_term_vec(:),porosity_vec(:),dispersion_vec(:),flux_vec(:),flux_coeffs(:)
    real(kind=8), allocatable :: d(:),e(:)
    integer(kind=4) :: parameters_flag,i,Num_cells,Num_time,info,n,adapt_ref_flag,scheme,int_method,r_flag,flux_ord,half_num_cells
    real(kind=8), parameter :: pi=4d0*atan(1d0), epsilon_x=1d-2, epsilon_t=1d-4
    character(len=200) :: filename
    logical :: evap,dimless
!****************************************************************************************************************************************************
! Dimensionless form flag
    dimless=.false.
    this%dimensionless=dimless
! Boundary conditions
    call my_BCs%read_BCs("BCs.dat")
    if (my_BCs%BCs_label(1)==1 .and. my_BCs%BCs_label(2)==1) then
        call my_BCs%read_Dirichlet_BCs("Dirichlet_BCs.dat")
        call my_BCs%read_flux_inf("flux_inflow.dat")
    else if (my_BCs%BCs_label(1)==3) then
        call my_BCs%read_Robin_BC_inflow("Robin_BC_inflow.dat")
    end if
    call this%set_BCs(my_BCs)
 ! Uniform mesh
    my_mesh=>my_homog_mesh
    call my_mesh%read_mesh("Delta_x_homog.dat")
    call this%set_spatial_discr(my_mesh)
    Num_cells=this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
! Uniform time discretisation
    my_time_discr=>my_homog_time_discr
    Delta_t=1d-1*my_homog_mesh%Delta_x**2
    Final_time=5d-1
    int_method=1                                                ! 1: Euler explicit
                                                                ! 2: Euler semi-implicit
                                                                ! 3: Euler fully implicit
                                                                ! 4: Crank-Nicolson
                                                                ! 5: RKF45
    call my_homog_time_discr%set_Delta_t_homog(Delta_t)
    call my_time_discr%set_int_method(int_method)
    call my_time_discr%set_Final_time(Final_time)
    call my_time_discr%compute_Num_time()
    call this%set_time_discr(my_time_discr)
!****************************************************************************************************************************************************
! Transport properties
    call my_props_tpt%read_props("tpt_props.dat",this%spatial_discr)
    if (my_props_tpt%source_term_order==0) then
        call my_props_tpt%compute_flux_lin(this%BCs%flux_inf,this%spatial_discr,this%BCs%flux_out)
    else
        ! esto es una chapuza pero bueno
        open(unit=2,file='flux_coeffs.dat',status='old',action='read')
        read(2,*) flux_ord
        allocate(flux_coeffs(flux_ord+1))
        read(2,*) flux_coeffs
        close(2)
        call my_props_tpt%set_source_term_order(flux_ord-1)
        call my_props_tpt%compute_flux_nonlin(flux_coeffs,this%spatial_discr,this%BCs%flux_out)
        call my_props_tpt%compute_source_term(this%spatial_discr,flux_coeffs)
    end if
    call my_props_tpt%set_source_term_flag(this%BCs)
    call this%set_tpt_props_heterog_obj(my_props_tpt)
!****************************************************************************************************************************************************
! Stability parameters
    call my_stab_params_tpt%compute_stab_params(this%tpt_props_heterog,my_homog_mesh%Delta_x,Delta_t)
    call this%set_stab_params_tpt(my_stab_params_tpt)
!****************************************************************************************************************************************************
! External concentration
    allocate(c_e(Num_cells))
    c_e=0d0
    call this%set_conc_ext(c_e)
    call this%set_conc_star_flag()
! Initial concentration
    allocate(c0(Num_cells))
    half_num_cells=nint(Num_cells/2d0)
    if (this%BCs%BCs_label(1)==1 .and. this%BCs%BCs_label(2)==1) then
        if (this%BCs%evap==.false.) then
            c0(1:half_num_cells)=1
            c0(half_num_cells+1:Num_cells)=0
        else
            x_1=0.749999
            x_2=3.499999
            x_3=3.749999
            Delta_x=this%spatial_discr%get_mesh_size()
            
            do i=1,Num_cells
                x=(2*i-1)*Delta_x/2d0
                if (x < 0.500001) then
                    c0(i)=0d0
                else if (x < x_1) then
                    c0(i)=4*x-2
                else if (x<x_2) then
                    c0(i)=1d0
                else if (x<x_3) then
                    c0(i)=-4*x+15
                else
                    c0(i)=0d0
                end if
            end do
        end if
    else
        c0=0d0
    end if
    if (this%spatial_discr%targets_flag==1) then
        if (this%BCs%BCs_label(1)==1) then ! Dirichlet inflow
            c0(1)=this%BCs%conc_inf
        end if
        if (this%BCs%BCs_label(2)==1) then ! Dirichlet outflow
            c0(Num_cells)=this%BCs%conc_out
        end if
    end if
    call this%set_conc_init(c0)
!****************************************************************************************************************************************************
! Matrices y vectors
    call this%allocate_trans_mat()
    call this%compute_trans_mat_PDE()
    call this%compute_source_term_PDE()
    call this%F_mat%allocate_matrix(Num_cells)
    call this%compute_F_mat_PDE()
!****************************************************************************************************************************************************
    nullify(my_mesh,my_time_discr)
end subroutine