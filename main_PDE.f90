subroutine main_PDE(this)
    use BCs_subroutines_m
    !use eigenvalues_eigenvectors_m
    !use prob_dens_fcts_m
    implicit none

    ! Variables
    class(PDE_1D_c) :: this
    !integer(kind=4), intent(in), optional :: opcion
    !class(props_c), pointer :: my_props=>null()
    !type(diff_props_c), target :: my_props_diff
    !class(spatial_discr_c), pointer :: my_mesh=>null()
    !type(spatial_discr_rad_c), target :: my_radial_mesh
    !class(time_discr_c), pointer :: my_time_discr=>null()
    !type(time_discr_homog_c), target :: my_homog_time_discr
    !type(time_discr_heterog_c), target :: my_heterog_time_discr
    !type(BCs_t) :: my_BCs
    !class(stab_params_c), pointer :: my_stab_params=>null()
    !type(stab_params_diff_c), target :: my_stab_params_diff
    
    integer(kind=4) :: N_t
    real(kind=8), allocatable :: Time_out(:),Delta_t(:),Delta_r(:),Delta_r_D(:)
    real(kind=8) :: Final_time, measure, Delta_x, length, L2_norm_vi, radius
    real(kind=8) :: a, Delta_r_0
    !real(kind=8) :: source_term,retardo
    !real(kind=8), allocatable :: porosity_vec(:), flux_vec(:), velocity_vec(:), dispersion_vec(:)
    !real(kind=8) :: theta,sigma
    !real(kind=8), allocatable :: A_sub(:),A_diag(:),A_super(:),A_mat_Dirichlet(:,:),T_sub(:),T_diag(:),T_super(:),lambda(:),lambda_star(:),eigenvectors(:,:),eigenvectors_star(:,:),z0(:),A_lambda(:,:,:),A_lambda_v(:)
    !real(kind=8) :: beta,courant,Gamma_sub,Gamma_diag,Gamma_super,radius
    integer(kind=4) :: scheme,int_method,Num_time,BCs,opcion_BCs,parameters_flag,i,j,Num_targets,eqn_flag,dim,niter_max,niter,Dirichlet_BC_location,targets_flag
    !real(kind=8), allocatable :: K(:,:)
    !real(kind=8), allocatable :: c0(:),c_e(:),Delta_x_vec(:),g(:),b(:),y0(:),y(:), Time_out(:)
    !real(kind=8), allocatable :: get_my_parameters(:), Dirichlet_BCs(:)
    !logical :: evap
    !real(kind=8), parameter :: tol=1d-5, pi=4d0*atan(1d0)
    !character(len=200) :: filename
!****************************************************************************************************************************************************
! Pre-process
    call this%initialise_PDE()
!****************************************************************************************************************************************************
! Process
    !print *, this%props%get_props()
    select type (this)
    class is (diffusion_1D_transient_c)
        open(unit=5,file='time_out.dat',status='old',action='read')
        read(5,*) N_t ! number of output times
        allocate(Time_out(N_t))
        read(5,*) Time_out ! output times
        close(5)
        !Time_out=[0d0,this%time_discr%Final_time]
        !select type (time=>this%time_discr)
        !type is (time_discr_heterog_c)
            !Delta_t=this%time_discr%get_Delta_t()
            !Time_out=[0d0,Delta_t]
            !Time_out=[0d0,this%time_discr%Final_time]
        !end select
    type is (diffusion_1D_transient_c)
        !Time_out=[1d-2,1d-1,3d-1,5d-1,1d0] ! dimensionless time
        !call this%solve_write_diffusion_transient(Time_out)
    !class is (PDE_1D_c)
        !Time_out=[]
    end select
    call this%solve_write_PDE_1D(Time_out)
!****************************************************************************************************************************************************
! Post-process
    !nullify(my_mesh,my_time_discr,my_props,my_stab_params)
end subroutine