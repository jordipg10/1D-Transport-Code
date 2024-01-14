subroutine solve_PDE_EI_Delta_t_homog(this,theta,Time_out,output)
    ! Solves 1D transient PDE with homogeneous time step using Lagr implicit method
    
    ! this: transient PDE object
    ! theta: time weighting factor
    ! Time_out: output time values
    ! output: concentration vs time output
    
    ! Results at all intermediate steps are written in binary mode in file conc_binary_EE.txt
    
    use BCs_subroutines_m
    use metodos_sist_lin_m
    implicit none
    
    ! Variables
    class(PDE_1D_transient_c) :: this
    real(kind=8), intent(in) :: theta
    real(kind=8), intent(in) :: Time_out(:)
    real(kind=8), intent(out) :: output(:,:)

    integer(kind=4) :: n,i,icol,k,out_freq,conc_r_flag,source_term_flag,Num_output
    real(kind=8) :: Time
    real(kind=8), parameter :: epsilon=1d-9
    real(kind=8), allocatable :: conc_old(:),conc_new(:),b(:)
    type(tridiag_matrix_c) :: E_mat,B_mat,A

    procedure(Dirichlet_BCs_PDE), pointer :: p_BCs=>null()
    
    n=this%spatial_discr%Num_targets
    select type (this)
    class is (diffusion_1D_transient_c)
        select type (time_discr=>this%time_discr)
        type is (time_discr_homog_c)
            conc_old=this%conc_init
            allocate(conc_new(n),b(n))
        ! BCs pointer
            if (this%BCs%BCs_label(1)==1 .and. this%BCs%BCs_label(2)==1) then
                !call Dirichlet_BCs_PDE(this)
                p_BCs=>Dirichlet_BCs_PDE
            else if (this%BCs%BCs_label(1)==2 .and. this%BCs%BCs_label(2)==2) then
                !call Neumann_homog_BCs(this)
                p_BCs=>Neumann_homog_BCs
            else
                error stop "Boundary conditions not implemented yet"
            end if
            call this%compute_A_mat_lin_syst(theta,A)
            Num_output=size(Time_out)
        ! Implicit Lagr
            open(unit=0,file="conc_binary_EI.txt",form="unformatted",access="sequential",status="unknown") 
            icol=1
            Time=0
            if (abs(Time-Time_out(icol))<epsilon) then
                output(:,icol)=conc_old
                icol=icol+1
            end if
            do k=1,time_discr%Num_time
                Time=k*time_discr%Delta_t
                write(0) Time, conc_old
            ! Linear system
                call this%compute_b_lin_syst(theta,conc_old,b,k)
                call Thomas(A,b,conc_new)
                if (abs(Time-Time_out(icol))<epsilon) then
                    output(:,icol)=conc_new
                    icol=icol+1
                    if (icol>Num_output) then
                        write(*,*) "Reached Num_output"
                        exit
                    end if
                end if
                conc_old=conc_new
            end do
            this%conc=conc_new
            deallocate(conc_old,conc_new)
            close(0)
        end select
    end select
end subroutine 