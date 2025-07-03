subroutine read_transport_data_WMA(this,unit,file_tpt)!,tpt_props,BCs,mesh,time_discr)
    use transport_transient_m
    class(transport_1D_transient_c) :: this
    integer(kind=4), intent(in) :: unit
    character(len=*), intent(in) :: file_tpt
    !type(matrix_real_c), intent(out) :: mixing_ratios
    !real(kind=8), intent(out), allocatable :: f_vec(:)
    !type(tpt_props_heterog_c), intent(out) :: tpt_props
    !type(BCs_t), intent(out) :: BCs
    !!real(kind=8), intent(in) :: Delta_x
    !type(mesh_1D_Euler_homog_c), intent(out) :: mesh ! homogeneous Euler mesh 1D
    !type(time_discr_homog_c), intent(out) :: time_discr ! homogeneous time discretisation
    
    character(len=256) :: str,label
    character(len=:), allocatable :: str_trim
    type(mesh_1D_Euler_homog_c) :: mesh
    type(time_discr_homog_c) :: time_discr
    type(BCs_t) :: BCs
    type(tpt_props_heterog_c) :: tpt_props
    integer(kind=4) :: i,j,BCs_label(2),comm_ind
    integer(kind=4), allocatable :: num_mix_waters(:),mix_wat_indices(:,:)
    real(kind=8) :: Delta_x,Delta_t,r,phi,D
    !real(kind=8), allocatable :: mixing_ratios(:,:)
    logical :: flag_props
    type(diag_matrix_c) :: F_mat
    
    !type(transport_1D_transient_c) :: this
    
    
    
    open(unit,file=file_tpt,status='old',action='read')
    do 
        read(unit,*) label
        if (label=='end') then
            rewind(unit)
            exit
        !else if (label=='MESH') then
        !    read(unit,*) Num_cells
        !    read(unit,*) Delta_x
        !    call mesh%set_targets(Num_cells,0)
        !    call mesh%set_Delta_x_homog(Delta_x)
        !    call mesh%compute_measure()
        !    call this%set_spatial_discr(mesh)
        !else if (label=='TIME') then
        !    read(unit,*) Num_time
        !    read(unit,*) Delta_t
        !    read(unit,*) int_method
        !    call time_discr%set_Num_time(Num_time)
        !    call time_discr%set_Delta_t_homog(Delta_t)
        !    call time_discr%set_int_method(int_method)
        !    call time_discr%compute_Final_time()
        !    call this%set_time_discr(time_discr)
        else if (label=='BOUNDARY CONDITIONS') then
            read(unit,*) BCs_label(1), BCs_label(2)
            read(unit,*) BCs%evap
            call BCs%set_BCs_label(BCs_label)
            call this%set_BCs(BCs)
        else if (label=='TRANSPORT PROPERTIES') then
            tpt_props%homog_flag=.true.
            !read(unit,*) flag_props, r
            !if (flag_props==.true.) then
            !    allocate(tpt_props%source_term(Num_cells))
            !    tpt_props%source_term=r
            !    call tpt_props%set_source_term_flag(BCs)
            !else
            !    tpt_props%homog_flag=.false.
            !    error stop
            !end if
            read(unit,*) flag_props, phi
            if (flag_props==.true.) then
                allocate(tpt_props%porosity(this%spatial_discr%Num_targets-this%spatial_discr%targets_flag),this%F_mat%diag(this%spatial_discr%Num_targets-this%spatial_discr%targets_flag))
                tpt_props%porosity=phi
                this%F_mat%diag=phi
            else
                tpt_props%homog_flag=.false.
                error stop
            end if
            read(unit,*) flag_props, D
            if (flag_props==.true.) then
                allocate(tpt_props%dispersion(this%spatial_discr%Num_targets-this%spatial_discr%targets_flag))
                tpt_props%dispersion=D
            else
                tpt_props%homog_flag=.false.
                error stop
            end if
            !read(unit,*) flag_props, q
            !if (flag_props==.true.) then
            !    allocate(tpt_props%porosity(this%spatial_discr%Num_targets-this%spatial_discr%targets_flag))
            !    tpt_props%porosity=phi
            !else
            !    error stop
            !end if
            call this%set_tpt_props_heterog_obj(tpt_props)
        else if (label=='MIXING RATIOS') then
            i=0 ! counter targets
            call this%allocate_mixing_ratios()
            !this%mixing_ratios(1,1)=0d0
            !read(unit,*) (this%mixing_ratios(i,j), j=2,4)
            do
                i=i+1
                if (i<=this%spatial_discr%Num_targets) then
                    read(unit,*) this%mixing_ratios%cols(i)%dim
                    allocate(this%mixing_ratios%cols(i)%col_1(this%mixing_ratios%cols(i)%dim))
                else
                    exit
                end if
            end do
            !read(unit,*) (this%mixing_ratios(this%spatial_discr%Num_targets,j), j=1,2), this%mixing_ratios(this%spatial_discr%Num_targets,4)!, f_vec(this%spatial_discr%Num_targets-this%spatial_discr%targets_flag)
            !this%mixing_ratios(this%spatial_discr%Num_targets,3)=0d0
        !else if (label=='MIXING WATERS') then
            !i=1 ! counter targets
            !allocate(num_mix_waters(this%spatial_discr%Num_targets-this%spatial_discr%targets_flag))
            !do
            !    read(unit,"(A10)") str
            !    if (str=='*') then
            !        if (i==this%spatial_discr%Num_targets-this%spatial_discr%targets_flag+1) then
            !            rewind(unit)
            !            exit
            !        else
            !            error stop
            !        end if
            !    else
            !        str_trim=trim(str)
            !        comm_ind=index(str_trim,'!')
            !        if (comm_ind>0) then
            !            str_trim=str_trim(1:comm_ind-1)
            !        else
            !            continue
            !        end if
            !        num_mix_waters(i)=nint((len(str_trim)+1)/2d0)
            !    end if
            !    i=i+1
            !end do
        else
            continue
        end if
    end do
    !do
        !read(unit,*) label
        !if (label=='end') then
        !    exit
        !else if (label=='MIXING RATIOS') then
        !    i=1 ! counter targets
        !    allocate(mix_wat_indices(this%spatial_discr%Num_targets-this%spatial_discr%targets_flag,maxval(num_mix_waters)))
        !    mix_wat_indices=0
        !    do
        !        read(unit,"(A10)") str
        !        if (str=='*') then
        !            if (i==this%spatial_discr%Num_targets-this%spatial_discr%targets_flag+1) then
        !                exit
        !            else
        !                error stop
        !            end if
        !        else
        !            read(unit,*) (mix_wat_indices(i,j), j=1,num_mix_waters(i))
        !        end if
        !        i=i+1
        !    end do
        !else
        !    continue
        !end if
    !end do
    !deallocate(mix_wat_indices)
    do 
        read(unit,*) label
        if (label=='end') then
            exit
        else if (label=='MIXING RATIOS') then
            i=0 ! counter targets
            do
                i=i+1
                if (i<=this%spatial_discr%Num_targets) then
                    read(unit,*) num_targets, (this%mixing_ratios%cols(i)%col_1(j), j=1,num_targets)
                else
                    exit
                end if
            end do
        else
            continue
        end if
    end do
    close(unit)
    do i=1,this%spatial_discr%Num_targets
        print *, this%mixing_ratios%cols(i)%dim
    end do
end subroutine