program main
    use MRMT_m
    use transport_transient_m
    use transport_m
    implicit none
! Variables
    class(PDE_1D_c), pointer :: my_PDE=>null()
    type(diffusion_1D_c), target :: my_diff
    type(transport_1D_c), target :: my_tpt
    type(diffusion_1D_transient_c), target :: my_diff_transient
    type(transport_1D_transient_c), target :: my_tpt_transient
    
    class(PDE_model_c), pointer :: my_model=>null()
    type(MRMT_c), target :: MRMT
    type(PDE_model_c), target :: st_model
        
    integer(kind=4) :: eqn_flag,opcion,info,method,model,i
    real(kind=8) :: res
    real(kind=8), allocatable :: MRMT_output(:,:),prueba(:)
    character(len=200) :: root
!****************************************************************************************************************************************************
    eqn_flag=4 ! 1: dif stat, 2: dif trans, 3: tpt stat, 4: tpt trans
    model=1 ! 1: traditional, 2: MRMT
    method=1 ! 1: numerical in space & time, 2: eigendecomposition
    root='C:\Users\user2319\OneDrive\Documentos\IDAEA\fortran\codigo\vscode\examples' !> root of the input files
    if (eqn_flag==1) then
        my_PDE=>my_diff
    else if (eqn_flag==2) then
        my_PDE=>my_diff_transient
    else if (eqn_flag==3) then
        my_PDE=>my_tpt
    else if (eqn_flag==4) then
        my_PDE=>my_tpt_transient
    end if
    call my_PDE%set_sol_method(method)
    
    if (model==1) then
        my_model=>st_model
    else if (model==2) then
        my_model=>MRMT
    end if
    call my_model%set_PDE(my_PDE)
    call my_model%PDE%main_PDE(root)
end program