function der_Gaussian_pdf_1D(sigma,mesh,Gaussian_pdf) result(der_pdf)
    use spatial_discr_1D_m
    use prob_dens_fcts_m
    implicit none
    real(kind=8), intent(in) :: sigma
    class(spatial_discr_c), intent(in) :: mesh
    real(kind=8), intent(in) :: Gaussian_pdf(:)
    real(kind=8), allocatable :: der_pdf(:)
    
    integer(kind=4) :: i,n
    real(kind=8) :: mu
    
    n=mesh%Num_targets
    if (size(Gaussian_pdf)/=n) error stop "Dimension error in Gaussian pdf"
    
    if (mod(mesh%Num_targets,2)==0) then
        mu=(n-1)/2d0
    else
        mu=floor(n/2d0)
    end if
    allocate(der_pdf(n))
    select type (mesh)
    type is (mesh_1D_Lagr_homog_c)
        forall (i=1:n)
            der_pdf(i)=Gaussian_pdf(i)*(-mesh%Delta_x*(i-1-mu)/sigma**2)
        end forall
    end select
end function