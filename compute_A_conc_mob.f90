! Computes matrix A in Euler implicit linear system for mobile zone concentrations in MRMT
! A*c_mob^(k+1)=b
subroutine compute_A_mat_conc_mob(this,theta,Delta_t,A_mat)
    use MRMT_m
    use matrices_m
    implicit none
    class(MRMT_c), intent(in) :: this
    real(kind=8), intent(in) :: theta ! time weighting factor
    !real(kind=8), intent(in) :: alpha(:) ! exchange rates
    !real(kind=8), intent(in) :: prob(:) ! probabilities
    real(kind=8), intent(in) :: Delta_t ! time step
    !real(kind=8), intent(in) :: phi_mob ! mobile porosity
    !real(kind=8), intent(in) :: phi_imm ! immmobile porosity
    class(tridiag_matrix_c), intent(out) :: A_mat
    
    integer(kind=8) :: i,j,n_imm
    real(kind=8) :: sum
    
    !n_imm=size(alpha) ! number of immobile zones
    if (theta<0d0 .or. theta>1d0) error stop "Theta must be between 0 and 1"
    A_mat%sub=(-theta*Delta_t/this%mob_zone%mob_por)*this%PDE%trans_mat%sub
    A_mat%super=(-theta*Delta_t/this%mob_zone%mob_por)*this%PDE%trans_mat%super
    sum=0d0
    do i=1,this%n_imm
        sum=sum+this%imm_zones(i)%imm_por*this%imm_zones(i)%prob*this%imm_zones(i)%exch_rate*(theta-(this%imm_zones(i)%exch_rate*Delta_t*theta**2)/(1+this%imm_zones(i)%exch_rate*Delta_t*theta))
    end do
    !do i=1,this%n_imm
        A_mat%diag=1d0-(theta*Delta_t/this%mob_zone%mob_por)*this%PDE%trans_mat%diag+(Delta_t/this%mob_zone%mob_por)*sum
    !end do
end subroutine