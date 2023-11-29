module eigenvectors_m
    use vectors_m
    implicit none
    save
    interface
        subroutine eigenvectors_tridiag_sym_matrix(a,b,lambda,v)
            implicit none
            real(kind=8), intent(in) :: a(:) ! diagonal elements
            real(kind=8), intent(in) :: b(:) ! non-diagonal elements
            real(kind=8), intent(in) :: lambda(:) ! eigenvalues
            real(kind=8), intent(out) :: v(:,:) ! eigenvectors
        end subroutine
    end interface
end module