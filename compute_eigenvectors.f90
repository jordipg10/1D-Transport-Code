subroutine compute_eigenvectors(this)
    use eigenvectors_m
    use matrices_m
    implicit none
    class(sq_matrix_c) :: this
    
    allocate(this%eigenvectors(size(this%eigenvalues),size(this%eigenvalues)))
    select type (this)
    type is (tridiag_sym_matrix_c)
        call eigenvectors_tridiag_sym_matrix(this%diag,this%sub,this%eigenvalues,this%eigenvectors)
    end select
end subroutine