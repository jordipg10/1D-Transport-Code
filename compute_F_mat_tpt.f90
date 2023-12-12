subroutine compute_F_mat_tpt(this)
! F_ii=phi_i
    use transport_transient_m
    use transport_properties_heterog_m
    implicit none
    class(transport_1D_transient_c) :: this
    
    integer(kind=4) :: n
    
    n=this%spatial_discr%Num_targets
     
    !select type (props=>this%props)
    !type is (tpt_props_homog_c)
    if (this%spatial_discr%adapt_ref==1) then
        deallocate(this%F_mat%diag)
        call this%F_mat%allocate_matrix(n)
    end if
    
    
    this%F_mat%diag=this%tpt_props_heterog%porosity!*this%tpt_props_heterog%retardo
    !type is (tpt_props_heterog_c)
        !this%F_mat%diag=props%porosity*props%retardo
    !type is (props_viena_c)
    !    call this%F_mat%allocate_matrix(1)
    !    this%F_mat%diag=props%porosity
    !end select
        
end subroutine