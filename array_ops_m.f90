! Array operations module
module array_ops_m
    implicit none
    save
    contains
        subroutine append_int_1D_array(array,new_elem) ! appends element in integer vector
            implicit none
            integer(kind=4), intent(inout), allocatable :: array(:)
            integer(kind=4), intent(in) :: new_elem
            
            integer(kind=4) :: i
            integer(kind=4), allocatable :: aux_array(:)
            
            print *, size(array)
            aux_array=array
            if (allocated(array)) then
                deallocate(array)
            end if
            allocate(array(size(aux_array)+1))
            do i=1,size(array)-1
                array(i)=aux_array(i)
            end do
            array(size(array))=new_elem
        end subroutine
        
        
end module