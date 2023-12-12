! Mobile zone module
module mob_zone_m
    implicit none
    save
    type, public :: mob_zone_c
        real(kind=8), allocatable :: conc(:)
        real(kind=8) :: mob_por ! mobile porosity
    end type
end module