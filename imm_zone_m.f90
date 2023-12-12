! Immobile zone module
module imm_zone_m
    implicit none
    save
    type, public :: imm_zone_c
        real(kind=8) :: conc
        real(kind=8) :: imm_por ! immobile porosity
        real(kind=8) :: exch_rate ! (alpha)
        real(kind=8) :: res_time ! residence time (tau)
        real(kind=8) :: prob ! probability
        real(kind=8) :: flux ! flux of immobile zone
    end type
end module