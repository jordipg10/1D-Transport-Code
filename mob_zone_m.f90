! Mobile zone module
module mob_zone_m
    implicit none
    save
    type, public :: mob_zone_c
        real(kind=8), allocatable :: conc(:)
        !integer(kind=4) :: n_imm ! number of immobile zones
        real(kind=8) :: mob_por ! mobile porosity
        !real(kind=8) :: imm_por ! immobile porosity
        !real(kind=8), allocatable :: exchange_rates(:) ! (alpha)
        !real(kind=8), allocatable :: res_times(:) ! residence times (tau)
        !real(kind=8), allocatable :: probs(:) ! probabilities
        !real(kind=8), allocatable :: flux(:) ! flux of immobile zones
    end type
end module