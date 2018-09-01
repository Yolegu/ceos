program ceos
    use pure_mod
    use math_mod
    use gnuplot_mod
    implicit none

    type(pure_type) :: pure
    integer :: eos_id, alpha_id
    real(8) :: tc, pc, acen
    real(8) :: v(3)
    integer :: n

    call read_options(eos_id, tc, pc, acen, alpha_id)
    call pure%init(eos_id, tc, pc, acen, alpha_id)
    call pure%show_isotherm([120.d0, 130.d0, 140.d0])

end program

subroutine read_options(eos_id, tc, pc, acen, alpha_id)
    implicit none

    integer, intent(out) :: eos_id
    real(8), intent(out) :: tc
    real(8), intent(out) :: pc
    real(8), intent(out) :: acen
    integer, intent(out) :: alpha_id

    open(20, file = 'options.txt', status = 'old')

    read(20,*)
    read(20,*)eos_id
    read(20,*)
    read(20,*)tc
    read(20,*)
    read(20,*)pc
    pc = pc * 1.d5
    read(20,*)
    read(20,*)acen
    read(20,*)
    read(20,*)alpha_id

    close(20)

end subroutine

