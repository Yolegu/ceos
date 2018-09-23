program ceos
    use read_file_mod
    use pure_mod
    use mix_mod
    use math_mod
    use gnuplot_mod
    implicit none

    type(pure_type), allocatable :: pure(:)
    type(options_type) :: options
    integer :: i
    real(8) :: psat

    call read_options(options)

    allocate(pure(options%nb_comp))

    do i = 1, options%nb_comp
        call pure(i)%init(options, i)
    end do

    call pure(1)%show_isotherm([0.95d0 * pure(1)%tc])

end program
