module gnuplot_mod

    implicit none
    private

    character(20), parameter :: instruction_file = "instructions.gnu"

    type, public :: serie_type
        real(8), allocatable :: x(:)
        real(8), allocatable :: y(:)
    contains
        procedure :: init => init_serie
    end type

    public :: plot

contains

    subroutine init_serie(self, x, y)

        class(serie_type) :: self
        real(8) :: x(:), y(:)

        allocate(self%x(size(x)))
        allocate(self%y(size(y)))

        self%x = x
        self%y = y

    end subroutine

    subroutine plot(serie, xlabel, ylabel)

        type(serie_type), intent(in) :: serie(:)
        character(*) :: xlabel, ylabel
        integer :: i, j

        open(20, file = 'data.gnu')

        do j = 1, size(serie)
            do i = 1, size(serie(j)%x)
                write(20,*)serie(j)%x(i), serie(j)%y(i)
            end do
            write(20,*)
            write(20,*)
        end do

        call write_header()
        open(30, file = instruction_file, position = 'append')

        write(30,*)"set xlabel '"//trim(xlabel)//"'"
        write(30,*)"set ylabel '"//trim(ylabel)//"'"
        write(30,*)"plot 'data.gnu' using 1:2 with lines"
        close(30)

        call execute_command_line('gnuplot -p '//trim(instruction_file))

        close(20)

    end subroutine

    subroutine write_header()

        open(30, file = instruction_file)
        write(30,*)"set terminal wxt size 800, 600"
        write(30,*)"set key off"
        write(30,*)"set xtics font 'Consolas,10'"
        write(30,*)"set ytics font 'Consolas,10'"
        close(30)

    end subroutine

end module
