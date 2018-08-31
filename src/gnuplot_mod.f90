module gnuplot_mod

    implicit none
    private

    character(20), parameter :: instruction_file = "instructions.gnu"

    public :: plot

contains

    subroutine plot(x, y)

        real(8), intent(in) :: x(:), y(:)
        integer :: i

        open(20, file = 'data.gnu')

        do i = 1, size(x)
            write(20,*)x(i), y(i)
        end do

        call write_header()
        open(30, file = instruction_file, position = 'append')
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
