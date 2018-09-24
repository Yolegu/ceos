module read_file_mod
implicit none

    type :: options_type
        character(100) :: txt
        character(100) :: eos_id
        integer :: nb_comp
        real(8) :: r1, r2
        character(100), allocatable :: name(:), alpha_id(:), c_id(:)
        real(8), allocatable :: tc(:), pc(:), acen(:), alpha_param(:,:), c_param(:,:)
    end type

contains

    subroutine read_options(options)

        type(options_type) :: options
        character(100) :: txt
        integer :: i, nb_comp

        open(20, file = 'options.txt', status = 'old')

        read(20,*) txt, options%eos_id

        if (trim(options%eos_id) == "USER") then
            backspace(20)
            read(20,*) txt, options%eos_id, options%r1, options%r2
        end if

        read(20,*) txt, nb_comp
        options%nb_comp = nb_comp

        allocate(options%name(nb_comp), options%tc(nb_comp), options%pc(nb_comp), options%acen(nb_comp)&
        &, options%alpha_id(nb_comp), options%c_id(nb_comp))
        allocate(options%alpha_param(nb_comp, 10), options%c_param(nb_comp, 10))

        do i = 1, nb_comp
            read(20,*)
            read(20,*)txt, options%name(i)
            read(20,*)txt, options%tc(i)
            read(20,*)txt, options%pc(i)
            read(20,*)txt, options%acen(i)

            read(20,*)txt, options%alpha_id(i)

            selectcase(trim(options%alpha_id(i)))
                case("TWU 88")
                    backspace(20)
                    read(20,*)txt, options%alpha_id(i), options%alpha_param(i, 1:2)
                case("TWU 91")
                    backspace(20)
                    read(20,*)txt, options%alpha_id(i), options%alpha_param(i, 1:3)
            end select

            read(20,*)txt, options%c_id(i)

            selectcase(trim(options%c_id(i)))
                case("CONSTANT")
                    backspace(20)
                    read(20,*)txt, options%c_id(i), options%c_param(i, 1)
            end select

        end do

        close(20)

    end subroutine

end module
