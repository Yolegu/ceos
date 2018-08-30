module math_mod
    use constants_mod
    implicit none

contains

    subroutine set_cubic_monic(coeff, coeff_monic)

        real(8), intent(in) :: coeff(:)
        real(8), intent(out) :: coeff_monic(size(coeff))

        coeff_monic = coeff / coeff(size(coeff))

    end subroutine

    subroutine cardan(coeff, x, n)

        real(8), intent(in) :: coeff(4)
        real(8), intent(out) :: x(3)
        integer, intent(out) :: n

        real(8) :: a, b, c, d
        real(8) :: q, r, delta, s, t, theta

        a = coeff(4)
        b = coeff(3)
        c = coeff(2)
        d = coeff(1)

        q = (3.d0 * c - b**2) / 9.d0
        r = (9.d0 * b * c - 27.d0 * d - 2.d0 * b**3) / 54.d0
        delta = q**3 + r**2

        if (delta > 0.d0) then
            s = cbrt(r + sqrt(delta))
            t = cbrt(r - sqrt(delta))
            x(1) = s + t - b / 3.d0
            x(2:3) = 0.d0
            n = 1
        else
            theta = acos(r / sqrt(-q**3))
            x(1) = 2.d0 * sqrt(-q) * cos(theta / 3.d0) - b / 3.d0
            x(2) = 2.d0 * sqrt(-q) * cos(theta / 3.d0 + 2.d0 * pi / 3.d0) - b / 3.d0
            x(3) = 2.d0 * sqrt(-q) * cos(theta / 3.d0 + 4.d0 * pi / 3.d0) - b / 3.d0
            n = 3
        end if

    end subroutine

    subroutine newton_after_cardan(coeff, x, n)

        ! Newton's method applied to the result of Cardano's method
        ! The function is f(x) = a x^3 + b x^2 + c x + d
        ! Function's derivative is f'(x) = 3a x^2 + 2b x + c

        real(8), intent(in) :: coeff(4)
        real(8), intent(inout) :: x(3)
        integer, intent(in) :: n
        real(8) :: f
        real(8) :: dfdx
        integer :: i
        real(8) :: a, b, c, d
        real(8), parameter :: eps = 1.d-15

        a = coeff(4)
        b = coeff(3)
        c = coeff(2)
        d = coeff(1)

        do i = 1, n
            f = x(i)**3 + b * x(i)**2 + c * x(i) + d
            if (f > eps) then
                dfdx = 3.d0 * x(i)**2 + 2.d0 * b * x(i) + c
                x(i) = x(i) - f / dfdx
            else
                cycle
            end if
        end do

    end subroutine

    subroutine swap(x, a, b)

        real(8) :: x(:)
        integer, intent(in) :: a, b
        real(8) :: temp

        temp = x(a)
        x(a) = x(b)
        x(b) = temp

    end subroutine

    subroutine sort_cardan_roots(x)
        ! Sort the three roots in decreasing order
        ! src: https://stackoverflow.com/questions/4793251/sorting-int-array-with-only-3-elements

        real(8), intent(inout) :: x(3)
        real(8) :: temp

        if (x(1) < x(2)) then
            call swap(x, 1, 2)
        end if

        if (x(2) < x(3)) then
            call swap(x, 2, 3)
        end if

        if (x(1) < x(2)) then
            call swap(x, 1, 2)
        end if

    end subroutine

    subroutine solve_cubic_polynomial(coeff, x, n)

        real(8), intent(in) :: coeff(4)
        real(8) :: coeff_monic(4)
        real(8), intent(out) :: x(3)
        integer :: n

        call set_cubic_monic(coeff, coeff_monic)
        call cardan(coeff_monic, x, n)
        call newton_after_cardan(coeff_monic, x, n)

        if (n == 3) then
            call sort_cardan_roots(x)
        end if

    end subroutine

    subroutine test_cardan()

        real(8) :: x(3)
        integer :: n

        write(*,*)'Solve x^3 + 6 x^2 + 11 x + 6 = 0'
        write(*,*)'Answer: x(1) = -2'
        write(*,*)'Answer: x(2) = -1'
        write(*,*)'Answer: x(3) = -3'

        call solve_cubic_polynomial([6.d0, 11.d0, 6.d0, 1.d0], x, n)
        write(*,*)'x_cardan = ', x

        write(*,*)
        write(*,*)'Solve 1000 x^3 - 1254 x^2 - 496 x + 191 = 0'
        write(*,*)'Answer: x(1) = 1.4997993055'
        write(*,*)'Answer: x(2) = -0.5003313644'
        write(*,*)'Answer: x(3) = 0.254532059'

        call solve_cubic_polynomial([191.d0, -496.d0, -1254.d0, 1000.d0], x, n)
        write(*,*)'x_cardan = ', x

    end subroutine

    function cbrt(x)

        real(8), intent(in) :: x
        real(8) :: cbrt

        cbrt = sign(abs(x)**(1.d0 / 3.d0), x)

    end function

end module
