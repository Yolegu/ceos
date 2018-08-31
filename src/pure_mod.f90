module pure_mod

    use constants_mod
    use math_mod
    use gnuplot_mod
    implicit none
    private

    type, public :: pure_type
        real(8) :: r1
        real(8) :: r2
        real(8) :: tc
        real(8) :: pc
        real(8) :: acen
        real(8) :: ac
        real(8) :: bc
        real(8) :: zc
        real(8) :: vc
        procedure(alpha_int), pointer :: alpha
    contains
        procedure :: init
        procedure :: a
        procedure :: p
        procedure :: solve_eos
        procedure :: show_isotherm
    end type

    abstract interface
        function alpha_int(self, tr)
            import pure_type
            class(pure_type) :: self
            real(8), intent(in) :: tr
            real(8) :: alpha_int
        end function
    end interface

contains

    subroutine init(self, eos_id, tc, pc, acen, alpha_id)

        class(pure_type) :: self
        integer, intent(in) :: eos_id
        real(8), intent(in) :: tc
        real(8), intent(in) :: pc
        real(8), intent(in) :: acen
        integer, intent(in) :: alpha_id
        real(8) :: r1
        real(8) :: r2
        real(8) :: etac
        real(8) :: omegaa
        real(8) :: omegab

        if (eos_id == 0) then
            r1 = 0.d0
            r2 = 0.d0
        elseif (eos_id == 1) then
            r1 = -1.d0
            r2 = 0.d0
        elseif (eos_id == 2) then
            r1 = -1.d0 - sqrt(2.d0)
            r2 = -1.d0 + sqrt(2.d0)
        end if

        self%r1 = r1
        self%r2 = r2
        self%tc = tc
        self%pc = pc
        self%acen = acen

        etac = 1.d0 / (((1.d0 - r1) * (1.d0 - r2))**(1.d0/3.d0) &
        & + ((1.d0 - r2) * (1.d0 - r1))**(1.d0/3.d0) + 1.d0)

        omegaa = (1.d0 - etac * r1) * (1.d0 - etac * r2) * (2.d0 - etac * (r1 + r2)) / &
        & ((1.d0 - etac) * (3.d0 - etac * (1.d0 + r1 + r2))**2)
        self%ac = omegaa * (rgp * tc)**2 / pc

        omegab = etac / (3.d0 - etac * (1.d0 + r1 + r2))
        self%bc = omegab * rgp * tc / pc

        self%vc = self%bc / etac

        self%zc = omegab / etac

        if (alpha_id == 0) then
            self%alpha => alpha_one
        end if

    end subroutine

    function p(self, t, v)

        class(pure_type) :: self
        real(8), intent(in) :: t
        real(8), intent(in) :: v
        real(8) :: p

        p = rgp * t / (v - self%bc) - self%a(t) / ((v - self%r1 * self%bc) * (v - self%r2 * self%bc))

    end function

    function alpha_one(self, tr)

        class(pure_type) :: self
        real(8), intent(in) :: tr
        real(8) :: alpha_one
        alpha_one = 1.d0

    end function

    function a(self, t)

        class(pure_type) :: self
        real(8) :: a
        real(8) :: t

        a = self%ac * self%alpha(t / self%tc)

    end function

    subroutine solve_eos(self, t, p, v, n)

        class(pure_type) :: self
        real(8), intent(in) :: t
        real(8), intent(in) :: p
        real(8), intent(out) :: v(3)
        integer, intent(out) :: n

        real(8) :: a, b, c, d
        real(8) :: coeff(4)

        a = 1.d0

        b = - (self%bc * (self%r1 + self%r2 + 1.d0) + rgp * t / p)

        c = self%bc**2 * (self%r1 * self%r2 + self%r1 + self%r2) &
            & + rgp * t * self%bc * (self%r1 + self%r2) / p &
            & + self%a(t) / p

        d = -self%bc * (self%r1 * self%r2 * self%bc**2 &
            & + self%r1 * self%r2 * self%bc * rgp * t / p + self%a(t) / p)

        coeff(1) = d
        coeff(2) = c
        coeff(3) = b
        coeff(4) = a

        call solve_cubic_polynomial(coeff, v, n)

    end subroutine

    subroutine show_isotherm(self, t)

        class(pure_type) :: self

        integer, parameter :: n_points = 1000

        real(8), intent(in) :: t
        real(8) :: p_arr(n_points), v_arr(n_points)
        integer :: i

        v_arr = linspace(1.01d0 * self%bc, 50.d0 * self%bc, n_points)

        do i = 1, n_points
            p_arr(i) = self%p(t, v_arr(i))
            if (p_arr(i) < 2.d0 * self%pc) then
                v_arr = linspace(v_arr(i), 50.d0 * self%bc, n_points)
                exit
            end if
        end do

        do i = 1, n_points
            p_arr(i) = self%p(t, v_arr(i))
        end do

        call plot(v_arr * 1.d6, p_arr * 1.d-5)

    end subroutine

end module
