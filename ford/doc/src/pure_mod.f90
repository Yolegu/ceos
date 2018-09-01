module pure_mod

    use constants_mod
    use math_mod
    use gnuplot_mod
    implicit none
    private

    type, public :: pure_type
        private
        real(8), public :: tc
        real(8), public :: pc
        real(8), public :: acen
        real(8) :: r1 !! ceos universal constant
        real(8) :: r2
        real(8) :: ac
        real(8) :: bc
        real(8) :: zc
        real(8) :: vc
        procedure(alpha_int), pointer, private :: alpha, dalpha_dtr, d2alpha_dtr2
        real(8) :: m_soave
    contains
        private
        procedure, public :: init
        procedure, public :: show_isotherm
        procedure :: a, da_dt, d2a_dt2
        procedure :: p, dp_dt, dp_dv
        procedure :: solve_eos
        procedure :: m_soave_rk, m_soave_pr
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

        !! Apply the critical specifications to estimate the pure compound parameters. Full details are given in the README section.

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
            self%dalpha_dtr => dalpha_one_dtr
            self%d2alpha_dtr2 => d2alpha_one_dtr2
        elseif (alpha_id == 1) then
            if (eos_id == 1) then
                self%m_soave = self%m_soave_rk()
            elseif (eos_id == 2) then
                self%m_soave = self%m_soave_rk()
            else
                stop 'No Soave function for this equation of state'
            end if
            self%alpha => alpha_soave
            self%dalpha_dtr => dalpha_soave_dtr
            self%d2alpha_dtr2 => d2alpha_soave_dtr2
        else
            stop 'The selected alpha function does not exists'
        end if

    end subroutine

    function p(self, t, v)
        !! Calculate the pressure in Pascal using the generalized cubic equation of state expression
        !! $$P = \frac{RT}{v-b} - \frac{a\left(T\right)}{(v-r_1 b)(v-r_2 b)}$$

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

    function dalpha_one_dtr(self, tr) result(dalpha_dtr)

        class(pure_type) :: self
        real(8), intent(in) :: tr
        real(8) :: dalpha_dtr

        dalpha_dtr = 0.d0

    end function

    function d2alpha_one_dtr2(self, tr) result(d2alpha_dtr2)

        class(pure_type) :: self
        real(8), intent(in) :: tr
        real(8) :: d2alpha_dtr2

        d2alpha_dtr2 = 0.d0

    end function

    function m_soave_rk(self) result(m)

        class(pure_type) :: self
        real(8) :: m

        m = 0.480d0 + 1.574d0 * self%acen - 0.176d0 * self%acen**2

    end function

    function m_soave_pr(self) result(m)

        class(pure_type) :: self
        real(8) :: m

        m = 0.37464d0 + 1.54226d0 * self%acen - 0.26992d0 * self%acen**2

    end function

    function alpha_soave(self, tr) result(alpha)

        class(pure_type) :: self
        real(8) :: alpha
        real(8), intent(in) :: tr

        alpha = (1.d0 + self%m_soave * (1.d0 - sqrt(tr)))**2

    end function

    function dalpha_soave_dtr(self, tr) result(dalpha_dtr)

        class(pure_type) :: self
        real(8), intent(in) :: tr
        real(8) :: dalpha_dtr

        dalpha_dtr = - self%m_soave * (1.d0 - self%m_soave * (sqrt(tr) - 1.d0)) / sqrt(tr)

    end function

    function d2alpha_soave_dtr2(self, tr) result(d2alpha_dtr2)

        class(pure_type) :: self
        real(8), intent(in) :: tr
        real(8) :: d2alpha_dtr2

        d2alpha_dtr2 = - self%m_soave * (self%m_soave + 1.d0) / (2.d0 * tr**(3.d0/2.d0))

    end function

    function a(self, t)

        class(pure_type) :: self
        real(8) :: a
        real(8) :: t

        a = self%ac * self%alpha(t / self%tc)

    end function

    function da_dt(self, t)

        class(pure_type) :: self
        real(8) :: da_dt
        real(8) :: t

        da_dt = self%ac / self%tc * self%dalpha_dtr(t / self%tc)

    end function

    function d2a_dt2(self, t)

        class(pure_type) :: self
        real(8) :: d2a_dt2
        real(8) :: t

        d2a_dt2 = self%ac / self%tc**2 * self%d2alpha_dtr2(t / self%tc)

    end function

    function dp_dt(self, t, v)

        class(pure_type) :: self
        real(8) :: dp_dt
        real(8), intent(in) :: t, v

        dp_dt = rgp / (v - self%bc) - self%da_dt(t) / &
        & ((v - self%r1 * self%bc) * (v - self%r2 * self%bc))

    end function

    function dp_dv(self, t, v)

        class(pure_type) :: self
        real(8) :: dp_dv
        real(8), intent(in) :: t, v

        dp_dv = - rgp * t / (v - self%bc)**2 &
        & + self%a(t) * (2.d0 * v - (self%r1 + self%r2) * self%bc) &
        & / ((v - self%bc * self%r1) * (v - self%bc * self%r2))**2

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

        real(8), intent(in) :: t(:)
        real(8) :: p_arr(n_points), v_arr(n_points)
        integer :: i, j
        type(serie_type), allocatable :: serie(:)

        allocate(serie(size(t)))

        do j = 1, size(t)

            v_arr = linspace(1.01d0 * self%bc, 50.d0 * self%bc, n_points)

            do i = 1, n_points
                p_arr(i) = self%p(t(j), v_arr(i))
                if (p_arr(i) < 2.d0 * self%pc) then
                    v_arr = linspace(v_arr(i), 50.d0 * self%bc, n_points)
                    exit
                end if
            end do

            do i = 1, n_points
                p_arr(i) = self%p(t(j), v_arr(i))
            end do

            call serie(j)%init(v_arr * 1.d6, p_arr * 1.d-5)

        end do

        call plot(serie, "v / [cm^3.mol^{-1}]", "p / bar")

    end subroutine

end module
