module pure_mod
    use iso_fortran_env
    use constants_mod
    use math_mod
    use gnuplot_mod
    use read_file_mod
    implicit none
    private

    type, public :: pure_type
        private
        character(100) :: eos_id
        real(8), public :: tc
        real(8), public :: pc
        real(8), public :: acen
        real(8) :: r1 !! ceos universal constant
        real(8) :: r2
        real(8) :: ac
        real(8), public :: bc
        real(8) :: zc
        real(8) :: vc
        real(8) :: cc
        procedure(alpha_int), pointer, private :: alpha, dalpha_dtr, d2alpha_dtr2
        procedure(ctrans_int), pointer, private :: c, dc_dtr, d2c_dtr2
        real(8) :: m_soave
        logical :: r1_neq_r2
    contains
        private
        procedure, public :: init
        procedure, public :: show_isotherm
        procedure, public :: show_psat
        procedure, public :: a, da_dt, d2a_dt2
        procedure :: b
        procedure, public :: dc_dt, d2c_dt2
        procedure :: p, dp_dt, dp_dv
        procedure :: s1 ,s2
        procedure, public :: psat
        procedure :: z
        procedure :: a_ec
        procedure :: ln_f, ln_phi
        procedure, public :: solve_eos
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

    abstract interface
        function ctrans_int(self, tr)
            import pure_type
            class(pure_type) :: self
            real(8), intent(in) :: tr
            real(8) :: ctrans_int
        end function
    end interface

contains

    subroutine init(self, options, i)

        !! Apply the critical specifications to estimate the pure compound parameters. Full details are given in the README section.

        class(pure_type) :: self
        type(options_type) :: options
        real(8) :: r1
        real(8) :: r2
        real(8) :: etac
        real(8) :: omegaa
        real(8) :: omegab
        integer :: i

        self%eos_id = options%eos_id

        selectcase(options%eos_id)
            case ("VAN DER WAALS")
                r1 = 0.d0
                r2 = 0.d0
                self%r1_neq_r2 = .false.
            case ("REDLICH-KWONG")
                r1 = -1.d0
                r2 = 0.d0
                self%r1_neq_r2 = .true.
            case ("PENG-ROBINSON")
                r1 = -1.d0 - sqrt(2.d0)
                r2 = -1.d0 + sqrt(2.d0)
                self%r1_neq_r2 = .true.
            case ("USER")
                r1 = options%r1
                r2 = options%r2
            case default
                stop "Unknown equation of state!"
        end select

        self%r1 = r1
        self%r2 = r2

        self%tc = options%tc(i)
        self%pc = options%pc(i) * 1.d5
        self%acen = options%acen(i)

        etac = 1.d0 / (((1.d0 - r1) * (1.d0 - r2)**2)**(1.d0/3.d0) &
        & + ((1.d0 - r2) * (1.d0 - r1)**2)**(1.d0/3.d0) + 1.d0)

        omegaa = (1.d0 - etac * r1) * (1.d0 - etac * r2) * (2.d0 - etac * (r1 + r2)) / &
        & ((1.d0 - etac) * (3.d0 - etac * (1.d0 + r1 + r2))**2)

        self%ac = omegaa * (rgp * self%tc)**2 / self%pc

        omegab = etac / (3.d0 - etac * (1.d0 + r1 + r2))

        self%bc = omegab * rgp * self%tc / self%pc

        self%vc = self%bc / etac

        self%zc = omegab / etac

        selectcase(options%alpha_id(i))
            case ("ONE")
                self%alpha => alpha_one
                self%dalpha_dtr => dalpha_one_dtr
                self%d2alpha_dtr2 => d2alpha_one_dtr2
            case ("INV SQRT")
                self%alpha => alpha_invsqrt
                self%dalpha_dtr => dalpha_invsqrt_dtr
                self%d2alpha_dtr2 => d2alpha_invsqrt_dtr2
            case ("SOAVE 72")
                selectcase(options%eos_id)
                    case ("REDLICH-KWONG")
                        self%m_soave = self%m_soave_rk()
                    case ("PENG-ROBINSON")
                        self%m_soave = self%m_soave_pr()
                    case default
                        stop 'No Soave function for this equation of state'
                end select
                self%alpha => alpha_soave
                self%dalpha_dtr => dalpha_soave_dtr
                self%d2alpha_dtr2 => d2alpha_soave_dtr2
            case default
                stop 'The selected alpha function does not exists'
        end select

        ! Volume translation parameter
        self%cc = options%c_param(i, 1)

        selectcase(options%c_id(i))
            case("CONSTANT")
                self%c => c_cst
                self%dc_dtr => dc_dtr_cst
                self%d2c_dtr2 => d2c_dtr2_cst
            case ("MAGOULAS-TASSIOS")
                self%c => c_magoulas_tassios
            case default
                stop 'The selected volume translation function does not exists'
        end select

        if (self%vc < 0.d0) then
            stop "Negative critical molar volume because of volume translation !"
        end if

    end subroutine

    function s1(self, t)

        class(pure_type) :: self
        real(8), intent(in) :: t
        real(8) :: tr, s1

        tr = t / self%tc
        s1 = self%r1 - self%c(tr) / self%b(t) * (1.d0 - self%r1)

    end function

    function s2(self, t)

        class(pure_type) :: self
        real(8), intent(in) :: t
        real(8) :: tr, s2

        tr = t / self%tc
        s2 = self%r2 - self%c(tr) / self%b(t) * (1.d0 - self%r2)

    end function

    function b(self, t)

        class(pure_type) :: self
        real(8), intent(in) :: t
        real(8) :: b

        b = self%bc - self%c(t / self%tc)

        if (b < 0.d0) then
            stop "Negative covolume because of volume translation !"
        end if

    end function

    function p(self, t, v)
        !! Calculate the pressure in Pascal using the generalized cubic equation of state expression
        !! $$P = \frac{RT}{v-b} - \frac{a\left(T\right)}{(v-r_1 b)(v-r_2 b)}$$

        class(pure_type) :: self
        real(8), intent(in) :: t
        real(8), intent(in) :: v
        real(8) :: p

        p = rgp * t / (v - self%b(t)) - self%a(t) / ((v - self%s1(t) * self%b(t)) * (v - self%s2(t) * self%b(t)))

    end function

    function z(self, t, v)

        class(pure_type) :: self
        real(8), intent(in) :: t, v
        real(8) :: z

        z = self%p(t,v) * v / (rgp * t)

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

    function alpha_invsqrt(self, tr)

        class(pure_type) :: self
        real(8), intent(in) :: tr
        real(8) :: alpha_invsqrt
        alpha_invsqrt = 1.d0 / sqrt(tr)

    end function

    function dalpha_invsqrt_dtr(self, tr) result(dalpha_dtr)

        class(pure_type) :: self
        real(8), intent(in) :: tr
        real(8) :: dalpha_dtr

        dalpha_dtr = -0.5d0 * tr**(-3.d0 / 2.d0)

    end function

    function d2alpha_invsqrt_dtr2(self, tr) result(d2alpha_dtr2)

        class(pure_type) :: self
        real(8), intent(in) :: tr
        real(8) :: d2alpha_dtr2

        d2alpha_dtr2 = 0.75d0 * tr**(-5.d0 / 2.d0)

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

    function c_cst(self, tr) result(c)

        class(pure_type) :: self
        real(8), intent(in) :: tr
        real(8) :: c

        c = self%cc

    end function

    function dc_dtr_cst(self, tr) result(dc_dt)

        class(pure_type) :: self
        real(8), intent(in) :: tr
        real(8) :: dc_dt

        dc_dt = 0.d0

    end function

    function d2c_dtr2_cst(self, tr) result(d2c_dt2)

        class(pure_type) :: self
        real(8), intent(in) :: tr
        real(8) :: d2c_dt2

        d2c_dt2 = 0.d0

    end function

    function c_magoulas_tassios(self, tr) result(c)

        class(pure_type) :: self
        real(8), intent(in) :: tr
        real(8) :: c
        real(8) :: d0, d1, d2, d3, d4, k0, k1, k2, k3, k4, l0, l1
        real(8) :: t0, beta, zc, zc0, tc

        selectcase(self%eos_id)
            case("VAN DER WAALS")
                d0 = 0.483798d0
                d1 = 1.643232d0
                d2 = -0.288718d0
                d3 = 0.066013d0
                d4 = 0.0d0
                k0 = 0.036722d0
                k1 = 0.063541d0
                k2 = -0.076221d0
                k3 = 0.060362d0
                k4 = -0.015772d0
                l0 = -7.099630d0
                l1 = -21.156900d0
                zc0 = 0.375d0
            case("PENG-ROBINSON")
                d0 = 0.384401d0
                d1 = 1.522760d0
                d2 = -0.213808d0
                d3 = 0.0346160
                d4 = -0.001976d0
                k0 = -0.014471d0
                k1 = 0.067498d0
                k2 = -0.084852d0
                k3 = 0.067298d0
                k4 = -0.017366d0
                l0 = -10.244700d0
                l1 = -28.631200d0
                zc0 = 0.3074d0
            case default
                stop "No Magoulas-Tassios volume translation correlation for the equation of state!"
        end select

        t0 = rgp * self%tc / self%pc * (k0 + k1 * self%acen + k2 * self%acen**2 + k3 * self%acen**3 + k4 * self%acen**4)
        beta = l0 + l1 * self%acen
        zc = 0.289d0 - 0.0701d0 * self%acen - 0.0207d0 * self%acen**2
        tc = rgp * self%tc / self%pc * (zc0 - zc)
        c = t0 + (tc - t0) * exp(beta * abs(1.d0 - tr))

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

    function dc_dt(self, t)

        class(pure_type) :: self
        real(8) :: dc_dt
        real(8) :: t

        dc_dt = 1.d0 / self%tc * self%dc_dtr(t / self%tc)

    end function

    function d2c_dt2(self, t)

        class(pure_type) :: self
        real(8) :: d2c_dt2
        real(8) :: t

        d2c_dt2 = 1.d0 / self%tc**2 * self%d2c_dtr2(t / self%tc)

    end function

    function dp_dt(self, t, v)

        class(pure_type) :: self
        real(8) :: dp_dt
        real(8), intent(in) :: t, v

        dp_dt = rgp / (v - self%b(t)) - self%da_dt(t) / &
        & ((v - self%s1(t) * self%b(t)) * (v - self%s2(t) * self%b(t)))

    end function

    function dp_dv(self, t, v)

        class(pure_type) :: self
        real(8) :: dp_dv
        real(8), intent(in) :: t, v

        dp_dv = - rgp * t / (v - self%b(t))**2 &
        & + self%a(t) * (2.d0 * v - (self%s1(t) + self%s2(t)) * self%b(t)) &
        & / ((v - self%b(t) * self%s1(t)) * (v - self%b(t) * self%s2(t)))**2

    end function

    function a_ec(self, t, v)

        class(pure_type) :: self
        real(8), intent(in) :: t, v
        real(8) :: a_ec

        if (self%r1_neq_r2) then
            a_ec = rgp * t * log(v / (v - self%b(t))) &
            & + self%a(t) / (self%b(t) * (self%s1(t) - self%s2(t))) &
            & * log((v - self%b(t) * self%s1(t)) / (v - self%b(t) * self%s2(t)))
        else
            a_ec = rgp * t * log(v / (v - self%b(t))) - self%a(t) / (v - self%s1(t))
        end if

    end function

    function ln_phi(self, t, v)

        class(pure_type) :: self
        real(8), intent(in) :: t, v
        real(8) :: ln_phi

        if (self%r1_neq_r2) then
            ln_phi = self%p(t, v) * v / (rgp * t) - 1.d0 &
            & - log(self%p(t, v) * (v - self%b(t)) / (rgp * t)) &
            & + self%a(t) / (rgp * t * self%b(t) * (self%s1(t) - self%s2(t))) &
            & * log((v - self%s1(t) * self%b(t))/(v - self%s2(t) * self%b(t)))
        else
            ln_phi = self%p(t, v) * v / (rgp * t) - 1.d0 - log(self%p(t, v) * (v - self%b(t)) / (rgp * t)) &
            & - self%a(t) / (rgp * t * v)
        end if

    end function


    function ln_f(self, t, v)

        class(pure_type) :: self
        real(8), intent(in) :: t, v
        real(8) :: ln_f

        if (self%r1_neq_r2) then
            ln_f = self%p(t, v) * v / (rgp * t) - 1.d0 &
            & - log((v - self%b(t)) / (rgp * t)) &
            & + self%a(t) / (rgp * t * self%b(t) * (self%s1(t) - self%s2(t))) &
            & * log((v - self%s1(t) * self%b(t))/(v - self%s2(t) * self%b(t)))
        else
            ln_f = self%p(t, v) * v / (rgp * t) - 1.d0 - log((v - self%b(t)) / (rgp * t)) - self%a(t) / (rgp * t * v)
        end if

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

        b = - (self%b(t) * (self%s1(t) + self%s2(t) + 1.d0) + rgp * t / p)

        c = self%b(t)**2 * (self%s1(t) * self%s2(t) + self%s1(t) + self%s2(t)) &
            & + rgp * t * self%b(t) * (self%s1(t) + self%s2(t)) / p &
            & + self%a(t) / p

        d = -self%b(t) * (self%s1(t) * self%s2(t) * self%b(t)**2 &
            & + self%s1(t) * self%s2(t) * self%b(t) * rgp * t / p + self%a(t) / p)

        coeff(1) = a
        coeff(2) = b
        coeff(3) = c
        coeff(4) = d

        call solve_cubic_polynomial(coeff, v, n)

        v = 1.d0 / v

    end subroutine

    subroutine show_isotherm(self, t)

        class(pure_type) :: self

        integer, parameter :: n_points = 1000

        real(8), intent(in) :: t(:)
        real(8) :: p_arr(n_points), v_arr(n_points), v_liq(n_points), v_int(n_points), v_gas(n_points)
        integer :: i, j
        type(serie_type), allocatable :: serie(:)
        integer :: n
        real(8) :: v(3)

        allocate(serie(size(t)))

        ! Calculation of the isotherm by changing the molar volume at fixed temperature and pressure

        do j = 1, size(t)

            v_arr = linspace(1.01d0 * self%b(t(j)), 50.d0 * self%b(t(j)), n_points)

            do i = 1, n_points
                p_arr(i) = self%p(t(j), v_arr(i))
                if (p_arr(i) < 2.d0 * self%pc) then
                    v_arr = linspace(v_arr(i), 50.d0 * self%b(t(j)), n_points)
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

    function psat(self, t)

        class(pure_type) :: self
        real(8), intent(in) :: t
        real(8) :: psat
        real(8) :: psat_old
        integer :: k
        real(8) :: v(3), v_liq, v_gas
        real(8) :: psi
        integer :: n

        k = 1
        do
            psat = self%p(t, (1.01d0**k) * self%vc)
            if (psat > 0.d0) then
                exit
            else
                k = k + 1
            end if
        end do

        do
            call self%solve_eos(t, psat, v, n)
            if (n /= 3) then
                stop 'Less than three volumes found in the psat calculation procedure!'
            end if

            v_liq = v(1)
            v_gas = v(3)
            psi = self%ln_f(t, v_gas) - self%ln_f(t, v_liq)

            psat_old = psat

            psat = psat_old - rgp * t * psi / (v_gas - v_liq)

            if (psat < 0.d0) then
                psat = psat_old * exp(-rgp * t / psat_old * psi / (v_gas - v_liq))
            end if

            if (abs(psi) < 1.d-10 .and. abs(psat_old - psat) / psat < 1.d-10) then
                exit
            end if

        end do

    end function

    subroutine show_psat(self)

        class(pure_type) :: self

        integer, parameter :: n_points = 1000

        real(8) :: p_arr(n_points), t_arr(n_points)
        integer :: i, j
        type(serie_type) :: serie

        t_arr = linspace(0.40d0 * self%tc, 0.999d0 * self%tc, n_points)

        do i = 1, n_points
            p_arr(i) = self%psat(t_arr(i))
        end do

        call serie%init(t_arr, p_arr * 1.d-5)

        call plot([serie], "T / K", "p / bar")

    end subroutine

end module
