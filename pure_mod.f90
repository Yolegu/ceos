module pure_mod

    use constants_mod
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
    contains
        procedure :: init
    end type

contains

    subroutine init(self, r1, r2, tc, pc)

        class(pure_type) :: self
        real(8), intent(in) :: r1
        real(8), intent(in) :: r2
        real(8), intent(in) :: tc
        real(8), intent(in) :: pc
        real(8) :: etac
        real(8) :: omegaa
        real(8) :: omegab

        self%r1 = r1
        self%r2 = r2
        self%tc = tc
        self%pc = pc

        etac = 1.d0 / (((1.d0 - r1) * (1.d0 - r2))**(1.d0/3.d0) &
        & + ((1.d0 - r2) * (1.d0 - r1))**(1.d0/3.d0) + 1.d0)

        omegaa = (1.d0 - etac * r1) * (1.d0 - etac * r2) * (2.d0 - etac * (r1 + r2)) / &
        & ((1.d0 - etac) * (3.d0 - etac * (1.d0 + r1 + r2))**2)
        self%ac = omegaa * (rgp * tc)**2 / pc

        omegab = etac / (3.d0 - etac * (1.d0 + r1 + r2))
        self%bc = omegab * rgp * tc / pc

        self%vc = self%bc / etac

        self%zc = omegab / etac

    end subroutine

end module
