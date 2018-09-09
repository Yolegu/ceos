module mix_mod
    use pure_mod
    implicit none
    private

    type, public :: mix_type
        type(pure_type), allocatable :: pure(:)
    end type

contains

    function a_quad(self, t, x) result(a)

        class(mix_type) :: self
        real(8), intent(in) :: t
        real(8), intent(in) :: x(:)
        real(8) :: a
        integer :: i, j

        a = 0.d0
        do i = 1, size(x)
            do j = 1, size(x)
                a = a + x(i) * x(j) * sqrt(self%pure(i)%a(t) * self%pure(j)%a(t))
            end do
        end do

    end function

end module
