program main
  implicit none

  real, allocatable, dimension(:,:) :: g_matrix
  real :: x0_cond, xn_cond, y0_cond, yn_cond
  integer :: n

  n = 10
  allocate(g_matrix(n, n))
  g_matrix = 0

  call apply_bound_conditions(n, g_matrix, 0.0, 0.0, 0.0, 0.0)
  call build_g_matrix(n, g_matrix)

  write (*,*) g_matrix

contains

  subroutine apply_bound_conditions(n, g_matrix, x0_cond, xn_cond, y0_cond, &
                                    yn_cond)
    implicit none

    integer, intent(in) :: n
    real, dimension(n, n), intent(inout) :: g_matrix
    real, intent(in) :: x0_cond, xn_cond, y0_cond, yn_cond

    g_matrix(:, 1) = x0_cond
    g_matrix(:, n) = xn_cond
    g_matrix(1, :) = y0_cond
    g_matrix(n, :) = yn_cond

    return
  end subroutine apply_bound_conditions

  subroutine build_g_matrix(n, g_matrix)
    integer :: n, i, j
    real, dimension(n, n) :: g_matrix
    real(8) :: x, y, xystep

    xystep = 1.0 / real(n - 1)
    do i = 2, n - 1
      do j = 2, n - 1
        g_matrix(i, j) = g(xystep*(i - 1), xystep*(j - 1))
      end do
    end do

  end subroutine build_g_matrix

  real function g(x, y)
    real(8) :: x, y, pi
    pi = 3.14
    g = (1/(5**3*(2*pi)**(1/3)))*exp((-(20*x-10)**2-(20*y-10)**2)/(2*5**2))
  end function g

end program main
