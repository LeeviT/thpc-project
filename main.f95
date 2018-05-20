program main
  implicit none

  real, allocatable, dimension(:,:) :: g_matrix
  real :: x0_cond, xn_cond, y0_cond, yn_cond
  integer :: n

  n = 10
  allocate(g_matrix(n, n))
  g_matrix = 0

  call apply_bound_conditions(n, g_matrix, 1.0, 1.0, 1.0, 1.0)

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
    integer :: n
    real, dimension(n, n) :: g_matrix

  end subroutine build_g_matrix

end program main
