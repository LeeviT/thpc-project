program main
  implicit none
  real, allocatable, dimension(:,:) :: g_matrix, f_matrix, red, black
  real :: x0_cond, xn_cond, y0_cond, yn_cond
  integer :: n, i

  n = 10
  do i=1, n-1
    write(*,*) i
  end do

  allocate(g_matrix(n, n), f_matrix(n, n))
  allocate(red((n-2)/2, (n-2)/2), black((n-2)/2, (n-2)/2))
  write (*,*) "done"
  g_matrix = 0
  !f_matrix = 0
  red = 0
  black = 0
  write (*,*) "done"
  call apply_bound_conditions(n, g_matrix, 0.0, 0.0, 0.0, 0.0)
  call build_g_matrix(n, g_matrix)
  call build_red(n, g_matrix, black, red, 0.5)

  write (*,*) red

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
    implicit none
    integer :: n, i, j
    real, dimension(n, n) :: g_matrix
    real :: x, y, xystep

    xystep = 1.0 / real(n - 1)
    do i = 2, n - 1
      do j = 2, n - 1
        g_matrix(i, j) = g(xystep*(i - 1), xystep*(j - 1))
      end do
    end do
  end subroutine build_g_matrix

  subroutine build_red(n, g_matrix, black, red, relax_par)
    implicit none
    integer, intent(in) :: n
    integer :: i, j, istart
    real, intent(in) :: relax_par
    real :: bound = 0
    real, dimension(n, n), intent(in) :: g_matrix
    real, dimension(n/2, n/2), intent(in) :: black
    real, dimension(n/2, n/2), intent(inout) :: red

    do j = 2, n-1
      if (mod(j, 2) == 0) then
        istart = 3
      else
        istart = 2
      end if
      do i = istart, n-1, 2
        write(*,*) i,j
        if (j == 2) then
          red(i/2, j/2) = (1 - relax_par)*red(i/2, j/2) + (relax_par/4)* &
          (black(i/2+1, i/2) + black(i/2-1, i/2) + black(i/2, i/2+1) + &
          bound)- (relax_par*g_matrix(i, j))/(4*N**2)
        elseif (j == n-1) then
          red(i/2, j/2) = (1 - relax_par)*red(i/2, j/2) + (relax_par/4)* &
          (black(i/2+1, i/2) + black(i/2-1, i/2) + bound + &
          black(i/2, i/2 - 1))- (relax_par*g_matrix(i, j))/(4*N**2)
        elseif (i == 2) then
          red(i/2, j/2) = (1 - relax_par)*red(i/2, j/2) + (relax_par/4)* &
          (black(i/2+1, i/2) + bound + black(i/2, i/2+1) + &
          black(i/2, i/2 - 1))- (relax_par*g_matrix(i, j))/(4*N**2)
        elseif (i == n-1) then
          red(i/2, j/2) = (1 - relax_par)*red(i/2, j/2) + (relax_par/4)* &
          (bound + black(i/2-1, i/2) + black(i/2, i/2+1) + &
          black(i/2, i/2 - 1))- (relax_par*g_matrix(i, j))/(4*N**2)
        else
          red(i/2, j/2) = (1 - relax_par)*red(i/2, j/2) + (relax_par/4)* &
          (black(i/2+1, i/2) + black(i/2-1, i/2) + black(i/2, i/2+1) + &
          black(i/2, i/2 - 1))- (relax_par*g_matrix(i, j))/(4*N**2)
        end if
      end do
    end do
    return
  end subroutine build_red

  real function g(x, y)
    real :: x, y, pi
    pi = 3.14
    g = (1/(5**3*(2*pi)**(1/3)))*exp((-(30*x-15)**2-(30*y-15)**2)/(2*5**2))
  end function g

end program main
