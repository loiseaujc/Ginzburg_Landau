module Ginzburg_Landau
  !> RKLIB module.
  use rklib_module
  !> Mesh/Physical parameters of the simulation.
  use params
  implicit none

  !> Number of function evaluations.
  integer :: fevals = 0

  private
  public :: rhs, adjoint_rhs, fevals

contains

  subroutine rhs(me, t, x, f)
    !> Time-integrator.
    class(rk_class), intent(inout)             :: me
    !> Current time.
    real(kind=wp)  , intent(in)                :: t
    !> State vector.
    real(kind=wp)  , dimension(:), intent(in)  :: x
    !> Time-derivative.
    real(kind=wp)  , dimension(:), intent(out) :: f

    !> Internal variables.
    integer :: i, j, k
    real(kind=wp), dimension(nx) :: u, du
    real(kind=wp), dimension(nx) :: v, dv
    real(kind=wp)                :: d2u, d2v, cu, cv

    !> Sets the internal variables.
    f = 0.0_wp
    u = x(1:nx)      ; du = f(1:nx)
    v = x(nx+1:2*nx) ; dv = f(nx+1:2*nx)

    !---------------------------------------------------
    !-----     Linear Ginzburg Landau Equation     -----
    !---------------------------------------------------

    !> Left most boundary points.
    cu = u(2) / (2*dx) ; cv = v(2) / (2*dx)
    du(1) = -(real(nu)*cu - aimag(nu)*cv) !> Convective term.
    dv(1) = -(aimag(nu)*cu + real(nu)*cv) !> Convective term.

    d2u = (u(2) - 2*u(1)) / dx**2 ; d2v = (v(2) - 2*v(1)) / dx**2
    du(1) = du(1) + real(gamma)*d2u - aimag(gamma)*d2v !> Diffusion term.
    dv(1) = dv(1) + aimag(gamma)*d2u + real(gamma)*d2v !> Diffusion term.

    du(1) = du(1) + mu(1)*u(1) !> Non-parallel term.
    dv(1) = dv(1) + mu(1)*v(1) !> Non-parallel term.

    !> Interior nodes.
    do i = 2, nx-1
       !> Convective term.
       cu = (u(i+1) - u(i-1)) / (2*dx)
       cv = (v(i+1) - v(i-1)) / (2*dx)
       du(i) = -(real(nu)*cu - aimag(nu)*cv)
       dv(i) = -(aimag(nu)*cu + real(nu)*cv)

       !> Diffusion term.
       d2u = (u(i+1) - 2*u(i) + u(i-1)) / dx**2
       d2v = (v(i+1) - 2*v(i) + v(i-1)) / dx**2
       du(i) = du(i) + real(gamma)*d2u - aimag(gamma)*d2v
       dv(i) = dv(i) + aimag(gamma)*d2u + real(gamma)*d2v

       !> Non-parallel term.
       du(i) = du(i) + mu(i)*u(i)
       dv(i) = dv(i) + mu(i)*v(i)
    enddo
    
    !> Right most boundary points.
    cu = -u(nx-1) / (2*dx) ; cv = -v(nx-1) / (2*dx)
    du(nx) = -(real(nu)*cu - aimag(nu)*cv) !> Convective term.
    dv(nx) = -(aimag(nu)*cu + real(nu)*cv) !> Convective term.

    d2u = (-2*u(nx) + u(nx-1)) / dx**2 ; d2v = (-2*v(nx) + v(nx-1)) / dx**2
    du(nx) = du(nx) + real(gamma)*d2u - aimag(gamma)*d2v !> Diffusion term.
    dv(nx) = dv(nx) + aimag(gamma)*d2u + real(gamma)*d2v !> Diffusion term.

    du(nx) = du(nx) + mu(nx)*u(nx) !> Non-parallel term.
    dv(nx) = dv(nx) + mu(nx)*v(nx) !> Non-parallel term.

    !> Copy results to the output array.
    f(1:nx) = du ; f(nx+1:2*nx) = dv

    fevals = fevals + 1
    return
  end subroutine rhs

  subroutine adjoint_rhs(me, t, x, f)
    !> Time-integrator.
    class(rk_class), intent(inout)             :: me
    !> Current time.
    real(kind=wp)  , intent(in)                :: t
    !> State vector.
    real(kind=wp)  , dimension(:), intent(in)  :: x
    !> Time-derivative.
    real(kind=wp)  , dimension(:), intent(out) :: f

    !> Internal variables.
    integer :: i, j, k
    real(kind=wp), dimension(nx) :: u, du
    real(kind=wp), dimension(nx) :: v, dv
    real(kind=wp)                :: d2u, d2v, cu, cv

    !> Sets the internal variables.
    f = 0.0_wp
    u = x(1:nx)      ; du = f(1:nx)
    v = x(nx+1:2*nx) ; dv = f(nx+1:2*nx)

    !---------------------------------------------------
    !-----     Linear Ginzburg Landau Equation     -----
    !---------------------------------------------------

    !> Left most boundary points.
    cu = u(2) / (2*dx) ; cv = v(2) / (2*dx)
    du(1) = (real(nu)*cu + aimag(nu)*cv) !> Convective term.
    dv(1) = (-aimag(nu)*cu + real(nu)*cv) !> Convective term.

    d2u = (u(2) - 2*u(1)) / dx**2 ; d2v = (v(2) - 2*v(1)) / dx**2
    du(1) = du(1) + real(gamma)*d2u + aimag(gamma)*d2v !> Diffusion term.
    dv(1) = dv(1) - aimag(gamma)*d2u + real(gamma)*d2v !> Diffusion term.

    du(1) = du(1) + mu(1)*u(1) !> Non-parallel term.
    dv(1) = dv(1) + mu(1)*v(1) !> Non-parallel term.

    !> Interior nodes.
    do i = 2, nx-1
       !> Convective term.
       cu = (u(i+1) - u(i-1)) / (2*dx)
       cv = (v(i+1) - v(i-1)) / (2*dx)
       du(i) = (real(nu)*cu + aimag(nu)*cv)
       dv(i) = (-aimag(nu)*cu + real(nu)*cv)

       !> Diffusion term.
       d2u = (u(i+1) - 2*u(i) + u(i-1)) / dx**2
       d2v = (v(i+1) - 2*v(i) + v(i-1)) / dx**2
       du(i) = du(i) + real(gamma)*d2u + aimag(gamma)*d2v
       dv(i) = dv(i) - aimag(gamma)*d2u + real(gamma)*d2v

       !> Non-parallel term.
       du(i) = du(i) + mu(i)*u(i)
       dv(i) = dv(i) + mu(i)*v(i)
    enddo

    !> Right most boundary points.
    cu = -u(nx-1) / (2*dx) ; cv = -v(nx-1) / (2*dx)
    du(nx) = (real(nu)*cu + aimag(nu)*cv) !> Convective term.
    dv(nx) = (-aimag(nu)*cu + real(nu)*cv) !> Convective term.

    d2u = (-2*u(nx) + u(nx-1)) / dx**2 ; d2v = (-2*v(nx) + v(nx-1)) / dx**2
    du(nx) = du(nx) + real(gamma)*d2u + aimag(gamma)*d2v !> Diffusion term.
    dv(nx) = dv(nx) - aimag(gamma)*d2u + real(gamma)*d2v !> Diffusion term.

    du(nx) = du(nx) + mu(nx)*u(nx) !> Non-parallel term.
    dv(nx) = dv(nx) + mu(nx)*v(nx) !> Non-parallel term.

    !> Copy results to the output array.
    f(1:nx) = du ; f(nx+1:2*nx) = dv

    fevals = fevals + 1
    return
  end subroutine adjoint_rhs

end module Ginzburg_Landau

! (a + ib) * (x + iy) = ax + i*ay + i*bx - by = (ax - by) + i * (bx + ay)
! (a - ib) * (x + iy) = ax + i*ay - i*bx + by = (ax + by) + i * (ay - bx)
