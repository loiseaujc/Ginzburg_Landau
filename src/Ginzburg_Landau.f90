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

  ! subroutine rhs(du, u, p)
  !   !> Time-derivative of the state vector.
  !   complex(kind=wp), intent(out)   :: du(:)
  !   !> State vector of the system.
  !   complex(kind=wp), intent(in)    :: u(:)
  !   !> Simulation parameters.
  !   type(parameters), intent(in)    :: p

  !   !> Miscellaneous.
  !   integer :: i, j, k

  !   ! --> Initialize time-derivative.
  !   du = (0.0_wp, 0.0_wp)

  !   ! --> Left most boundary point.
  !   du(1) = -p%nu * u(2) / (2.0_wp * p%dx)
  !   du(1) = du(1) + p%gamma*(u(2) - 2.0_wp*u(1)) / p%dx**2
  !   du(1) = du(1) + p%mu(1) * u(1)

  !   ! --> Interior nodes.
  !   do i = 2, p%nx-1
  !      ! --> Convective term.
  !      du(i) = -p%nu * (u(i+1) - u(i-1)) / (2.0_wp * p%dx)
  !      ! --> Diffusion.
  !      du(i) = du(i) + p%gamma * (u(i+1) - 2.0_wp*u(i) + u(i-1)) / p%dx**2
  !      ! --> Driving term.
  !      du(i) = du(i) + p%mu(i) * u(i)
  !   enddo

  !   ! --> Right most boundary point.
  !   du(p%nx) = p%nu * u(p%nx-1) / (2.0_wp * p%dx)
  !   du(p%nx) = du(p%nx) + p%gamma * (-2.0_wp * u(p%nx) - u(p%nx-1)) / p%dx**2
  !   du(p%nx) = du(p%nx) + p%mu(p%nx) * u(p%nx)

  !   return
  ! end subroutine rhs

  ! subroutine forcing(du, f, t, p)
  !   !> Time-derivative of the state vector.
  !   complex(kind=wp), intent(inout) :: du(:)
  !   !> State vector of the forcing term.
  !   complex(kind=wp), intent(in)    :: f(:)
  !   !> Current simulation time.
  !   real(kind=wp), intent(in)       :: t
  !   !> Simulation parameters.
  !   type(parameters), intent(in)    :: p
  !   !> Miscellaneous.
  !   integer :: i, j, k

  !   do i = 1, p%nx-1
  !      du(i) = du(i) + f(i) * cmplx(cos(p%omega*t), sin(p%omega*t))
  !   enddo

  !   return
  ! end subroutine forcing

  ! subroutine adjoint_forcing(du, f, t, p)
  !   !> Time-derivative of the state vector.
  !   complex(kind=wp), intent(inout) :: du(:)
  !   !> State vector of the forcing term.
  !   complex(kind=wp), intent(in)    :: f(:)
  !   !> Current simulation time.
  !   real(kind=wp), intent(in)       :: t
  !   !> Simulation parameters.
  !   type(parameters), intent(in)    :: p
  !   !> Miscellaneous.
  !   integer :: i, j, k

  !   do i = 1, p%nx-1
  !      du(i) = du(i) + f(i) * cmplx(cos(p%omega*t), -sin(p%omega*t))
  !   enddo

  !   return
  ! end subroutine adjoint_forcing

  ! subroutine forced_step(u, f, t, p)
  !   !> State vector of the system.
  !   complex(kind=wp), intent(inout) :: u(:)
  !   !> Forcing term.
  !   complex(kind=wp), intent(in)    :: f(:)
  !   !> Current simulation time.
  !   real(kind=wp)   , intent(in)    :: t
  !   !> Simulation parameters.
  !   type(parameters), intent(in)    :: p

  !   !> RK substeps.
  !   complex(kind=wp) :: k(size(u), 4)
  !   real(kind=wp)    :: coefs(4) = [1.0_wp, 2.0_wp, 2.0_wp, 1.0_wp] / 6.0_wp

  !   !> Runge-Kutta integration.
  !   call rhs(k(:, 1), u, p) ; call forcing(k(:, 1), f, t, p)
  !   call rhs(k(:, 2), u + p%dt/2.0_wp * k(:, 1), p) ; call forcing(k(:, 2), f, t+p%dt/2.0_wp, p)
  !   call rhs(k(:, 3), u + p%dt/2.0_wp * k(:, 2), p) ; call forcing(k(:, 3), f, t+p%dt/2.0_wp, p)
  !   call rhs(k(:, 4), u + p%dt * k(:, 3), p) ; call forcing(k(:, 4), f, t+p%dt, p)

  !   !> New state.
  !   u(:) = u(:) + p%dt * matmul(k, coefs)

  !   return
  ! end subroutine forced_step

  ! !----------------------------------------------------------
  ! !-----                                                -----
  ! !-----     DEFINITION OF THE ADJOINT TIME-STEPPER     -----
  ! !-----                                                -----
  ! !----------------------------------------------------------

  ! subroutine adjoint_rhs(du, u, p)
  !   !> Time-derivative of the state vector.
  !   complex(kind=wp), intent(out) :: du(:)
  !   !> State vector of the system.
  !   complex(kind=wp), intent(in)  :: u(:)
  !   !> Simulation parameters.
  !   type(parameters), intent(in)  :: p

  !   !> Miscellaneous.
  !   integer :: i, j, k

  !   ! --> Initialize time-derivative.
  !   du = (0.0_wp, 0.0_wp)

  !   ! --> Left-most boundary point.
  !   du(1) = conjg(p%nu) * u(2) / (2.0_wp*p%dx)
  !   du(1) = du(1) + conjg(p%gamma) * (u(2) - 2.0_wp *u(1)) / p%dx**2
  !   du(1) = du(1) + p%mu(1) * u(1)

  !   ! --> Inner points.
  !   do i = 2, p%nx-1
  !      ! --> Convective term.
  !      du(i) = conjg(p%nu) * (u(i+1) - u(i-1)) / (2.0_wp * p%dx)
  !      ! --> Diffusion.
  !      du(i) = du(i) + conjg(p%gamma) * (u(i+1) - 2.0_wp * u(i) + u(i-1)) / p%dx**2
  !      ! --> Driving term.
  !      du(i) = du(i) + p%mu(i) * u(i)
  !   enddo

  !   ! --> Right-most boundary point.
  !   du(p%nx) = conjg(p%nu) * (-u(p%nx-1)) / (2.0_wp * p%dx)
  !   du(p%nx) = du(p%nx) + conjg(p%gamma) * (-2.0_wp * u(p%nx) + u(p%nx-1)) / p%dx**2
  !   du(p%nx) = du(p%nx) + p%mu(p%nx) * u(p%nx)

  !   return
  ! end subroutine adjoint_rhs

end module Ginzburg_Landau

! (a + ib) * (x + iy) = ax + i*ay + i*bx - by = (ax - by) + i * (bx + ay)
! (a - ib) * (x + iy) = ax + i*ay - i*bx + by = (ax + by) + i * (ay - bx)
