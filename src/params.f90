module params
  use LightKrylov
  use stdlib_math , only : linspace
  implicit none

  !------------------------------
  !-----     PARAMETERS     -----
  !------------------------------

  ! --> Mesh related parameters.
  real(kind=wp), parameter :: L = 200.0_wp !> Domain length
  integer      , parameter :: nx = 512     !> Number of grid points (excluding boundaries).
  real(kind=wp), parameter :: dx = L / nx  !> Grid size.

  ! --> Physical parameters.
  complex(kind=wp), parameter :: nu    = cmplx(2.0_wp, 0.2_wp, kind=wp)
  complex(kind=wp), parameter :: gamma = cmplx(1.0_wp, -1.0_wp, kind=wp)
  real(kind=wp)   , parameter :: mu_0  = 0.38_wp
  real(kind=wp)   , parameter :: c_mu  = 0.2_wp
  real(kind=wp)   , parameter :: mu_2  = -0.01_wp
  real(kind=wp)               :: mu(1:nx)

contains

  subroutine initialize_parameters()
    implicit none
    !> Mesh array.
    real(kind=wp), allocatable :: x(:)

    !> Construct mesh.
    x = linspace(-L/2, L/2, nx+2)

    !> Construct mu(x)
    mu(:) = (mu_0 - c_mu**2) + (mu_2 / 2.0_wp) * x(2:nx+1)**2

    return
  end subroutine initialize_parameters

end module params
