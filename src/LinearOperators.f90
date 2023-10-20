module LinearOperators
  use iso_fortran_env
  !> Definition of the direct and adjoint linearized GL equations.
  use Ginzburg_Landau
  !> Krylov methods.
  use LightKrylov
  !> Physical parameters.
  use params
  !> Time integrators.
  use rklib_module
  implicit none

  private

  !-------------------------------------------
  !-----     LIGHTKRYLOV VECTOR TYPE     -----
  !-------------------------------------------

  type, extends(abstract_vector), public :: vector
     real(kind=wp) :: x(2*nx) = 0.0_wp
   contains
     private
     procedure, pass(self), public :: zero
     procedure, pass(self), public :: dot
     procedure, pass(self), public :: scal
     procedure, pass(self), public :: axpby
  end type vector

  !------------------------------------------
  !-----     EXPONENTIAL PROPAGATOR     -----
  !------------------------------------------

  type, extends(abstract_linop), public :: exponential_prop
     real(kind=wp) :: t
   contains
     private
     procedure, pass(self), public :: matvec  => direct_solver
     procedure, pass(self), public :: rmatvec => adjoint_solver
  end type exponential_prop

  !--------------------------------------
  !-----     RESOLVENT OPERATOR     -----
  !--------------------------------------

  type, extends(abstract_linop), public :: resolvent_op
     real(kind=wp) :: omega
   contains
     private
     procedure, pass(self), public :: matvec  => direct_map
     procedure, pass(self), public :: rmatvec => adjoint_map
  end type resolvent_op

contains

  !----------------------------------------------------
  !-----     TYPE-BOUND PROCEDURE FOR VECTORS     -----
  !----------------------------------------------------

  subroutine zero(self)
    class(vector), intent(inout) :: self
    self%x = 0.0_wp
    return
  end subroutine zero
  
  real(kind=wp) function dot(self, vec) result(alpha)
    class(vector), intent(in)          :: self
    class(abstract_vector), intent(in) :: vec

    select type(vec)
    type is(vector)
       alpha = dot_product(self%x, vec%x)
    end select
    return
  end function dot

  subroutine scal(self, alpha)
    class(vector), intent(inout) :: self
    real(kind=wp), intent(in)    :: alpha
    self%x = self%x * alpha
    return
  end subroutine scal

  subroutine axpby(self, alpha, vec, beta)
    class(vector)         , intent(inout) :: self
    class(abstract_vector), intent(in)    :: vec
    real(kind=wp)         , intent(in)    :: alpha, beta

    select type(vec)
    type is(vector)
       self%x = alpha*self%x + beta*vec%x
    end select
    return
  end subroutine axpby

  !--------------------------------------------------------------
  !-----     TYPE-BOUND PROCEDURE FOR EXPONENTIAL PROP.     -----
  !--------------------------------------------------------------

  subroutine direct_solver(self, vec_in, vec_out)
    class(exponential_prop), intent(in)  :: self
    class(abstract_vector) , intent(in)  :: vec_in
    class(abstract_vector) , intent(out) :: vec_out

    !> Time-integrator.
    type(rkt54_class) :: prop
    !> Tentative time-step.
    real(kind=wp)     :: dt = 1.0_wp

    select type(vec_in)
    type is(vector)
       select type(vec_out)
       type is(vector)

          !> Initialize propagator.
          call prop%initialize(n=2*nx, f=rhs, rtol=[rtol], atol=[atol])
          !> Integrate ODE.
          call prop%integrate(0.0_wp, vec_in%x, dt, self%t, vec_out%x)

       end select
    end select
    return
  end subroutine direct_solver

  subroutine adjoint_solver(self, vec_in, vec_out)
    class(exponential_prop), intent(in)  :: self
    class(abstract_vector) , intent(in)  :: vec_in
    class(abstract_vector) , intent(out) :: vec_out

    !> Time-integrator.
    type(rkt54_class) :: prop
    !> Tentative time-step.
    real(kind=wp)     :: dt = 1.0_wp

    select type(vec_in)
    type is(vector)
       select type(vec_out)
       type is(vector)

          !> Initialize propagator.
          call prop%initialize(n=2*nx, f=adjoint_rhs, rtol=[rtol], atol=[atol])
          !> Integrate ODE.
          call prop%integrate(0.0_wp, vec_in%x, dt, self%t, vec_out%x)

       end select
    end select
    return
  end subroutine adjoint_solver

  !-------------------------------------------------------------------
  !-----     TYPE-BOUND PROCEDURE FOR THE RESOLVENT OPERATOR     -----
  !-------------------------------------------------------------------

  subroutine direct_map(self, vec_in, vec_out)
    class(resolvent_op)   , intent(in)  :: self
    class(abstract_vector), intent(in)  :: vec_in
    class(abstract_vector), intent(out) :: vec_out

    !> Exponential propagator.
    class(exponential_prop), allocatable :: A
    class(identity_linop)  , allocatable :: Id

    !> Internal variables.
    real(kind=wp)             :: x0(2*nx) = 0.0_wp, xf(2*nx) = 0.0_wp
    real(kind=wp)             :: forcing(2*nx)
    type(rkt54_class)         :: prop
    real(kind=wp)             :: dt = 1.0_wp
    character(len=:), allocatable :: message

    class(axpby_linop), allocatable :: S
    class(vector), allocatable :: b
    type(gmres_opts)          :: opts
    integer                   :: info

    !> Miscellaneous.
    real(kind=wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    real(kind=wp)            :: tau

    !> Integration time.
    if (self%omega == 0) then
       tau = 100.0_wp
    else
       tau = 2.0_wp * pi / abs(self%omega)
    endif
       
    !> Initialize needed operators.
    Id = identity_linop() ; A = exponential_prop(tau)

    !> Compute int_{0}^t exp( (t-s)*A ) * f(s) ds
    select type(vec_in)
    type is(vector)
       allocate(b, source=vec_in) ; call b%zero()
       !> Sets the forcing spatial support.
       forcing(:) = vec_in%x(:)
    end select
    
    !> Compute the right-hand side vector for gmres.
    call prop%initialize(n=2*nx, f=forced_rhs, rtol=[rtol], atol=[atol])
    call prop%integrate(0.0_wp, x0, dt, tau, b%x)
    ! call prop%status(message=message)
    ! write(output_unit, '(A)') message

    !> GMRES solve to compute the post-transient response.
    opts = gmres_opts(verbose=.false., atol=atol, rtol=rtol)
    S = axpby_linop(Id, A, 1.0_wp, -1.0_wp)
    call gmres(S, b, vec_out, info, options=opts)

    return

  contains

    subroutine forced_rhs(me, t, x, f)
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

      !> Add the forcing term.
      du(:) = du(:) + forcing(1:nx)*cos(self%omega*t) - forcing(nx+1:2*nx)*sin(self%omega*t)
      dv(:) = dv(:) + forcing(1:nx)*sin(self%omega*t) + forcing(nx+1:2*nx)*cos(self%omega*t)
      
      !> Copy results to the output array.
      f(1:nx) = du ; f(nx+1:2*nx) = dv

      return
    end subroutine forced_rhs
    
  end subroutine direct_map

  subroutine adjoint_map(self, vec_in, vec_out)
    class(resolvent_op)   , intent(in)  :: self
    class(abstract_vector), intent(in)  :: vec_in
    class(abstract_vector), intent(out) :: vec_out

    !> Exponential propagator.
    class(exponential_prop), allocatable :: A
    class(identity_linop)  , allocatable :: Id

    !> Internal variables.
    real(kind=wp)             :: x0(2*nx) = 0.0_wp, xf(2*nx) = 0.0_wp
    real(kind=wp)             :: forcing(2*nx)
    type(rkt54_class)         :: prop
    real(kind=wp)             :: dt = 1.0_wp
    character(len=:), allocatable :: message

    class(axpby_linop), allocatable :: S
    class(vector), allocatable :: b
    type(gmres_opts)          :: opts
    integer                   :: info

    !> Miscellaneous.
    real(kind=wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    real(kind=wp)            :: tau

    !> Integration time.
    if (self%omega == 0) then
       tau = 100.0_wp
    else
       tau = 2.0_wp * pi / abs(self%omega)
    endif

    !> Initialize needed operators.
    Id = identity_linop() ; A = exponential_prop(tau)

    !> Compute int_{0}^t exp( (t-s)*A ) * f(s) ds
    select type(vec_in)
    type is(vector)
       allocate(b, source=vec_in) ; call b%zero()
       !> Sets the forcing spatial support.
       forcing(:) = vec_in%x(:)
    end select

    !> Compute the right-hand side vector for gmres.
    call prop%initialize(n=2*nx, f=forced_adjoint_rhs, rtol=[rtol], atol=[atol])
    call prop%integrate(0.0_wp, x0, dt, tau, b%x)
    ! call prop%status(message=message)
    ! write(output_unit, '(A)') message

    !> GMRES solve to compute the post-transient response.
    opts = gmres_opts(verbose=.false., atol=atol, rtol=rtol)
    S = axpby_linop(Id, A, 1.0_wp, -1.0_wp, .false., .true.)
    call gmres(S, b, vec_out, info, options=opts, transpose=.false.)

    return

  contains

    subroutine forced_adjoint_rhs(me, t, x, f)
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

      !> Add the forcing term.
      du(:) = du(:) + forcing(1:nx)*cos(self%omega*t) + forcing(nx+1:2*nx)*sin(self%omega*t)
      dv(:) = dv(:) - forcing(1:nx)*sin(self%omega*t) + forcing(nx+1:2*nx)*cos(self%omega*t)
     
      !> Copy results to the output array.
      f(1:nx) = du ; f(nx+1:2*nx) = dv
      
      return
    end subroutine forced_adjoint_rhs
    
  end subroutine adjoint_map

end module LinearOperators
