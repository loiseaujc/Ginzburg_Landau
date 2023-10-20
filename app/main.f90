program main
  use iso_fortran_env
  !> Fortran Standard Library.
  use stdlib_io_npy  , only : save_npy
  use stdlib_math    , only : linspace, logspace
  !> Runge-Kutta time integrators library.
  use rklib_module
  !> Mesh/Physical parameters of the simulation.
  use params
  !> Ginzburg-Landau equations.
  use Ginzburg_Landau
  !>
  use LinearOperators
  !>
  use stdlib_io_npy, only : save_npy
  implicit none

  !> Range of frequencies investigated.
  real(kind=wp), allocatable :: omega(:)
  !> Elapsed time.
  real(kind=wp), allocatable :: elapsed_time(:)
  integer      , allocatable :: func_evals(:)
  !> Miscellaneous.
  integer :: i, j, k

  !> Initialize physical parameters.
  call initialize_parameters()

  !>
  ! call linear_stability_analysis()
  ! call transient_growth(linspace(1.0_wp, 60.0_wp, 60))
  ! call resolvent_analysis(linspace(-5.0_wp, 5.0_wp, 101))

  !> Performance timings.
  omega = linspace(-5.0_wp, 5.0_wp, 101)
  allocate(elapsed_time(1:size(omega))) ; elapsed_time = 0.0_wp
  allocate(func_evals(1:size(omega)))   ; func_evals   = 0

  do i = 1, size(omega)
     write(output_unit, *) "--: Circular frequency             :", omega(i)
     call resolvent_single_freq(omega(i), elapsed_time(i))
     func_evals(i) = fevals ; fevals = 0
     write(output_unit, *) "    Elapsed time                   :", elapsed_time(i), "seconds."
     write(output_unit, *) "    Number of function evaluations :", func_evals(i)
     write(output_unit, *)
  enddo

  call save_npy("Elapsed_time_nx128_rkt54.npy", elapsed_time)
  call save_npy("Function_evaluations_nx128_rkt54.npy", func_evals)

contains

  subroutine linear_stability_analysis()
    !> Exponential Propagator.
    class(exponential_prop), allocatable :: A
    !> Krylov subspace.
    integer, parameter :: kdim = 2*nx
    class(vector), allocatable :: X(:)
    !> Coordinates of the eigenvectors in the Krylov basis.
    complex(kind=wp) :: v(kdim, kdim)
    !> Eigenvalues.
    complex(kind=wp) :: lambda(kdim)
    !> Residual.
    real(kind=wp)    :: residuals(kdim)
    !> Information flag.
    integer          :: info

    !> Miscellaneous.
    integer :: i, j, k
    real(kind=wp) :: alpha
    class(abstract_vector), allocatable :: wrk
    complex(kind=wp) :: eigenvector(2*nx)

    !>
    A = exponential_prop(1.0_wp)

    ! --> Initialize Krylov subspace.
    allocate(X(1:kdim+1)) ; call random_number(X(1)%x)
    alpha = X(1)%norm() ; call X(1)%scal(1.0_wp / alpha)

    ! --> Eigenvalue analysis.
    call eigs(A, X, v, lambda, residuals, info, nev=128, transpose=.true.)

    ! --> Transform eigenspectrum.
    lambda = log(lambda) / A%t

    ! --> Save the eigenspectrum.
    open(unit=1234, file="adjoint_eigenspectrum.dat")
    do i = 1, size(lambda)
       write(1234, *) real(lambda(i)), aimag(lambda(i)), residuals(i)
    enddo
    close(1234)
    ! call save_eigenspectrum(real(lambda), aimag(lambda), residuals, "Direct_Eigenspectrum.npy")
    ! --> Save leading eigenvector.
    call get_vec(wrk, X(1:kdim), real(v(:, 1)))
    select type(wrk)
    type is(vector)
       call save_npy("Leading_eigenvector.npy", wrk%x)
    end select

    return
  end subroutine linear_stability_analysis

  subroutine transient_growth(times)
    !> Time instants at which to compute the optimal growth.
    real(kind=wp), intent(in) :: times(:)
    !> Exponential propagator.
    type(exponential_prop), allocatable :: A
    !> Krylov subspaces.
    integer, parameter :: kdim = 2*nx
    type(vector), allocatable :: U(:), V(:)
    !> Singular triplets.
    real(kind=wp) :: sigma(kdim), uvecs(kdim, kdim), vvecs(kdim, kdim)
    real(kind=wp) :: residuals(kdim)
    !> Information flag.
    integer :: info
    !> Gains.
    real(kind=wp) :: gains(size(times), 5)
    !> Miscellaneous.
    integer :: i, j, k
    real(kind=wp) :: alpha

    ! --> Initialize variables.
    allocate(U(kdim+1)) ; allocate(V(kdim+1))
    A = exponential_prop(0.0_wp)
    gains = 0.0_wp ; sigma = 0.0_wp ; uvecs = 0.0_wp ; vvecs = 0.0_wp ; residuals = 0.0_wp

    do i = 1, size(times)

       write(*, *) "--: Integration time : ", times(i)

       !> Set integration time for the exponential propagator.
       A%t = times(i)

       !> Initialize Krylov subspaces.
       do j = 1, size(U)
          call U(j)%zero() ; call V(j)%zero()
       enddo
       call random_number(U(1)%x) ; alpha = U(1)%norm() ; call U(1)%scal(1.0_wp / alpha)

       !> Singular value computation.
       call svds(A, U, V, uvecs, vvecs, sigma, residuals, info, nev=10, tolerance=rtol)

       !> Store computed gains.
       do k = 1, size(gains, 2)
          gains(i, k) = sigma(2*k-1)**2
       enddo
       write(*, *) "    Optimal gains :", gains(i, :), sigma(1:2*size(gains, 2))**2
       write(*, *)
    enddo

    call save_npy("Optimal_gains.npy", gains)

    return
  end subroutine transient_growth

  subroutine resolvent_analysis(omegas)
    !> Time instants at which to compute the optimal growth.
    real(kind=wp), intent(in) :: omegas(:)
    !> Exponential propagator.
    type(resolvent_op), allocatable :: R
    !> Krylov subspaces.
    integer, parameter :: kdim = 2*nx
    type(vector), allocatable :: U(:), V(:)
    !> Singular triplets.
    real(kind=wp) :: sigma(kdim), uvecs(kdim, kdim), vvecs(kdim, kdim)
    real(kind=wp) :: residuals(kdim)
    !> Information flag.
    integer :: info
    !> Gains.
    real(kind=wp) :: gains(size(omegas), 5)
    !> Miscellaneous.
    integer :: i, j, k
    real(kind=wp) :: alpha

    ! --> Initialize variables.
    allocate(U(kdim+1)) ; allocate(V(kdim+1))
    R = resolvent_op(0.0_wp)
    gains = 0.0_wp ; sigma = 0.0_wp ; uvecs = 0.0_wp ; vvecs = 0.0_wp ; residuals = 0.0_wp

    do i = 1, size(omegas)

       write(*, *) "--: Circular frequency : ", omegas(i)

       !> Set the forcing frequency for the Resolvent.
       R%omega = omegas(i)

       !> Initialize Krylov subspaces.
       do j = 1, size(U)
          call U(j)%zero() ; call V(j)%zero()
       enddo
       call random_number(U(1)%x) ; alpha = U(1)%norm() ; call U(1)%scal(1.0_wp / alpha)

       !> Singular value computation.
       call svds(R, U, V, uvecs, vvecs, sigma, residuals, info, nev=10, tolerance=rtol)

       !> Store computed gains.
       do k = 1, size(gains, 2)
          gains(i, k) = sigma(2*k-1)
       enddo
       write(*, *) "    Optimal gains :", gains(i, :)
       write(*, *)
    enddo

    call save_npy("Resolvent_gains.npy", gains)

    return
  end subroutine resolvent_analysis

  subroutine resolvent_single_freq(omega, elapsed_time)
    !> Time instants at which to compute the optimal growth.
    real(kind=wp), intent(in) :: omega
    !> Elapsed time.
    real(kind=wp), intent(out) :: elapsed_time
    real(kind=wp)              :: start_time = 0.0_wp, end_time = 0.0_wp
    !> Exponential propagator.
    type(resolvent_op), allocatable :: R
    !> Krylov subspaces.
    integer, parameter :: kdim = 2*nx
    type(vector), allocatable :: U(:), V(:)
    !> Singular triplets.
    real(kind=wp) :: sigma(kdim), uvecs(kdim, kdim), vvecs(kdim, kdim)
    real(kind=wp) :: residuals(kdim)
    !> Information flag.
    integer :: info
    !> Miscellaneous.
    integer :: i, j, k
    real(kind=wp) :: alpha

    ! --> Initialize variables.
    allocate(U(1:kdim+1)) ; allocate(V(1:kdim+1))
    R = resolvent_op(omega)
    sigma = 0.0_wp ; uvecs = 0.0_wp ; vvecs = 0.0_wp ; residuals = 0.0_wp

    !> Initialize Krylov subspaces.
    do j = 1, size(U)
       call U(j)%zero() ; call V(j)%zero()
    enddo
    call random_number(U(1)%x) ; alpha = U(1)%norm() ; call U(1)%scal(1.0_wp / alpha)

    !> Singular value computation.
    call cpu_time(start_time)
    call svds(R, U, V, uvecs, vvecs, sigma, residuals, info, nev=10, tolerance=rtol)
    call cpu_time(end_time)

    !> Elapsed time.
    elapsed_time = end_time - start_time

    return
  end subroutine resolvent_single_freq

end program main
