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

  ! !> Initialize physical parameters.
  ! call initialize_parameters() ; call extract_matrix()

  ! !> Eigenspectrum of the direct and adjoint operator + leading eigenvectors.
  ! call linear_stability_analysis()

  ! !> Optimal perturbation analysis.
  ! call transient_growth(linspace(0.5_wp, 50.0_wp, 100))

  ! !> Resolvent analysis.
  ! call resolvent_analysis(linspace(-5.0_wp, 5.0_wp, 201))

  call resolvent_computational_statistics()

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
    call eigs(A, X, v, lambda, residuals, info, nev=128, transpose=.false.)

    ! --> Transform eigenspectrum.
    lambda = log(lambda) / A%t

    ! --> Save the eigenspectrum.
    call save_eigenspectrum(real(lambda), aimag(lambda), residuals, "Direct_Eigenspectrum.npy")
    ! --> Save leading eigenvector.
    call get_vec(wrk, X(1:kdim), real(v(:, 1)))
    select type(wrk)
    type is(vector)
       call save_npy("Leading_direct_eigenvector.npy", wrk%x)
    end select

    ! --> Initialize Krylov subspace.
    call random_number(X(1)%x)
    alpha = X(1)%norm() ; call X(1)%scal(1.0_wp / alpha)

    ! --> Eigenvalue analysis.
    call eigs(A, X, v, lambda, residuals, info, nev=128, transpose=.true.)

    ! --> Transform eigenspectrum.
    lambda = log(lambda) / A%t

    ! --> Save the eigenspectrum.
    call save_eigenspectrum(real(lambda), aimag(lambda), residuals, "Adjoint_Eigenspectrum.npy")
    ! --> Save leading eigenvector.
    call get_vec(wrk, X(1:kdim), real(v(:, 1)))
    select type(wrk)
    type is(vector)
       call save_npy("Leading_adjoint_eigenvector.npy", wrk%x)
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

    write(output_unit, *) "--: Circular frequency : ", omega


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

    write(output_unit, *) "Elapsed time :", elapsed_time

    return
  end subroutine resolvent_single_freq

  subroutine resolvent_computational_statistics()
    !> Range of frequencies.
    real(kind=wp), allocatable :: omegas(:)
    integer                    :: nfreqs = 51
    !> Timings.
    real(kind=wp), allocatable :: timings(:)
    !> Functions evaluations.
    integer, allocatable       :: func_evals(:)
    !> Miscellaneous.
    integer :: i, j, k
    character(len=100) :: filename

    !----------------------------------
    !-----     INITIALIZATION     -----
    !----------------------------------

    omegas = linspace(-5.0_wp, 5.0_wp, nfreqs)
    allocate(timings(nfreqs)) ; allocate(func_evals(nfreqs))
    timings = 0.0_wp ; func_evals = 0

    do i = 1, nfreqs
       call resolvent_single_freq(omegas(i), timings(i))
       func_evals(i) = fevals ; fevals = 0
    enddo

    write(filename, '(a, i0.4, a)') 'Elapsed_time_nx=', nx, '.npy'
    call save_npy(filename, timings)
    write(filename, '(a, i0.4, a)') 'Func_evals_nx=', nx, '.npy'
    call save_npy(filename, func_evals)

    return
  end subroutine resolvent_computational_statistics

end program main
