module neoclassical_nlopt_ver2

    use iso_Fortran_env, only: rk => real64, ik => int64
    use numerics
    use io
    use ogpf
    use nlopt_wrap
    use nlopt_enum
    implicit none

    ! ---------- !
    ! parameters !
    ! ---------- !
    real(rk), parameter :: zss = 1
    real(rk), parameter :: tol = 1.0D-6           ! D is e in matlab, double precision
    real(rk), parameter :: zeta = 1.38_rk
    real(rk), parameter :: beta = 0.96_rk         ! yearly interest rate 4%
    real(rk), parameter :: delta = 0.065_rk       ! I/K ~ 10%
    real(rk), parameter :: psi = 2.15_rk          ! hours of worked = 1/3
    real(rk), parameter :: alpha = 0.27_rk        ! K/Y ~ 2.5
    real(rk), parameter :: nu = 0.6_rk            ! labor share: 60% in US
    real(rk), parameter :: rho_z = 0.909_rk       ! Autocor of AR(1) agg. shock
    real(rk), parameter :: sigma_z = 0.014_rk     ! std of AR(1) agg. shock
    real(rk), parameter :: gamma = 0.0260_rk      ! capital reallocation cost
    real(rk), parameter :: eta = 0.75_rk          ! I-tech: new v.s. used ratio
    real(rk), parameter :: s = 5.0_rk             ! I-tech: CES coefficient
    real(rk), parameter :: rho_e = 0.687_rk       ! Autocor of AR(1) idio. shock
    real(rk), parameter :: sigma_e = 0.117_rk     ! std of AR(1) idio. shock
    integer(ik), parameter :: knum = 10_ik
    integer(ik), parameter :: enum = 7_ik

    ! ----- !
    ! grids !
    ! ----- !
    real(rk), dimension(2) :: kbounds = (/0.05_rk, 6.0_rk/)
    real(rk), dimension(knum) :: kgrid
    real(rk), dimension(enum) :: egrid
    real(rk), dimension(enum) :: piess
    real(rk), dimension(enum, enum) :: pie

    ! --------- !
    ! solutions !
    ! --------- !
    real(rk), dimension(:, :), allocatable :: v, gk

    ! ------------------- !
    ! func_data for nlopt !
    ! ------------------- !
    type value_data
        real(rk) :: evvec(knum)
    end type


contains

    subroutine valueiter()

        integer(ik) :: indk, inde, indef
        integer(ik) :: iter, maxiter, iunit
        real(rk), dimension(knum, enum) :: tv, tgk, ev
        real(rk), dimension(knum) :: evvec
        real(rk) :: dist, vTol, distv, distgk
        real(rk) :: epsval, kval, tempval, nval, yval

        iter = 0_ik
        maxiter = 1000_ik

        tv = 0.0_rk
        tgk = 0.0_rk
        ev = 0.0_rk
        evvec = 0.0_rk

        write(*, *) 'solve for value function use pp-linear'

        vTol = tol*10.0_rk
        dist = 2.0*vTol

        do while (dist > vTol .AND. iter <= maxiter)

            ! -------------- !
            ! DEBUG
            ! -------------- !
            iter = maxiter !
            ! -------------- !

            iter = iter + 1

            ! --------------------------------- !
            ! calculate conditional expectation !
            ! --------------------------------- !

            do indk = 1, knum, 1
                do inde = 1, enum, 1
                    tempval = 0.0_rk
                    do indef = 1, enum, 1
                        tempval = tempval + pie(inde, indef)*v(indk, indef)
                    enddo
                    ev(indk, inde) = tempval
                enddo
            enddo

            do inde = 1, enum, 1
                epsval = egrid(inde)
                evvec = ev(:, inde)
                call kupTarget(evvec, epsval)
            enddo




        enddo

    end subroutine valueiter

    subroutine kupTarget(evvec, epsval)
        real(rk) :: epsval, fmax
        real(rk), dimension(:) :: evvec
        real(rk), dimension(knum) :: guess
        integer :: stat
        intent(in) :: evvec, epsval
        type(nlopt_opt) :: opt

        guess = kgrid(knum/2)

        ! write(*, *) knum/2

        ! call pprint(guess)

        call create(opt, algorithm_from_string("LN_SBPLX"), int(knum, 4))
        call opt%set_lower_bounds(kgrid(1))
        call opt%set_upper_bounds(kgrid(knum))
        call opt%set_xtol_rel(tol)
        call opt%set_max_objective(nlopt_func(kupObj, evvec))
        call opt%optimize(guess, fmax, stat)
        call destroy(opt)

    end subroutine kuptarget

    function kupObj(guess, gradient, func_data) result(f)
        real(rk), dimension(:), intent(in) :: guess
        real(rk), dimension(:), intent(inout), optional :: gradient
        class(*), intent(in), optional :: func_data
        type(value_data) :: vdata

        integer(ik) :: kidx
        real(rk) :: f, kw, evval
        real(rk), dimension(knum) :: evvec

        select type(func_data)
            type is (value_data)
            vdata = func_data
        end select

        IF (PRESENT(gradient)) THEN
            ! Why?!
        END IF

        call linear_interpolation(kgrid, knum, guess(1), kidx, kw)
        evval = vdata%evvec(kidx)*kw + vdata%evvec(kidx+1_ik)*(1.0_rk - kw)

        f = -1.0_rk*guess(1) + beta*evval

    end function kupObj

    subroutine initGrids()
        integer(ik) :: i
        real(rk) :: tempterm(enum, enum)

        call tauchen(0.0_rk, sigma_e, rho_e, 2.575_rk, enum, egrid, pie)
        egrid = exp(egrid)
        kgrid = linspace(kbounds(1), kbounds(2), knum)

        tempterm = eye(enum)
        do i = 1, 1000, 1
            tempterm = matmul(tempterm, pie)
        enddo
        piess = tempterm(1, :)

    end subroutine initGrids

    subroutine initSol()

        if (allocated(v)) deallocate(v)
        if (allocated(gk)) deallocate(gk)

        allocate(v(knum, enum), source=0.0_rk)
        allocate(gk(knum, enum), source=0.0_rk)

    end subroutine initSol


end module neoclassical_nlopt_ver2
