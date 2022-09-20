module ReplicateLanteri2018

    use iso_Fortran_env, only: rk => real64, ik => int32
    use numerics
    use io
    use ogpf
    use nlopt_wrap
    use nlopt_enum
    implicit none

    ! ------ !
    ! output !
    ! ------ !
    character(len=*), parameter :: resDir = "./results/"
    character(len=*), parameter :: figDir = "./figures/"
    character(len=128) :: sep
    character(len=1), parameter :: tab = char(9)
    logical, parameter :: drawfig = .true.
    logical, parameter :: show_valueiter = .true.


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
    integer(ik), parameter :: knum = 1000_ik
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
    type solutions
        real(rk), dimension(:, :), allocatable :: v, gk
        real(rk), dimension(:, :), allocatable :: kupmat, kdnmat
    end type

    ! ------------------- !
    ! func_data for nlopt !
    ! ------------------- !
    type configurations
        real(rk) :: zval
        real(rk) :: wval
        real(rk) :: Qbuy
        real(rk) :: qsell
        real(rk) :: evvec(knum)
    end type


!     ! integer(ik) :: objcount = 0

!     ! ------------------- !
!     ! func_data for nlopt !
!     ! ------------------- !
!     type obj_data
!         real(rk) :: evvec(knum)
!         real(rk) :: Qbuy
!         real(rk) :: qsell
!     end type


contains

    subroutine valueiter(sol, conf)

        type(solutions), intent(inout) :: sol
        type(configurations), intent(inout) :: conf
        integer(ik) :: indk, inde, indef
        integer(ik) :: iter, maxiter
        real(rk), dimension(knum, enum) :: v, gk, tv, tgk, ev
        real(rk), dimension(knum) :: evvec
        real(rk) :: dist, vTol, distv, distgk
        real(rk) :: epsval, kval, tempval, nval, yval, zval, wval, vval, kstay
        real(rk) :: kupstar, evupstar, kdnstar, evdnstar, vstar, kstar
        ! type(obj_data) :: fdata

        iter = 0_ik
        maxiter = 1000_ik

        v = 0.0_rk
        gk = 0.0_rk
        tv = 0.0_rk
        tgk = 0.0_rk
        ev = 0.0_rk

        ! initialize fdata
        ! fdata%Qbuy = conf%Qbuy
        ! fdata%qsell = conf%qsell

        write(*, *) 'solve for value function use pp-linear'

        vTol = tol*10.0_rk
        dist = 2.0*vTol

        do while (dist > vTol .AND. iter <= maxiter)

            ! -------------- !
            ! DEBUG
            ! -------------- !
            ! iter = maxiter !
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
                conf%evvec = evvec
                call kupTarget(kupstar, evupstar, conf)
                call kdnTarget(kdnstar, evdnstar, conf)
                ! call kuptarget_orig(kupstar, evupstar, evvec, conf)
                ! call kdntarget_orig(kdnstar, evdnstar, evvec, conf)
                do indk = 1, knum, 1
                    kval = kgrid(indk)
                    kstay = (1.0_rk - delta)*kval
                    nval = ( ( nu * conf%zval * epsval * kval**alpha ) / conf%wval )**(1/(1-nu))
                    yval = conf%zval * epsval * kval**alpha * nval**nu
                    call krule(vstar, kstar, evvec, yval, nval, wval, kstay, kupstar, evupstar, kdnstar, evdnstar, conf)

                    ! vstar = yval + (1-delta)*conf%qsell*kval + evdnstar
                    ! kstar = kdnstar
                    ! vstar = yval + (1-delta)*conf%Qbuy*kval + evupstar
                    ! kstar = kupstar
                    tv(indk, inde) = vstar
                    tgk(indk, inde) = kstar
                    sol%kupmat(indk, inde) = kupstar
                    sol%kdnmat(indk, inde) = kdnstar
                    ! if (indk == 8_ik) then
                    !     write(*, *) evupstar, evdnstar, vstar, kstar
                    ! endif
                enddo
            enddo

            distv = maxval(dabs(v - tv))
            distgk = maxval(dabs(gk - tgk))
            dist = dmax1(distv, distgk)

            v = tv
            gk = tgk

            if (iter == 1 .and. show_valueiter) then
                ! print on terminal
                write(*, '(a4, 4(a20))') 'iter', 'distv', 'distgk', 'kmin', 'kmax'
                ! write the same content into the variable sep
                write(sep, '(a4, 4(a20))') 'iter', 'distv', 'distgk', 'kmin', 'kmax'
                ! repeatively print = with the length of sep
                write(*, '(a)') repeat('=', len_trim(sep))
            endif


            if (mod(iter, 25) == 0 .and. show_valueiter) then
            ! if (show_valueiter) then
                write(*, '(I4, 4(ES20.6))') iter, distv, distgk, minval(gk), maxval(gk)
            endif

        enddo

        ! call pprint(ev)

        sol%v = v
        sol%gk = gk


    end subroutine valueiter

    subroutine krule(vstar, kstar, evvec, yval, nval, wval, kstay, kupstar, evupstar, kdnstar, evdnstar, conf)
        real(rk), intent(in) :: yval, nval, wval, kstay, kupstar, evupstar, kdnstar, evdnstar
        real(rk), dimension(knum) :: evvec
        type(configurations), intent(inout) :: conf
        real(rk), intent(out) :: vstar, kstar
        real(rk) :: kup, evup, kdn, evdn, vup, vdn, evval, kw
        integer(ik) :: kidx

        if (kupstar >= kstay) then
            kup = kupstar; evup = evupstar
        else
            kup = kstay
            evup = evaluate(evvec, kup, conf%Qbuy)
        endif

        if (kdnstar <= kstay) then
            kdn = kdnstar; evdn = evdnstar
        else
            kdn = kstay
            evdn = evaluate(evvec, kdn, conf%qsell)
        endif

        vup = yval - conf%wval*nval + conf%Qbuy*kstay + evup
        vdn = yval - conf%wval*nval + conf%qsell*kstay + evdn

        if (vup .ge. vdn) then
            vstar = vup
            kstar = kup
        else
            vstar = vdn
            kstar = kdn
        endif

    end subroutine krule

    subroutine kuptarget_orig(kupstar, evupstar, evvec, conf)
        type(configurations), intent(in) :: conf
        real(rk), dimension(knum), intent(in) :: evvec
        real(rk), intent(out) :: kupstar, evupstar
        integer(ik) :: kidx, iter
        real(rk) :: a, b, c, d, z, fval, fc, fd, rg
        real(rk) :: Jval, kfval, evval, kw, dist, gssTol

        iter = 0

        Jval = conf%Qbuy

        rg = (3.0_rk - dsqrt(5.0_rk)) / 2.0_rk
        a = kgrid(1)
        b = kgrid(knum)
        c = a + rg*(b-a);
        d = a + (1-rg)*(b-a);

        ! initial c evaluation
        kfval = c
        fval = evaluate(evvec, kfval, Jval)
        fc = -fval

        ! initial d evaluation
        kfval = d
        fval = evaluate(evvec, kfval, Jval)
        fd = -fval

        gssTol = tol
        dist = 2.0_rk*gssTol

        do
            if (dist < gssTol) exit

            iter = iter + 1

            if (fc .ge. fd) then

                z = c + (1.0_rk-rg)*(b-c)
                ! case 1 [a c d b] <--- [c d z b]
                a = c
                c = d
                fc = fd
                d = z

                kfval = d
                fval = evaluate(evvec, kfval, Jval)
                fd = -fval
            else
                z = a + rg*(d-a)
                ! case 2 [a c d b] <--- [a z c d]
                b = d
                d = c
                fd = fc
                c = z

                kfval = c
                fval = evaluate(evvec, kfval, Jval)
                fc = -fval
            endif

            dist = b - a

            ! write (*, '((a, i3), 3(a, f9.6))') 'iter = ', iter, ', dist = ', dist, ', kfval = ', kfval, ', fval = ', fval

        enddo

        kupstar = kfval
        evupstar = fval

    end subroutine kuptarget_orig

    subroutine kdntarget_orig(kdnstar, evdnstar, evvec, conf)
        type(configurations), intent(in) :: conf
        real(rk), dimension(knum), intent(in) :: evvec
        real(rk), intent(out) :: kdnstar, evdnstar
        integer(ik) :: kidx, iter
        real(rk) :: a, b, c, d, z, fval, fc, fd, rg
        real(rk) :: Jval, kfval, evval, kw, dist, gssTol

        iter = 0

        Jval = conf%qsell

        rg = (3.0_rk - dsqrt(5.0_rk)) / 2.0_rk
        a = kgrid(1)
        b = kgrid(knum)
        c = a + rg*(b-a);
        d = a + (1-rg)*(b-a);

        ! initial c evaluation
        kfval = c
        fval = evaluate(evvec, kfval, Jval)
        fc = -fval

        ! initial d evaluation
        kfval = d
        fval = evaluate(evvec, kfval, Jval)
        fd = -fval

        gssTol = tol
        dist = 2.0_rk*gssTol

        do
            if (dist < gssTol) exit

            iter = iter + 1

            if (fc .ge. fd) then

                z = c + (1.0_rk-rg)*(b-c)
                ! case 1 [a c d b] <--- [c d z b]
                a = c
                c = d
                fc = fd
                d = z

                kfval = d
                fval = evaluate(evvec, kfval, Jval)
                fd = -fval
            else
                z = a + rg*(d-a)
                ! case 2 [a c d b] <--- [a z c d]
                b = d
                d = c
                fd = fc
                c = z

                kfval = c
                fval = evaluate(evvec, kfval, Jval)
                fc = -fval
            endif

            dist = b - a

            ! write (*, '((a, i3), 3(a, f9.6))') 'iter = ', iter, ', dist = ', dist, ', kfval = ', kfval, ', fval = ', fval

        enddo

        kdnstar = kfval
        if (kfval .ge. 0.0_rk) then
            evdnstar = fval
        else
            kfval = 0.0_rk
            fval = evaluate(evvec, kfval, Jval)
            kdnstar = kfval
            evdnstar = fval
        endif


    end subroutine kdntarget_orig

    subroutine kdnTarget(kdnstar, evdnstar, conf)
        type(configurations), intent(inout) :: conf
        real(rk), intent(out) :: kdnstar, evdnstar
        type(nlopt_opt) :: opt
        real(rk), dimension(1) :: guess
        real(rk) :: fmax
        integer :: stat

        guess = kgrid(knum/2)

        ! remember that the last number is in int4 but not 8
        call create(opt, algorithm_from_string("LN_SBPLX"), 1)
        call opt%set_lower_bounds(kgrid(1))
        call opt%set_upper_bounds(kgrid(knum))
        associate( f => nlopt_func(kdnObj, conf) )
            call opt%set_xtol_rel(tol)
            call opt%set_max_objective(f)
            call opt%optimize(guess, fmax, stat)
        end associate
        call destroy(opt)

        if (stat < NLOPT_SUCCESS) then
            write (*, "(A,I5)") "NLopt failed with code ", stat
        end if

        kdnstar = guess(1)
        evdnstar = fmax

    end subroutine kdnTarget

    function kdnObj(guess, gradient, func_data) result(f)
        real(rk), dimension(:), intent(in) :: guess
        real(rk), dimension(:), intent(inout), optional :: gradient
        class(*), intent(in), optional :: func_data
        type(configurations) :: conf
        integer(ik) :: kidx
        real(rk) :: f, kw, evval, kfval
        real(rk), dimension(knum) :: evvec

        select type(func_data)
            type is (configurations)
            conf = func_data
        end select

        IF (PRESENT(gradient)) THEN
            ! Why?!
        END IF

        kfval = guess(1)

        ! call linear_interpolation(kgrid, knum, kfval, kidx, kw)
        kidx = gridlookup(kgrid, knum, kfval)
        kw = gridweight(kgrid, knum, kfval, kidx)
        evval = conf%evvec(kidx)*kw + conf%evvec(kidx+1_ik)*(1.0_rk - kw)

        f = -1.0_rk*conf%qsell*kfval + beta*evval

    end function kdnObj

    subroutine kupTarget(kupstar, evupstar, conf)
        type(configurations), intent(inout) :: conf
        real(rk), intent(out) :: kupstar, evupstar
        type(nlopt_opt) :: opt
        real(rk), dimension(1) :: guess
        real(rk) :: fmax
        integer :: stat

        guess = kgrid(knum/2)

        ! call create(opt, algorithm_from_string("LN_SBPLX"), int(knum, 4))

        ! remember that the last number is in int4 but not 8
        call create(opt, algorithm_from_string("LN_SBPLX"), 1)
        call opt%set_lower_bounds(kgrid(1))
        call opt%set_upper_bounds(kgrid(knum))
        associate( f => nlopt_func(kupObj, conf) )
            call opt%set_xtol_rel(tol)
            call opt%set_max_objective(f)
            call opt%optimize(guess, fmax, stat)
        end associate
        call destroy(opt)

        if (stat < NLOPT_SUCCESS) then
            write (*, "(A,I5)") "NLopt failed with code ", stat
        end if

        kupstar = guess(1)
        evupstar = fmax

    end subroutine kupTarget

    function kupObj(guess, gradient, func_data) result(f)
        real(rk), dimension(:), intent(in) :: guess
        real(rk), dimension(:), intent(inout), optional :: gradient
        class(*), intent(in), optional :: func_data
        type(configurations) :: conf
        integer(ik) :: kidx
        real(rk) :: f, kw, evval, kfval
        real(rk), dimension(knum) :: evvec

        select type(func_data)
            type is (configurations)
            conf = func_data
        end select

        IF (PRESENT(gradient)) THEN
            ! Why?!
        END IF

        kfval = guess(1)

        ! call linear_interpolation(kgrid, knum, kfval, kidx, kw)
        kidx = gridlookup(kgrid, knum, kfval)
        kw = gridweight(kgrid, knum, kfval, kidx)
        evval = conf%evvec(kidx)*kw + conf%evvec(kidx+1_ik)*(1.0_rk - kw)

        f = -1.0_rk*conf%Qbuy*kfval + beta*evval

        ! call pprint(guess)
        ! write(*, *) "dim of guess is ", size(guess)

        ! write(*, '(2(a, f9.6), a, I4)') "guess(1) = ", guess, ", evval = ", evval

    end function kupObj

    subroutine initGrids()
        integer(ik) :: i
        real(rk) :: tempterm(enum, enum)

        call tauchen(0.0_rk, sigma_e, rho_e, 2.575_rk, enum, egrid, pie)
        ! call rouwenhorst(rho_e, sigma_e, enum, egrid, pie)
        egrid = exp(egrid)
        kgrid = linspace(kbounds(1), kbounds(2), knum)

        tempterm = eye(enum)
        do i = 1, 1000, 1
            tempterm = matmul(tempterm, pie)
        enddo
        piess = tempterm(1, :)

    end subroutine initGrids

    subroutine initSol(sol)

        type(solutions), intent(inout) :: sol

        if (allocated(sol%v)) deallocate(sol%v)
        if (allocated(sol%gk)) deallocate(sol%gk)
        if (allocated(sol%kupmat)) deallocate(sol%kupmat)
        if (allocated(sol%kdnmat)) deallocate(sol%kdnmat)


        allocate(sol%v(knum, enum), source=0.0_rk)
        allocate(sol%gk(knum, enum), source=0.0_rk)
        allocate(sol%kupmat(knum, enum), source = 0.0_rk)
        allocate(sol%kdnmat(knum, enum), source = 0.0_rk)

    end subroutine initSol

    function evaluate(evvec, kfval, Jval) result(fval)
        integer(ik) :: kidx
        real(rk) :: evvec(knum), kfval, fval, kw, evval, Jval
        intent(in) :: evvec, kfval

        kidx = gridlookup(kgrid, knum, kfval)
        kw = gridweight(kgrid, knum, kfval, kidx)
        ! call linear_interpolation(kgrid, knum, kfval, kidx, kw)
        evval = evvec(kidx)*kw + evvec(kidx+1_ik)*(1.0_rk - kw)
        fval = -Jval*kfval + beta*evval

    end function evaluate


end module ReplicateLanteri2018
