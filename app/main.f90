program main

    use ReplicateLanteri2018
    implicit none

    type(gpf) :: gp
    type(solutions) :: sol
    type(configurations) :: conf

    if (.not. .d. resDir) call execute_command_line('mkdir ' // resDir);
    if (.not. .d. figDir) call execute_command_line('mkdir ' // figDir);

    call initGrids()
    call initSol(sol)

    conf%zval = zss
    conf%qsell = 0.918790087890625_rk
    conf%wval = 1.087975585937500_rk
    conf%Qbuy = ( eta + (1 - eta)*(conf%qsell + gamma)**(1-s) )**( 1/(1-s) )

    call valueiter(sol, conf)

    ! call pprint(sol%kupmat)
    ! write(*, *) ''
    ! call pprint(sol%kdnmat)

    ! write(*, *) 'v:'
    ! call pprint(sol%v)
    write(*, *) 'gk:'
    call pprint(sol%gk)

    if (drawfig) then
        call gp%options('set terminal pdfcairo enhanced color dashed font "Alegreya, 14" rounded size 16 cm, 9.6 cm')
        call gp%options('set output "' // figDir // 'gk.pdf"')
        call gp%preset(.false.)
        call gp%xlabel('k')
        call gp%ylabel('gk')
        call gp%title('capital decision rule')
        call gp%options('set key outside')
        call gp%plot( &
            kgrid, &
            reshape( (/(1.0_rk-delta)*kgrid, sol%gk(:, 2), sol%gk(:, 5), sol%gk(:, 7)/), (/knum, 4/) ), &
            lspec =' &
                    title "(1-delta)*k" dt 2 lc -1; &
                    title "e = ' // num2str(real(egrid(2), 4)) // '" ls 1; &
                    title "e = ' // num2str(real(egrid(5), 4)) // '" ls 2; &
                    title "e = ' // num2str(real(egrid(7), 4)) // '" ls 3; &
            ')
        ! call gp%plot( &
        !     x1 = kgrid, y1 = (1.0_rk-delta)*kgrid, ls1 = 't "(1-delta)*k" dt 2 lc -1', &
        !     ! x1 = kgrid, y1 = sol%gk(:, 1), ls1 = 't "e = ' // num2str(real(egrid(1), 4)) // '" ls 1', &
        !     x2 = kgrid, y2 = sol%gk(:, 2), ls2 = 't "e = ' // num2str(real(egrid(2), 4)) // '" ls 2', &
        !     x3 = kgrid, y3 = sol%gk(:, enum/2), ls3 = 't "e = ' // num2str(real(egrid(enum/2), 4)) // '" ls 3', &
        !     x4 = kgrid, y4 = sol%gk(:, enum-2), ls4 = 't "e = ' // num2str(real(egrid(enum-2), 4)) // '" ls 4'&
        !     )
        call gp%reset() ! use reset to be able to save next figure to pdf
    endif



end program main
