module io

    use iso_Fortran_env, only: rk => real64, ik => int32
    implicit none
    private
    public operator(.f.), operator(.d.), pprint, head, tail

    interface operator(.f.)
        module procedure file_exists
    end interface

    interface operator(.d.)
        module procedure directory_exists
    end interface

    ! print the whole matrix
    interface pprint
        module procedure pprint_1d_int
        module procedure pprint_1d_real
        module procedure pprint_2d_int
        module procedure pprint_2d_real
    end interface

    ! print the head of matrix, default 10 rows
    interface head
        module procedure head_int
        module procedure head_real
    end interface head

    ! print the tail of matrix, default 10 rows
    interface tail
        module procedure tail_int
        module procedure tail_real
    end interface tail

contains

    function file_exists(filename) result(res)
        implicit none
        character(len=*),intent(in) :: filename
        logical                     :: res
        ! Check if the file exists
        inquire( file=trim(filename), exist=res )
    end function file_exists

    function directory_exists(dirname) result(res)
        implicit none
        character(len=*),intent(in) :: dirname
        logical                     :: res
        ! Check if the file exists
        inquire( file=trim(dirname // '.' ), exist=res )
    end function directory_exists

    ! subroutine pprint_3d_real(mat, 1dbegin, 1dend, 2dbegin, 2dend, 3dbegin, 3dend)
    !     real(rk), intent(in) :: mat(:, :, :)
    !     integer(ik), optional :: 1dbegin, 1dend, 2dbegin, 2dend, 3dbegin, 3dend
    !     integer(ik) :: 1b, 1e, 2b, 2e, 3b, 3e, res(3), i, j, k
    !     character(len=30) :: f, str(3)

    !     res = shape(mat)
    !     1b = 1; if (present(1dbegin)) 1b = 1dbegin
    !     1e = res(1); if (present(1dend)) 1e = 1dend
    !     2b = 2; if (present(2dbegin)) 2b = 2dbegin
    !     2e = res(2); if (present(2dend)) 2e = 2dend
    !     3b = 3; if (present(3dbegin)) 3b = 3dbegin
    !     3e = res(3); if (present(3dend)) 3e = 3dend

    !     write(unit=str, fmt='(I0)') res
    !     f = trim('(' // trim(str(2)) // 'g13.5)')

    ! end subroutine pprint_3d_real

    subroutine pprint_2d_real(mat, rowbegin, rowend, colbegin, colend)
        real(rk), intent(in) :: mat(:, :)
        integer(ik), optional :: rowbegin, rowend, colbegin, colend
        integer(ik) :: res(2), i, j, rb, re, cb, ce
        character(len=30) :: f, str(2)
        res = shape(mat)

        rb = 1; if (present(rowbegin)) rb = rowbegin
        re = res(1); if (present(rowend)) re = rowend
        cb = 1; if (present(colbegin)) cb = colbegin
        ce = res(2); if (present(colend)) ce = colend

        write(unit=str, fmt='(I0)') res
        f = trim('(' // trim(str(2)) // 'g13.5)')
        do i = rb, re, 1
            write(*, f) ( mat(i,j), j = cb, ce, 1 )
        enddo
    end subroutine pprint_2d_real

    subroutine pprint_2d_int(mat, rowbegin, rowend, colbegin, colend)
        integer(ik), intent(in) :: mat(:, :)
        integer(ik), optional :: rowbegin, rowend, colbegin, colend
        integer(ik) :: res(2), i, j, rb, re, cb, ce
        character(len=30) :: f, s, str(2)
        res = shape(mat)

        rb = 1; if (present(rowbegin)) rb = rowbegin
        re = res(1); if (present(rowend)) re = rowend
        cb = 1; if (present(colbegin)) cb = colbegin
        ce = res(2); if (present(colend)) ce = colend

        write(unit=str, fmt='(I0)') res
        write(unit=s, fmt='(I0)') maxval(mat)
        write(unit=s, fmt='(I0)') len_trim(s) + 1
        f = trim('(' // trim(str(2)) // 'i' // trim(s) // ')')
        do i = rb, re, 1
            write(*, f) ( mat(i,j), j = cb, ce, 1 )
        enddo
    end subroutine pprint_2d_int

    subroutine pprint_1d_real(mat, rowbegin, rowend)
        real(rk), intent(in) :: mat(:)
        integer(ik), optional :: rowbegin, rowend
        integer(ik) :: res, i, rb, re
        character(len=30) :: f, s, str

        res = size(mat)
        rb = 1; if (present(rowbegin)) rb = rowbegin
        re = res; if (present(rowend)) re = rowend

        write(unit=str, fmt='(I0)') res
        f = trim('(' // trim(str) // 'g13.5)')
        do i = rb, re, 1
            write(*, f) ( mat(i) )
        enddo

    end subroutine pprint_1d_real

    subroutine pprint_1d_int(mat, rowbegin, rowend)
        integer(ik), intent(in) :: mat(:)
        integer(ik), optional :: rowbegin, rowend
        integer(ik) :: res, i, rb, re
        character(len=30) :: f, s, str

        res = size(mat)
        rb = 1; if (present(rowbegin)) rb = rowbegin
        re = res; if (present(rowend)) re = rowend

        write(unit=str, fmt='(I0)') res
        write(unit=s, fmt='(I0)') maxval(mat)
        write(unit=s, fmt='(I0)') len_trim(s) + 1
        f = trim('(' // trim(str) // 'i' // trim(s) // ')')
        do i = rb, re, 1
            write(*, f) ( mat(i) )
        enddo

    end subroutine pprint_1d_int

    subroutine head_real(mat, num)
        real(rk), intent(in) :: mat(:, :)
        integer(ik), optional :: num
        integer(ik) :: n, res(2), i, j
        character(len=30) :: f, str(2)
        res = shape(mat)
        write(unit=str, fmt='(I0)') res
        n = 10_ik
        if (present(num)) n = num
        f = trim('(' // trim(str(2)) // 'g13.5)')
        do i = 1, n, 1
            write(*, f) ( mat(i,j), j=1, res(2), 1 )
        enddo
    end subroutine head_real

    subroutine head_int(mat, num)
        integer(ik), intent(in) :: mat(:, :)
        integer(ik), optional :: num
        integer(ik) :: n, res(2), i, j
        character(len=30) :: f, s, str(2)
        res = shape(mat)
        write(unit=str, fmt='(I0)') res
        write(unit=s, fmt='(I0)') maxval(mat)
        write(unit=s, fmt='(I0)') len_trim(s) + 1
        n = 10_ik
        if (present(num)) n = num
        f = trim('(' // trim(str(2)) // 'i' // trim(s) // ')')
        do i = 1, n, 1
            write(*, f) ( mat(i,j), j=1, res(2), 1 )
        enddo
    end subroutine head_int

    subroutine tail_real(mat, num)
        real(rk), intent(in) :: mat(:, :)
        integer(ik), optional :: num
        integer(ik) :: n, res(2), i, j
        character(len=30) :: f, str(2)
        res = shape(mat)
        write(unit=str, fmt='(I0)') res
        n = 10_ik
        if (present(num)) n = num
        f = trim('(' // trim(str(2)) // 'g13.5)')
        do i = res(1) - n + 1, res(1), 1
            write(*, f) ( mat(i,j), j=1, res(2), 1 )
        enddo
    end subroutine tail_real

    subroutine tail_int(mat, num)
        integer(ik), intent(in) :: mat(:, :)
        integer(ik), optional :: num
        integer(ik) :: n, res(2), i, j
        character(len=30) :: f, s, str(2)
        res = shape(mat)
        write(unit=str, fmt='(I0)') res
        write(unit=s, fmt='(I0)') maxval(mat)
        write(unit=s, fmt='(I0)') len_trim(s) + 1
        n = 10_ik
        if (present(num)) n = num
        f = trim('(' // trim(str(2)) // 'i' // trim(s) // ')')
        do i = res(1) - n + 1, res(1), 1
            write(*, f) ( mat(i,j), j=1, res(2), 1 )
        enddo
    end subroutine tail_int

end module io
