! gfortran -o local_Chern.out local_Chern.f90 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
program local_Chern
    use,intrinsic :: iso_fortran_env
    implicit none
    !parameter
    ! 前回のループのloopmax
    integer, parameter :: loop_0 = 0
    ! xを進める回数
    integer, parameter :: loopmax = 10
    ! xの初期値
    double precision, parameter :: pos_xmin = -1.
    ! xの進み幅
    double precision, parameter :: dx = 0.01
    integer, parameter :: N = 2869
    integer, parameter :: allhop = 5620
    double precision, parameter :: pos_y = -1.707106781
    double precision, parameter :: kappa = 0.1
    double precision, parameter :: U = 1.5
    double precision, parameter :: tildemu = 0.0
    double precision, parameter :: hsq = 0.1
    double precision, parameter :: lambda_R = 0.5
    double precision, parameter :: lambda_D = 0.0
    double precision, parameter :: pi = 4*atan(1.)
    character(7), parameter :: hop_file = "hop.txt"
    character(19), parameter :: Delta_file = "pair_potential.txt", PN_file = "particle_number.txt"
    character(24), parameter :: pos_file = "coordinates.txt", C_file = "Cherns.txt"
    character(10), parameter :: result_file = "result.txt"
    !variable
    integer :: unit_write_result, loop
    double precision :: mu, h_z = sqrt(hsq)
    double precision :: PN(2*N)
    complex(kind(0d0)) :: Delta(N), Localizer(8 * N, 8 * N)
    double precision :: pos_x
    double precision :: Cherns(loopmax) = 0
    ! for ZHEEVD of lapack
    integer :: INFO, IWORK(1)
    double precision :: W(8 * N), RWORK(8 * N)
    complex(kind(0d0)) :: WORK(8 * N + 1)
    ! to clock time
    integer(int64) :: time_begin_c,time_end_c, CountPerSec, CountMax

    call system_clock(time_begin_c, CountPerSec, CountMax)

    call read_quantities()
    mu = tildemu - U*sum(PN) / (2 * N)

    do loop = 1, loopmax
        pos_x = pos_xmin + dx * (loop - 1)
        call make_Localizer()
        call ZHEEVD('N', 'U', 8 * N, Localizer, 8 * N, W, WORK, 8 * N + 1,&
        RWORK, 8 * N, IWORK, 1, INFO)
        call input_Chern()
        call write_files_E()
        print *, loop
    end do

    call system_clock(time_end_c)
    call write_files_result()
contains
    subroutine make_Localizer()
        integer :: i, j, k, l, m
        double precision :: x, y
        complex(kind(0d0)) :: L_ij

        Localizer = 0

        ! diagonal elements
        do k = 1, N
            ! spin up
            i = 2 * k - 1
            j = i
            L_ij = -mu - U * PN(2 * k) + h_z

            Localizer(i, j) = L_ij
            Localizer(i + 2 * N, j + 2 * N) = -L_ij

            Localizer(i + 4 * N, j + 4 * N) = -L_ij
            Localizer(i + 6 * N, j + 6 * N) = L_ij

            ! spin down
            i = 2 * k
            j = i
            L_ij = -mu - U * PN(2 * k - 1) - h_z !

            Localizer(i, j) = L_ij
            Localizer(i + 2 * N, j + 2 * N) = -L_ij

            Localizer(i + 4 * N, j + 4 * N) = -L_ij
            Localizer(i + 6 * N, j + 6 * N) = L_ij
        end do

        ! hopping elements
        open(newunit = unit_write_result, file = hop_file)
        do k = 1, allhop
            read(unit_write_result, *) l, m
            l = l + 1
            m = m + 1

            ! spin up
            i = 2 * l - 1
            j = 2 * m - 1
            L_ij = -1.

            Localizer(i, j) = L_ij
            Localizer(j, i) = L_ij
            Localizer(i + 2 * N, j + 2 * N) = -L_ij
            Localizer(j + 2 * N, i + 2 * N) = -L_ij

            Localizer(i + 4 * N, j + 4 * N) = -L_ij
            Localizer(j + 4 * N, i + 4 * N) = -L_ij
            Localizer(i + 6 * N, j + 6 * N) = L_ij
            Localizer(j + 6 * N, i + 6 * N) = L_ij

            ! spin down
            i = 2 * l
            j = 2 * m
            L_ij = -1.

            Localizer(i, j) = L_ij
            Localizer(j, i) = L_ij
            Localizer(i + 2 * N, j + 2 * N) = -L_ij
            Localizer(j + 2 * N, i + 2 * N) = -L_ij

            Localizer(i + 4 * N, j + 4 * N) = -L_ij
            Localizer(j + 4 * N, i + 4 * N) = -L_ij
            Localizer(i + 6 * N, j + 6 * N) = L_ij
            Localizer(j + 6 * N, i + 6 * N) = L_ij
        end do
        close (unit_write_result)

        !spin-orbit coupling
        open (newunit = unit_write_result, file=hop_file)
        do k = 1, allhop
            read (unit_write_result, *) l, m, x, y
            l = l + 1
            m = m + 1

            !iup, jdown
            i = 2*l - 1
            j = 2*m
            L_ij = lambda_R*cmplx(-x, y, kind(0d0)) + lambda_D*cmplx(-y, x, kind(0d0))

            Localizer(i, j) = L_ij
            Localizer(j, i) = conjg(L_ij)
            Localizer(i + 2 * N, j + 2 * N) = -conjg(L_ij)
            Localizer(j + 2 * N, i + 2 * N) = -L_ij

            Localizer(i + 4 * N, j + 4 * N) = -L_ij
            Localizer(j + 4 * N, i + 4 * N) = -conjg(L_ij)
            Localizer(i + 6 * N, j + 6 * N) = conjg(L_ij)
            Localizer(j + 6 * N, i + 6 * N) = L_ij

            !idown, jup
            i = 2*l
            j = 2*m - 1
            L_ij = lambda_R*cmplx(x, y, kind(0d0)) + lambda_D*cmplx(y, x, kind(0d0))

            Localizer(i, j) = L_ij
            Localizer(j, i) = conjg(L_ij)
            Localizer(i + 2 * N, j + 2 * N) = -conjg(L_ij)
            Localizer(j + 2 * N, i + 2 * N) = -L_ij

            Localizer(i + 4 * N, j + 4 * N) = -L_ij
            Localizer(j + 4 * N, i + 4 * N) = -conjg(L_ij)
            Localizer(i + 6 * N, j + 6 * N) = conjg(L_ij)
            Localizer(j + 6 * N, i + 6 * N) = L_ij
        end do
        close(unit_write_result)

        !pair potential
        do k = 1, N
            !(L_ij)c_{idown}^dag c_{iup}^dag
            i = 2*k
            j = 2*k - 1 + 2*N
            L_ij = -Delta(k)

            Localizer(i, j) = L_ij
            Localizer(j, i) = conjg(L_ij)

            Localizer(i + 4 * N, j + 4 * N) = -L_ij
            Localizer(j + 4 * N, i + 4 * N) = -conjg(L_ij)

            !(L_ij)c_{iup}^dag c_{idown}^dag
            i = 2*k - 1
            j = 2*k + 2*N
            L_ij = Delta(k)

            Localizer(i, j) = L_ij
            Localizer(j, i) = conjg(L_ij)

            Localizer(i + 4 * N, j + 4 * N) = -L_ij
            Localizer(j + 4 * N, i + 4 * N) = -conjg(L_ij)
        end do

        ! position operator
        open(newunit = unit_write_result, file = pos_file)
        do k = 1, N
            read (unit_write_result, *) x, y
            ! spin up
            i = 2 * k - 1
            j = 2 * k - 1 + 4 * N
            L_ij = kappa * cmplx(x - pos_x, y - pos_y, kind(0d0))

            Localizer(i, j) = L_ij
            Localizer(j, i) = conjg(L_ij)
            Localizer(i + 2 * N, j + 2 * N) = L_ij
            Localizer(j + 2 * N, i + 2 * N) = conjg(L_ij)

            ! spin down
            i = 2 * k
            j = 2 * k + 4 * N
            L_ij = kappa * cmplx(x - pos_x, y - pos_y, kind(0d0))

            Localizer(i, j) = L_ij
            Localizer(j, i) = conjg(L_ij)
            Localizer(i + 2 * N, j + 2 * N) = L_ij
            Localizer(j + 2 * N, i + 2 * N) = conjg(L_ij)
        end do
        close(unit_write_result)
    end subroutine make_Localizer

    subroutine input_Chern()
        integer :: i
        double precision :: Chern
        Chern = 0

        do i = 1, 8 * N
            if ( W(i) > 0 ) then
                Chern = Chern + 1
            else if ( W(i) < 0 ) then
                Chern = Chern - 1
            end if
        end do

        Cherns(loop) = Chern / 2
    end subroutine input_Chern

    pure function int32_to_string(num) result(str)
        use, intrinsic :: iso_fortran_env
        implicit none

        integer(int32), intent(in) :: num
        character(:), allocatable :: str  ! retval

        integer(int32), parameter :: Extra_Digits = 0 ! 追加する桁数
        integer(int32) :: num_digits                 ! 整数numの桁数
        integer(int32) :: dgt_digits                 ! 整数numの桁数num_digitsの桁数
        character(:), allocatable :: str_digits
        character(:), allocatable :: fmt

        num_digits = get_digits_of(num) + Extra_Digits ! 整数numの桁数を取得
        dgt_digits = get_digits_of(num_digits)         ! 整数numの桁数num_digitsの桁数を取得

        ! 整数の桁数を文字列に変換
        allocate(character(dgt_digits)::str_digits)
        write(str_digits,'(I0)') num_digits

        ! 書式を作成
        fmt = "(I"//str_digits//"."//str_digits//")"

        ! 整数numを文字列に変換
        allocate(character(num_digits)::str)
        write(str,fmt) num
    end function int32_to_string

    pure function get_digits_of(num) result(num_digit)
        implicit none
        integer(int32),intent(in) :: num
        integer(int32) :: num_digit

        num_digit = int(log10(dble(num)))+1
    end function get_digits_of

    subroutine write_files_E()
        integer :: i
        open(newunit=unit_write_result, file="eigvals/x_" // int32_to_string(loop + loop_0) // ".txt")
            do i = 1, 8 * N
                write(unit_write_result, *) W(i)
            end do
        close(unit_write_result)
    end subroutine write_files_E

    subroutine write_files_result()
        integer :: i
        
        open(newunit=unit_write_result, file=C_file)
            do i = 1, loopmax
                write(unit_write_result, *) Cherns(i)
            end do
        close(unit_write_result)

        open(newunit=unit_write_result, file=result_file)
            write(unit_write_result, *) "Run time=", real(time_end_c - time_begin_c, kind(0d0)) / CountPerSec, "sec"
            write (unit_write_result, *) "N =", N
            write (unit_write_result, *) "U =", U
            write (unit_write_result, *) "tildemu =", tildemu
            write (unit_write_result, *) "h_z^2 =", hsq
            write (unit_write_result, *) "lambda_R =", lambda_R
            write (unit_write_result, *) "lambda_D =", lambda_D
            write (unit_write_result, *) "pos_y =", pos_y
        close(unit_write_result)
    end subroutine write_files_result

    subroutine read_quantities()
        integer :: i
        double precision :: s, t

        open(newunit = unit_write_result, file = Delta_file)
        do i = 1, N
            read(unit_write_result, *) s, t
            Delta(i) = cmplx(s * cos(t), s * sin(t), kind(0d0))
        end do
        close(unit_write_result)

        open(newunit = unit_write_result, file = PN_file)
        do i = 1, N
            read(unit_write_result, *) s, t
            PN(2 * i - 1) = s
            PN(2 * i) = t
        end do
        close(unit_write_result)
    end subroutine read_quantities
end program local_Chern