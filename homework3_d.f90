program homework3_d

    real :: mu 
    real :: xin
    ! Xvalues for averaging
    real :: x, dx
    real :: dxaverage
    character(30) :: filename

    integer :: i

    print *, 'Enter mu:'
    read *, mu
    write(filename,"('dxaverage',i0,'.txt')") nint(mu)

    open(42,file=filename)
    open(43,file="dxon.txt")

    dxaverage = 0

    do i=0,500
        xin = 0.2 + (i * .0001)
        x = (xin * mu) * (1 - xin)
        dx = x - xin 

        if (dx .LT. 0) then
            dx = (-1 * dx)
        end if

        dxaverage = (dx + dxaverage)/(i + 1)
        write(43,*) i, dx
    end do 

    write(42,*) dxaverage
    ! Call on gnu
    CALL SYSTEM('gnuplot -p plotcoded.plt')

end program homework3_d