program hw5final
    
    real, allocatable, dimension(:,:,:) :: cube
    real, allocatable, dimension(:,:,:) :: cubenext
    real, allocatable, dimension(:,:,:) :: cdensity
    real, allocatable, dimension(:,:,:) :: e
    real, allocatable, dimension(:,:,:) :: r
    
    integer :: i, j, k, intstep, imaging, reader, ct, cs, halfs, tolchoice
    real :: stepsize
    real :: lim, dv
    real :: iteration

    !---------------------------------------------------------Begin
    ! Let's first create a grid
    ! Note that fortran starts arrays at slot 1, not 0
    print *, "Small cube (0) or big cube (1)?"
    read *, ct

    print *, "Big tolerance (1) or low tolerance (0)?"
    read *, tolchoice

    if (ct .EQ. 0) then
        cs = 101
        halfs = 51

        allocate(cube(cs,cs,cs))
        allocate(cdensity(cs,cs,cs))
        intstep = cs
        stepsize = 0.2

        ! Set all values to 0
        cube = 0
        cdensity = 0

        ! Introduce the quadripole
        ! Each charge is a finite
        cdensity(46,46,51) = 1 / stepsize**3
        cdensity(47,46,51) = 1 / stepsize**3
        cdensity(47,47,51) = 1 / stepsize**3
        cdensity(46,47,51) = 1 / stepsize**3
        cdensity(46,46,52) = 1 / stepsize**3
        cdensity(47,46,52) = 1 / stepsize**3
        cdensity(47,47,52) = 1 / stepsize**3
        cdensity(46,47,52) = 1 / stepsize**3

        cdensity(46,56,51) = -1 / stepsize**3
        cdensity(46,55,51) = -1 / stepsize**3
        cdensity(47,56,51) = -1 / stepsize**3
        cdensity(47,55,51) = -1 / stepsize**3
        cdensity(46,56,52) = -1 / stepsize**3
        cdensity(46,55,52) = -1 / stepsize**3
        cdensity(47,56,52) = -1 / stepsize**3
        cdensity(47,55,52) = -1 / stepsize**3
        
        cdensity(56,56,51) = 1 / stepsize**3
        cdensity(56,55,51) = 1 / stepsize**3
        cdensity(55,56,51) = 1 / stepsize**3
        cdensity(55,55,51) = 1 / stepsize**3
        cdensity(56,56,52) = 1 / stepsize**3
        cdensity(56,55,52) = 1 / stepsize**3
        cdensity(55,56,52) = 1 / stepsize**3
        cdensity(55,55,52) = 1 / stepsize**3
        
        cdensity(56,46,51) = -1 / stepsize**3
        cdensity(56,47,51) = -1 / stepsize**3
        cdensity(55,46,51) = -1 / stepsize**3
        cdensity(55,47,51) = -1 / stepsize**3
        cdensity(56,46,52) = -1 / stepsize**3
        cdensity(56,47,52) = -1 / stepsize**3
        cdensity(55,46,52) = -1 / stepsize**3
        cdensity(55,47,52) = -1 / stepsize**3

    else if (ct .EQ. 1) then
        cs = 201
        halfs = 101

        allocate(cube(cs,cs,cs))
        allocate(cdensity(cs,cs,cs))
        intstep = cs
        stepsize = 0.1

        ! Set all values to 0
        cube = 0
        cdensity = 0

        ! Introduce the quadripole
        ! Each charge is a finite
        cdensity(91,91,halfs) = 1 / stepsize**3
        cdensity(92,91,halfs) = 1 / stepsize**3
        cdensity(92,92,halfs) = 1 / stepsize**3
        cdensity(91,92,halfs) = 1 / stepsize**3
        cdensity(91,91,102) = 1 / stepsize**3
        cdensity(92,91,102) = 1 / stepsize**3
        cdensity(92,92,102) = 1 / stepsize**3
        cdensity(91,92,102) = 1 / stepsize**3

        cdensity(91,111,halfs) = -1 / stepsize**3
        cdensity(91,110,halfs) = -1 / stepsize**3
        cdensity(92,111,halfs) = -1 / stepsize**3
        cdensity(92,110,halfs) = -1 / stepsize**3
        cdensity(91,111,102) = -1 / stepsize**3
        cdensity(91,110,102) = -1 / stepsize**3
        cdensity(92,111,102) = -1 / stepsize**3
        cdensity(92,110,102) = -1 / stepsize**3
        
        cdensity(111,111,halfs) = 1 / stepsize**3
        cdensity(111,110,halfs) = 1 / stepsize**3
        cdensity(110,111,halfs) = 1 / stepsize**3
        cdensity(110,110,halfs) = 1 / stepsize**3
        cdensity(111,111,102) = 1 / stepsize**3
        cdensity(111,110,102) = 1 / stepsize**3
        cdensity(110,111,102) = 1 / stepsize**3
        cdensity(110,110,102) = 1 / stepsize**3
        
        cdensity(111,91,halfs) = -1 / stepsize**3
        cdensity(111,92,halfs) = -1 / stepsize**3
        cdensity(110,91,halfs) = -1 / stepsize**3
        cdensity(110,92,halfs) = -1 / stepsize**3
        cdensity(111,91,102) = -1 / stepsize**3
        cdensity(111,92,102) = -1 / stepsize**3
        cdensity(110,91,102) = -1 / stepsize**3
        cdensity(110,92,102) = -1 / stepsize**3

    end if

    !---------------------------------------------------------LOOPING
    open(30,file="2iterationchange.txt")
    allocate(cubenext(cs,cs,cs))
    
    if (tolchoice .EQ. 0) then
        lim = 0.00001 * 101**3
    else
        lim = 0.00001 * 201**3
    end if

    iteration = 1
    dv = 100

    do while (dv.GT. lim)
        call SYSTEM('cls')
        print *, "Iteration", iteration
        print *, "Previous dv: ", dv
        dv = 0

        do i=2,(intstep-1)
            do j=2,(intstep-1)
                do k=2,(intstep-1)
                    ! These additions are separated in order to keep from having one massive line of code
                    cubenext(i,j,k) = (1.0/6.0) * (cube(i-1,j,k) + cube(i+1,j,k))
                    cubenext(i,j,k) = cubenext(i,j,k) + (1.0/6.0) * (cube(i,j-1,k)+cube(i,j+1,k))
                    cubenext(i,j,k) = cubenext(i,j,k) + (1.0/6.0) * (cube(i,j,k-1)+cube(i,j,k+1))
                    cubenext(i,j,k) = cubenext(i,j,k) + ((1.0/6.0) * cdensity(i,j,k) * stepsize)
                    !cubenext(i,j,k) = cubenext(i,j,k) + ((1.0/6.0) * cdensity(i,j,k))
                    dv = dv + abs(cubenext(i,j,k) - cube(i,j,k))
                end do
            end do
        end do
        cube = cubenext
        write(30,*) iteration, dv
        iteration = iteration + 1
    end do
    CALL SYSTEM('gnuplot -p 1iterationchange.plt')
    close(30)
    deallocate(cubenext)

    !---------------------------------------------------------Graphing/Output

    !---------------------------------------------------------part a

    ! Open the output files
    open(42,file="2allpots.txt")
    open(43,file="2xypots.txt")
    open(44,file="2xzpots.txt")

    print *, "Plot density (2, 5, and 10 are recommended): "
    read *, imaging

    ! Plotting cube 
    do i=1,intstep
        call SYSTEM('cls')
        call loadingbar(i, intstep)
        do j=1,intstep
            do k=1,intstep
                !3D Image
                if (((mod(i,imaging) - 1) .EQ. 0) .AND. ((mod(j,imaging) - 1) .EQ. 0) .AND. ((mod(k,imaging) - 1) .EQ. 0)) then
                    if (cube(i,j,k) .NE. 0) then
                        write(42,*) ((i-1)*stepsize),((j-1)*stepsize),((k-1)*stepsize),cube(i,j,k)
                    end if
                end if
            end do
        end do
    end do
    close(42)

    ! Plotting xy
    do i=1,intstep
        call SYSTEM('cls')
        call loadingbar(i, intstep)
        do j=1,intstep
            write(43,*) ((i-1)*stepsize), ((j-1)*stepsize), cube(i,j,halfs)
        end do
    end do
    close(43)

    ! Plotting xz
    do i=1,intstep
        call SYSTEM('cls')
        call loadingbar(i, intstep)
        do k=1,intstep
            write(44,*) ((i-1)*stepsize), ((k-1)*stepsize), cube(i,halfs,k) !* 1000)
        end do
    end do
    close(44)

    CALL SYSTEM('gnuplot -p 1volumeplot.plt')
    CALL SYSTEM('gnuplot -p 1sliceplot.plt')

    print *, "Part a complete. Continue?"
    read *, reader

    !---------------------------------------------------------part b

    ! Create r values to plot
    allocate(r(cs,cs,cs))
    open(42,file="2rx.txt")
    open(43,file="2rdiag.txt")
    do i=1,intstep
        do j=1,intstep
            do k=1,intstep
                r(i,j,k) = sqrt((i-(halfs*stepsize))**2+(j-(halfs*stepsize))**2+(k-(halfs*stepsize))**2)
                if ((i .GE. halfs) .AND. (i .EQ. j) .AND. (j .EQ. k)) then
                    write(43,*) r(i,j,k), cube(i,j,k)
                end if
                if ((i .GE. halfs) .AND. (j .EQ. halfs) .AND. (k .EQ. halfs)) then
                    write(42,*) r(i,j,k), cube(i,j,k)
                end if
            end do
        end do
    end do

    close(42)
    close(43)

    CALL SYSTEM('gnuplot -p 1rvplot.plt')

    print *, "Part b complete. Continue?"
    read *, reader

    !---------------------------------------------------------part b

    ! Find the whole electric field
    allocate(e(cs,cs,cs))
    e = 0
    open(42,file="2ex.txt")
    open(43,file="2ediag.txt")

    do i=2,(intstep-1)
        call SYSTEM('cls')
        call loadingbar(i,(intstep-1))
        do j=2,(intstep-1)
            do k=2,(intstep-1)
                ! Another separation to keep the lines shorter
                e(i,j,k) = 0.5*(((cube(i-1,j,k)-cube(i,j,k))/stepsize)+((cube(i,j,k)-cube(i+1,j,k))/stepsize))
                e(i,j,k) = e(i,j,k) + 0.5*(((cube(i,j-1,k)-cube(i,j,k))/stepsize)+((cube(i,j,k)-cube(i,j+1,k))/stepsize))
                e(i,j,k) = e(i,j,k) + 0.5*(((cube(i,j,k-1)-cube(i,j,k))/stepsize)+((cube(i,j,k)-cube(i,j,k+1))/stepsize))
                if ((i .GE. halfs) .AND. (i .EQ. j) .AND. (j .EQ. k)) then
                    write(43,*) r(i,j,k), e(i,j,k)
                end if
                if ((i .GE. halfs) .AND. (j .EQ. halfs) .AND. (k .EQ. halfs)) then
                    write(42,*) r(i,j,k), e(i,j,k)
                end if
            end do
        end do
    end do

    close(42)
    close(43)

    CALL SYSTEM('gnuplot -p 1replot.plt')

    print *, "Part c complete. Continue?"
    read *, reader

    deallocate(e)
    deallocate(cube)
    deallocate(r)

    print *, "Done..."
    read *, reader

end program hw5final

subroutine loadingbar(i, intstep)
    integer :: i, intstep
    real :: percent
    percent = ((real(i)/real(intstep)) * 100.0)
    print *, percent
    if (percent .EQ. 100) then 
        print *, '.____________________.  complete'
    else if (percent .GT. 90.0) then
        print *, '.__________________  .'
    else if (percent .GT. 80.0) then
        print *, '.________________    .'
    else if (percent .GT. 70.0) then
        print *, '.______________      .'
    else if (percent .GT. 60.0) then
        print *, '.____________        .'
    else if (percent .GT. 50.0) then
        print *, '.__________          .'
    else if (percent .GT. 40.0) then
        print *, '.________            .'
    else if (percent .GT. 30.0) then
        print *, '.______              .'
    else if (percent .GT. 20.0) then
        print *, '.____                .'
    else if (percent .GT. 10.0) then
        print *, '.__                  .'
    else
        print *, '.                    .'
    end if

end subroutine