program hw6final
    double precision, allocatable, dimension(:) :: string, stringlast, stringnext !String y values associated with data location (x)
    double precision, allocatable, dimension(:) :: ysingle !String's fourier lists
    integer :: xpoints, tpoints, ftpoints, psplit !The sushi of the problem (and the plot split)
    integer :: ppoint !Pluck starting point (index location)
    integer :: i, t !Arbitrary integers for loops
    double precision :: pspot !starting location (in meters) for the pluck
    double precision :: r, c, k, yo, expo !numerical speed, speed, width, and y coefficients
    double precision :: dx, dt, time, x, ynow !x and time steps, dummy y value, time
    double precision :: freq, w, finalfreq, yfour, pi, printer, freqprint, wprint !frequency, final frequency, weight, pi, arbitrary printables

    !------------------------------------------INITIALIZATION-----------------------------------------------
    !Create the string
    print *, "Numerical coefficient (r):"
        read *, r
    call sushi(xpoints, tpoints, dx, dt, r, k, pspot, ppoint, yo, ftpoints, c) !This function just initiates the initial values
    allocate(string(xpoints+1))
    allocate(stringlast(xpoints+1))
    allocate(stringnext(xpoints+1))
    allocate(ysingle(ftpoints))
    string = 0

    pi = 4.0 * atan(1.0) !pi
    
    !Set up the guassian pluck
    do i = 2,(xpoints)
        x = i * dx
        expo = -k * (x - pspot)**2
        string(i) = yo * exp(expo)
    end do
    stringlast = string

    open(42,file="stringwiggle.txt") !Prepare the position output file
    open(43,file="stringfreqs.txt") !Prepare the fourier output file  
    
    !------------------------------------------PRIMARY LOOP-----------------------------------------------
    !For each step in time, calculate a y value for each step along the string
    print *, "Plot split: "
    read *, psplit

    do t = 0,tpoints !Time loop
        call SYSTEM('cls')
        call loadingbar(t, tpoints)
        time = t * dt !This is the time we output

        do i = 1,(xpoints+1) !X loop
            x = i * dx !This is the x we output
            ynow = string(i) !dummy y value for clarity

            if ((i .EQ. 1) .OR. (i .EQ. (xpoints+1))) then
                continue
            end if !Gotta keep the ends at zero

            stringnext(i) = (r**2 * (string(i-1) + string(i+1))) + (2 * (1 - r**2) * string(i)) - stringlast(i)

            if ((i .EQ. nint(xpoints/2.0)) .AND. (t .LE. ftpoints))then
                ysingle(t) = stringnext(i)
            end if

            if (mod(t,psplit) .EQ. 0) then
                write(42,*) time, x, string(i), stringlast(i), (r**2 * (string(i-1) + string(i+1))) 
            end if

        end do

        stringlast = string
        string = stringnext

    end do

    close(42)

    !------------------------------------------FOURIER LOOP-----------------------------------------------
    freq = 0
    print *, "Max frequency to check:"
    read *, finalfreq
    !open(44,file="troubleshoot.txt")
    do while (freq .LE. finalfreq) !For each frequency in this range...
        call SYSTEM('cls')
        print *, freq
        yfour=0
        do t = 1,(ftpoints) !... find the summation of the fourier transform with these time and y values.
            yfour = yfour + (dt * ysingle(t) * cos(2.0 * pi * freq * real(t) * dt)) !note that yfour is the total value of the summation
        end do
        printer = exp(10000.0 * yfour**2) !Power scale is arbitrary, so this gets us up to all positive values and just under 1
        w = 2 * pi * freq
        !Thin out the peaks
        if (printer .GT. 1) then
            freqprint = freq - 0.1
            wprint = w - 0.1
            write(43,*) freqprint, wprint, 1.0000

            write(43,*) freq, w, printer

            freqprint = freq + 0.1
            wprint = w + 0.1
            write(43,*) freqprint, wprint, 1.0000
        else 
            write(43,*) freq, w, printer
        end if
        freq = freq + (1/((ftpoints+1) * dt)) !Set to test the next frequency (j=f N dt)
    end do
    !close(44)
    close(43)

    !------------------------------------------CALL GNUPLOT-----------------------------------------------
    CALL SYSTEM('gnuplot -p stringwiggle.plt')
    CALL SYSTEM('gnuplot -p stringfreqs.plt')

    !------------------------------------------CONFIRMATIONS-----------------------------------------------
    print *, "dt:", dt
    print *, "r:", r
    print *, "dx:", dx
    print *, "c:", c
    printer = 1/(2 * dt)
    print *, "Nyquist freq:", printer

    read *, r

end program hw6final

subroutine sushi(xpoints, tpoints, dx, dt, r, k, pspot, ppoint, yo, ftpoints, c)
    integer :: xpoints, tpoints, ppoint
    integer :: echoice, ftpoints
    double precision :: dx, dt, totalt, r
    double precision :: ft, mass, c, length
    double precision :: fwhm, k, pspot, yo

    print *, "Manual entry (1)?"
    read *, echoice

    if (echoice .EQ. 1) then
        print *, "All manual entries should be integers."
        print *, "String Length (meters):"
        read *, length
        print *, "Points to simulate:"
        read *, points
        print *, "Force of tension (N): "
        read *, ft
        print *, "Mass of string (g):"
        read *, mass
        print *, "Full width half max: "
        read *, fwhm
        print *, "Pluck location (meters): "
        read *, pspot
        print *, "Pluck height divided by length of string: "
        read *, yo
        c = sqrt(ft * length / mass)
        totalt = length/c
    else
        length = 1.0
        points = 1000
        ft = 850.0
        mass = .025
        fwhm = 0.1
        pspot = 0.4
        yo = 0.25
        c = sqrt(ft * length / mass)
        totalt = length/c
    end if

    dx = length/points
    dt = r * dx / c
    tpoints = nint(2.5 * totalt / dt)
    ftpoints = nint(2.0 * totalt / dt) - 1
    xpoints = nint(length / dx)
    k = 1/(2 * (fwhm/2.36)**2)
    ppoint = nint(pspot/dx)

end subroutine sushi

subroutine loadingbar(t, tpoints)
    integer :: t, tpoints
    real :: percent
    percent = ((real(t)/real(tpoints)) * 100.0)
    print *, percent, "%"
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