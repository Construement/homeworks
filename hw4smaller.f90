program hw4smaller

    ! Introduce constants
    real, dimension(9) :: const

    ! Introduce values for analytical solution
    real, dimension(5) :: bog

    ! Other values
    real :: dt, steps, periods

    call valueprep(const, bog, dt, steps, periods)
    call functo(const, bog, dt, steps)

end program hw4smaller

!######################################################################### initialize

subroutine valueprep(const, bog, dt, steps, periods)
    ! Preserve arguments
    real, dimension(9) :: const
    real, dimension(5) :: bog
    real :: dt, steps, periods
    integer :: systype

    ! Set the constants
    const(1) = 1.495978707E11 !meters per AU 
    const(9) = 3.154E7 !seconds per year
    const(2) = 6.67428E-11 * const(9)**2 / const(1)**3 !G in AU^3 kg^-1 year^-2
    const(3) = 4.0 * atan(1.0) !pi
    const(4) = 0.25 !eccentricity
    const(5) = 40.0 !Semimajor axis
    const(6) = const(5) * sqrt(1 - const(4)**2) !semiminor axis
    
    print *, "Pluto & Sun system (1) or custom system (2)?"
    read *, systype
    if (systype .EQ. 1) then
        const(7) = 1.9891E30 !msun
        const(8) = 3.0E-6 * const(7) !mpluto
    else
        print *, "Large mass (kg): "
        read *, const(7)
        print *, "Smaller mass (kg): "
        read *, const(8)
    end if

    ! Let's find the important analytical values
    bog(1) = const(5) * (1 - const(4)) !rp (AU)
    bog(2) = const(5) * (1 + const(4)) !ra (AU)
    bog(5) = sqrt(4 * const(3)**2 * const(5)**3 / (const(2) * (const(7) +  const(8)))) !period (years)
    !bog(3) = sqrt(const(2) * const(7) * ((2 / bog(1)) - (1 / const(5)))) !vperi (AU / year)
    bog(3) = 0.5
    bog(4) = sqrt(const(2) * const(7) * ((2 / bog(2)) - (1 / const(5)))) !vaph (AU / year)

    ! Output values for part a
    open(40,file="partavalues.txt")
    write(40,*) "Semimajor Axis: ", const(5)
    write(40,*) "Semiminor Axis: ", const(6)
    write(40,*) "Perihelion: ", bog(1)
    write(40,*) "Velocity at Perihelion: ", bog(3)
    write(40,*) "Aphelion: ", bog(2)
    write(40,*) "Velocity at Aphelion: ", bog(4)
    write(40,*) "Period of Orbit: ", bog(5)
    close(40)

    ! Now lets set up and run the loop
    print*, "How many periods to run?"
    read *, periods
    print*, "Enter the time step (years):"
    read *, dt
    steps = bog(5) * periods / dt
    steps = nint(steps)

end subroutine valueprep

!######################################################################### selection loop

subroutine functo(const, bog, dt, steps)
    ! Preserve arguments
    real, dimension(9) :: const
    real, dimension(5) :: bog
    real :: dt, steps, periods
    integer :: functionreq

    open(41,file="precoverbeta.txt")

    do
        print *, "Function request: "
        print *, "Change the amount of periods (0)"
        print *, "Analytic (1)"
        print *, "Euler (2)"
        print *, "Euler-Cromer (3)"
        print *, "---Euler-Cromer includes dA/dt analysis"
        print *, "---Euler-Cromer includes variable beta"
        print *, "Personal Method (4)"
        print *, "Precession (5)"
        print *, "End (any other value)"
        read *, functionreq

        if (functionreq .EQ. 0) then
            print*, "How many periods to run?"
            read *, periods
            print*, "Enter the time step (years):"
            read *, dt
            steps = bog(5) * periods / dt
            steps = nint(steps)
        else if (functionreq .EQ. 1) then
            call analytic(const)
        else if (functionreq .EQ. 2) then
            call euler(const, bog, dt, steps)
        else if (functionreq .EQ. 3) then
            call eulercromer(const, bog, dt, steps)
        else if (functionreq .EQ. 4) then
            call mymethod(const, bog, dt, steps)
        else if (functionreq .EQ. 5) then
            call precession(const, bog, dt, steps)
        else 
            close(41)
            CALL SYSTEM('gnuplot -p precoverbeta.plt')
            exit
        end if
    end do

end subroutine functo

!######################################################################### analytic solution

subroutine analytic(const)
    ! Preserve argument
    real, dimension(8) :: const

    ! Values for solution
    real :: x, y, stepsize, sunspot
    integer :: xstep, xstepfin, xn, yn

    stepsize = 0.25
    xstepfin = nint(const(5) / stepsize)

    ! Prepare the output files
    open(42,file="analyticsolution.txt")
    open(43,file="sunspot.txt")

    do yn=0,1
        do xn=0,1
            do xstep=0,xstepfin
                x = (-1)**xn * xstep * stepsize
                y = (-1)**yn * const(6) * sqrt(1 - (x / const(5))**2)
                write(42,*) x, y
            end do
        end do
    end do
    
    sunspot = const(5) * const(4)
    write(43,*) sunspot, 0.00000

    close(42)
    close(43)

    ! Call on gnu
    CALL SYSTEM('gnuplot -p analytic.plt')

    print *, 'Analysis complete beep boop'

end subroutine analytic

!######################################################################### euler approximation

subroutine euler(const, bog, dt, steps)
    ! Preserve arguments
    real, dimension(9) :: const
    real, dimension(5) :: bog
    real :: dt, steps

    ! Values for solution
    real, dimension(2) :: plutopos, sunpos, plutov, acceleration, r
    real, dimension(2) :: lastpos, lastv
    real :: rmag, vmag, amag, angle, t
    real :: potential, kinetic, etotal
    integer :: i

    ! Initial conditions
    plutopos(1) = const(5)
    plutopos(2) = 0
    plutov(1) = 0
    plutov(2) = bog(3)
    sunpos(1) = const(5) * const(4)
    r(1) = bog(1)
    r(2) = 0
    lastpos = plutopos
    lastv = plutov
    t = 0

    ! Prepare outputs
    open(42,file="euler.txt")
    open(43,file="sunspot.txt")
    write(42,*) plutopos(1), plutopos(2)

    do i=1,nint(steps)
        t = i * dt
        ! Find force, then divide by pluto's mass to get acceleration magnitude
        r(1) = plutopos(1) - sunpos(1)
        r(2) = plutopos(2)
        angle = atan(r(2)/r(1)) ! LOOKOUT - angle is in radians
        if (r(1) .LT. 0) then
            angle = const(3) + angle
        end if
        rmag = sqrt(r(1)**2 + r(2)**2)
        amag = const(2) * const(7) / (rmag**2)
        acceleration(1) = (-1) * amag * cos(angle)
        acceleration(2) = (-1) * amag * sin(angle)
        plutov(1) = lastv(1) + (acceleration(1) * dt)
        plutov(2) = lastv(2) + (acceleration(2) * dt)
        plutopos(1) = lastpos(1) + (lastv(1) * dt)
        plutopos(2) = lastpos(2) + (lastv(2) * dt)

        vmag = sqrt(plutov(1)**2 + plutov(2)**2)
        potential = (-1) * const(2) * const(7) * const(8) / rmag
        kinetic = 0.5 * const(8) * vmag**2
        etotal = kinetic + potential

        write(42,*) plutopos(1), plutopos(2), potential, kinetic, etotal, t
        lastpos = plutopos
        lastv = plutov

    end do

    sunspot = const(5) * const(4)
    write(43,*) sunspot, 0.00000

    close(42)
    close(43)

    ! Call on gnu
    CALL SYSTEM('gnuplot -p euler.plt')
    CALL SYSTEM('gnuplot -p eulerenergy.plt')

end subroutine euler

!######################################################################### euler-cromer approximation

subroutine eulercromer(const, bog, dt, steps)
    ! Preserve arguments
    real, dimension(9) :: const
    real, dimension(5) :: bog
    real :: dt, steps

    ! Values for solution
    real, dimension(2) :: plutopos, sunpos, plutov, acceleration, r
    real, dimension(2) :: lastpos, lastv, dr
    real :: rmag, vmag, amag, angle, t, area
    real :: potential, kinetic, etotal
    integer :: i, image

    ! Initial conditions
    plutopos(1) = const(5)
    plutopos(2) = 0
    plutov(1) = 0
    plutov(2) = bog(3)
    sunpos(1) = const(5) * const(4)
    r(1) = bog(1)
    r(2) = 0
    lastpos = plutopos
    lastv = plutov

    ! Prepare outputs
    open(42,file="eulercromer.txt")
    open(43,file="sunspot.txt")
    open(44,file="areastuffs.txt")
    write(42,*) plutopos(1), plutopos(2)

    do i=1,nint(steps)
        t = i * dt
        ! Find force, then divide by pluto's mass to get acceleration magnitude
        r(1) = plutopos(1) - sunpos(1)
        r(2) = plutopos(2)
        angle = atan(r(2)/r(1)) ! LOOKOUT - angle is in radians
        if (r(1) .LT. 0) then
            angle = const(3) + angle
        end if
        rmag = sqrt(r(1)**2 + r(2)**2)
        amag = const(2) * const(7) / (rmag**2)
        acceleration(1) = (-1) * amag * cos(angle)
        acceleration(2) = (-1) * amag * sin(angle)
        plutov(1) = lastv(1) + (acceleration(1) * dt)
        plutov(2) = lastv(2) + (acceleration(2) * dt)
        dr(1) = plutov(1) * dt
        dr(2) = plutov(2) * dt
        plutopos(1) = lastpos(1) + dr(1)
        plutopos(2) = lastpos(2) + dr(2)

        ! Energy stuff
        vmag = sqrt(plutov(1)**2 + plutov(2)**2)
        potential = (-1) * const(2) * const(7) * const(8) / rmag
        kinetic = 0.5 * const(8) * vmag**2
        etotal = kinetic + potential

        ! Area stuff
        area = 0.5 * ((r(1) * dr(2)) - (dr(1) * r(2)))

        ! Output and reset
        write(42,*) plutopos(1), plutopos(2), potential, kinetic, etotal, t
        write(44,*) t, area
        lastpos = plutopos
        lastv = plutov

    end do

    sunspot = const(5) * const(4)
    write(43,*) sunspot, 0.00000

    close(42)
    close(43)
    close(44)

    print *, "Print euler-cromer alone (1), with euler (2), or both (3)?"
    read *, image

    ! Call on gnu
    if (image .EQ. 1) then
        CALL SYSTEM('gnuplot -p eulercromeralone.plt')
    else if (image .EQ. 2) then
        CALL SYSTEM('gnuplot -p eulercromer.plt')
    else if (image .EQ. 2) then
        CALL SYSTEM('gnuplot -p eulercromer.plt')
        CALL SYSTEM('gnuplot -p eulercromeralone.plt')
    end if
    
    CALL SYSTEM('gnuplot -p eulercromerenergy.plt')
    CALL SYSTEM('gnuplot -p areastuffs.plt')

end subroutine eulercromer

!######################################################################### Personal method

subroutine mymethod(const, bog, dt, steps)
    ! Preserve arguments
    real, dimension(9) :: const
    real, dimension(5) :: bog
    real :: dt, steps

    ! Values for solution
    real, dimension(2) :: plutopos, sunpos, plutov, acceleration, r
    real, dimension(2) :: lastpos, lastv, dv
    real :: rmag, vmag, amag, angle, t
    real :: potential, kinetic, etotal
    integer :: i

    ! Initial conditions
    plutopos(1) = const(5)
    plutopos(2) = 0
    plutov(1) = 0
    plutov(2) = bog(3)
    sunpos(1) = const(5) * const(4)
    r(1) = bog(1)
    r(2) = 0
    lastpos = plutopos
    lastv = plutov

    ! Prepare outputs
    open(42,file="mymethod.txt")
    open(43,file="sunspot.txt")
    write(42,*) plutopos(1), plutopos(2)

    do i=1,nint(steps)
        t = i * dt
        ! Find force, then divide by pluto's mass to get acceleration magnitude
        r(1) = plutopos(1) - sunpos(1)
        r(2) = plutopos(2)
        angle = atan(r(2)/r(1)) ! LOOKOUT - angle is in radians
        if (r(1) .LT. 0) then
            angle = const(3) + angle
        end if
        rmag = sqrt(r(1)**2 + r(2)**2)
        amag = const(2) * const(7) / (rmag**2)
        acceleration(1) = (-1) * amag * cos(angle)
        acceleration(2) = (-1) * amag * sin(angle)
        dv(1) = acceleration(1) * dt
        dv(2) = acceleration(2) * dt
        plutov(1) = lastv(1) + dv(1)
        plutov(2) = lastv(2) + dv(2)
        plutopos(1) = lastpos(1) + (lastv(1) * dt) + (0.5 * dv(1) * dt)
        plutopos(2) = lastpos(2) + (lastv(2) * dt) + (0.5 * dv(2) * dt)

        vmag = sqrt(plutov(1)**2 + plutov(2)**2)
        potential = (-1) * const(2) * const(7) * const(8) / rmag
        kinetic = 0.5 * const(8) * vmag**2
        etotal = kinetic + potential

        write(42,*) plutopos(1), plutopos(2), potential, kinetic, etotal, t
        lastpos = plutopos
        lastv = plutov

    end do

    sunspot = const(5) * const(4)
    write(43,*) sunspot, 0.00000

    close(42)
    close(43)

    ! Call on gnu
    CALL SYSTEM('gnuplot -p mymethod.plt')
    CALL SYSTEM('gnuplot -p mymethodenergy.plt')

end subroutine mymethod

!######################################################################### Precession

subroutine precession(const, bog, dt, steps)
    ! Preserve arguments
    real, dimension(9) :: const
    real, dimension(5) :: bog
    real :: dt, steps

    ! Values for solution
    real, dimension(2) :: plutopos, sunpos, plutov, acceleration, r
    real, dimension(2) :: lastpos, lastv, dr
    real :: rmag, angle, t, beta, tlast, precangle, dbeta, depo
    real :: vmaglast, vmag, degree, degreelast, precspeed
    !real :: potential, kinetic, etotal
    integer :: i, vtracker, count

    print *, "Enter beta value: "
    read *, beta

    ! Initial conditions
    plutopos(1) = const(5)
    plutopos(2) = 0
    plutov(1) = 0
    !plutov(2) = bog(3)
    !sunpos(1) = 10
    plutov(2) = bog(4)
    sunpos(1) = -10
    r(1) = plutopos(1) - sunpos(1)
    r(2) = 0
    lastpos = plutopos
    lastv = plutov

    ! Prepare outputs
    open(42,file="precession.txt")
    open(43,file="sunspot.txt")
    open(44,file="precessioninfo.txt")
    write(42,*) plutopos(1), plutopos(2)

    vtracker = 1
    tlast = 0
    degreelast = 180
    count = 0
    depo = beta + 1.0

    do i=1,nint(steps)
        t = i * dt
        ! Find force, then divide by pluto's mass to get acceleration magnitude
        r(1) = plutopos(1) - sunpos(1)
        r(2) = plutopos(2)
        angle = atan(r(2)/r(1)) ! LOOKOUT - angle is in radians
        if (r(1) .LT. 0) then
            angle = const(3) + angle
        end if

        rmag = sqrt(r(1)**2 + r(2)**2)
        acceleration(1) = (-1) * const(2) * const(7) * r(1) / ((rmag**beta) * rmag) ! const(2) = G, const(7) = msun, and these don't cause any data issues
        acceleration(2) = (-1) * const(2) * const(7) * r(2) / ((rmag**beta) * rmag)
        plutov(1) = lastv(1) + (acceleration(1) * dt)
        plutov(2) = lastv(2) + (acceleration(2) * dt)
        dr(1) = plutov(1) * dt
        dr(2) = plutov(2) * dt
        plutopos(1) = lastpos(1) + dr(1)
        plutopos(2) = lastpos(2) + dr(2)

        if (t .GT. 5) then
            ! Precession calculations
            vmaglast = sqrt(lastv(1)**2 + lastv(2)**2)
            vmag = sqrt(plutov(1)**2 + plutov(2)**2)
            if (vtracker .EQ. 1) then
                if (vmag .GT. vmaglast) then

                    ! If at aphelion, get the degree value
                    degree = atan(r(2) / r(1)) * 180 / const(3)
                    if (r(1) .LT. 0) then
                        degree = 180 + degree
                    else if (r(2) .LT. 0) then
                        degree = 360 + degree
                    end if             

                    precangle = degree - degreelast
                    
                    ! Subtract over degree=0
                    if ((precangle) .LT. 0) then
                        degreelast = degreelast - 360
                        precangle = degree - degreelast
                    end if

                    ! Get the speed of precession
                    if ((t-tlast) .GT. 3) then
                        precspeed = (degree-degreelast) / (t - tlast)
                        if (count .GT. 0) then
                            write(44,*) t, degree, precangle, precspeed
                            print *, "The speed of precession is ", precspeed, " degrees per year."
                        end if
                        vtracker = 0
                        tlast = t
                        degreelast = degree
                        count = count + 1
                    end if

                end if
            end if
            if (vtracker .EQ. 0) then
                if (vmag .LT. vmaglast) then
                    if ((t-tlast) .GT. 5) then
                        vtracker = 1
                    end if
                end if
            end if
        end if

        ! Output and reset
        write(42,*) plutopos(1), plutopos(2), t
        lastpos = plutopos
        lastv = plutov

    end do

    sunspot = (-1) * const(5) * const(4)
    write(43,*) sunspot, 0.00000

    dbeta = beta - 2
    write(41,*) dbeta, precspeed


    close(42)
    close(43)
    close(44)

    ! Call on gnu
    CALL SYSTEM('gnuplot -p precession.plt')
    CALL SYSTEM('gnuplot -p paot.plt')

end subroutine precession