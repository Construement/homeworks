program hw2better
    
    ! Define Variables
    real :: t, vx, vy, vxo, vyo, x, y, xo, yo
    real :: fdrag, fxd, fyd
    real :: m, C, p, A, vo, t_in, B, G, r, pi, wind
    real :: angledeg, anglerad
    real :: dt, dvx, dvy, i
    real :: xex, yex, v_in, vx_in, vy_in, x_in, y_in
    real :: vtotal, ymax, ycheck, pchoice
    
    ! Create the filename variables
    !character(20) :: filename
    
    ! Enter initial values
    x_in=0
    y_in=0
    t_in=0
    m=140
    G=9.81
    C=0.15
    r=0.09
    print *, 'Enter initial velocity [m/s]'
    read *, v_in
    print *, 'Enter wind speed (negative if against flight) [m/s]'
    read *, wind
    print *, 'Is air density constant? (1 for yes)'
    read *, pchoice
    
    pi = 3.14159265
    dt = 0.1
    A = pi * r**2
    
    ! Opening files to write to
    open(42,file='posdata.txt')
    write(42,*) '     angle              x               y             x_exact           y_exact            vx               vy'
    open(43,file='shpeed.txt')
    write(43,*) '     angle              t               v'
    open(44,file='maxes.txt')
    write(44,*) '     angle            range            y_max          tflight'
    
    ! Let's create the loop over all angles
    do angledeg=25,50
        
        ! Code-exclusive variables
        i = 1
        vo = v_in
        xex=x_in
        yex=y_in
        
        anglerad = pi * angledeg/180.0
        
        !print *, angledeg
        !print *, anglerad
        
        ! Initial x, y, and directional v and f
        vxo = vo * cos(anglerad)
        vyo = vo * sin(anglerad)
        !fdrag = B * vo**2
        !fxd = -fdrag * (vxo + wind)
        !fyd = -fdrag * vyo
        vx_in = vxo
        vy_in = vyo
        xo = x_in
        yo = y_in
        xex = x_in 
        yex = y_in
        
        write(42,*) angledeg, x_in, y_in, xex, yex, vxo, vyo
        write(43,*) angledeg, t_in, vo
        
        ! Main function
        ! do is a function in fortran. If I meant subroutine, I would say subroutine.
        ! The main function of this program is to simulate the projectile's motion.
        
        do
            t = t_in + (i * dt)
            
            ! In case of changing air density
            p = 1.29
            if (pchoice .NE. 1) then
                p = 1.29 * (1 - ((0.0065 * yo)/280.0))**2.5
                end if
            B = 0.5 * C * p * A
            
            ! Find drag from initial v, then find dv to get the new v
            fdrag = B * sqrt((vxo + wind)**2 + vyo**2)
            fxd = -fdrag * (vxo + wind)
            fyd = -fdrag * vyo
            dvx = fxd/m * dt
            dvy = (fyd/m - G) * dt
            vx = vxo + dvx
            vy = vyo + dvy
            vtotal=sqrt((vxo + wind)**2 + vyo**2)
            
            ! Basic Kinematic equation: x=xo+vo*t+.5*a*t^2 where a is constant dv/dt
            ! but dv/dt * dt^2 = dv * dt
            x = xo + (vxo * dt) + (0.5 * dvx * dt)
            y = yo + (vyo * dt) + (0.5 * dvt * dt)
            ycheck = y-yo
            if (ycheck .GT. 0) then
                ymax = y
                end if
            
            ! Exact solution
            xex = x_in + (vx_in * t)
            yex = y_in + (vy_in * t) - (0.5 * G * t**2)
            
            ! Correct any negative y_exact values
            if (yex .LE. 0) then
                yex = 0
                end if
            
            ! End the loop when the projectile lands
            if (y .LE. 0) then
                write(42,*) angledeg, x, y, xex, yex, vxo, vyo
                write(43,*) angledeg, t, vtotal
                write(44,*) angledeg, x, ymax, t
                exit
                end if
            
            write(42,*) angledeg, x, y, xex, yex, vxo, vyo
            write(43,*) angledeg, t, vtotal
            
            ! Reset values for the loop
            xo = x
            yo = y
            vxo = vx
            vyo = vy
            i = i + 1
            end do
    
        end do
    
    end program hw2better