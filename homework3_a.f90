program homework3_a
    
    ! Define Variables
    real :: mu, xnow, xlast
    real :: stepmu, stepsizemu
    real :: stepx, initstep, useful
    integer :: i, imu, ix, filevalue
    
    ! User Input
    
    
    ! Initialize Others
    
    initstep = 0.2
    stepsizemu = .00001
    stepmu = 4.0/stepsizemu
    stepx = 200
    mu = 1.5
    
    ! Open files
    open (42,file="basicoutput.txt")
    open (43,file="xmoutput1.txt")
    open (44,file='xmoutput2.txt')
    open (45,file='xmoutput3.txt')
    open (46,file='xmoutput4.txt')
    
    do ix=1,4
        
        xinit = initstep * ix
        print *, xinit
        filevalue = 42 + ix
        print *, filevalue
        
        do imu=1,stepmu
            mu = stepsizemu * imu
            xlast = xinit
            useful = 200 - ((0.8) * exp(mu))
        
            ! Loop for x stability
            do i=0,stepx
                
                xnow = mu * xlast * (1 - xlast)
                !write (42,*) mu, i , xnow
                
                if (i .GT. useful) then
                    write (filevalue,*) mu, xnow
                    end if
                    
                xlast = xnow
                    
                end do
        
            end do
        
        end do
    
    ! Call on gnu
    CALL SYSTEM('gnuplot -p plotcodea.plt')
    
    
    end program homework3_a