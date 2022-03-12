program homework3_b

    ! Define Variables
    real :: mu, mint
    real :: xcnow, xclast
    real :: xinow, xilast
    real :: webxo, webx
    integer :: i, reach

    ! Open files for output
    open(42,file="curveout.txt")
    open(43,file="idout.txt")
    open(44,file="web.txt")

    ! Initialize some variables
    print *, 'mu value:'
    read *, mu

    i = 0

    ! Curve and id portion
    do
        xclast = i * 0.001
        xilast = i * 0.001
        if (xclast .GT. 1) then
            exit
        end if
        xcnow = (mu * xclast) * (1 - xclast)
        xinow = xilast
        write(42,*) xclast, xcnow
        write(43,*) xilast, xinow
        i = i + 1
    end do

    ! Cobweb portion
    ! ID starting point
    webxo = 0.2
    write(44,*) webxo, webxo
    mint = mu * 1.5
    reach = nint(exp(mint))

    do i=0,reach
        ! Curve, then ID, then repeat
        webx = (mu * webxo) * (1 - webxo)
        write(44,*) webxo, webx
        write(44,*) webx, webx
        webxo = webx
    end do


    ! Call on gnu
    CALL SYSTEM('gnuplot -p plotcodeb.plt')

end program homework3_b