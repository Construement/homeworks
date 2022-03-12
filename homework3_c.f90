program homework3_c
    
    real :: mu
    real :: x
    real :: f1, f2, f3, fid
    real :: f3_123, f3_45, f3_678
    integer :: i

    print *, 'Enter mu:'
    read *, mu

    open(42,file="allfoutput.txt")

    do i=0,100
        x = i * .01
        ! f functions
        f1 = (mu * x) - (mu * x**2)
        f2 = (x * mu**2) - (x**2 * (mu**3 + mu**2)) + (2 * x**3 * mu**3) - (x**4 * mu**3)
        ! Split f3 into 3 parts for clarity, then just added them back together
        f3_123 = (x * mu**3) - (x**2 * (mu**3 + mu**4 + mu**5)) + (2 * x**3 * (mu**4 + mu**5 + mu**6))
        f3_45 = (x**5 * ((4 * mu**7) + (6 * mu**6))) - (x**4 * (mu**4 + mu**5 + (6 * mu**6) + mu**7))
        f3_678 = (4 * x**7 * mu**7) - (x**6 * ((2 * mu**6) + (6 * mu**7))) - (x**8 * mu**7)
        f3 = f3_123 + f3_45 + f3_678
        ! ID function
        fid = x
        write(42,*) x, fid, f1, f2, f3
    end do

    ! Call on gnu
    CALL SYSTEM('gnuplot -p plotcodec.plt')

end program homework3_c