program homework1general

    implicit none

    ! inital variables
    real :: N_in, t_in, t_fin, steps, stepsize, grat, drat
    ! dt = t_next-t_now (it's the step size)
    ! steps = (t_fin-t_in)/stepsize
    ! drat is our death rate
    ! grat is our growth rate
    
    ! Things to solve for
    real :: N_now, N_next, t_now, c, N_nowex, Rt, Rtx
    ! Rt is the slope of our approximation
    ! Rt is the slope of the exact solution
    
    ! Other variables
    integer :: i
    
    ! Open documents to write on
    open(42,file='exact.txt')
    open(43,file='euler.txt')
    open(44,file='alldata.txt')
    
    ! Let's enter some values
    ! Created a subroutine so I can collapse the print and read lines
    call enters(N_in, t_in, t_fin, steps, stepsize, grat, drat)
    
    ! Initial datapoint
    t_now = t_in
    N_now = N_in
    N_nowex = N_in
    Rt = (grat*N_now)-(drat*N_now**2)
    Rtx = Rt
    
    ! Write the initial datapoint and find c
    write(42,*) '        t             N              dN/dt'
    write(42,*) t_now, N_now, Rtx
    write(43,*) '        t             N              dN/dt'
    write(43,*) t_now, N_now, Rt
    write(44,*) '  t                N_euler          N_exact          dN/dt_euler      dN/dt_exact'
    write(44,*) t_now, N_now, N_now, Rt, Rtx
    
    c = (1/grat)*log(1/((grat/N_now)-drat))
    
    ! Main function
    do i=1,steps
            
        t_now = t_in + (stepsize * i)
            
        ! Exact solution
        N_nowex = (grat*exp((grat*c)+(grat*t_now)))/(drat*exp((grat*c)+(grat*t_now))+1)
        Rtx = (grat*N_nowex)-(drat*N_nowex**2)
        
        ! Euler method
        N_next = N_now + (grat * N_now * stepsize) - (drat * N_now**2 * stepsize)
        N_now = N_next
        Rt = (grat*N_now)-(drat*N_now**2)
        write(43,*) t_now, N_now, Rt
        write(42,*) t_now, N_nowex, Rtx
        write(44,*) t_now, N_now, N_nowex, Rt, Rtx
        
        end do
    
    
    end program homework1general
    
subroutine enters(N1, t1, t2, s, ss, g, d)
    
    implicit none
    real :: N1, t1, t2, s, ss, g, d
    integer :: choice
    
        ! Just comment out unused values(this is a general template)
        print *, 'Enter initial population'
        read *, N1
        print *, 'Enter starting time [years]'
        read *, t1
        print *, 'Do we know the final time?'
        print *, "Respond with 1 for yes or 2 for no"
        read *, choice
        if (choice==1) then
            print *, 'Enter the final time [years]'
            read *, t2
        else if (choice==2) then
            print *, 'You must enter both steps and step size later.'
        else
            t2 = choice
            end if
            
        print*, 'Enter the growth rate [decimal]'
        read *, g
        print *, 'Enter the death rate [decimal]'
        read*, d
        print *, 'Do we know the step size or the amount of steps?'
        print *, "Respond with 1 for size, 2 for steps, or 3 for both."
        read *, choice
        if (choice==1) then
            print *, 'Enter the step size'
            read *, ss
            s = (t2-t1)/ss
        else if (choice==2) then
            print *, 'Enter the amount of steps'
            read *, s
            ss = (t2-t1)/s
        else if (choice==3) then
            print *, 'Enter the amount of steps'
            read *, s
            print *, 'Enter the step size'
            read *, ss
            end if

            
    end subroutine enters