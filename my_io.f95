module my_io
    implicit none


    contains
    
    subroutine warr(a, ofile)
        use initial_donditions
        real(mp)  a(0:n_steps-1, 0:n_steps-1)
        character(*) ofile
        integer i
        open(1, file=ofile, position="append")
        do i = 0, n_steps-1
            write(1,*) a(i,:)
        enddo
        write(1,*) '#' 
        !close(1)
    end 
    
    subroutine warri(a, ofile)
        use initial_donditions
        integer a(0:n_steps-1, 0:n_steps-1)
        character(*) ofile
        integer i
        open(1, file=ofile)
        do i = 0, n_steps-1
            write(1,*) a(i,:)
        enddo
        close(1)
    end 
end module
