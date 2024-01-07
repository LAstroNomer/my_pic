module particles
    use initial_conditions
    implicit none

    subroutine get_particles_destribution(X, Y, N, )
        real(mp) particles(0:n_particles-1, 3), h_tmp
        integer nmax

        nmax = int(sqrt(float(n_particles))
        h_tmp = 
        do i =0, nmax - 1
            do j = 0, nmax - 1
                
                particles(i*100 + j, 1) =  j*h_tmp
                particles(i*100 + j, 2) =  i*h_tmp

            enddo
        enddo 

    end subroutine

    !-----------------------------------------------------------

    subroutine count_part_in_cell(X, Y, N)
        real(mp) X(:), Y(:), N(:, :)
        integer i, j
        


    end subroutine

end module
