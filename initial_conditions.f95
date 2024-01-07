    module initial_donditions
    implicit none

    ! float precision
    integer, parameter :: mp = 8
    
    ! particels parameters
    integer, parameter :: n_init = 9
    integer, parameter :: n_particles_types = 2

    ! grid sizes
    integer, parameter :: n_steps = 100
    integer, parameter :: n_particles = n_init * n_steps * n_steps

    real(mp), parameter :: limit = 0.5d0
    real(mp), parameter :: g = 10.0d0
    contains

    function phi0(x, y)
        real(mp) x, y, phi0
        phi0 = g * y
    end function

    function rho0(x, y)
        real(mp) x, y, rho0, r

        if (y >= limit) then
            rho0 = 1.2920d0
        else
            rho0 = 1.1455d0
        endif
    end function rho0

    function u0(x, y)
        real(mp) x, y, u0
        u0 = 0.0d0 
    end function u0

    function v0(x, y)
        real(mp) x, y, v0
        v0 = 0d0 
    end function v0
    
    function eps0(x, y)
        real(mp) x, y, eps0, eps1, eps2
        eps2 = 156457.0d0 !T = 0 C
        eps1 = 176515.0d0 !T = +35 C

        eps0 = eps1 + y*(eps2 - eps1)
    end function eps0

    !subroutine borders(u, v, e)
    !    real(mp) u(0:n_steps-1, 0:n_steps-1), v(0:n_steps-1,0:n_steps-1), e(0:n_steps-1, 0:n_steps-1, 0:n_particles_types)
    !    real(mp) eps1, eps2
    !    eps2 = 156457.0d0 !T = 0 C
    !    eps1 = eps2 !176515.0d0 !T = +35 C

    !    v(:, 0) = 0
    !    v(:, n_steps-1) = 0

    !    e(:, 0, :) = eps1
    !    e(:, n_steps-1, :) = eps2

    !end subroutine

    function equation_of_state(M, eps, h)
        real(mp) M(0:n_steps-1, 0:n_steps-1, 0:n_particles_types-1)
        real(mp) eps(0:n_steps-1, 0:n_steps-1, 0:n_particles_types-1)
        real(mp) val(0:n_steps-1, 0:n_steps-1)
        real(mp) equation_of_state(0:n_steps-1, 0:n_steps-1)
        real(mp) h, coef

        val = (M(:, :, 0)* eps(:,:,0) + M(:, :, 1)* eps(:,:,1))/h**2
        coef = 4d0/3d0
        equation_of_state = (coef - 1) * val
    end 


    subroutine particles_destribution(x, y, m_p, M, h)
        real(mp) M(0:n_steps-1, 0:n_steps-1, 0:n_particles_types-1)   
        real(mp) x(1: n_particles, 0:n_particles_types-1)           
        real(mp) y(1: n_particles, 0:n_particles_types-1)           
        real(mp) m_p(1: n_particles, 0:n_particles_types-1)         
        real(mp) hx(9), hy(9)
        
        real(mp) xk, yk, m_one, h, theta_x, theta_y
        integer k, l, i, j, alpha 

        hx(1:3) = [ -h/4d0, 0d0, h/4d0 ]
        hx(4:6) = [ -h/4d0, 0d0, h/4d0 ]
        hx(7:9) = [ -h/4d0, 0d0, h/4d0 ]


        hy(1:3) = [ -h/4d0, -h/4d0, -h/4d0 ]
        hy(4:6) = [ 0d0, 0d0, 0d0 ]
        hy(7:9) = [ h/4d0, h/4d0, h/4d0 ]
        ! Распределение частиц по ячейкам
        ! В каждой ячейке n_init частиц
        ! alpha = 1 более тяжелые частицы
        k = 1
        l = 0
        do i = 0, n_steps-1
            do j = 0, n_steps-1
                do while (k <= size(hx)) 

                    !call random_number(theta_x) ! Рандомное число (0, 1)
                    !theta_x = theta_x - 0.5d0
                    xk = (i + 0.5d0)*h +hx(k) !+ theta_x*h
                    !call random_number(theta_y)
                    !theta_y = theta_y - 0.5d0
                    yk = (j + 0.5d0)*h + hy(k) !+ theta_y*h

                    if (yk >= limit) then
                        alpha = 1
                    else 
                        alpha = 0
                    endif

                    m_one = rho0(xk, yk)*h**2/n_init 

                    x(k + l*n_init, alpha) = xk 
                    y(k + l*n_init, alpha) = yk
                    m_p(k + l*n_init, alpha) = m_one
                    M(i, j, alpha) = M(i, j, alpha) + m_one
                    k = k + 1
                enddo
                k = 1
                l = l + 1 
            enddo
        enddo

        !write(*, *) maxval(x), maxval(y)
        !write(*, *) minval(abs(x)), minval(abs(y))

        !do k = 1, n_particles
        !    if ((x(k, 0) < 0) .and. (x(k,1) < 0)) then
        !        write(*,*) k
        !    endif
        !enddo
    end subroutine


end module
