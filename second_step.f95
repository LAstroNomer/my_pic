module second_step
    use initial_donditions
    use my_io
    implicit none

    contains

    subroutine second(x, y, m_p,  M, e, u, v, tau, h, phi)
    real(mp) x(1:n_particles, 0:n_particles_types-1)
    real(mp) y(1:n_particles, 0:n_particles_types-1) 
    real(mp) m_p(1:n_particles, 0:n_particles_types-1)
    real(mp) x_new(1:n_particles, 0:n_particles_types-1)
    real(mp) y_new(1:n_particles, 0:n_particles_types-1) 
    real(mp) M(0:n_steps-1, 0:n_steps-1, 0:n_particles_types-1)
    real(mp) dw(0:n_steps-1, 0:n_steps-1, 0:n_particles_types-1)
    real(mp) M_all(0:n_steps-1, 0:n_steps-1)
    real(mp) e(0:n_steps-1, 0:n_steps-1, 0:n_particles_types-1)
    real(mp) w(0:n_steps-1, 0:n_steps-1, 0:n_particles_types-1)
    real(mp) u(0:n_steps-1, 0:n_steps-1)
    real(mp) v(0:n_steps-1, 0:n_steps-1)
    real(mp) phi(0:n_steps-1, 0:n_steps-1)
    real(mp) du(0:n_steps-1, 0:n_steps-1)
    real(mp) dv(0:n_steps-1, 0:n_steps-1)
    real(mp) u_new(0:n_steps-1, 0:n_steps-1)
    real(mp) v_new(0:n_steps-1, 0:n_steps-1)
    real(mp) e_new(0:n_steps-1, 0:n_steps-1, 0:n_particles_types-1)
    real(mp) w_new(0:n_steps-1, 0:n_steps-1, 0:n_particles_types-1)
    real(mp) N(0:n_steps-1, 0:n_steps-1)
    real(mp) M_new(0:n_steps-1, 0:n_steps-1, 0:n_particles_types-1)
    real(mp) tau, h, ek, min_mass
    integer i, j, k, l, i_old, i_new, j_old, j_new, alpha


    write(*, *) '====SECOND===='
    call move(x, y, x_new, y_new, u, v, tau, h)

    M_new = 0d0
    N     = 0d0

    do k = 1, n_particles 
        alpha = 0

        if (x_new(k, alpha) < 0) then
            alpha = 1
        endif

        i = int(x_new(k, alpha)/h)
        j = int(y_new(k, alpha)/h)
        
        N(i, j) = N(i, j) +  1

        M_new(i, j, alpha) = M_new(i, j, alpha) + m_p(k, alpha) 
    enddo

    !open(1, file='xy.out')
    !do k =1, size(x)
    !    write(1, *) x(k), x_new(k),  y(k), y_new(k)
    !enddo
    !call warr(N, 'N.out')

    write(*,*) 'N sum', sum(N)
    write(*,*) 'old mass', sum(M)
    write(*,*) 'new mass', sum(M_new)


    du = 0d0
    dv = 0d0
    dw = 0d0

    do k = 1, n_particles

        alpha = 0
        if (m_p(k, alpha) < 0) then
            alpha = 1
        endif

        i_old = int(x(k, alpha)/h)
        j_old = int(y(k, alpha)/h)
                
        i_new = int(x_new(k, alpha)/h)
        j_new = int(y_new(k, alpha)/h)

        if ((abs(i_new - i_old) + abs(j_new - j_old)) == 0) then
            !write(*,*) 'Pass'
            continue
        endif


        du(i_old, j_old) = du(i_old, j_old) - m_p(k, alpha) * u(i_old, j_old) 
        du(i_new, j_new) = du(i_new, j_new) + m_p(k, alpha) * u(i_old, j_old) 

        dv(i_old, j_old) = dv(i_old, j_old) - m_p(k, alpha) * v(i_old, j_old) 
        dv(i_new, j_new) = dv(i_new, j_new) + m_p(k, alpha) * v(i_old, j_old) 
        
        
        ek = m_p(k, alpha) * (e(i_old, j_old, alpha) + phi(i_old, j_old) + u(i_old, j_old)**2/2d0 + v(i_old, j_old)**2/2d0)
        dw(i_old, j_old, alpha) =  dw(i_old, j_old, alpha) - ek/h**2
        dw(i_new, j_new, alpha) =  dw(i_new, j_new, alpha) + ek/h**2

    enddo

    write(*, *) 'summ du', sum(du)
    write(*, *) 'summ dv', sum(dv)

    M_all = M(:, :, 0) + M(:, :, 1)
    u_new = (M_all * u + du)/(M_new(:, :, 0) + M_new(:, :, 1))
    v_new = (M_all * v + dv)/(M_new(:, :, 0) + M_new(:, :, 1))
    !call borders(u_new, v_new, e_new)

    w(:, :, 0) = M(:, :, 0) * (e(:, :, 0) + phi +  (u**2 + v**2)/2d0)/h**2
    w(:, :, 1) = M(:, :, 1) * (e(:, :, 1) + phi + (u**2 + v**2)/2d0)/h**2

    write(*,*) 'Energy in', sum(w)*h**2
    
    
    w_new = w + dw   
    !call warr(dw(:, :, 1), 'dw.out')
    write(*, *) 'summ dw', (sum(w_new) -sum(w))*h**2
    min_mass = minval(abs(m_p)) 
    e_new = e
    do alpha = 0, 1
        do i = 0, n_steps-1
            do j = 0, n_steps-1
                if (M_new(i, j, alpha) < min_mass) then
                    e_new(i, j, alpha) = e(i, j, alpha)
                else
                    e_new(i, j, alpha) = h*h*w_new(i, j, alpha)/M_new(i, j, alpha)
                    e_new(i, j, alpha) = e_new(i, j, alpha) - (u_new(i, j)**2 + v_new(i, j)**2)/2d0 - phi(i, j)
                endif
            enddo
        enddo
    enddo


    u = u_new
    v = v_new
    e = e_new
    x = x_new
    y = y_new
    M = M_new
    write(*, *) 'summ e' , sum(u), sum(v), sum(M)
    
    w_new(:, :, 0) = M(:, :, 0) * (e(:, :, 0) + phi +  (u**2 + v**2)/2d0)/h**2
    w_new(:, :, 1) = M(:, :, 1) * (e(:, :, 1) + phi + (u**2 + v**2)/2d0)/h**2
    write(*, *) 'summ dw', (sum(w_new) -sum(w))*h**2

    end subroutine 

    !-------------------------------------------------------------------------
    subroutine move(x, y, x_new, y_new, u_arr, v_arr, tau, h)
    real(mp) x(1:n_particles, 0:n_particles_types-1), y(1:n_particles, 0:n_particles_types-1)
    real(mp) x_new(1:n_particles, 0:n_particles_types-1), y_new(1:n_particles, 0:n_particles_types-1)
    real(mp) u_arr(0:n_steps-1, 0:n_steps-1), v_arr(0:n_steps-1, 0:n_steps-1)
    real(mp) tmp(0:n_steps-1, 0:n_steps-1)
    real(mp) uk, vk, u00, u10, u01, u11, v00, v10, v01, v11, tau, h
    integer i, j, i_left, i_right, j_up, j_down,  k, alpha

    x_new = -1d0
    y_new = -1d0

    
    do k = 1, n_particles
        alpha = 0

        if (x(k, alpha) < 0) then
            alpha = 1
        endif    
        
        i = int(x(k, alpha)/h)
        j = int(y(k, alpha)/h)

        if (y(k, alpha) >= (j + 0.5d0) * h) then
            j_down  = j
            j_up = j+1
        else
            j_down  = j-1
            j_up = j
        endif

        if (x(k, alpha) >= (i + 0.5d0) * h) then
            i_left  = i
            i_right = i+1
        else
            i_left  = i-1
            i_right = i
        endif


        ! coord fix
        if (i_right > (n_steps-1)) then
            i_right = 0
        endif
        
        if (i_left < 0) then
            i_left = n_steps - 1
        endif

        if (j_down < 0) then
            j_down = 0
        endif

        if (j_up > (n_steps-1)) then
            j_up = n_steps - 1
        endif

        u00 = u_arr(i_left, j_down) 
        u10 = u_arr(i_right, j_down) 
        u01 = u_arr(i_left, j_up) 
        u11 = u_arr(i_right, j_up) 


        v00 = v_arr(i_left, j_down) 
        v10 = v_arr(i_right, j_down) 
        v01 = v_arr(i_left, j_up) 
        v11 = v_arr(i_right, j_up) 
        
        if (i_left == n_steps-1) then 
            if (x(k, alpha) < 0.5d0) then
                uk =  interp(u00, u10, u01, u11, x(k, alpha) + h/2d0, y(k, alpha) - j_down*h)
                vk =  interp(v00, v10, v01, v11, x(k, alpha) + h/2d0, y(k, alpha) - j_down*h)
            else
                uk =  interp(u00, u10, u01, u11, x(k, alpha) - i_left*h, y(k, alpha) - j_down*h)
                vk =  interp(v00, v10, v01, v11, x(k, alpha) - i_left*h, y(k, alpha) - j_down*h)
            endif
        else
            uk =  interp(u00, u10, u01, u11, x(k, alpha) - i_left*h, y(k, alpha) - j_down*h)
            vk =  interp(v00, v10, v01, v11, x(k, alpha) - i_left*h, y(k, alpha) - j_down*h)
        endif
        !if ((tau*uk > h) .or. (tau*vk > h)) then
        !    write(*,*) 'BAD TAU'
            !write(*,*) uk, u00, u01, u10, u11, i_left, i_right, j_down, j_up
            !write(*,*) vk, v00, v01, v10, v11, i_left, i_right, j_down, j_up
        !endif

        x_new(k, alpha) = x(k, alpha) + tau*uk
        y_new(k, alpha) = y(k, alpha) + tau*vk

        
        if (x_new(k, alpha) < 0) then 
            x_new(k, alpha)  = 1 + x_new(k, alpha) 
        endif

        if (y_new(k, alpha) < 0) then 
            y_new(k, alpha) = -y_new(k, alpha)
        endif

        if (x_new(k, alpha) > 1) then
            x_new(k, alpha) = x_new(k, alpha) - 1
        endif

        if (y_new(k, alpha) > 1) then
            y_new(k, alpha) = 2 - y_new(k, alpha)
        endif
    enddo
        
    end subroutine

    !-------------------------------------------------------------------------
    function interp(z0, z1, z2, z3, x, y)
        real(mp) interp, z0, z1, z2, z3, x, y
        real(mp) a00, a10, a01, a11

        if ((x > 1) .or. (y > 1)) then
            write(*, *) 'BAD INTERP'
        endif

        a00 = z0
        a10 = z1 - z0
        a01 = z2 - z0
        a11 = z3 - z2 - z1 + z0

        interp = a00 + a10 * x + a01 *y + a11*x*y

    end function 
end module
