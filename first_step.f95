module first_step
    use initial_donditions
    use my_io
    implicit none

    contains

    subroutine first(p, e, u, v, tau, h, M, m_p,  phi)

        real(mp) M(0:n_steps-1, 0:n_steps-1, 0:n_particles_types-1)
        real(mp) e(0:n_steps-1, 0:n_steps-1, 0:n_particles_types-1)
        real(mp) w(0:n_steps-1, 0:n_steps-1)
        real(mp) m_p(1:n_particles, 0:n_particles_types-1)
        real(mp) phi(0:n_steps-1, 0:n_steps-1)
        real(mp) p(0:n_steps-1, 0:n_steps-1)
        real(mp) u(0:n_steps-1,0:n_steps-1)
        real(mp) v(0:n_steps-1,0:n_steps-1)
        real(mp) dE(0:n_steps-1,0:n_steps-1)
        real(mp) du(0:n_steps-1, 0:n_steps-1)
        real(mp) dv(0:n_steps-1,0:n_steps-1)
        real(mp) dw(0:n_steps-1, 0:n_steps-1)
        real(mp) rho_all(0:n_steps-1, 0:n_steps-1)
        real(mp) u_new(0:n_steps-1, 0:n_steps-1)
        real(mp) v_new(0:n_steps-1, 0:n_steps-1)
        real(mp) w_new(0:n_steps-1, 0:n_steps-1)
        real(mp) e_new(0:n_steps-1, 0:n_steps-1, 0:n_particles_types-1)
        real(mp) A, B, tau, h, pp, pm, hp, hm, min_mass
        integer i, j

        ! begin
        write(*, *) '==========FIRST========'

        rho_all = (M(:, :, 0) + M(:, :, 1))/h**2
            
        du = 0d0
        dv = 0d0

        do i = 0, n_steps-1
            do j = 1, n_steps-2
                ! x
                if (i == 0) then
                    pm = (p(n_steps-1, j) + p(0, j))/2d0
                else
                    pm = (p(i-1, j) + p(i, j))/2d0  
                endif
                
                if (i == n_steps-1) then
                    pp = (p(n_steps-1, j) + p(0, j))/2d0
                else
                    pp = (p(i+1, j) + p(i, j))/2d0 
                endif 
                
                du(i, j) = (pp - pm)/rho_all(i, j)
                

                ! y 
                pm = (p(i, j-1) + p(i, j))/2d0  
                
                
                pp = (p(i, j+1) + p(i, j))/2d0 
                
                dv(i, j) = (pp - pm)/rho_all(i, j) + g*h
            enddo
        enddo

        dv(:, 0) = 0
        dv(:, n_steps-1) = 0

        write(*, *) 'summ du', sum(du)
        write(*, *) 'summ dv', sum(du)
        u_new = u - tau/h * du
        v_new = v - tau/h * dv
        
        !call borders(u_new, v_new, e_new)

        ! calculate w_i_j
        w = rho_all * ((u**2 + v**2)/2d0 + phi) + M(:, :, 0) * e(:, :, 0)/h**2 + M(:, :, 1) * e(:, :, 1)/h**2
        write(*,*) 'Energy in', sum(w)*h**2
        dw = 0d0

        do i = 0, n_steps-1
            do j = 1, n_steps-2
                
                ! x
                if (i == n_steps -1) then
                    pp = (p(n_steps-1, j) + p(0, j))/2d0
                    A = pp *(u_new(n_steps-1, j) + u(n_steps-1, j) + u_new(0, j) + u(0, j))/4d0 
                else
                    pp = (p(i+1, j) + p(i, j))/2d0 
                    A = pp *(u_new(i, j) + u(i, j) + u_new(i+1, j) + u(i+1, j))/4d0 
                endif
                
                if (i == 0) then
                    pm = (p(n_steps-1, j) + p(0, j))/2d0
                    A = A - pm * (u_new(0, j) + u(0, j) + u_new(n_steps-1, j) + u(n_steps-1, j))/4d0    
                else
                    pm = (p(i-1, j) + p(i, j))/2d0  
                    A = A - pm * (u_new(i, j) + u(i, j) + u_new(i-1, j) + u(i-1, j))/4d0
                endif
                

                pm = (p(i, j-1) + p(i, j))/2d0  
                pp = (p(i, j+1) + p(i, j))/2d0 

                B = pp * (v_new(i, j) + v(i, j) + v_new(i, j+1) + v(i, j+1))/4d0 
                B = B - pm * (v_new(i, j) + v(i, j) + v_new(i, j-1) + v(i, j-1))/4d0    

                dw(i, j) = (A + B)
            enddo
        enddo

        w_new = w - tau/h*dw
        write(*,*) 'Energy diff', (sum(w_new) - sum(w))*h**2

        e_new(:, :, 0) = w_new/rho_all - (phi + (u_new**2 + v_new**2)/2d0) 
        dE = 0d0
        dE = ((M(:, :, 0) + M(:, :, 1)) * e_new(:, :, 0) - M(:, :, 0) * e(:, :, 0) - M(:, :, 1) * e(:, :, 1))/(M(:,:,0) + M(:, :,1))
        !call warr(dE, 'de.out')

        e_new = e
        min_mass = minval(abs(m_p))
        do i = 0, n_steps-1
            do j = 1, n_steps-2
                
                !if ((dE(i, j) + e_new(i, j, 0)) > 0) then
                !    e_new(i, j, 0) = e_new(i, j, 0) + dE(i, j)
                !endif

                !if ((dE(i, j) + e_new(i, j, 1)) > 0) then
                !    e_new(i, j, 1) = e_new(i, j, 1) + dE(i, j)
                !endif
                if (M(i, j, 0) < min_mass) then 
                    e_new(i, j, 1) = e(i, j, 1) + dE(i, j)
                else 
                    if (M(i, j, 1) < min_mass) then 
                        e_new(i, j, 0) = e(i, j, 0) + dE(i, j)
                    else
                        e_new(i, j, 1) = e(i, j, 1) + dE(i, j)
                        e_new(i, j, 0) = e(i, j, 0) + dE(i, j)
                    endif
                endif
            enddo
        enddo

        e = e_new

        w_new = rho_all * ((u_new**2 + v_new**2)/2d0 + phi) + M(:, :, 0) * e(:, :, 0)/h**2 + M(:, :, 1) * e(:, :, 1)/h**2

        write(*,*) 'Energy diff', (sum(w_new) - sum(w) ) *h**2
        call warr( (w - w_new) *h**2, 'de.out')


        u = u_new
        v = v_new 
        
    end subroutine

end module
