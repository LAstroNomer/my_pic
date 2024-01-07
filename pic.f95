    program pic

    ! Подгрузка модулей
    ! Начальные условия
    use initial_donditions
    ! Модуль вывода данных в файл
    use my_io
    ! Первый шаг 
    use first_step
    ! Второй шаг 
    use second_step
    implicit none

    ! Необходимые переменные

    ! 3D
    real(mp)  M(0:n_steps-1, 0:n_steps-1, 0:n_particles_types-1)                   ! Расперделение
                                                                               ! массы           
    real(mp) e_arr(0:n_steps-1,0: n_steps-1, 0:n_particles_types-1)                ! Сетка внутренней энергии
    real(mp) w_arr(0:n_steps-1,0: n_steps-1)                ! Сетка полной энергии
    
    ! 2D
    real(mp) rho_arr(0:n_steps-1, 0:n_steps-1)             ! Сетка плотности 
    real(mp) M_all(0:n_steps-1, 0:n_steps-1)               ! Сетка общей массы
    real(mp) u_arr(0:n_steps-1,0: n_steps-1)               ! Сетка скорости вдоль оси х
    real(mp) v_arr(0:n_steps-1,0: n_steps-1)               ! Сетка скорости вдоль оси y
    real(mp) p_arr(0:n_steps-1,0: n_steps-1)               ! Сетка давления
    real(mp) phi_arr(0:n_steps-1,0: n_steps-1)             ! Сетка потенциальной энергии

    
    real(mp) x(1: n_particles, 0:n_particles_types-1)           ! Координата Х частиц
    real(mp) y(1: n_particles, 0:n_particles_types-1)           ! Координата Y частиц
    real(mp) m_p(1: n_particles, 0:n_particles_types-1)         ! Массы частиц
    
    ! 1D
    real(mp) h, tau, total_mass, tmp_rho, m_one, tmp_max, theta_x, theta_y
    real(mp) t, tmax, xk, yk, tmp
    integer i, j, k, l, f, alpha


    ! Зануляем все массивы
    rho_arr = 0d0
    u_arr = 0d0
    v_arr = 0d0
    e_arr = 0d0
    p_arr = 0d0
    phi_arr = 0d0
    w_arr = 0d0
    x     = -1d0
    y     = -1d0
    m_p   = -1d0
    M     = 0d0
    M_all = 0d0

    ! Шаг по сетке
    h = 1d0/n_steps
    write(*, *) 'h = ', h


    ! Распределение частиц по ячейкам
    ! В каждой ячейке n_init частиц
    ! В каждой ячейке частицы распределены равномерно
    ! alpha = 1 более тяжелые частицы
    call particles_destribution(x, y, m_p, M, h) ! from initial conditions

    ! Распределение на сетке
    M_all = M(:, :, 0) + M(:, :, 1)
    do i = 0, n_steps-1
        do j = 0, n_steps-1
            u_arr(i, j)   =  u0((i+0.5)*h, (j+0.5)*h) 
            v_arr(i, j)   =  v0((i+0.5)*h, (j+0.5)*h) 
            phi_arr(i, j) =  phi0((i+0.5)*h, (j+0.5)*h) 
            !if ((i+0.5)*h >= limit) then
            !    alpha = 1
            !else
            !    alpha = 0
            !endif
            ! Одно физ. в-во, удельная энергия одинакова
            e_arr(i, j, 0)   =  eps0((i+0.5)*h, (j+0.5)*h)
            e_arr(i, j, 1)   =  eps0((i+0.5)*h, (j+0.5)*h)
        enddo
    enddo

    tmax = 1
    t = 0
    do while (t <= tmax)
        write(*, *) '======STEP====='
        
        ! Вычисление давления
        p_arr      = equation_of_state(M, e_arr, h)
        tmp_max = sqrt(maxval(u_arr**2) + maxval(v_arr**2))

        tau = min(0.01d0*h, h/tmp_max)
        write(*, *) 'xy', minval(abs(x)), minval(abs(y)), maxval(x), maxval(y)
        write(*, *) 'tau=', tau
        t = t + tau

        M_all = (M(:, :, 0) + M(:, :, 1))
        w_arr = M_all *(phi_arr +  (u_arr**2 + v_arr**2)/2d0)/h**2 &
            + M(:, :, 0) * e_arr(:, :, 0)/h**2 + M(:, :, 1) * e_arr(:, :, 1)/h**2
        write(*,*) 'Energy', sum(w_arr)*h**2
        write(*, *) 'Mass', sum(M)
        
        ! Первый шаг
        call first(p_arr, e_arr, u_arr, v_arr, tau, h, M, m_p,  phi_arr)
        
        ! Второй шаг
        call second(x, y, m_p,  M, e_arr, u_arr, v_arr, tau, h, phi_arr)
        
        !open(1, file='xy1.out', position="append")
        !open(2, file='xy2.out',position="append")
        !do k = 1, n_particles
        !    alpha = 0
        !    if (m_p(k, alpha) < 0) then
        !        alpha = 1
        !    endif
        !    write(alpha+1, *) x(k, alpha), y(k, alpha)
        !enddo
        !write(1, *) '#'
        !write(2, *) '#'
        !call warr(p_arr, 'p.out')
        !call warr(u_arr, 'u.out')
        !call warr(v_arr, 'v.out')
        !call warr(e_arr(:, :, 0), 'e0.out')
        !call warr(e_arr(:, :, 1), 'e1.out')
        !call warr(w_arr, 'w.out')
        !call warr((M(:,:,0))/h**2, 'm.out')

   enddo
end
