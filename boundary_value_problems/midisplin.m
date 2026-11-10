function [t, x] = midisplin(metodo, p, q, r, C1, C2, interval, x0, N) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% INPUT: p, q, r funciones de la EDO 
% C1=0 Dirichlet en t_0 (x(t0) = a)
% C1=1 Neumann en t_0   (x'(t0) = a)
% C2=0 Dirichlet en t_N (x(T) = b)
% C2=1 Neumann en t_N   (x'(T) = b)
% metodo = 'euler', 'rk4' u 'ode45'

funcdisplin1 = @(t,y) [y(2); p(t).*y(2) + q(t).*y(1) + r(t)];
funcdisplin2 = @(t,y) [y(2); p(t).*y(2) + q(t).*y(1)];

a = x0(1); 
b = x0(2); 

% Condición inicial (C1)
if C1 == 0 % Dirichlet en t0
    u0 = [a; 0];
    v0 = [0; 1];
elseif C1 == 1 % Neumann en t0
    u0 = [0; a];
    v0 = [1; 0];
else
    error('C1 debe ser 0 o 1');
end

% Dependiendo del método
if strcmp(metodo, 'euler')
    [t, u] = mieuler(funcdisplin1, interval, u0, N);
    [t, v] = mieuler(funcdisplin2, interval, v0, N);

elseif strcmp(metodo, 'rk4')
    [t, u] = mirk4(funcdisplin1, interval, u0, N);
    [t, v] = mirk4(funcdisplin2, interval, v0, N);

elseif strcmp(metodo, 'ode45')
    options = odeset('RelTol', 1e-4, 'AbsTol', 1e-5);
    t_vec = linspace(interval(1), interval(2), N);

    [t, u] = ode45(funcdisplin1, t_vec, u0, options);
    [t, v] = ode45(funcdisplin2, t_vec, v0, options);

    t = t'; 
    u = u'; 
    v = v';
else
    error('Método no reconocido. Usa ''euler'', ''rk4'' u ''ode45''.');
end

% Valores finales (C2)
if C2 == 0 % Dirichlet en T
    u_T = u(1, end);
    v_T = v(1, end);
elseif C2 == 1 % Neumann en T
    u_T = u(2, end);
    v_T = v(2, end);
else
    error('C2 debe ser 0 o 1');
end

% Cálculo de s y de x
if v_T ~= 0 
    s = (b - u_T) / v_T;
    x = u + s .* v;
    disp('La solucion es única y viene dada por:'); 
    x(1, :)
else 
    if abs(u_T - b) < 1e-10 
        x = u;
        disp('Hay infinitas soluciones.')
    else 
        x = []; 
        disp('El problema no tiene solución'); 
    end 
end