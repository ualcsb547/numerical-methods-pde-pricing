function [t_final, y_final, s_optimo, iter] = midisp_nolineal(f, intv, alpha, beta, s_a, s_b, Tol, Nmax, N_pasos)
% Problemas no lineales mediante bisección
% y'' = f(t, y, y'), y(a) = alpha, y(b) = beta.

a = intv(1);
b = intv(2);

options = odeset('RelTol', 1e-4, 'AbsTol', 1e-5);
t_vec = linspace(a, b, N_pasos);

% F(s) = H(s) - beta
function err = evaluar_error(s_val)
    y0 = [alpha; s_val];
    [~, y] = ode45(f, t_vec, y0, options);
    err = y(end, 1) - beta; 
end

% Bisección
s0 = s_a;
s1 = s_b;

F0 = evaluar_error(s0);
F1 = evaluar_error(s1);

if F0 * F1 > 0
    error('F(s) no cambia de signo en el intervalo [s_a, s_b].');
end

iter = 0;

while iter < Nmax
    iter = iter + 1;
    s_c = (s0 + s1)/2;
    Fc = evaluar_error(s_c);

    if abs(Fc) <= Tol       
        s_optimo = s_c;      
        break;
    end

    if F0 * Fc < 0 
        s1 = s_c; 
        F1 = Fc;
    else    
        s0 = s_c; 
        F0 = Fc;
    end
end

if iter >= Nmax  
    warning('Se alcanzó Nmax sin llegar a la tolerancia requerida.');       
    s_optimo = s_c;   
end

% Resultados
y0_opt = [alpha; s_optimo];
[t_final, y_sol] = ode45(f, t_vec, y0_opt, options);

t_final = t_final';
y_final = y_sol(:, 1)'; 
end