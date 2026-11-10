function [t, y] = mdf_dirichlet(p, q, r, intervalo, y0, N)

a = intervalo(1);
b = intervalo(2);
alpha = y0(1);
beta = y0(2);

% Comprobar si N es adecuada para ver si es estrictamente diagonal
% dominante:
j = 0.01; 
n0 = floor((b-a)/j);
t0 = a + j*(0:n0);
cota = (b-a) * max(abs(p(t0))) / 2; % Sacado de mis apuntes de grado

if N > cota
    
    h = (b - a) / (N + 1);
    t = linspace(a, b, N + 2)'; % Vector columna de nodos (incluyendo fronteras)
    
    % Nodos interiores
    ti = t(2:end-1);
    
    % p, q, r en los nodos interiores
    pt = p(ti);
    qt = q(ti);
    rt = r(ti);
    
    % Matriz tridiagonal
    diag_abajo = 1 - (h/2) * pt(2:end);
    diag_principal = -2 + h^2 * qt;
    diag_arriba = 1 + (h/2) * pt(1:end-1);
    
    A = diag(diag_principal) + diag(diag_arriba, 1) + diag(diag_abajo, -1);
    
    % B
    B = h^2 * rt;
    B(1) = B(1) - (1 - (h/2) * pt(1)) * alpha;
    B(end) = B(end) - (1 + (h/2) * pt(end)) * beta;
    
    % Sistema (nodos interiores)
    yi = A \ B;

    y = [alpha; yi; beta];
else
    % Si N no supera la cota, devolvemos vectores vacíos
    t = []; 
    y = []; 
    disp('N es demasiado pequeño, tome uno adecuado');
end
end