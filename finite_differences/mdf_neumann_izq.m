function [t, y] = mdf_neumann_izq(p, q, r, intervalo, y0, N)

% Neumann a la izquierda:

a = intervalo(1);
b = intervalo(2);
alpha = y0(1); % y'(a) = alpha
beta = y0(2);  % y(b) = beta

% Comprobar si N es adecuada para ver si es estrictamente diagonal
% dominante
j = 0.01; 
n0 = floor((b-a)/j);
t0 = a + j*(0:n0);
cota = (b-a) * max(abs(p(t0))) / 2;

if N > cota
    h = (b - a) / (N + 1);
    t = linspace(a, b, N + 2)';
    
    % p, q, r en los nodos interiores
    ti = t(2:end-1);
    pt = p(ti);
    qt = q(ti);
    rt = r(ti);
    
    % A 
    A = zeros(N + 1, N + 1);
    
    % x0 (Neumann)
    A(1, 1) = 1 - q(a) * (h^2 / 2);
    A(1, 2) = -1;
    
    % x1 - xN
    diag_abajo = 1 - (h/2) * pt(2:end);
    diag_principal = -2 + h^2 * qt;
    diag_arriba = 1 + (h/2) * pt(1:end-1);
    
    A(2:end, 2:end) = diag(diag_principal) + diag(diag_arriba, 1) + diag(diag_abajo, -1);
    
    % x1 con y0
    A(2, 1) = 1 - (h/2) * pt(1);
    
    % Términos independientes
    B = zeros(N + 1, 1);
    
    % x0
    B(1) = -alpha * h + p(a) * alpha * (h^2 / 2) - r(a) * (h^2 / 2);
    
    % x1 - xN
    B(2:end) = h^2 * rt;
    
    % Dirichlet a la derecha
    B(end) = B(end) - (1 + (h/2) * pt(end)) * beta;
    
    % Resolución:
    yi = A \ B;
    
    % Con extremos
    y = [yi; beta];
else
    % Si N no supera la cota, devolvemos vectores vacíos
    t = []; 
    y = []; 
    disp('N es demasiado pequeño, tome uno adecuado');
end
end