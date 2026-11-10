function [t, y] = mdf_neumann_der(p, q, r, intervalo, y0, N)

% Neumann a la derecha:

a = intervalo(1);
b = intervalo(2);
alpha = y0(1); % y(a) = alpha
beta = y0(2);  % y'(b) = beta

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

    % Coeficientes de las yi
    A = zeros(N + 1, N + 1);

    % xN (Neumann)
    A(N+1, N+1) = 1 - q(b) * (h^2 / 2);
    A(N+1, N) = -1;

    % Nodos interiores
    diag_abajo = 1 - (h/2) * pt(2:end);
    diag_principal = -2 + h^2 * qt;
    diag_arriba = 1 + (h/2) * pt(1:end-1);

    A(1:N, 1:N) = diag(diag_principal) + diag(diag_arriba, 1) + diag(diag_abajo, -1);

    % xN con y{N+1}
    A(N, N+1) = 1 + (h/2) * pt(end);

    % Términos independientes
    B = zeros(N + 1, 1);

    % x{N+1}
    B(end) = beta * h + p(b) * beta * (h^2 / 2) - r(b) * (h^2 / 2);

    % x1 - xN
    B(1:N) = h^2 * rt;

    % Dirichlet a la izq 
    B(1) = B(1) - (1 + (h/2) * pt(1)) * alpha;

    % Resolucíon:
    yi = A \ B;

    y = [alpha; yi];
else
    % Si N no supera la cota, devolvemos vectores vacíos
    t = []; 
    y = []; 
    disp('N es demasiado pequeño, tome uno adecuado');
end
end