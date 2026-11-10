clc
clear all
close all
format long

a=0; b=1; c=0; T=1; 

M = 1000;
dt = (T-c)/(M+1); % partición temporal
t1 = c : dt : T; 


N = 20;
h = (b-a)/(N+1); % partición espacial
x1 = a : h : b;

k = 1;
gamma = dt*k/(h.^2);

% Como ahora tenemos estabilidad siempre, no hacemos warning

[t, x] = meshgrid(t1, x1); % cuadrícula
z = zeros(length(x1), length(t1)); 

%  u(x,0)
z(:, 1) = 4*x1' - 4.*(x1'.^2); 

% Condiciones de frontera
z(1, :) = 0;   
z(end, :) = 0; 

% Construcción de matrices A (futuro) y B (presente)

N = length(x1) - 2; % Nodos espaciales interiores

% Matriz A 
diag_princ_A = (1 + gamma) * ones(N, 1);
diag_sec_A = (-gamma/2) * ones(N - 1, 1);
A = diag(diag_princ_A) + diag(diag_sec_A, 1) + diag(diag_sec_A, -1);

% Matriz B
diag_princ_B = (1 - gamma) * ones(N, 1);
diag_sec_B = (gamma/2) * ones(N - 1, 1);
B = diag(diag_princ_B) + diag(diag_sec_B, 1) + diag(diag_sec_B, -1);

for i = 2:length(t1)
    % Extraemos los nodos interiores de i-1
    U_presente = z(2:end-1, i-1);

    % Resolución
    U_futuro = A \ (B * U_presente);

    % Guardamos el resultado para i
    z(2:end-1, i) = U_futuro;
end

% i) Solución en 2D para T = 1
figure
plot(x1, z(:, end), 'r', 'LineWidth', 2)
title('Solución en T=1 (Crank-Nicholson)')
xlabel('x'); ylabel('u(x, 1)')

% ii) Solución en 3D
figure
surf(t, x, z)
shading interp 
title('Evolución temporal del calor (Crank-Nicholson)')
xlabel('t'); ylabel('x'); zlabel('u(x,t)')