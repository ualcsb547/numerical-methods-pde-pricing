clc
clear all
close all
format long


a=0; b=2; c=0; T=2; 

M = 1000;
dt = (T-c)/(M+1); % partición temporal
t1 = c : dt : T; 


N = 20;
h = (b-a)/(N+1); % partición espacial
x1 = a : h : b;

k = 1;
gamma = dt*k/h;  % Número de Courant


% Comprobación de estabilidad de ondas
if gamma > 1
    warning(['Este esquema no es estable. La velocidad de propagación de la ' ...
        'información en la malla es mayor que en el espacio físico.'])
end

[t, x] = meshgrid(t1, x1); 
z = zeros(length(x1), length(t1));


% Fronteras espaciales
z(1, :) = 0; 
z(end, :) = 0;

% Posición inicial en t=0 
z(:, 1) = sin(pi .* x1); 

% Velocidad inicial
g = pi .* sin(pi .* x1'); 

% Cálculo del instante t=2 usando la velocidad g
for j = 2:length(x1)-1
    z(j, 2) = (1 - gamma^2)*z(j, 1) + (gamma^2/2)*(z(j+1, 1) + z(j-1, 1)) + dt * g(j);
end 

% Interior
for i = 3:length(t1)
    for j = 2:length(x1)-1
        z(j, i) = 2*(1 - gamma^2)*z(j, i-1) + gamma^2*(z(j+1, i-1) + z(j-1, i-1)) - z(j, i-2);
    end
end

% Solución exacta
u_exacta = (cos(pi.*t) + sin(pi.*t)) .* sin(pi.*x);

% Gráficas 
% i) Solución real y aproximada en T=2 (última columna)
figure
plot(x1, z(:, end), 'r', 'LineWidth', 2)
hold on
plot(x1, u_exacta(:, end), 'b--', 'LineWidth', 2)
grid on
title('Solución en T=2 (Ecuación de Ondas)')
legend('Aproximada', 'Exacta')
xlabel('x'); ylabel('u(x, 2)')

% ii) Gráfica 3D de la solución aproximada
figure
surf(t, x, z)
shading interp
title('Evolución temporal de la onda ')
xlabel('t'); ylabel('x'); zlabel('u(x,t)')