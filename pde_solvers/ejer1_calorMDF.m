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

% Hagamos un warning que nos avise de que si gamma > 1/2, entonces no podemos
% seguir porque el esquema no sería estable.
if gamma > 1/2
    warning(['Este esquema no es estable, la velocidad de propagación de la ' ...
        'información en la malla es mayor que en el espacio físco'])
end

[t, x] = meshgrid (t1,x1); % cuadrícula

% Inicializamos la matriz no cuadrada por como hemos cogido la partición 
z = zeros (length (x1), length (t1));

z (: ,1) = 4*x1' - 4.*x1'.^2; % primera columna (condición inicial)
z (1 ,:) = 0; % primera fila (condiciones de frontera)
z(end,:) = 0; % última

% Resto de filas y columnas
for i = 2:length(t1)
    for j = 2:length(x1)-1
        z(j,i)= z(j,i-1) +k * gamma*(z(j+1,i-1)-2.*z(j,i-1)+z(j-1 ,i-1));
    end
end

% i) Solución en 2D en T = 1
figure
plot(x1, z(:, end), 'r', 'LineWidth', 2)
title('Solución a tiempo T=1 (Ecuación del calor)')
xlabel('x'); ylabel('u(x, 1)')

% ii) Solución en 3D
figure
surf (t,x,z)
shading interp
title('Evolución temporal de la ecuación del calor')
xlabel('t'); ylabel('x'); zlabel('u(x,t)')

% EJERCICIO 2: hecho en el archivo: ejercicio2.pdf
