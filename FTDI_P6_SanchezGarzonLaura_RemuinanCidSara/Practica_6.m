%% Ejercicio 1: Diseño a partir de filtros ideales centrados
%% FPB ideal
clc
clear all
close all

% Generación del filtro frecuencial con los parámetros deseados
fs=400; % Frecuencia de muestreo en el retículo ortogonal (400 Hz)
v1=1/fs; v2=1/fs; % Vectores retículo ortogonal; número de muestras sea impar en ambas dimensiones
x=[0:v1:1]; % Coordenadas del retículo (normalizamos para que todas las frecuencias se representen entre 0 y 1)
y=[0:v2:1];
[X,Y]=meshgrid(x,y); % Genera matrices bidimensionales

f0x=0.5; % Establece las frecuencias centrales del filtro en las direcciones horizontal y vertical.
f0y=0.5;
Xc = (X - f0x);Yc = (Y - f0y); % Posición central del filtro (en el centro de la imagen)

D0=fs/4; % Frecuencia de corte (horizontal y vertical)
f_filter =double(abs(Xc) <= D0/fs & abs(Yc) <= D0/fs ); % Creamos filtro en forma de ventana cuadrada centrada en (0.5, 0.5) y con una frecuencia de corte=10

% Obtención de su respuesta impulsiva
f_filter_d = ifftshift(f_filter); % Mover primero el centro del filtro (frecuencia cero) al origen de la matriz 
h_d = real(ifft2(f_filter_d)); % Deshacer esta operación en la transformada inversa de Fourier 
h = fftshift(h_d); % Se vuelve a desplazar la respuesta al impulso a la posición (0,0)

% Selección de rango –anchura de la máscara- y obtención de valores enteros en el rango
% deseado cuyo error de redondeo de lugar a mínimo error cuadrático medio respecto a los
% valores racionales originales
cx=1+fs/2; cy=1+fs/2; % (cx,cy) es la posición del centro del filtro
w=1;
filter_mask=h(cy-w:cy+w,cx-w:cx+w); % máscara de filtrado 3x3

[m,j] = min(filter_mask(:)); % Mínimo valor en filter_mask
filter_mask = filter_mask./m; % El menor valor será uno

R = ceil(127./max(abs(filter_mask(:)))); % Máximo factor de ponderación utilizable (rango de la máscara: [−128,128))

% búsqueda de la versión entera que minimiza el MSE
mse = Inf.*ones(1,R);
for j =1:R,
dif = (filter_mask - double(round(filter_mask.*j)./j)).^2;
mse(j) = mean(dif(:));
end
[~,j]=find(mse==min(mse),1);
filter_mask = round(filter_mask.*j);

% Ajuste de la respuesta de la máscara a una señal constante
C=sum(sum(filter_mask));
filter_mask=filter_mask/C;

figure('Name', 'FPB ideal');
subplot(1,3,1); imshow(f_filter); title(sprintf('Filtro paso bajo; E=%g', calcular_energia(f_filter))); colorbar;
subplot(1,3,2); imshow(h); title(sprintf('Impulso; E=%g', calcular_energia(h))); colorbar;
subplot(1,3,3); plot(h); title(sprintf('Impulso; E=%g', calcular_energia(h))); colorbar;

figure('Name', '3D');
subplot(1,2,1); surf(filter_mask); title(sprintf('filter_mask; E=%g', calcular_energia(filter_mask))); colorbar;
subplot(1,2,2); surf(f_filter, "EdgeColor","None"); title(sprintf('Máscara 3x3; E=%g', calcular_energia(f_filter))); colorbar;


[ima,map]=imread("Xray_gray.jpg");
ima=double(ima);
ima_filtered=imfilter(ima,filter_mask);
ima_diferencia=((ima)-(ima_filtered)).^2;

figure('Name', 'FPB ideal');
subplot(1,3,1); imshow(uint8(ima)); title(sprintf('ima original; E=%g', calcular_energia(ima))); colorbar;
subplot(1,3,2); imshow(uint8(ima_filtered)); title(sprintf('ima filtrada; E=%g', calcular_energia(ima_filtered))); colorbar;
subplot(1,3,3); imshow(uint8(ima_diferencia)); title(sprintf('imagen diferencia; E=%g', calcular_energia(ima_diferencia))); colorbar;

%% Ejercicio 2: Diseño a partir de filtros ideales circulares centrados
clc
clear all
close all
% Generación del filtro frecuencial con los parámetros deseados
fs=400; % Frecuencia de muestreo en el retículo ortogonal (400 Hz)
v1=1/fs; v2=1/fs; % Vectores retículo ortogonal; número de muestras sea impar en ambas dimensiones
x=[0:v1:1]; % Coordenadas del retículo (normalizamos para que todas las frecuencias se representen entre 0 y 1)
y=[0:v2:1];
[X,Y]=meshgrid(x,y); % Genera matrices bidimensionales

f0x=0.5; % Establece las frecuencias centrales del filtro en las direcciones horizontal y vertical.
f0y=0.5;
Xc = (X - f0x);Yc = (Y - f0y); % Posición central del filtro (en el centro de la imagen)

D0=fs/4; % Frecuencia de corte (horizontal y vertical)
f_filter =double((Xc).^2 + (Yc).^2) <= (D0/fs).^2; % Creamos filtro en forma de ventana cuadrada centrada en (0.5, 0.5) y con una frecuencia de corte=10

% Obtención de su respuesta impulsiva
f_filter_d = ifftshift(f_filter); % Mover primero el centro del filtro (frecuencia cero) al origen de la matriz 
h_d = real(ifft2(f_filter_d)); % Deshacer esta operación en la transformada inversa de Fourier 
h = fftshift(h_d); % Se vuelve a desplazar la respuesta al impulso a la posición (0,0)

% Selección de rango –anchura de la máscara- y obtención de valores enteros en el rango
% deseado cuyo error de redondeo de lugar a mínimo error cuadrático medio respecto a los
% valores racionales originales
cx=1+fs/2; cy=1+fs/2; % (cx,cy) es la posición del centro del filtro
w=5; % kernel 11x11 --> (11-1)/2 =5
filter_mask=h(cy-w:cy+w,cx-w:cx+w); % máscara de filtrado 3x3

[m,j] = min(filter_mask(:)); % Mínimo valor en filter_mask
filter_mask = filter_mask./m; % El menor valor será uno

R = ceil(127./max(abs(filter_mask(:)))); % Máximo factor de ponderación utilizable (rango de la máscara: [−128,128))

% búsqueda de la versión entera que minimiza el MSE
mse = Inf.*ones(1,R);
for j =1:R,
dif = (filter_mask - double(round(filter_mask.*j)./j)).^2;
mse(j) = mean(dif(:));
end
[~,j]=find(mse==min(mse),1);
filter_mask = round(filter_mask.*j);

% Ajuste de la respuesta de la máscara a una señal constante
C=sum(sum(filter_mask));
h_f = zeros(size(h));
h_f(cy-w:cy+w,cx-w:cx+w) = filter_mask./C;
h_f_d=fftshift(h_f);
f_filter_t_d=fft2(h_f_d);
f_filter_t=abs(ifftshift(f_filter_t_d));
f_filter_t=f_filter_t.*(1/max(max(f_filter_t)));

figure('Name', 'FPB ideal');
subplot(1,3,1); imshow(f_filter); title(sprintf('Filtro paso bajo; E=%g', calcular_energia(f_filter))); colorbar;
subplot(1,3,2); imshow(h_f); title(sprintf('Impulso; E=%g', calcular_energia(h))); colorbar;
subplot(1,3,3); plot(h_f); title(sprintf('Impulso; E=%g', calcular_energia(h))); colorbar;

figure('Name', '3D');
subplot(1,2,1); surf(filter_mask); title(sprintf('filter_mask; E=%g', calcular_energia(filter_mask))); colorbar;
subplot(1,2,2); surf(f_filter,"EdgeColor","None"); title(sprintf('Máscara 3x3; E=%g', calcular_energia(f_filter))); colorbar;

[ima,map]=imread("Xray_gray.jpg");
ima=double(ima);
ima_filtered=imfilter(ima,filter_mask);
ima_diferencia=((ima)-(ima_filtered)).^2;

figure('Name', 'FPB ideal');
subplot(1,3,1); imshow(uint8(ima)); title(sprintf('ima original; E=%g', calcular_energia(ima))); colorbar;
subplot(1,3,2); imshow(uint8(ima_filtered)); title(sprintf('ima filtrada; E=%g', calcular_energia(ima_filtered))); colorbar;
subplot(1,3,3); imshow(uint8(ima_diferencia)); title(sprintf('imagen diferencia; E=%g', calcular_energia(ima_diferencia))); colorbar;





function energia = calcular_energia(imagen)

imagen=double(imagen); % para evitar desbordamientos en caso de unit, logical, ...
energia = sum(sum(imagen .* imagen));

% Si la imagen tiene múltiples canales (por ejemplo, RGB), sumar las energías de cada canal
if size(imagen, 3) > 1
    energia = sum(energia(:));
end

end