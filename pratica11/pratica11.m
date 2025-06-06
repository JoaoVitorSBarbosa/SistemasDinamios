%% Laboratório de Sistemas Dinâmicos
% Prática 11
% Data: 25/08/2024
% Autores: Ana Clara Gomes & João Vitor Barbosa

% Limpar Workspace
clear all;
close all;
clc;

%%
% Script Prática 11

%% Exercicio 1 
dados = load('ensaio_prbs.txt');
t = dados(:,1);
u = dados(:,2);
y = dados(:,3);

yd = y(1:4:end);
ud = u(1:4:end);

% Estimação para 1 atraso
Y = yd(2:end);
X = [yd(1:end-1), ud(1:end-1)];

teta = pinv(X) * Y;
A = teta(1);
B = teta(2);

ym1 = zeros(size(y));  % Inicializa vetor de saída modelada
ym1(1) = y(1);
for k = 2:numel(y)
    ym1(k) = A * ym1(k-1) + B * u(k-1);
end

% Gráfico separado para 1 atraso
figure;
plot(t, y, 'b'); grid on; hold on;
plot(t, ym1, 'r'); grid on;
title('Saída Real vs Modelada (1 atraso)');
xlabel('Tempo');
ylabel('Saída');
legend('Saída Real', 'Modelo 1 atraso');
err1 = immse(y, ym1);
fprintf('Erro Médio Quadrático para 1 atraso: %f\n', err1);

% Estimação para 2 atrasos
Y = yd(3:end);
X = [yd(1:end-2), yd(2:end-1), ud(1:end-2), ud(2:end-1)];

teta2 = pinv(X) * Y;
A2 = teta2(1);
B2 = teta2(2);
C2 = teta2(3);
D2 = teta2(4);

ym2 = zeros(size(y));  % Inicializa vetor de saída modelada
ym2(1) = y(1);
ym2(2) = y(2);
for k = 3:numel(y)
    ym2(k) = A2 * ym2(k-1) + B2 * ym2(k-2) + C2 * u(k-1) + D2 * u(k-2);
end

% Gráfico separado para 2 atrasos
figure;
plot(t, y, 'b'); grid on; hold on;
plot(t, ym2, 'g'); grid on;
title('Saída Real vs Modelada (2 atrasos)');
xlabel('Tempo');
ylabel('Saída');
legend('Saída Real', 'Modelo 2 atrasos');
err2 = immse(y, ym2);
fprintf('Erro Médio Quadrático para 2 atrasos: %f\n', err2);

% Gráfico de comparação
figure;
plot(t, y, 'b'); grid on; hold on;
plot(t, ym1, 'r'); 
plot(t, ym2, 'g'); 
title('Comparação entre Saída Real e Modelos');
xlabel('Tempo');
ylabel('Saída');
legend('Saída Real', 'Modelo 1 atraso', 'Modelo 2 atrasos');
