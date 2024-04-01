%% Laboratório de Sistemas dinâmicos
% Prática 02
% Data: 25/03/2024
% Autores: Ana Clara Gomes & João Vitor Barbosa


%%Limpar Workspace
clear all;
close all;
clc;

%%
% Script Prática 2

%% Ex. 01

syms a b c d;
M = [a b; c d];

disp('DET:');
pretty(det(M));
disp(' ');
disp('INV:');
pretty(inv(M));
disp(' ');
disp('TRACO:');
pretty(trace(M));
disp(' ');

%% Ex. 02 (a)

n = 0:10; 
x_n = (-1).^n;

subplot(3,1,1);
stem(n,x_n);
xlabel('n');
ylabel('x[n]');

subplot(3,1,2);
stairs(n,x_n);
xlabel('n');
ylabel('x[n]');

subplot(3,1,3);
bar(n,x_n);
xlabel('n');
ylabel('x[n]');

%% Ex. 02 (b)
n = 0:10; 
x_n = cos(((pi / 12).*n) + (pi/4));

subplot(3,1,1);
stem(n,x_n);
xlabel('n');
ylabel('x[n]');

subplot(3,1,2);
stairs(n,x_n);
xlabel('n');
ylabel('x[n]');

subplot(3,1,3);
bar(n,x_n);
xlabel('n');
ylabel('x[n]');


%% Ex. 03

y = zeros(1,10);
y(1) = 10;
y(2) = 22;

for i = 3:10
   y(i) = (2 * y(i -1)) + 2;  
end

disp(y);


%% Ex. 04

y = zeros(1,10);
n = -2:1:10;
y(1) = 1;
y(2) = 2;

for i = 3:10
   y(i) = y(i -1) - (0.24 * y(i - 2)) + (n(i)) - (2 * (n(i) - 1));  
end

disp(y);
