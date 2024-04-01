%% Laboratório de Sistemas dinâmicos
% Prática 03
% Data: 01/04/2024
% Autores: Ana Clara Gomes & João Vitor Barbosa


%%Limpar Workspace
clear all;
close all;
clc;

%%
% Script Prática 3

%% Ex. 01 (a)

s = tf('s');
q_s = s + 1;
p_s = (s^2) + (2 * s) + 1;

disp('p(s) * q(s): ');
res = p_s * q_s;
printsys(res.num{1}, res.den{1});
disp(' ');

%% Ex. 01 (b)
s = tf('s');
q_s = s + 1;
p_s = (s^2) + (2 * s) + 1;

G_s = q_s / p_s;

disp('Polos: ');
disp(pole(G_s))
disp(' ');
disp('Zeros: ');
disp(zero(G_s))
disp(' ');

%% Ex. 01 (c)

s = tf('s');

p_s = (s^2) + (2 * s) + 1;


disp(evalfr(p_s,-1));


%% Ex. 01 (d)

G_s = tf([1 1],[1 2 1]);

pzmap(G_s);
grid on;

%% Ex. 02 (a)

s = tf('s');


c_s = 1 / (s+1);
g_s = 1 / (s+3);

y_s = series(c_s,g_s);



printsys(y_s.num{1}, y_s.den{1});
