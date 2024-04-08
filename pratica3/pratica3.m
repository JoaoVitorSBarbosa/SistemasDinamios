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

%% Ex. 02 (b)

s = tf('s');


c_s = 1 / (s+1);
g_s = 1 / (s+3);

y_s = series(c_s,g_s);


step(y_s);

%% Ex. 02 (c)

s = tf('s');


c_s = 1 / (s+1);
g_s = 1 / (s+3);

y_s = series(c_s,g_s);

subplot(2,1,1);
step(c_s);
subplot(2,1,2);
step(g_s);

%% Ex. 02 (d)

s = tf('s');


c_s = 1 / (s+1);
g_s = 1 / (s+3);

y_s = series(c_s,g_s);


impulse(y_s);

%% Ex. 03 (a)

s = tf('s');


c_s = 1 / (s+1);
g_s = 1 / (s+3);

[R, P, K] = tf2zp(c_s.Numerator{1}, c_s.Denominator{1});
[RR, PP, KK] = tf2zp(g_s.Numerator{1}, g_s.Denominator{1});


disp('C(s) Params:');
disp('Ganho ');
disp(K);
disp('Polos');
disp(c_s.Denominator{1});
disp('Zeros');
disp(c_s.Numerator{1});

disp('G(s) Params:');
disp('Ganho ');
disp(KK);
disp('Polos');
disp(g_s.Denominator{1});
disp('Zeros');
disp(g_s.Numerator{1});

%% Ex. 03 (b)

s = tf('s');


c_s = 1 / (s+1);
g_s = 1 / (s+3);

y_s = series(c_s,g_s);


c = step(c_s);
g = step(g_s);

save step.mat;

figure(1);
subplot(5,1,1);
t1 = [0:0.1:1];
y1= step(c_s,t1);
plot(t1,y1);

xlabel('Tempo(s)');
ylabel('Amplitude');
title('Resposta ao degrau'); 
grid('on');

subplot(5,1,2);
t2 = [0:0.1:2];
y2= step(c_s,t2);
plot(t2,y2);

xlabel('Tempo(s)');
ylabel('Amplitude');
title('Resposta ao degrau'); 
grid('on');

subplot(5,1,3);
t3 = [0:0.1:3];
y3= step(c_s,t3);
plot(t3,y3);

xlabel('Tempo(s)');
ylabel('Amplitude');
title('Resposta ao degrau'); 
grid('on');

subplot(5,1,4);
t4 = [0:0.1:4];
y4= step(c_s,t4);
plot(t4,y4);

xlabel('Tempo(s)');
ylabel('Amplitude');
title('Resposta ao degrau'); 
grid('on');

subplot(5,1,5);
t5 = [0:0.1:5];
y5= step(c_s,t5);
plot(t5,y5);

xlabel('Tempo(s)');
ylabel('Amplitude');
title('Resposta ao degrau'); 
grid('on');

%G(s)

figure(2);

subplot(5,1,1);

y6= step(g_s,t1);
plot(t1,y6);

xlabel('Tempo(s)');
ylabel('Amplitude');
title('Resposta ao degrau'); 
grid('on');

subplot(5,1,2);
t2 = [0:0.1:2];
y7= step(g_s,t2);
plot(t2,y7);

xlabel('Tempo(s)');
ylabel('Amplitude');
title('Resposta ao degrau'); 
grid('on');

subplot(5,1,3);
t3 = [0:0.1:3];
y8= step(g_s,t3);
plot(t3,y8);

xlabel('Tempo(s)');
ylabel('Amplitude');
title('Resposta ao degrau'); 
grid('on');

subplot(5,1,4);
t4 = [0:0.1:4];
y9= step(g_s,t4);
plot(t4,y9);

xlabel('Tempo(s)');
ylabel('Amplitude');
title('Resposta ao degrau'); 
grid('on');

subplot(5,1,5);
t5 = [0:0.1:5];
y10= step(g_s,t5);
plot(t5,y10);

xlabel('Tempo(s)');
ylabel('Amplitude');
title('Resposta ao degrau'); 
grid('on');


%% Ex. 04 (a)

s = tf('s');

u_s = 2 / (s + 0.5);
d_s = 2.5 / (s + 0.5);

y_s = u_s - d_s;

printsys(y_s.Numerator{1}, y_s.Denominator{1});

%% Ex. 04 (b)

s = tf('s');

u_s = 2 / (s + 0.5);
d_s = 2.5 / (s + 0.5);

y_s = u_s - d_s;

for i = 3:10
   y(i) = (2 * y(i -1)) + 2;  
end
