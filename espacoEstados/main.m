%% Limpa Workspace

clear;
clc;

%% Prática

A = [ 
    0 1 0;
    0 0 1;
    0 -1 -2
    ]
B = [0; 0; 1];

C = [1 0 0 ];

D = [0];

% Requisitos

UP = 0.05;
Tp = 1; % <1s
e = 0; % rampa e degrau

% Definicação

Gs = ss(A,B,C,D);

zeta = -log(UP)/sqrt(pi()^2 + log(UP)^2);

wn = pi() / (Tp * sqrt(1-zeta^2));

% Polos dominantes
sd1 = -zeta*wn+j*wn*sqrt(1-zeta^2);
sd2 = -zeta*wn-j*wn*sqrt(1-zeta^2);

% Polos adicionais
%sd3 = zero(Gs);
sd3 = -5*zeta*wn;

sd = [sd1 sd2 sd3]

Ab = [ A, zeros(length(A),1); -C, 0];
Bb = [B; 0];

kb = place(A,B,sd);

ks = [kb(1), kb(2), kb(3)]

Gc = ss(A-B * ks, B, C, D);

% Correção do erro em regime permanente
Kd = 1/dcgain(Gc);

Gcs = ss(A-B * ks, B * Kd, C, D);

Gfs = feedback(Gs,1);
step(Gfs,Gcs);
disp(stepinfo(Gcs));
I = [ 1 0 0; 0 1 0; 0 0 1];