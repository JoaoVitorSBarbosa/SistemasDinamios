%% Laboratório de Sistemas dinâmicos
% Prática 01
% Data: 11/03/2024
% Autores: Ana Clara Gomes & João Vitor Barbosa


%%Limpar Workspace
clear all;
close all;
clc;

%%
% Script Prática 1

%% Ex. 01 (a)

arr1 = [1 1 6;5 -2 1;-8 2 -3];
arr2 = [2 9; -5 -1; 9 2];

sizeArr1 = size(arr1);
sizeArr2 = size(arr2);

if (sizeArr1(1) == sizeArr1(2))
    disp('Matriz 1 eh quadrada'); 
else
    disp('Matriz 1 não eh quadrada');
end

if (sizeArr2(1) == sizeArr2(2))
    disp('Matriz 2 eh quadrada'); 
else
    disp('Matriz 2 não eh quadrada');
end
%% Ex. 01 (b)

arr1 = [1 1 6;5 -2 1;-8 2 -3];
arr2 = [2 9; -5 -1; 9 2];

sizeArr1 = size(arr1);
sizeArr2 = size(arr2);

disp('Matriz 1');
disp(' ');
for linha = 1:sizeArr1(1)
    for coluna = 1:sizeArr1(2)
        if (arr1(linha,coluna) == 2)
            string = ['Achou numero 2 na linha ' int2str(linha) ' coluna ' int2str(coluna)];
            disp(string);
        end
    end 
end

disp(' ');
disp('Matriz 2');
disp(' ');
for linha = 1:sizeArr2(1)
    for coluna = 1:sizeArr2(2)
        if (arr2(linha,coluna) == 2)
            string = ['Achou numero 2 na linha ' int2str(linha) ' coluna ' int2str(coluna)];
            disp(string);
        end
    end 
end
%% Ex. 01 (c)
arr1 = [1 1 6;5 -2 1;-8 2 -3];
arr2 = [2 9; -5 -1; 9 2];

sizeArr1 = size(arr1);
sizeArr2 = size(arr2);

disp('Matriz 1');
disp(' ');
for linha = 1:sizeArr1(1)
    for coluna = 1:sizeArr1(2)
        if (arr1(linha,coluna) < 0)
            string = ['Achou numero negativo na linha ' int2str(linha) ' coluna ' int2str(coluna)];
            disp(string);
        end
    end 
end

disp(' ');
disp('Matriz 2');
disp(' ');
for linha = 1:sizeArr2(1)
    for coluna = 1:sizeArr2(2)
        if (arr2(linha,coluna) < 0)
            string = ['Achou numero negativo na linha ' int2str(linha) ' coluna ' int2str(coluna)];
            disp(string);
        end
    end 
end

%% Ex. 02
t=-10:0.01:10;
y = exp(-2 * t).* sin(2 * pi * 3 * t);
plot(t,y);

xlabel('Tempo (s)');
ylabel('Amplitude');

%% Ex. 03
t=-10:0.01:10;
x = cos(t).* sin(20*t);
plot(t,x);

xlabel('Tempo (s)');
ylabel('Amplitude');

%% Ex. 04 (a)

[R, P, K] = residue([6 6],[1 4.59 0.58 0]);

disp('R');
disp(R);
disp(' ');
disp('P');
disp(P);
disp(' ');
disp('K');
disp(K);
disp(' ');


    

%% Ex. 04 (b)

[R, P, K] = residue([1 2 3],[1 3 3 1]);

disp('R');
disp(R);
disp(' ');
disp('P');
disp(P);
disp(' ');
disp('K');
disp(K);
disp(' ');
