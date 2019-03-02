clear; clc; close all;

%% input params
Nx = 101;
Ny = 101;
Lx = 10;
Ly = 10;
T = 1;
%% read data
filePath = 'data.txt';
fileID = fopen(filePath,'r');
fgets(fileID);
U = zeros(Ny,Nx);
for j=1:Ny
    tline = fgets(fileID);
    U(j,:) = cell2mat(textscan(tline,'%f'))';
end
V = zeros(Ny,Nx);
fgets(fileID);
for j=1:Ny
    tline = fgets(fileID);
    V(j,:) = cell2mat(textscan(tline,'%f'))';
end
fclose(fileID);

%% plot 2D surface

x = linspace(-Lx/2,Lx/2, Nx);
y = linspace(-Ly/2,Ly/2, Ny);
[X,Y] = meshgrid(x,y);

% U
figure;
mesh(X, Y, U);
title('U velocity field at T=1');

% V
figure;
mesh(X, Y, V);
title('V velocity field at T=1');