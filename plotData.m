clear; clc; close all;

%% input params
Nt = 101;
Nx = 101;
Ny = 101;
Lx = 10;
Ly = 10;
T = 1;
%% read data
filePath = 'data.txt';
fileID = fopen(filePath,'r');
fgets(fileID);
U = zeros(Nt,Ny,Nx);
for k=1:Nt
    fgets(fileID);
    for j=1:Ny
        tline = fgets(fileID);
        U(k,j,:) = cell2mat(textscan(tline,'%f'))';
    end
end
fgets(fileID);
V = zeros(Nt,Ny,Nx);
for k=1:Nt
    tline = fgets(fileID);
    for j=1:Ny
        tline = fgets(fileID);
        V(k,j,:) = cell2mat(textscan(tline,'%f'))';
    end
end
fclose(fileID);

%% plot 2D surface

x = linspace(-Lx/2,Lx/2, Nx);
y = linspace(-Ly/2,Ly/2, Ny);
[X,Y] = meshgrid(x,y);
figure;
zlim manual;
for k=1:Nt
    mesh(X, Y, squeeze(U(k,:,:)));
    zlim([0 2]);
    pause(0.01);
end
figure;
zlim manual
for k=1:Nt
    mesh(X, Y, squeeze(V(k,:,:)));
    zlim([0 2]);
    pause(0.01);
end