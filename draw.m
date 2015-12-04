clc; clear; close all;

nsrc = 281;

dt=0.008;
dx=0.5;
t=(0:8001)*dt;
x=(0:nsrc)*dx;
f1=fopen('initw.dat')
uz=fread(f1,[8001,nsrc],'float32');
subplot(1,2,1)
imagesc(uz)

f2=fopen('initu.dat')
ux=fread(f2,[8001,nsrc],'float32');
subplot(1,2,2)
imagesc(ux)


