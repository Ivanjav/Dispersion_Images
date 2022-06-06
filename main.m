clc; clear all; close all

%% Cargar Datos
load vertical
nd=20;  %Diezmado
Data=vertical(1:nd:end,10:end);
dt=(5e-5)*nd;
fs=1/dt;
dx=1;
fmin=15; fmax=60;
N=4001;      % Número de puntos para la fft2D
v=100:800;
x=0:dx:(size(Data,2)-1)*dx;
t=0:dt:(size(Data,1)-1)*dt;

%% Shot gather
figure
imagesc(x,t,Data), colormap(gray), set(gca,'fontsize',14), 
title('Shot gather (Original)','FontSize',22,'Interpreter','Latex')
xlabel('Offset (m)','FontSize',18,'Interpreter','Latex'), ylabel('Time (s)','FontSize',18,'Interpreter','Latex')

%% FV domain by FK transform and mapping
[FKdata,FVdata,frequency,wavenumber,velocity]=FVFKdomain(Data,dx,dt,N,fmin,fmax,v(1),v(end),v(2)-v(1));

figure
imagesc(wavenumber,frequency,(abs(FKdata')./max(abs(FKdata')))',[0 1]), set(gca,'fontsize',18,'TickLabelInterpreter','latex'), colormap jet,
title('$f-k$ domain','FontSize',22,'Interpreter','Latex')
ylabel('Frequency (Hz)','FontSize',22,'Interpreter','Latex'), xlabel('Wavenumber (m$^{-1}$)','FontSize',22,'Interpreter','Latex')
set(gca,'YDir','normal')

%% FV domain by Radon transform
[f,v,FVdata_LRT]=fv_domain_LRT(Data,dt,x,v,fmin,fmax);

%%
figure
subplot(121), imagesc(frequency,velocity,abs(FVdata)./max(abs(FVdata))), set(gca,'fontsize',18,'TickLabelInterpreter','latex'), colormap jet,
title('$f-v$ domain','FontSize',22,'Interpreter','Latex')
xlabel('Frequency (Hz)','FontSize',22,'Interpreter','Latex'), ylabel('Velocity (m/s)','FontSize',22,'Interpreter','Latex')
set(gca,'YDir','normal')
subplot(122), imagesc(f,v,abs(FVdata_LRT)./max(abs(FVdata_LRT))), set(gca,'fontsize',18,'TickLabelInterpreter','latex'), colormap jet,
title('$f-v$ domain LRT','FontSize',22,'Interpreter','Latex')
xlabel('Frequency (Hz)','FontSize',22,'Interpreter','Latex'), ylabel('Velocity (m/s)','FontSize',22,'Interpreter','Latex')
set(gca,'YDir','normal')
