%% IQ to vUS data processing for ex vivo data using the basic model
clear all;
addpath('./Functions/');
%% load data 
disp(['Loading data...']);
% load ('./DATA/exvivoData5a.mat');  % angled flow, preset speed 5 mm/s
load ('./DATA/exvivoData15a.mat');  % angled flow, preset speed 15 mm/s
% load ('./DATA/exvivoData9t.mat');  % transverse flow, preset speed 9 mm/s
% load ('./DATA/exvivoData25t.mat');  % transverse flow, preset speed 25 mm/s
% IQ: beamformed complex quadratue data
[nz,nx,nt]=size(IQ);
tic
%% DAQ infomation and DATA processing parameter
DAQinfo.C=1540;                    % sound speed, m/s
DAQinfo.FWHM=[125 90]*1e-6;        % (X,Z) spatial resolution, Full Width at Half Maximum of point spread function, m
DAQinfo.rFrame=5000;               % sIQ frame rate, Hz
DAQinfo.f0=16.625E6;               % Transducer center frequency, Hz
PRSSinfo.g1nT=nt;                  % g1 calculation sample number
PRSSinfo.g1nTau=100;               % maximum number of time lag
PRSSinfo.SVDrank=[3 nt];           % SVD rank [low high]
PRSSinfo.HPfC=10;                  % High pass filtering cutoff frequency, Hz
PRSSinfo.NEQ=1;                    % 0: no noise equalization; 1: apply noise equalization
%% Clutter rejection
disp('SVD Processing ... ');
[sIQ]=IQ2sIQ(IQ,DAQinfo,PRSSinfo); 
% sIQ: SVD filtered data
%% vUS data processing
clear IQ
PRSSinfo.rfnScale=1; % spatial resize scale
disp('vUS Processing ...(NOTE: it taks around 6-10 mins depending on computing power)');
[F, V, Vz, SigmaVz, Vcz, R]=sIQ2vUS_SV(sIQ, DAQinfo,PRSSinfo);
toc
%% Results plot
[VzCmap,VzCmapDn]=Colormaps_fUS;
Coor.x=[1:nx]*0.05/PRSSinfo.rfnScale;
Coor.z=[1:nz]*0.05/PRSSinfo.rfnScale;
Fig=figure;
set(Fig,'Position',[400 400 1700 350])
subplot(1,3,1)
h1=imagesc(Coor.x,Coor.z,abs(V)); 
colormap(VzCmapDn);
caxis([0 30]);
colorbar
axis equal tight;
xlabel('x [mm]')
ylabel('z [mm]')
title('vUS-V [mm/s]')

subplot(1,3,2)
h2=imagesc(Coor.x,Coor.z,abs(Vz)); 
colormap(VzCmapDn);
caxis([0 30]);
colorbar
axis equal tight;
xlabel('x [mm]')
ylabel('z [mm]')
title('vUS-Vz [mm/s]')

subplot(1,3,3)
h3=imagesc(Coor.x,Coor.z,abs(Vcz)); 
colormap(VzCmapDn);
caxis([0 30]);
colorbar
axis equal tight;
xlabel('x [mm]')
ylabel('z [mm]')
title('Color Doppler-Vz [mm/s]')
