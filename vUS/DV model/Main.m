%% IQ to vUS data processing for in vivo experiment
clear all;
addpath('./Functions/');

SpatialMsk = questdlg('Use ULM mask for vUS spatial constrain?', ...
    'Select', ...
    'YES', 'NO','Cancel', 'Cancel');
%% load data 
disp(['Loading data...']);
load ('./DATA/invivoData.mat'); 
% IQ: beamformed complex quadratue data
% BB: microbubble accumuation map (resize to 25 um pixel size),
    % BB(:,:,1): up flow
    % BB(:,:,2): down flow
    % BB(:,:,3): all flow
% BBV: microbubble velocity map (resize to 25 um pixel size),
    % BBV(:,:,1): up flow
    % BBV(:,:,2): down flow
% BBVz: microbubble axial velocity map (resize to 25 um pixel size),
    % BBVz(:,:,1): up flow
    % BBVz(:,:,2): down flow
[nz,nx,nt]=size(IQ);
%% DAQ infomation and DATA processing parameter
DAQinfo.C=1540;                    % sound speed, m/s
DAQinfo.FWHM=[125 100]*1e-6;        % (X, Z) spatial resolution, Full Width at Half Maximum of point spread function, m
DAQinfo.rFrame=5000;               % sIQ frame rate, Hz
DAQinfo.f0=16.625E6;               % Transducer center frequency, Hz
PRSSinfo.g1nT=nt;                  % g1 calculation sample number
PRSSinfo.g1nTau=100;               % maximum number of time lag
PRSSinfo.SVDrank=[25 nt];          % SVD rank [low high]
PRSSinfo.HPfC=25;                  % High pass filtering cutoff frequency, Hz
PRSSinfo.NEQ=1;                    % 0: no noise equalization; 1: apply noise equalization
%% Clutter rejection
disp('SVD Processing ...');
[sIQ, sIQHP, sIQHHP]=IQ2sIQ(IQ,DAQinfo,PRSSinfo); 
% sIQ: Singular Value Decomposition (SVD) spatiaotemporal filtered signal
% sIQHP: SVD+High Pass filtering
% sIQHHP: SVD+High Pass filtering with 70 HZ cutoff frequency
% eqNoise: obtained noise map
%% vUS data processing
switch SpatialMsk
    case 'YES'
        %% vUS data processing using the ULM mask (BB)
        clear IQ sIQ sIQHHP
        PRSSinfo.rfnScale=2; % spatial resize scale
        PRSSinfo.MskType=1; % use ULM spatial mask
        disp('vUS Processing ...(NOTE: it taks 90-150 mins depending on computing power)');
        tic
        [F, V, Vz, SigmaVz, R]=sIQ2vUS_NP_DV(sIQHP, BB, DAQinfo,PRSSinfo);   
        toc
        %% Results plot
        [VzCmap]=Colormaps_fUS;
        xCoor=[1:nx]*0.05/PRSSinfo.rfnScale;
        zCoor=[1:nz]*0.05/PRSSinfo.rfnScale;

        Fig=figure;
        set(Fig,'Position',[400 400 1300 450])
        subplot(1,2,1)
        h1=imagesc(xCoor,zCoor,V(:,:,1)); % up flow
        colormap(VzCmap);
        caxis([-30 30]);
        colorbar
        axis equal tight;
        hold on;
        h2=imagesc(xCoor,zCoor,V(:,:,2)); % down flow
        AlphaMsk=(abs(V(:,:,2))/3).^4;
        AlphaMsk(AlphaMsk>1)=1;
        set(h2,'AlphaData',AlphaMsk)
        colormap(VzCmap);
        caxis([-30 30]);
        colorbar
        hold off
        axis equal tight;
        xlabel('x [mm]')
        ylabel('z [mm]')
        title('vUS-V [mm/s]')
        
        subplot(1,2,2)
        h1=imagesc(xCoor,zCoor,Vz(:,:,1)); % up flow
        colormap(VzCmap);
        caxis([-30 30]);
        colorbar
        axis equal tight;
        hold on;
        h2=imagesc(xCoor,zCoor,Vz(:,:,2)); % down flow
        AlphaMsk=(abs(V(:,:,2))/3).^4;
        AlphaMsk(AlphaMsk>1)=1;
        set(h2,'AlphaData',AlphaMsk)
        colormap(VzCmap);
        caxis([-25 25]);
        colorbar
        hold off
        axis equal tight;
        xlabel('x [mm]')
        ylabel('z [mm]')
        title('vUS-Vz [mm/s]')
        %% save data
        save('./vUSBB.mat','-v7.3','F','V','Vz','SigmaVz','R','BB','BBV','BBVz','P');
    case 'NO'
        %% PDI processing
        disp('PDI Processing ...');
        [PDIHHP]=sIQ2PDI(sIQHHP);
        clear IQ sIQ sIQHHP
        PRSSinfo.rfnScale=1; % spatial resize scale
        PRSSinfo.MskType=0; % not use ULM spatial mask
        disp('vUS Processing ... (NOTE: it taks 60-90 mins depending on computing power)');
        tic
        [F, V, Vz, SigmaVz, R]=sIQ2vUS_NP_DV(sIQHP,(PDIHHP).^0.5, DAQinfo, PRSSinfo); 
        toc
        %% Results plot
        [VzCmap]=Colormaps_fUS;
        Coor.x=[1:nx]*0.05/PRSSinfo.rfnScale;
        Coor.z=[1:nz]*0.05/PRSSinfo.rfnScale;
        
        Fig=figure;
        set(Fig,'Position',[400 400 1300 450])
        subplot(1,2,1)
        h1=imagesc(Coor.x,Coor.z,V(:,:,1)); % up flow
        colormap(VzCmap);
        caxis([-30 30]);
        colorbar
        axis equal tight;
        hold on;
        h2=imagesc(Coor.x,Coor.z,V(:,:,2)); % down flow
        AlphaMsk=(abs(V(:,:,2))/5).^4;
        AlphaMsk(AlphaMsk>1)=1;
        set(h2,'AlphaData',AlphaMsk)
        colormap(VzCmap);
        caxis([-30 30]);
        colorbar
        hold off
        axis equal tight;
        xlabel('x [mm]')
        ylabel('z [mm]')
        title('vUS-V [mm/s]')
        
        subplot(1,2,2)
        h1=imagesc(Coor.x,Coor.z,Vz(:,:,1)); % up flow
        colormap(VzCmap);
        caxis([-30 30]);
        colorbar
        axis equal tight;
        hold on;
        h2=imagesc(Coor.x,Coor.z,Vz(:,:,2)); % down flow
        AlphaMsk=(abs(Vz(:,:,2))/5).^4;
        AlphaMsk(AlphaMsk>1)=1;
        set(h2,'AlphaData',AlphaMsk)
        colormap(VzCmap);
        caxis([-25 25]);
        colorbar
        hold off
        axis equal tight;
        xlabel('x [mm]')
        ylabel('z [mm]')
        title('vUS-Vz [mm/s]')
        %% plot PDI-based HSV velocity map  
        figure;
        PLOTwtV(V,(PDIHHP).^0.5,Coor,[-30 30])
        title('vUS-V [mm/s]');
        %% save data
        save('./vUS.mat','-v7.3','F','V','Vz','SigmaVz','R','PDIHHP','P');
    case 'Cancel'
end
