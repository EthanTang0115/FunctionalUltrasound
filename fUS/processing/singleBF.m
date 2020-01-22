% example of beamforming method with plane wave emissions    
sampleDepth = 3456;
% load 1908264  % RF data of a ultrasound imaging phantom
RFmat = RcvBuffer_001;
RFframe = RFmat(1+2*sampleDepth:3*sampleDepth,:,1);
RFframe(1:200,:) = 0;

% beamforming parameters
p.c = 1.49;             % sound speed (mm/µs)
p.Bdx = 0.1979;           % pitch of the array (mm)
p.Rfech= 65;            % sampling frequency (Mhz)
p.Rret=2.72;          % delay of first sampling data after emission (µs)
p.DF =1;                % numerical apperture of the array. 
% p.angle=-0.14;              % emission angle in radians
p.angle=(-6)*3.14/180;              % emission angle in radians

sampleSpacing = 1/(p.Rfech*1e6)*p.c*1e6/2;

% axis of the pixels to image
x=[1:128].*p.Bdx;       % x axis for beamforming mm
% y=[10:1/50*1.5/2:37];   % z axis for beamforming mm
y=[1:sampleDepth]*sampleSpacing;   % z axis for beamforming mm

D=double(RFframe);
%tic;
bf=beamforming(D,p,x,y);  % beamforming procedure
disp("bf running ...")
%toc;

% visualization in log scale
b=abs(bf)./max(abs(bf(:)));
figure(3);
imagesc(x,y,log10(b)*20); caxis([-35 -5]); colormap gray;
%finalmatrix(:,:,1) = finalmatrix(:,:,1)+bf;
       
%disp("finished")
%filename = 'CompoundImage.mat';
% save(filename)
