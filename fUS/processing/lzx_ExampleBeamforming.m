% example of beamforming method with plane wave emissions
sampleDepth = 3456;
finalmatrix = double(zeros(sampleDepth,128,210));

signthird = 1;
n = zeros(sampleDepth*17,128,30);
for signthird = 1:7
    if signthird == 1
        n = RcvBuffer_001;
    elseif signthird == 2
        n = RcvBuffer_002;
    elseif signthird == 3
        n = RcvBuffer_003;
    elseif signthird == 4
        n = RcvBuffer_004;
    elseif signthird == 5
        n = RcvBuffer_005;
    elseif signthird == 6
        n = RcvBuffer_006;
    elseif signthird == 7
        n = RcvBuffer_007;
    else
        finish;
    end
    for signsecond = 1:30
        for signfirst = 0:16 % Through this loop, you can get the image compounded by 17 angle.
            % load 1908264  % RF data of a ultrasound imaging phantom
            RFmat = n;
            RFframe = RFmat(signfirst*sampleDepth+1:(signfirst + 1)*sampleDepth,:,signsecond);
            RFframe(1:500,:) = 0;
%             size(RFframe)

            % beamforming parameters
            p.c = 1.49;             % sound speed (mm/µs)
            p.Bdx = 0.1979;           % pitch of the array (mm)
            p.Rfech= 65;            % sampling frequency (Mhz)
            p.Rret=2.72;          % delay of first sampling data after emission (µs)
            p.DF =1;                % numerical apperture of the array. 
            % p.angle=-0.14;              % emission angle in radians
            p.angle=(signfirst-8)*3.14/180;              % emission angle in radians

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
%             figure(3);
%             imagesc(x,y,log10(b)*20); caxis([-35 -5]); colormap gray;
            finalmatrix(:,:,signthird*30-30+signsecond) = finalmatrix(:,:,signthird*30-30+signsecond)+bf;
        end
    end
end
disp("finished")
filename = 'CompoundImage.mat';
save(filename)
