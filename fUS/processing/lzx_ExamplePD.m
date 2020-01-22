clear;
load('CompoundImage.mat');  % read data of a small part of the brain cortex (IQR matrix)  

% IQR is a 3D matrix of 52x49x200 elements
% the first dimension is the depth (52 poins)
% the second dimension is the length (49 poins)
% the third dimension is the time (200 frames)
depth = 3456/60e6*1470*1000/2;
for i = 1:size(finalmatrix,3)
    finalmatrix2(:,:,i) = abs(hilbert(finalmatrix(:,:,i)));
    finalmatrix3(:,:,i) = imresize(finalmatrix2(:,:,i),[depth*10 128*0.2*10]);
end
% example of one ultrasound image
figure(1)
% imagesc(abs(finalmatrix3(:,:,1)));

imagesc(db(finalmatrix3(:,:,1)/max(max(finalmatrix3(:,:,1))))); 
caxis([-35 0]);

colormap gray;
title('Ultrasound image');

% computation of the power Doppler image
[B,A]=butter(4,70/1000*2,'high');    %coefficients for the high pass filter

% sustraction of the first image
% the signal stat at 0 and minimises filter oscilatons
for i=1:size(finalmatrix3,1)
    for j=1:size(finalmatrix3,2)
        finalmatrix3(i,j,:)=finalmatrix3(i,j,:)-finalmatrix3(i,j,1); 
    end
end

sb=filter(B,A,finalmatrix3,[],3);    % blood signal (filtering in the time dimension)
sb=sb(:,:,4:end);           % the first 4 temporal samples are eliminates (filter oscilations)
PDI=mean(abs(sb).^2,3);     % computing the intensity of the blood signal the 
                             
% display the power doppler image
figure(2);
imagesc(10*log10(PDI./max(PDI(:)))); 
caxis([-35 0]);
colormap gray;
title('PDI image');
axis image