
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beamforming method for plane wave emissions
% D(nch,nsamples) matrix of RF data
%    of size nch (number of channels) and nsamples (number of samples)
% p parameter structure see ExampleBeamforming 
% x,z coordinates of the points to be beamformed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bf=beamforming(D,p,x,z);

nx=length(x);
nz=length(z);
bf=zeros(nz,nx);      % init beamforming matrix                  

p.nch=size(D,2);
p.nsamples=size(D,1);

for ix=1:nx            % main loop, coordinates of the x points
%     ix
    for iz=1:nz
        [r,ap]=retard(p,x(ix),z(iz));      % calculate delays and apodization             %
        for j=1:p.nch
            bf(iz,ix)=bf(iz,ix)+D(r(j),j)*ap(j);   % delay and add the RF data 
        end
    end
end

bf=hilbert(bf);    %analitical signal



function [ri,ap]=retard(p,x0,z0)

x = (0:p.nch-1)*p.Bdx;                                        % array axis 
r= (cos(p.angle)*z0+x0*sin(p.angle)+sqrt(z0*z0+(x-x0).*(x-x0)))./p.c;  % delay in time
ri= round((r-p.Rret)*p.Rfech);                                % index of the data
ap= abs((x-x0)/z0) < p.DF ;                                   % square window apodisation
ri(find(ri<=0 | ri>=p.nsamples))=1;                           % limits (r=1 if data does not exist)

