% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL12_5_38mmFlashAngles_saveRFallFrames.m
% Example of plane wave imaging with steering angle transmits and saving RF data in realtime.
% This script shows how to use external function to save the RF data of all frames in a receive buffer
% in realtime. Here, only non-zeros channel data will be saved to reduce the
% saving time. The default frame limit is set to 100 RcvBuffers. Please change
% to a bigger number for more buffers.
%
% Description:
%   Sequence programming file for L12-5v 38mm Linear array, using plane wave
%   transmits with multiple steering angles. All 128 transmit and receive
%   channels are active for each acquisition. Processing is asynchronous
%   with respect to acquisition.
%
% For convenience, this script is currently set to launch VSX
% automatically. In order to save each frame in a correct order, a "sync"
% command is required for the hardware to wait for the software to finish
% saving RF data, image reconstruction and image display.
% Therefore, a warning message "timeToNextAcq duration too short" might occur
% if the interval (SeqControl(3).argument in this script) between two frames is
% not long enough.
%
% Last update:
% 05/23/2016 - test with SW 3.0.7

clear all
P.startDepth = 2;   % Acquisition depth in wavelengths
P.endDepth = 192;   % This should preferrably be a multiple of 128 samples.
NumBuffer = 7;

% Set the number of scanning angles and scanning interval
na = 17;      % Set na = number of angles.
if (na > 1)
    dtheta = (16*pi/180)/(na-1); % from -8 to +8 degrees
    P.startAngle = -16*pi/180/2;
else
    dtheta = 0;
    P.startAngle=0;
end % set dtheta to range over +/- 18 degrees.

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L12-5 38mm';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans.frequency = 15.000;   % The center frequency for the A/D 4xFc sampling.
Trans = computeTrans(Trans);  % L11-5v transducer is 'known' transducer so we can use computeTrans.
% Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.

% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, 0.5];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
for N = 1:NumBuffer
    Resource.RcvBuffer(N).datatype = 'int16';
    Resource.RcvBuffer(N).rowsPerFrame = na*4096; % this size allows for maximum range(1 aquisition for one angle)
    Resource.RcvBuffer(N).colsPerFrame = Resource.Parameters.numRcvChannels;
    Resource.RcvBuffer(N).numFrames = 30;    % 30 frames stored in RcvBuffer.
end

Resource.InterBuffer(1).numFrames = 1;   % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L12-5_38mmFlashAngles';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', kaiser(Resource.Parameters.numTransmit,1)', ...
                   'aperture', 1, ...%%%%%%
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Trans.numelements)), 1, na);
               
% - Set event specific TX attributes.
if fix(na/2) == na/2       % if na even
    P.startAngle = (-(fix(na/2) - 1) - 0.5)*dtheta;
else
    P.startAngle = -fix(na/2)*dtheta;
end
for n = 1:na   % na transmit events
    TX(n).Steer = [(P.startAngle+(n-1)*dtheta),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,141,275,404,510,603,702,782];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays.
% - We need na Receives for every frame.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
% Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
%                         'startDepth', P.startDepth, ...
%                         'endDepth', maxAcqLength,...
%                         'TGC', 1, ...
%                         'bufnum', 1, ...
%                         'framenum', 1, ...
%                         'acqNum', 1, ...
%                         'sampleMode', 'NS200BW', ...
%                         'mode', 0, ...
%                         'callMediaFunc', 0), 1, (na*Resource.RcvBuffer(1).numFrames)*NumBuffer);
Receive = repmat(struct('Apod', ones(1,Resource.Parameters.numTransmit), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength,...
                        'aperture', 1, ...%%%%%%
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, (na*Resource.RcvBuffer(1).numFrames)*NumBuffer);

% - Set event specific Receive attributes for each frame.
for N = 1:NumBuffer
    ind = Resource.RcvBuffer(N).numFrames*na*(N-1);
    for i = 1:Resource.RcvBuffer(N).numFrames
        Receive(ind+na*(i-1)+1).callMediaFunc = 1;
        for j = 1:na
            Receive(ind+na*(i-1)+j).bufnum = N;
            Receive(ind+na*(i-1)+j).framenum = i;
            Receive(ind+na*(i-1)+j).acqNum = j;
        end
    end
end

% Specify Recon structure arrays.
% - We need one Recon structures which will be used for each frame.
Recon = repmat(struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums', 1:na),1,NumBuffer);
for N = 2:NumBuffer
    Recon(N).RINums = (1+na*(N-1):na*N);
end

% Define ReconInfo structures.
% We need na ReconInfo structures for na steering angles.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, na*NumBuffer);
% - Set specific ReconInfo attributes.
for N = 1:NumBuffer
    if na>1
        ReconInfo(1+na*(N-1)).mode = 'replaceIQ'; % replace IQ data
        for j = 1:na  % For each row in the column
            ReconInfo(na*(N-1)+j).txnum = j;
            ReconInfo(na*(N-1)+j).rcvnum = Resource.RcvBuffer(N).numFrames*na*(N-1)+j;
        end
        ReconInfo(na*N).mode = 'accumIQ_replaceIntensity'; % accum and detect
    else
        ReconInfo(N).mode = 'replaceIntensity';
    end
end

% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',0,...      % display image after processing, 1
                         'displayWindow',1};

for N = 2:NumBuffer+1 %2
    % Save Data Process
    Process(N).classname = 'External';
    Process(N).method = 'saveRFallFrames'; % calls the 'saveRFallFrames' function
    Process(N).Parameters = {'srcbuffer','receive',... % name of buffer to process.
        'srcbufnum', N-1,...
        'srcframenum',0,... % process the all frames in a RcvBuffer
        'dstbuffer','none'};
end

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between s
% synthetic aperture acquisitions
SeqControl(2).argument = 58;  % 58 usec = 1ms/17 * 1000
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 5000 - (na-1)*58;  % 25 msec
SeqControl(4).command = 'returnToMatlab';
SeqControl(5).command = 'sync';
SeqControl(5).argument = 1e6;  % 1s

nsc = 6; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;
% tic;
for N = 1:NumBuffer

    ind = Resource.RcvBuffer(N).numFrames*na*(N-1);
    
    for i = 1:Resource.RcvBuffer(N).numFrames
        for j = 1:na                      % Acquire frame
            Event(n).info = 'Full aperture.';
            Event(n).tx = j;
            Event(n).rcv = ind+na*(i-1)+j;
            Event(n).recon = 0;
            Event(n).process = 0;
            Event(n).seqControl = 2;
            n = n+1;
        end
        Event(n-1).seqControl = [3,nsc]; % modify last acquisition Event's seqControl
        SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
        nsc = nsc+1;
    end
    %disp("collection ends");

    Event(n).info = 'Sync';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 5;
    n = n+1;
    
    tic;
    Event(n).info = 'Save RF data';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = N+1;%N+1
    Event(n).seqControl = 0;
    n=n+1;
    toc;
    %disp("save ends");

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;%N
    Event(n).process = 0;%1
    Event(n).seqControl = 0;
    if floor(i/5) == i/5     % Exit to Matlab every 5th frame
        Event(n).seqControl = 4;
    end
    n = n+1;
%     disp("recon ends");

end
% toc;
% No need to go back to the start, 
% Back to the start of the script
% Event(n).info = 'Jump back';
% Event(n).tx = 0;
% Event(n).rcv = 0;
% Event(n).recon = 0;
% Event(n).process = 0;
% Event(n).seqControl = 1;

% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoffCallback');

% - Range Change
MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');

EF(1).Function = text2cell('%saveRFallFrames');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 680;

% Save all the structures to a .mat file.
filename = 'L12-5_38mmFlashAngles_saveRFallFrame';
save(['MatFiles/',filename]);
VSX;

return

% **** Callback routines to be converted by text2cell function. ****
%SensCutoffCallback - Sensitivity cutoff change
ReconL = evalin('base', 'Recon');
for i = 1:size(ReconL,2)
    ReconL(i).senscutoff = UIValue;
end
assignin('base','Recon',ReconL);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Recon'};
assignin('base','Control', Control);
return
%SensCutoffCallback

%RangeChangeCallback - Range change
simMode = evalin('base','Resource.Parameters.simulateMode');
% No range change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.endDepth'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
P.endDepth = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        P.endDepth = UIValue*scaleToWvl;
    end
end
assignin('base','P',P);

evalin('base','PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));');
evalin('base','PData(1).Region = computeRegions(PData(1));');
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback

%saveRFallFrames - save RF
saveRFallFrames(RData)

persistent RcvBufferNum RFfilename

% file size can be reduced by elimating all zeros
TXApod = evalin('base','TX.Apod');
endSample = evalin('base','Receive(end).endSample');

RcvNumLimit = 100;
numLength = max(ceil(log10(abs(RcvNumLimit))),1)+1;

if isempty(RcvBufferNum)
    RcvBufferNum = 1;
    RFfilename = ['RFdata_',datestr(now,'dd-mmmm-yyyy_HH-MM-SS')];
else
    RcvBufferNum = RcvBufferNum + 1;
    if RcvBufferNum > RcvNumLimit, RcvBufferNum = 1; end % set a limit for testing
end

fname = ['RcvBuffer_',num2str(RcvBufferNum,['%0',num2str(numLength),'.0f'])];
newRData = RData(1:endSample,find(TXApod),:);
eval([fname,' = newRData;']);

tic
if isequal(RcvBufferNum,1)
    save(RFfilename,fname,'-v6');
else
    save(RFfilename,fname,'-v6','-append');
end
fprintf('saving time for buffer %g = %g s \n',RcvBufferNum,toc);
return
%saveRFallFrames

