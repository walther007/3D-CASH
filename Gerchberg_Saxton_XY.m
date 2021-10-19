% GERCHBERG SAXTON 1D phase retrieval;
% WALTHER April 2014
% Gerberg_Saxton_XY_vers_6

currentFolder = cd;

if ~exist(fullfile(currentFolder, 'Hologram'), 'dir'), mkdir(fullfile(currentFolder, 'Hologram')); end
if ~exist(fullfile(currentFolder, 'Functions'), 'dir'), error('Functions folder missing in the active directory'); end

path(path, fullfile(currentFolder, 'Functions'));
    
% **** Select Target pattern *******************************************
if ~exist('patternSELECT', 'var'); patternSELECT = 3; end
if ~any(ismember(patternSELECT, 1:4)); error('Pattern select out of range'); end 
% Available patterns are dfined in the case structure further below
% ***********************************************************************

% ******************** GENERAL PARAMETERS *******************************
lambda = 940e-6;  % wavelength in [mm]

rObj = 15.12/2;  % [mm] radius of objective BFP
fObj = 8; %[mm] objective focal length  OLYMPUS 25x
% OLYMPUS tube lemgth = 180 nm   f = 180/25 = 8 mm

rAOD = 15/2;  % in [mm] 
rBeam = 13/2;   % radius where the amplitude of a Gaussian beam decayed to 1/e, measured in AOD coordinates

% Take angular magnification by tube/scan telescope into account
f1 = 300;
f2 = 200;

Mlat = f2/f1;
Mang = 1/Mlat;

xtot = 30;
nx = 2^nextpow2(10000);

ytot = 30;
ny = 2^nextpow2(10000);

vAOD = 650; % AOD acoustique velocity [mm/ms]
tAOD = 2*rAOD/vAOD;  % in [ms]
writeAOD = 14000; % [kHz=1/ms] % AOD write clock [kHz]
nAOD = ceil(tAOD * writeAOD); % number of DDS discretisation points

% ---- Initialisation of flags and parameters ----

ampliMod_X = false;
ampliMod_Y = false;

dynTarget_X = true;
dynTarget_Y = true;

initFocal = true;

fixInit = false;
valInit_X = 100;
valInit_Y = 100;
yPhaseInverse = false;

xBWlimit = 8; % [MHz]
yBWlimit = 8; % [MHz]

X_smoothKernel = floor(2*rAOD/xtot * nx/nAOD /2)*2+1;
Y_smoothKernel = floor(2*rAOD/ytot * ny/nAOD /2)*2+1;

% -----------------------------------------------------------
% Definition of target pattern
xTargetAmpl = [];
yTargetAmpl = [];

switch patternSELECT
    case 1
        xTarget = 0;
        yTarget = 0;
    
    case 2 % 3x3 pattern
        xTarget = 0.008 * [-1,0,1]/2;
        yTarget = 0.008 * [-1,0,1]/2;
        %yTarget = 0;
        
        ampliMod_X = true;
        ampliMod_Y = true;

        dynTarget_X = false;
        dynTarget_Y = false;

        initFocal = true;
        fixInit = true;
        valInit_X = 1000;
        valInit_Y = 1000;
        yPhaseInverse = true;
        smoothON = false;
    
    case 3 % 5x5 (16x16 um) pattern 
        xTarget = 0.016 * linspace(-1,1,5)/2;
        yTarget = 0.016 * linspace(-1,1,5)/2;
        
        ampliMod_X = true;
        ampliMod_Y = true;

        dynTarget_X = false;
        dynTarget_Y = false;
        
        initFocal = true;
        fixInit = true;
        valInit_X = 10000;
        valInit_Y = 10000;
        yPhaseInverse = true;        
        smoothON = false;
        
    case 4 % 2x5 pattern
        xTarget = 0.024 * [-1,1]/2;
        yTarget = 0.024 * linspace(-1,1,5)/2;
        
        ampliMod_X = false;
        ampliMod_Y = true;
        
        dynTarget_X = false;
        dynTarget_Y = false;
        
        initFocal = true;
        fixInit = true;
        valInit_X = 1000;
        valInit_Y = 1000;
        yPhaseInverse = false;
        smoothON = false;
end
        
% ----------------------------------------------------------
% *************************************************************************
% ************** Definition of optical fields ******************
% 1 D fields in holographic plane 
% hx(:,1) = 1D spatial x-coordinates
% hx(:,2) = amplitude function
% hx(:,3) = 1D phase function
% HX = Fourier transform of hx

% 1 D fields in Fourier plane 
% fx(:,1) = 1D spatial x-coordinates
% fx(:,2) = 1D amplitude function
% fx(:,3) = 1D phase function
% FX = fourrier transform of fx

% hy(:,1) = 1D spatial y-coordinates
% hy(:,2) = amplitude function
% hy(:,3) = 1D phase function
% HY = Fourier transform of hy

% 1 D fields in Fourier plane 
% fy(:,1) = 1D spatial y-coordinates
% fy(:,2) = 1D amplitude function
% fy(:,3) = 1D phase function
% FY = fourrier transform of fy

% spatial coordinates in ao plane 
xao = xtot/Mlat * linspace(-0.5,0.5,nx)';  % in [mm]
yao = ytot/Mlat * linspace(-0.5,0.5,ny)';  % in [mm]

% spatial coordinates in the holographic objective plane 
xh = xtot * linspace(-0.5,0.5,nx)';  % in [mm]
yh = ytot * linspace(-0.5,0.5,ny)';  % in [mm]

% **** offset for amplitude signal in the holographic plane ****

% Spatial frequencies in holographic plane
XH = nx/xtot * linspace(-0.5,0.5,nx)';  % in [1/mm]
YH = ny/ytot * linspace(-0.5,0.5,ny)';  % in [1/mm]

% spatial coordinates in front focal plane
xf = nx*fObj*lambda/xtot * linspace(-0.5,0.5,nx)';  % in [mm]
yf = ny*fObj*lambda/ytot * linspace(-0.5,0.5,ny)';  % in [mm]

% spatial frequencies in front focal plane
XF = xtot/(fObj*lambda) * linspace(-0.5,0.5,nx)';  % in [1/mm]
YF = ytot/(fObj*lambda) * linspace(-0.5,0.5,ny)';  % in [1/mm]

% field in holographic plane
hx = [xh, ones(length(xh),1), zeros(length(xh),1)];
hy = [yh, ones(length(yh),1), zeros(length(yh),1)];

% Fourier spectrum of field in holograpic plane
%HX = [XH, ones(length(XH),1), zeros(length(XH),1)];
%HY = [YH, ones(length(YH),1), zeros(length(YH),1)];

% field in Fourier plane = focal plane
%fx = [xf, ones(length(xf),1), zeros(length(xf),1)];
%fy = [yf, ones(length(yf),1), zeros(length(yf),1)];

% Fourier spectrum of field in focal plane
%FX = [XF, ones(length(XF),1), zeros(length(XF),1)];
%FY = [YF, ones(length(YF),1), zeros(length(YF),1)];

%% Find X hologram

% *** Initialization
% define apertures in the OBA plane
insideObj = abs(xh) <= rObj;
insideBeam = abs(xao) <= rBeam;
insideAOD = abs(xao) <= rAOD;

% Define start holographic input for X
% Standard initialisation
hxAmplIN = ones(length(xh),1);
%hxAmplIN = exp(-(xh.^2)/(Mlat*rBeam)^2); % Gauss beam
hxAmplIN(~insideObj) = 0; 

% Define focal target pattern for X
if max(abs(xTarget)) > max(abs(xf)), error('Target values out of range'); end

realField = ones(length(xh),1);
%realField(~insideObj) = 0;
field = complex(realField, zeros(length(xh),1));

spec = fftshift(fft(ifftshift(field),nx));
psf = abs(spec) / max(abs(spec));

% ****** Gerchberg Saxton iteration ******
if smoothON
    [ hxAmplOUT, hxPhaseOUT, fxAmplOUT, fxStartPhase, fxTarget ]  = GerchbergSaxton4(hxAmplIN, xf, xTarget, xTargetAmpl, psf, dynTarget_X, 100, 0.98, X_smoothKernel, initFocal, fixInit, valInit_X);
else
    [ hxAmplOUT, hxPhaseOUT, fxAmplOUT, fxStartPhase, fxTarget ] = GerchbergSaxton4(hxAmplIN, xf, xTarget, xTargetAmpl, psf, dynTarget_X, 100, 0.98, [], initFocal, fixInit, valInit_X);
end
      
% *** result of phase optimization
% hxOUT(:,1) = 1D spatial x-coordinates in the AOD optical plane [mm]
% hxOUT(:,2) = x-amplitude function
% hxOUT(:,3) = x-phase function
% hxOUT(:,4) = x-amplitude function after smoothing
% hxOUT(:,5) = x-frequency function [MHz]
% hxOUT(:,6) = x-frequency function [MHz] after smoothing

hxPhaseOUT = hxPhaseOUT - 0.5 *(max(hxPhaseOUT(insideObj))+min(hxPhaseOUT(insideObj)));

if ampliMod_X
    hxAmplOUT = hxAmplOUT/max(hxAmplOUT(insideAOD));
    hxAmplOUT(hxAmplOUT>1) = 1;
    hxOUT = [xao, hxAmplOUT, hxPhaseOUT, zeros(length(xao),3)];
else
    hxOUT = [xao, hxAmplIN, hxPhaseOUT, zeros(length(xao),3)];
    field = hxAmplIN .* exp(1i * 2*pi* hxPhaseOUT);
    spec = fftshift(fft(ifftshift(field), nx))/nx;
    fxAmplOUT = abs(spec);
end

fxAmplOUT = fxAmplOUT / max(fxAmplOUT);

hxOUT(:,2) = hxOUT(:,2) / max(exp(-xao(insideBeam).^2/rBeam^2).*hxOUT(insideBeam,2));
hxOUT(:,4) = smooth(hxOUT(:,2), X_smoothKernel);
hxOUT(:,5) = phase2freq2(hxOUT(:,1), hxOUT(:,3));
hxOUT(:,5) = hxOUT(:,5) - mean(hxOUT(insideBeam,5));
hxOUT(:,6) = smooth(hxOUT(:,5), X_smoothKernel);

% ****************************

xBW = max(hxOUT(insideAOD,6)) - min(hxOUT(insideAOD,6));

hxOUT(hxOUT(:,6) > xBWlimit/2,6) = xBWlimit/2; 
hxOUT(hxOUT(:,6) < -xBWlimit/2,6) = -xBWlimit/2; 

% Calculate AM-induced power drop
ampliGauss = hxOUT(:,4) .* exp(-hxOUT(:,1).^2/rBeam^2);
hxPowDrop = 1 - ( sum(ampliGauss.^2) / sum(exp(-2*hxOUT(:,1).^2/rBeam^2)) );
hxPowTrans = 1 - hxPowDrop;

%% Find Y hologram

% *** Initialization
% define apertures
insideObj = abs(yh) <= rObj;
insideBeam = abs(yao) <= rBeam;
insideAOD = abs(yao) <= rAOD;

% Define start holographic input for Y
% Standard initialisation
hyAmplIN = ones(length(yh),1);
%hyAmplIN = exp(-yao.^2/rBeam^2); % Gauss beam
hyAmplIN(~insideObj) = 0;

% Define focal target pattern for X
if max(abs(yTarget)) > max(abs(yf)), error('Target values out of range'); end

realField = ones(length(yh),1);
%realField(~insideObj) = 0; 
field = complex(realField, zeros(length(yh),1));
spec = fftshift(fft(ifftshift(field),ny));
psf = abs(spec);

% ****** Gerchberg Saxton iteration ******
if smoothON
    [ hyAmplOUT, hyPhaseOUT, fyAmplOUT, fyStartPhase, fyTarget ] = GerchbergSaxton4(hyAmplIN, yf, yTarget, yTargetAmpl, psf, dynTarget_Y, 100, 0.99, Y_smoothKernel, initFocal, fixInit, valInit_Y); 
else
    [ hyAmplOUT, hyPhaseOUT, fyAmplOUT, fyStartPhase, fyTarget ] = GerchbergSaxton4(hyAmplIN, yf, yTarget, yTargetAmpl, psf, dynTarget_Y, 100, 0.99, [], initFocal, fixInit, valInit_Y);
end
% *** result of phase optimization
% hyOUT(:,1) = 1D spatial y-coordinates [mm] in the AOD optical plane
% hyOUT(:,2) = y-amplitude function
% hyOUT(:,3) = y-phase function
% hyOUT(:,4) = y-amplitude function after smoothing
% hyOUT(:,5) = y-frequency function [MHz]
% hyOUT(:,6) = y-frequency function [MHz] after smoothing

hyPhaseOUT = hyPhaseOUT - 0.5 *(max(hyPhaseOUT(insideObj))+min(hyPhaseOUT(insideObj)));
if yPhaseInverse, hyPhaseOUT = (-1) * hyPhaseOUT; end

if ampliMod_Y
    hyAmplOUT = hyAmplOUT/max(hyAmplOUT(insideAOD));
    hyAmplOUT(hyAmplOUT>1) = 1;
    hyOUT = [yao, hyAmplOUT, hyPhaseOUT, zeros(length(yao),3)];
else
    hyOUT = [yao, hyAmplIN, hyPhaseOUT, zeros(length(yao),3)];
    field = hyAmplIN .* exp(1i * 2*pi* hyPhaseOUT);
    spec = fftshift(fft(ifftshift(field), ny))/ny;
    fyAmplOUT = abs(spec);
end

fyAmplOUT = fyAmplOUT / max(fyAmplOUT);

hyOUT(:,2) = hyOUT(:,2) / max(exp(-yh(insideBeam).^2/rBeam^2) .* hyOUT(insideBeam,2));
hyOUT(:,4) = smooth(hyOUT(:,2), Y_smoothKernel);
hyOUT(:,5) = phase2freq2(hyOUT(:,1), hyOUT(:,3));
hyOUT(:,5) = hyOUT(:,5) - mean(hyOUT(insideBeam,5));
hyOUT(:,6) = smooth(hyOUT(:,5), Y_smoothKernel);

yBW = max(hyOUT(insideAOD,6)) - min(hyOUT(insideAOD,6));

hyOUT(hyOUT(:,6) > yBWlimit/2,6) = yBWlimit/2; 
hyOUT(hyOUT(:,6) < -yBWlimit/2,6) = -yBWlimit/2; 

% Calculate AM-induced power drop
ampliGauss = hyOUT(:,4) .* exp(-hyOUT(:,1).^2/rBeam^2);
hyPowDrop = 1 - ( sum(ampliGauss.^2) / sum(exp(-2*hyOUT(:,1).^2/rBeam^2)) );
hyPowTrans = 1 - hyPowDrop;

%% figures 
fig1 = figure;
set(fig1, 'Name', 'X-Hologramm', 'units', 'normalized', 'outerposition', [0 0 1 1]);

xhRange = find(insideAOD);
xfRange = find(abs(xf) < 0.05);

subplot(2,3,1)
plot(1000*xf(xfRange), fxTarget(xfRange), '-k'); axis tight;
title('Target Pattern');
xlabel('\mum');
subplot(2,3,2)
plot(hxOUT(xhRange,1), hxOUT(xhRange,2),'-k'); hold on; 
plot(hxOUT(xhRange,1), hxOUT(xhRange,4),'-r'); hold off;
axis tight;
xlabel('mm');
ylim([0,1.05]);
title('Holographic Amplitude');
subplot(2,3,3)
plot(hxOUT(xhRange,1), hxOUT(xhRange,6),'-r'); axis tight;
xlabel('mm');
title('AOD frequency (MHz)');
subplot(2,3,4)
plot(1000*xf(xfRange), fxStartPhase(xfRange), '-k'); axis tight;
ylim([-0.5,0.5]);
title('Start Phase (Focal Plane)');
xlabel('\mum');
subplot(2,3,5)
plot(hxOUT(xhRange,1), hxOUT(xhRange,3),'-k'); axis tight;
xlabel('mm');
ylim([-0.5,0.5]);
title('Holographic Phase (lambda)');
subplot(2,3,6)
plot(1000*xf(xfRange), fxAmplOUT(xfRange)/max(fxAmplOUT(xfRange)), '-k'); hold on;
%plot(1000*xf(xfRange), (fxAmplOUT(xfRange)/max(fxAmplOUT(xfRange))).^2, '-r');
%plot(1000*xf(xfRange), (fxAmplOUT(xfRange)/max(fxAmplOUT(xfRange))).^4, 'color', [0,0.7,0]);
axis tight;
title('Focal Amplitude');
xlabel('\mum')
%legend({'Amplitude'}, 'FontSize', 12, 'FontWeight', 'bold');
%legend('Location', 'northwest');
%legend('boxoff');

fig2 = figure;
set(fig2, 'Name', 'Y-Hologramm', 'units','normalized','outerposition',[0 0 1 1]);

yhRange = find(insideAOD);
yfRange = find(abs(yf) < 0.05);

subplot(2,3,1)
plot(1000*yf(yfRange), fyTarget(yfRange), '-k'); axis tight;
title('Target Pattern');
xlabel('\mum');
subplot(2,3,2)
plot(hyOUT(yhRange,1), hyOUT(yhRange,2),'-k'); hold on; 
plot(hyOUT(yhRange,1), hyOUT(yhRange,4),'-r'); hold off;
axis tight;
xlabel('mm');
ylim([0,1.05]);
title('Holographic Amplitude');
subplot(2,3,3)
plot(hyOUT(yhRange,1), hyOUT(yhRange,6),'-r'); axis tight;
xlabel('mm');
title('AOD frequency (MHz)');
subplot(2,3,4)
plot(1000*yf(yfRange), fyStartPhase(yfRange), '-k'); axis tight;
ylim([-0.5,0.5]);
title('Start Phase (Focal Plane)');
xlabel('\mum');
subplot(2,3,5)
plot(hyOUT(yhRange,1), hyOUT(yhRange,3),'-k'); axis tight;
ylim([-0.5,0.5]);
xlabel('mm');
title('Holographic Phase (lambda)');
subplot(2,3,6)
plot(1000*yf(yfRange), fyAmplOUT(yfRange), '-k'); hold on;
%plot(1000*yf(yfRange), (fyAmplOUT(yfRange)/max(fyAmplOUT(yfRange))).^2, '-r'); 
%plot(1000*yf(yfRange), (fyAmplOUT(yfRange)/max(fyAmplOUT(yfRange))).^4, 'Color', [0,0.7,0]);
axis tight;
title('Focal Amplitude');
xlabel('\mum')
%legend({'Amplitude'}, 'FontSize', 12, 'FontWeight', 'bold'); 
%legend('Location', 'northwest');
%legend('boxoff');

%%
save(fullfile(currentFolder, 'Hologram', 'GS holo plane X_data.mat'), 'hxOUT', '-mat');
outText = table({'xTarget'; 'xBW'; 'xBWlimit'; 'power drop'; 'lambda'}, {xTarget; xBW; xBWlimit; hxPowDrop; lambda*1e6});
writetable(outText, fullfile(currentFolder, 'Hologram', 'GS holo plane X_parameters.txt'), 'Delimiter', '\t', 'WriteVariableNames', false);

save(fullfile(currentFolder, 'Hologram', 'GS holo plane Y_data.mat'), 'hyOUT', '-mat');
outText = table({'yTarget'; 'yBW'; 'yBWlimit'; 'power drop'; 'lambda' }, {yTarget; yBW; yBWlimit; hyPowDrop; lambda*1e6});
writetable(outText, fullfile(currentFolder, 'Hologram', 'GS holo plane Y_parameters.txt'), 'Delimiter', '\t', 'WriteVariableNames', false);


