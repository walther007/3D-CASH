function [ hxAmplOUT, hxPhaseOUT, fxAmplOUT, fxStartPhase, fxTarget] = GerchbergSaxton4(hxAmplIN, xf, xTarget, xTargetAmpl, psf, dynTarget, iterMax, FoMtarget, smoothKernel, backwards, fixInit, valInit)
% Calculation of the phase + amplitude distribution in the holographic plane by iterative FFT (Gerchberg-Saxton)
% WALTHER, Apr 2014

switch nargin
    case {1,2,3,4}
        error('GERCHBERSAXTON: input error!');
    case 5
        dynTarget = false;
        iterMax = 100;
        FoMtarget = 0.99;
        smoothKernel = [];
        backwards = true;
        fixInit = true;
        valInit = 1000;
    case 6
        iterMax = 100;
        FoMtarget = 0.99;
        smoothKernel = [];
        backwards = true;
        fixInit = true;
        valInit = 1000;
    case 7
        FoMtarget = 0.99;
        smoothKernel = [];
        backwards = true;
        fixInit = true;
        valInit = 1000;
    case 8
        smoothKernel = [];
        backwards = true;
        fixInit = true;
        valInit = 1000;
    case 9
        backwards = true;
        fixInit = true;
        valInit = 1000;
    case 10
        fixInit = true;
        valInit = 1000;
    case 11
        valInit = 1000;
end

if ~isequal(length(hxAmplIN), length(xf), length(psf)), error('Wrong Input'); end
    
if isequal(mod(smoothKernel,2),0), smoothKernel = smoothKernel + 1; end

nx = length(hxAmplIN);

if fixInit
    rng(valInit);   % *** fix rnd seed
else
    rng('shuffle');
end
rngInit = rng;
disp(['rng = ' num2str(rngInit.Seed)]);

fxStartPhase = -0.1+0.2*rand(length(xf),1);

nSeed = 2*ceil(sum(xf>=0 & xf <=1.2e-3)/2);

fxTargetInd = zeros(1,length(xTarget));
fxTarget = zeros(length(xf),1);
for i = 1:length(xTarget)
    fxTargetInd(i) = find(xf>=xTarget(i),1,'first');
    if isempty(xTargetAmpl)
        fxTarget(fxTargetInd(i)) = 1;
    else
        fxTarget(fxTargetInd(i)) = xTargetAmpl(i);
    end
    fxStartPhase(fxTargetInd(i)-nSeed/2+(0:nSeed+1)) = -0.25+0.5*rand;
end

%fxTargetINIT = conv(fxTarget, psf, 'same');
fxTargetINIT = fxTarget;

fxPhase = fxStartPhase * 2*pi;

fxTarget = conv(fxTarget, psf, 'same');
fxTarget = fxTarget / max(fxTarget);

if backwards
    for i=1:iterMax
        % backpropagate to holographic plane
        spec = fxTarget .* exp(1i * fxPhase);
        field = fftshift(ifft(ifftshift(spec), nx));
        hxPhase = angle(field) - mean(angle(field));
        %if addPhase,  hxPhase = hxPhase + suppPhase; end
        hxAmpl = abs(field);

        if ~isempty(smoothKernel), hxPhase = smooth(hxPhase, smoothKernel); end
        % propagate to Fourier plane
        field = hxAmplIN .* exp(1i * hxPhase);
        spec = fftshift(fft(ifftshift(field),nx))/nx;
        fxPhase = angle(spec);
        
        fxAmpl = abs(spec)/max(abs(spec));

        % compare with Target pattern
        FoM = xcorr(fxAmpl, fxTargetINIT, 0) / sqrt(xcorr(fxAmpl,0) * xcorr(fxTargetINIT,0));
        if i==1, disp('**** GS phase optimization ****'); 
            disp(['Start target Corr. = ', num2str(FoM)]);
        elseif i<20 || ~ mod(i,20)
            disp(['Iteration # ', num2str(i), ' Target Corr. = ', num2str(FoM)]);
        end

        if FoM >= FoMtarget, break; end
        
        if dynTarget
            %fxTarget(fxTargetInd) = fxTarget(fxTargetInd) - 0.5* (fxAmpl(fxTargetInd) - mean(fxAmpl(fxTargetInd)));
            fxTarget(fxTargetInd) = fxTarget(fxTargetInd) - 0.5* (fxAmpl(fxTargetInd) - fxTargetINIT(fxTargetInd));
            fxTarget = fxTarget/max(fxTarget);
        end
       
        %plot(fxTarget); pause(1);
         %if addPhase,  fxPhase = fxPhase + suppPhase; end
    end
else
   hxPhase = zeros(nx,1);
   for i=1:iterMax
        %if addPhase,  hxPhase = hxPhase + suppPhase; end
        if ~isempty(smoothKernel), hxPhase = smooth(hxPhase, smoothKernel); end
        % propagate to Fourier plane
        
        field = hxAmplIN .* exp(1i * hxPhase);
        spec = fftshift(fft(ifftshift(field),nx))/nx;
        fxPhase = angle(spec);
        
        fxAmpl = abs(spec) / max(abs(spec));
        
        if dynTarget
            fxTarget(fxTargetInd) = fxTarget(fxTargetInd) - 0.5* (fxAmpl(fxTargetInd) - mean(fxAmpl(fxTargetInd)));
        end
        
        % compare with Target pattern
        FoM = xcorr(fxAmpl, fxTargetINIT, 0) / sqrt(xcorr(fxAmpl,0) * xcorr(fxTargetINIT,0));
        if i==1, disp('**** GS phase optimization ****'); 
            disp(['Start target Corr. = ', num2str(FoM)]);
        elseif i<20 || ~ mod(i,20)
            disp(['Iteration # ', num2str(i), ' Target Corr. = ', num2str(FoM)]);
        end
        
        % backpropagate to holographic plane
        spec = fxTarget .* exp(1i * fxPhase);
        field = fftshift(ifft(ifftshift(spec), nx));
        hxPhase = angle(field) - mean(angle(field));
        hxAmpl = abs(field);
        
        if FoM >= FoMtarget, break; end
        %if addPhase,  fxPhase = fxPhase + suppPhase; end
   end
end

hxPhaseOUT = hxPhase / (2*pi);
hxAmplOUT = hxAmpl;

% final approximation error
field = hxAmplOUT .* exp(1i * 2*pi* hxPhaseOUT);
spec = fftshift(fft(ifftshift(field), nx))/nx;
FoM = xcorr(abs(spec)/max(abs(spec)), fxTarget, 0) / sqrt(xcorr(abs(spec)/max(abs(spec)),0) * xcorr(fxTarget,0)); 
disp(['Final Corr Error = ', num2str(1-FoM)]);
fxAmplOUT = abs(spec);
end

