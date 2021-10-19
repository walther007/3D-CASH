function freq = phase2freq( x, phase )
% Convert phase into a AOD radio-frequency signal
% WALTHER Apr 2014
% freq [MHz]

vAOD = 650; % in [m/s] 

threshold = 0.5;

if ~isvector(x) || ~isvector(phase), error('unwrapPhase: WRONG INPUT'); end
if ~iscolumn(x), x = x'; end
if ~iscolumn(phase), phase = phase'; end

% define phase step
phaseStep = heaviside(x);
phaseStepDer = diff(phaseStep)/(x(end)-x(1)) * (length(x)-1);

phase = phase - median(phase);
phaseDer= diff(phase)/(x(end)-x(1)) * (length(x)-1);

% find positive phase jumps
stepInd = find(phaseDer > threshold * max(phaseStepDer));
phaseJump = false(length(x),1);
i=1;
while i <= length(stepInd)
    ii = 1;
    if i+ii < length(stepInd)
        while stepInd(i+ii) == stepInd(i)+ii
            ii = ii +1;
            if i+ii >= length(stepInd), break; end
        end
    end
    phaseJump(stepInd(floor(i+ii/2))) = true;
    i = i+ii;
end
posPhaseJumpInd = find(phaseJump);
    
% find negative phase jumps
stepInd = find(phaseDer < -threshold * max(phaseStepDer));
phaseJump = false(length(x),1);
i=1;
while i <= length(stepInd)
    ii = 1;
    if i +ii < length(stepInd)
        while stepInd(i+ii) == stepInd(i)+ii
            ii = ii +1;
            if i+ii == length(stepInd), break; end
        end
    end
    
    phaseJump(stepInd(floor(i+ii/2))) = true;
    i = i+ii;
end
negPhaseJumpInd = find(phaseJump);

% Correct phase derivative
for i = 1:length(posPhaseJumpInd)
    ind1 = max(1, round(posPhaseJumpInd(i)-3));
    ind2 = min(length(phase), round(posPhaseJumpInd(i)+3));
    phaseDiff = abs(phase(ind1) - phase(ind2));
    if phaseDiff >= 0.9
        %phase(ind2:end) = phase(ind2:end) - phaseDiff;
        phaseDer(ind1:ind2) = NaN;
    end       
end

for i = 1:length(negPhaseJumpInd)
    ind1 = max(1, round(negPhaseJumpInd(i)-1));
    ind2 = min(length(phase), round(negPhaseJumpInd(i)+1));
    phaseDiff = abs(phase(ind1) - phase(ind2));
    if phaseDiff >= 0.9
        %phase(ind2:end) = phase(ind2:end) + phaseDiff;
        phaseDer(ind1:ind2) = NaN;
    end
end   
% Continuation of the phase
phaseDer = interp1(x(isfinite(phaseDer)), phaseDer(isfinite(phaseDer)), x, 'pchip'); 
freq = phaseDer * vAOD / 1000;
end

