function result = RegimeDet(rdPars, times, vars, DEqns_comp, DEpars, varBounds, isCycle, ...
    algOptions, verboseTog, silentRegimes, figHandle)
% Reduced regime determination using epochs and scale thresholds.
% Version 1.18  (c) Robert Clewley, Department of Mathetmatics, Cornell
%   University, February 2006
% For use by DSSRT system.

%%%%%%%%%%%%%%   DEVELOPMENT NOTES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 1.18
%   Bug fix for inverted fast/slow labelling in verbose output
%
% Version 1.17
%   Minor bug fixes, esp. dealing with cycles consisting of a single epoch
%    (which turn out not to have any speed changes, so that the epoch becomes a
%    whole regime).
%   Now compatible with gap junctions in ODE system.
%
% Version 1.16
%   Added speedFussyLevel parameter to give two levels of fussiness
%   Made epoch time scale and derivative checks "look-ahead" to *next* tPos so that any new epoch is
%    started at the current tPos.
%   New epochs due to the above checked changes now force a new regime at
%   these points w/ fussy level 2.
%
% Version 1.15
%   option added to prevent long epoch durations from changing `slow` vars to `normal`
%   speedFussy option properly implemented, to search all data points for time scale status changes
%     (useful for orbits in which active set does not change but time scales do e.g. Michaelis-Menten)
%
%%%%%%%%%%%%%%   OTHER NOTES    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% silentRegimes is for automatic searching of regimes (like AttEst) ... not yet implemented
%
%% Content fields of `rdPars`:
% TS, Eseq, startTime, focusSet, tScaleThresh, varChangeThresh, varnames, DEixMap, caIxMap, inputsIx,
%   numExt, numInt   (NB caIxMap is candActsIxMapAbs)
%
%% Content fields of `algOptions`:
% speedFussy, fastLeaveFussy, longFussy
%
%% Content fields of `regimes` (the returned structure):
% epochs, timeInt, dynVars, dynCross, nonDynVars, fastVars, slowVars, dimensionA, dimensionB,
%   dimensionC, qsPars
%
% Use of speedChange.flag assumes that only one variable will change timescale status at a time,
%   otherwise only the last variable to be checked will be correctly updated. Therefore, a small
%   enough time-resolution of data is assumed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = {0, []}; % default (error) value

if nargin ~= 11
    disp('RegimeDet:  Wrong number of arguments passed. Expected 11')
    return
end

%% Load global constants
GlobalConstants % this line will be replaced by the contents of this script file when compiled (see DSSRT.m)

numTot = rdPars.numExt + rdPars.numInt;
numSamples = 5; % for near-constant variable status checking over a regime

% Time step dt is according to assumed FIXED step size from integrated data
dt = (round(times(2)/TIMERES)-round(times(1)/TIMERES))*TIMERES; % accurate to TIMERES resolution

if silentRegimes
    verboseTog=0; % force this off
end

%% Decide passive variables (ones without associated differential equations)
passiveVars = find(rdPars.DEixMap==0); % in range 0..numTot ( = numExt + numInt )

if isempty(rdPars.focusSet)
    disp('RegimeDet:  Error - transition sequence not focused on an observable variable')
    return
else
    primaryFocVarIx = rdPars.focusSet(1);
end
numVarsFocused = length(rdPars.focusSet);

primaryFocEqnIx = rdPars.DEixMap(primaryFocVarIx);
if primaryFocEqnIx == 0
    disp('RegimeDet:  Internal error - primary focused variable does not have an associated differential equation')
    return
end
primaryFocDE = DEqns_comp{primaryFocEqnIx};
Gamma1       = primaryFocDE{DE_GAMMA1i};

%% Set of ix's for search and DE list for those
ixSet = setdiff(rdPars.caIxMap{primaryFocVarIx},passiveVars);
DEset = cell(1,length(ixSet));
for ix=ixSet
    thisDEix  = rdPars.DEixMap(ix);
    DEset{ix} = DEqns_comp{thisDEix};
end

% Prepare derivative thresholds for deciding resolution of bolster epochs
varMaxIntervals = zeros(1,numTot);
for varIx = 1:numTot
    varMaxIntervals(varIx) = varBounds(varIx,2) - varBounds(varIx,1);
end

drawnow

numEps_orig = length(rdPars.TS{TS_TSEQ});

%% Bolster regimes with additional epochs for changing timescales, if speedFussy option set
if algOptions.speedFussy
    loopTS = rdPars.TS; % model it on the original!
    newTS = {};
    disp('RegimeDet:  Beginning exhaustive search for time scale status changes')
    disp('             to create `bolster` epochs')
    if isCycle & numEps_orig ~= 1
        % for single epoch in a cycle the end TS is not repeated!
        endEpIx = length(rdPars.TS{TS_TSEQ})-1;
    else
        endEpIx = length(rdPars.TS{TS_TSEQ});
    end
    newEpoch = false; % initial value
    newTSlist = []; % initial value for list of TS #'s that are added because of time scale or derivative changes
    for epIx = 1:endEpIx
        % get time pos limits, etc. for this original epoch
        thisEpoch   = rdPars.TS{TS_TSEQ}{epIx};
        tlo         = thisEpoch{TSEQ_TIME} + rdPars.startTime;
        if (isCycle || epIx < endEpIx) && numEps_orig ~= 1
            thi     = rdPars.TS{TS_TSEQ}{epIx+1}{TSEQ_TIME} + rdPars.startTime; % epIx+1 is always well-defined here
        else % epIx == endEpIx & ~isCycle | numEps_orig == 1
            thi     = times(rdPars.TS{TS_TSEQ}{1}{TSEQ_POSN} + rdPars.TS{TS_ENDP});
        end
        [isp poslo] = ismember(1,times <= tlo);
        [isp poshi] = ismember(1,times <= thi);
        newTS       = {newTS{:}, { thisEpoch{TSEQ_ACTS}, tlo-rdPars.startTime, poslo } }; % initialize this epoch
        tData       = times(poslo:poshi);
        vData       = vars(poslo:poshi,:);
        acts        = rdPars.caIxMap{primaryFocVarIx}( thisEpoch{TSEQ_ACTS}{primaryFocVarIx} ); % in absolute indices
        pots        = rdPars.caIxMap{primaryFocVarIx}( rdPars.Eseq{poslo}{ESEQ_POTENTS_BAR}{primaryFocVarIx} ); % in absolute indices
        nonDynVars  = intersect(acts,passiveVars);
        dynVars     = setdiff(acts,nonDynVars);
        dynCross    = []; % initial value
        
%        poslist = poslo; % no longer used

        % see if there are any non-primary (internal) variables in active set making up non-simple input terms to the ODE
        for g1t = 1:length(Gamma1)
            intVarIx = Gamma1{g1t}{DE_INTVARNAMEi};
            extVarIx = Gamma1{g1t}{DE_GAMVARNAMEi};
            if ismember(extVarIx,dynVars) && (intVarIx > 0 & ~ismember(intVarIx,passiveVars)) % then non-passive internal cross-multiplying variable specified
                dynCross = [dynCross, intVarIx];
            end
        end
        
        % go through set of vars in these limits
        velocityStatus.prev = 0;
        velocityStatus.curr = 0;
        varDerivIx = -1; % dummy value
        derivDir = 0; % dummy value
        for tPos = poslo:poshi-1
            if tPos > poslo
                speeds.curr = speeds.next; % from previous loop step
                vLine.curr  = vLine.next; % from previous loop step
            else
                speeds.curr  = zeros(1,numTot);
                vLine.curr   = vars( tPos, : );
                for ix = ixSet
                    speedResult = isAbnormalSpeed(DEset{ix},primaryFocDE,vLine.curr,rdPars.tScaleThresh);
                    speeds.curr(ix)  = speedResult(1);
                end
            end
            speeds.next  = zeros(1,numTot);            
            vLine.next   = vars( tPos+1, : );
            for ix = ixSet
                speedResult = isAbnormalSpeed(DEset{ix},primaryFocDE,vLine.next,rdPars.tScaleThresh);
                speeds.next(ix)  = speedResult(1);
            end

            for var=dynCross % go through internal variables cross-multiplying primary variables
                thisDE      = DEqns_comp{ rdPars.DEixMap(var) };
                if tPos == poslo
                    speedResult = isAbnormalSpeed(thisDE,primaryFocDE,vLine.curr,rdPars.tScaleThresh);
                    speeds.curr(var) = speedResult(1);
                end
                speedResult = isAbnormalSpeed(thisDE,primaryFocDE,vLine.next,rdPars.tScaleThresh);
                speeds.next(var) = speedResult(1);
            end
            % later, implement this check for Gamma2 set, when they are allowed cross-mult variables too

            % change in approximate time derivatives
            dVarRel = zeros(1,numTot);
            varIxList = [dynVars,primaryFocVarIx];
            for varIx = varIxList
                dVarRel(varIx) = abs((vLine.curr(varIx) - vLine.next(varIx))/varMaxIntervals(varIx));
            end
            % test size of derivatives
            [maxval maxix] = max(dVarRel);
            if maxval >= dt*rdPars.derivThresh
                velocityStatus.curr = 1;
                if velocityStatus.prev ~= 1
                    newEpoch = true;
                    varDerivIx = maxix;
                    derivDir = 1; % speeding up
                    if verboseTog
                        fprintf('RegimeDet:  Added epoch at var %s`s change to HIGH velocity at time %.3f\n',rdPars.varnames{varDerivIx},times(tPos))
                    end
                end
            else
                velocityStatus.curr = -1;
                if velocityStatus.prev == 1
                    newEpoch = true;
                    derivDir = -1; % slowing down
                    if verboseTog
                        fprintf('RegimeDet:  Added epoch at a var`s velocity change to LOW at time %.3f\n',times(tPos))
                    end
                end
            end
            
            % change in time scale status
            if sum(speeds.curr == speeds.next) ~= numTot % then there exists mismatch
                newEpoch = true;
                if verboseTog
                    fprintf('RegimeDet:  Added epoch due to time scale change at time %.3f\n',times(tPos))
                end
            end
            
            if newEpoch
                if newTS{length(newTS)}{3} == tPos % so that don't add a new epoch onto the same, already existent, one
                    newTS   = {newTS{1:length(newTS)-1}, { thisEpoch{TSEQ_ACTS}, times(tPos)-rdPars.startTime, tPos } }; % add new epoch
                else
                    newTS   = {newTS{:}, { thisEpoch{TSEQ_ACTS}, times(tPos)-rdPars.startTime, tPos } }; % add new epoch
%                     poslist = [poslist, tPos]; % interval start positions (in Eseq times)
                end
                newTSlist = [newTSlist; [length(newTS),varDerivIx,derivDir]]; % record of the TS # for which there has been a time scale or derivative change
                newEpoch = false; % reset
                if mod(length(newTSlist),10) == 1
                    drawnow % refresh screen at beginning and then if taking a long time
                    if ~verboseTog && length(newTSlist) > 1 % let user know that program is working!
                        disp('            ... still calculating `bolster` epochs')
                    end
                end
            end
            velocityStatus.prev = velocityStatus.curr;
        end
    end
    if verboseTog
        disp('RegimeDet:  Search completed')
    end
    % insert new set of epochs with events at time-speed status changes
    if isCycle && numEps_orig ~= 1
        newTS = { newTS{:}, rdPars.TS{TS_TSEQ}{endEpIx+1} };
    end
    loopTS{TS_TSEQ} = newTS;
else
    loopTS = rdPars.TS;
end

%% Prepare for main loop
if isCycle & numEps_orig ~= 1
    % if numEps_orig ~= 1 then assume we were given a single epoch for a limit
    % cycle, so that the TS would not have the repeated end entry at the
    % end of the period
    numEps = length(loopTS{TS_TSEQ})-1; % discount final one for isCycle (the one which is there to show the end time, etc.)
else
    numEps = length(loopTS{TS_TSEQ});
end
regimeNum = 0;
firstEpochFlag = false;
lastEpochFlag = false;
doneEpochs = 0;
speedChange.flag = false;
const.prev = zeros(1,numTot);

if verboseTog
    disp(' ') % leave extra space
    disp(['RegimeDet:  Time scales shown relative to that of primary focus variable ' rdPars.varnames{primaryFocVarIx}]);
    fprintf('              using time scale threshold 1/%.4f\n',rdPars.tScaleThresh)
    drawnow
end
interruptCounter = 0;

if isCycle
    endpos = numEps*2-1; % go an extra time around in case need to cover part-started regimes at beginning
else
    endpos = numEps;
    global_period = 0; % won't use this, just set it to zero for the regime structure
end

disp('RegimeDet:  Determining reduced dynamical regimes')
%% Step through epochs (more than once if start mid-way through a regime)
for pos = 1:endpos
    newRegime = false; % reset this for this epoch (initialize for pos=1)
    if verboseTog | ~silentRegimes
        interruptCounter = interruptCounter + 1; % check interrupt status every 4 epoch steps
        if interruptCounter == 4
            interruptCounter = 0;
            drawnow % only refresh screen if we'll be slowing down this loop with console display
            if get(figHandle,'CurrentCharacter') == 27 % escape key
                beep
                disp('RegimeDet:  Reduced regime determination cancelled by user')
                return
            end
        end
    end
    epochPos = mod(pos-1,numEps)+1;
    if verboseTog
        if pos > numEps % then isCycle and going round again
            fprintf('RegimeDet:  Epoch %i/%i second time around cycle ',epochPos,numEps)
        else
            fprintf('RegimeDet:  Epoch %i/%i ',epochPos,numEps)
        end
    end
    if isCycle
        prevEpochPos = mod(epochPos-2,numEps)+1;
        nextEpochPos = mod(epochPos,numEps)+1;
    else
        prevEpochPos = epochPos-1;
        nextEpochPos = epochPos + 1;
    end
    epochTS.curr = loopTS{TS_TSEQ}{epochPos};
    if ~isCycle && prevEpochPos == 0 || numEps == 1
        firstEpochFlag = true; % the current epoch is the first in the set
    else
        firstEpochFlag = false;
    end
    if ~isCycle && nextEpochPos == numEps+1
        lastEpochFlag = true; % the current epoch is the last in the set
    end % since this would be set in the last epoch there's no need to reset this on other times round loop

    %%%%% Info for the current epoch
    tlo.curr         = epochTS.curr{TSEQ_TIME}+rdPars.startTime;
    if (isCycle || epochPos < numEps) && numEps_orig ~= 1
        thi.curr     = loopTS{TS_TSEQ}{epochPos+1}{TSEQ_TIME}+rdPars.startTime; % epochPos+1 is always well-defined here
    else % (epochPos == numEps & ~isCycle) | numEps == 1
        thi.curr     = times(loopTS{TS_TSEQ}{1}{TSEQ_POSN} + loopTS{TS_ENDP} - 1);
    end
    [isp poslo.curr] = ismember(1,times <= tlo.curr);
    [isp poshi.curr] = ismember(1,times <= thi.curr);
    tData.curr       = times(poslo.curr:poshi.curr);
    vData.curr       = vars(poslo.curr:poshi.curr,:);
    acts.curr        = rdPars.caIxMap{primaryFocVarIx}( epochTS.curr{TSEQ_ACTS}{primaryFocVarIx} ); % in absolute indices
    pots.curr        = rdPars.caIxMap{primaryFocVarIx}( rdPars.Eseq{poslo.curr}{ESEQ_POTENTS_BAR}{primaryFocVarIx} ); % in absolute indices
    nonDynVars.curr  = intersect(acts.curr,passiveVars);
    dynVars.curr     = setdiff(acts.curr,nonDynVars.curr);
    dynCross.curr    = []; % initial value
    % see if there are any non-primary (internal) variables in active set making up non-simple input terms to the ODE
    for g1t = 1:length(Gamma1)
        intVarIx = Gamma1{g1t}{DE_INTVARNAMEi};
        extVarIx = Gamma1{g1t}{DE_GAMVARNAMEi};
        if ismember(extVarIx,dynVars.curr) && (intVarIx > 0 & ~ismember(intVarIx,passiveVars)) % then non-passive internal cross-multiplying variable specified
            dynCross.curr = [dynCross.curr, intVarIx];
        end
    end
    % later, implement this check for Gamma2 set, when they are allowed cross-mult variables too

    if verboseTog
        dispStr = '';
        for var=acts.curr
            dispStr = [dispStr ' ' rdPars.varnames{var}];
        end
        disp(['-- all actives: [' dispStr ' ]'])
    end

    if verboseTog && algOptions.speedFussy
        fprintf('            t = [ %.3f, %.3f ),  + [ %.3f, %.3f ) rel.\n',tlo.curr,thi.curr,tlo.curr-rdPars.startTime,thi.curr-rdPars.startTime)
    end

    % Work out relative timescales:
    % With option speedFussy = false, we are making the assumption that the timescale of a variable does not change
    %   significantly over an epoch. For each dynamic var, check its relative speed to primaryFocVar in 'middle' of epoch.
    % However, if speedFussy is selected, we've already separated epochs to check for speed status changes. After that, checking
    %   in middle is also safe...
    midPoint.curr    = round((poshi.curr+poslo.curr)/2);
    if midPoint.curr >= poshi.curr | midPoint.curr < poslo.curr
        midPoint.curr = poslo.curr;
    end
    vLine.curr       = vars( midPoint.curr, : ); % representative variable data line in 'middle' of epoch
    speeds.curr      = zeros(1,numTot);
    if verboseTog
        disp('            Time scale ratios of focused variable timescale for dynamic variables:')
    end
    for ix=ixSet
        speedResult = isAbnormalSpeed(DEset{ix},primaryFocDE,vLine.curr,rdPars.tScaleThresh);
        speeds.curr(ix) = speedResult(1);
        if verboseTog
            switch speedResult(1)
                case -1
                    dispStr  = '  ';
                    speedStr = 'slow';
                case 0
                    if speedResult(3) % faster
                        dispStr  = '1/';
                    else
                        dispStr  = '  ';
                    end
                    speedStr = 'normal speed';
                case 1
                    dispStr  = '1/';
                    speedStr = 'fast';
            end
            fprintf('               %s: %s%.4f ( %s )\n',rdPars.varnames{ix},dispStr,speedResult(2),speedStr)
        end
    end
    for var=dynCross.curr % go through internal variables cross-multiplying primary variables
        thisDE      = DEqns_comp{ rdPars.DEixMap(var) };
        speedResult = isAbnormalSpeed(thisDE,primaryFocDE,vLine.curr,rdPars.tScaleThresh);
        speeds.curr(var) = speedResult(1);
        if verboseTog
            switch speedResult(1)
                case -1
                    dispStr  = '  ';
                    speedStr = 'slow';
                case 0
                    if speedResult(3) % faster
                        dispStr  = '1/';
                    else
                        dispStr  = '  ';
                    end
                    speedStr = 'normal speed';
                case 1
                    dispStr  = '1/';
                    speedStr = 'fast';
            end
            fprintf('               %s: %s%.4f ( %s )\n',rdPars.varnames{var},dispStr,speedResult(2),speedStr)
        end
    end
    [p i]    = find(speeds.curr==-1); % p unused
    slowVars = intersect( union(dynVars.curr,dynCross.curr), i);
    [p i]    = find(speeds.curr== 1); % p unused
    fastVars = intersect( union(dynVars.curr,dynCross.curr), i);

    % only calc global_period first time around, so use pos not epochPos
    if isCycle && pos == 1
        tlo_global    = tlo.curr;
    end
    if isCycle && pos == numEps
        thi_global    = thi.curr;
        global_period = thi_global - tlo_global;
    end

    %%%%% Info for previous epoch, if available
    if ~firstEpochFlag
        epochTS.prev     = loopTS{TS_TSEQ}{prevEpochPos};
        tlo.prev         = epochTS.prev{TSEQ_TIME}+rdPars.startTime;
        thi.prev         = tlo.curr;
        [isp poslo.prev] = ismember(1,times <= tlo.prev);
        [isp poshi.prev] = ismember(1,times <= thi.prev);
        tData.prev       = times(poslo.prev:poshi.prev);
        vData.prev       = vars(poslo.prev:poshi.prev,:);
        acts.prev        = rdPars.caIxMap{primaryFocVarIx}( epochTS.prev{TSEQ_ACTS}{primaryFocVarIx} );
        if poslo.curr == 1
            poslo_prev = poshi.prev;
        else
            poslo_prev = poslo.curr - 1;
        end
        pots.prev        = rdPars.caIxMap{primaryFocVarIx}( rdPars.Eseq{poslo_prev}{ESEQ_POTENTS_BAR}{primaryFocVarIx} );
        nonDynVars.prev  = intersect(acts.prev,passiveVars);
        dynVars.prev     = setdiff(acts.prev,nonDynVars.prev);
        dynCross.prev    = []; % initial value
        % see if there are any non-primary (internal) variables in active set making up non-simple input terms to the ODE
        for g1t = 1:length(Gamma1)
            intVarIx = Gamma1{g1t}{DE_INTVARNAMEi};
            extVarIx = Gamma1{g1t}{DE_GAMVARNAMEi};
            if ismember(extVarIx,dynVars.prev) && (intVarIx > 0 & ~ismember(intVarIx,passiveVars)) % then non-passive internal cross-multiplying variable specified
                dynCross.prev = [dynCross.prev, intVarIx];
            end
            % later, implement this check for Gamma2 set, when they are allowed cross-mult variables too
        end
        % Work out relative timescales
        speeds.prev      = zeros(1,numTot);
        midPoint.prev    = round((poshi.prev+poslo.prev)/2);
        if midPoint.prev >= poshi.prev | midPoint.prev < poslo.prev
            midPoint.prev = poslo.prev;
        end
        vLine.prev       = vars( midPoint.prev, :); % representative variable data line in 'middle' of epoch
        for ix=setdiff(rdPars.caIxMap{primaryFocVarIx},passiveVars)
            thisDEix     = rdPars.DEixMap(ix);
            thisDE       = DEqns_comp{thisDEix};
            speedResult  = isAbnormalSpeed(thisDE,primaryFocDE,vLine.prev,rdPars.tScaleThresh);
            speeds.prev(ix) = speedResult(1);
        end
        for var=dynCross.prev % go through internal variables cross-multiplying primary variables
            thisDEix     = rdPars.DEixMap(var);
            thisDE       = DEqns_comp{thisDEix};
            speedResult  = isAbnormalSpeed(thisDE,primaryFocDE,vLine.prev,rdPars.tScaleThresh);
            speeds.prev(var) = speedResult(1);
        end
        if speedChange.flag && ~speedChange.newreg % then on prev epoch, a speed was changed, and other speedChange fields are defined
            speeds.prev(speedChange.var) = speedChange.dest; % must update this to be consistent
        end
    end

    %%%%% Info for next epoch, if available
    if ~lastEpochFlag
        epochTS.next     = loopTS{TS_TSEQ}{nextEpochPos};
        acts.next        = rdPars.caIxMap{primaryFocVarIx}( epochTS.next{TSEQ_ACTS}{primaryFocVarIx} );
        nonDynVars.next  = intersect(acts.next,passiveVars);
        dynVars.next     = setdiff(acts.next,nonDynVars.next);
        qsParsCand       = setdiff(dynVars.next,dynVars.curr); % forwards in time
        qsPars           = [];
        for var=pots.curr
            if ismember(var,qsParsCand)
                qsPars = [qsPars, var];
            end
        end
    else
        qsPars = pots.curr; % have no more information to help us here
    end

    constVars  = []; % initialize these here for both cycle and non-cycle types of orbit
    const.curr = zeros(1,numTot);
    speedChange.flag = false;
    speedChange.num  = 0;
    speedChange.var  = [];
    speedChange.source = [];
    speedChange.dest   = [];
    speedChange.newreg = false;
    leavers = [];
    joiners = [];
    stayers = [];

    %%%%% Do 'new regime?' tests
    if ~firstEpochFlag
        leavers = setdiff(dynVars.prev,dynVars.curr); % forwards in time, previous->current
        joiners = setdiff(dynVars.curr,dynVars.prev); % forwards in time, previous->current
        stayers = union( intersect(dynVars.curr,dynVars.prev) , intersect(dynCross.curr,dynCross.prev) ); % includes assoc'd cross-multiplying non-passive internals

        % test 1 -- if a non-passive variable joins
        if ~isempty(joiners)
            newRegime = true;
            if verboseTog
                dispStr = '';
                for var=joiners
                    dispStr = [dispStr rdPars.varnames{var} ' '];
                end
                disp(['            * new regime due to joiner(s) ' dispStr])
            end
        end

        % test 2 -- if a non-passive variable changes its speed status
        for var=stayers % includes checking time scales of cross-multiplying vars assoc'd with active observables
            speeds_temp = [speeds.prev(var), speeds.curr(var)];
            if speeds_temp == [-2, 0] % then 'near-constant' -> `undetermined` (not necessarily `normal`, because haven't determined these for this epoch yet), so skip for now
                continue
            end % only other alternatives are [-2, -1] or [-2, 1] for previous near-constant variable ...
            if speeds_temp(1) ~= speeds_temp(2)  &&  ( (newRegime & verboseTog) | ~newRegime | algOptions.speedFussy )
                if algOptions.speedFussy && algOptions.speedFussyLevel >= 1
                    newRegime = true;
                    if verboseTog
                        disp('            * new regime because speed change and speed fussiness option ON ...')
                    end
                end
                if verboseTog
                    switch speeds_temp(1)
                        case -2
                            speedStr.prev = 'near-constant';
                        case -1
                            speedStr.prev = 'slow';
                        case 0
                            speedStr.prev = 'normal speed';
                        case 1
                            speedStr.prev = 'fast';
                    end
                    switch speeds_temp(2)
                        case -1
                            speedStr.curr = 'slow';
                        case 0
                            speedStr.curr = 'normal speed';
                        case 1
                            speedStr.curr = 'fast';
                    end
                end
                if ~newRegime % o/w already changing regime so allow variable speed to change
                    if verboseTog & speeds_temp(2) ~= 0
                        speedStr.curr = [speedStr.curr ' (but forcing to normal speed)'];
                    end
                    speedChange.flag   = true;
                    speedChange.num    = speedChange.num + 1;
                    speedChange.var(speedChange.num)    = var;
                    speedChange.source(speedChange.num) = speeds_temp(1);
                    speedChange.dest(speedChange.num)   = 0; % force to normal speed
                    speedChange.newreg = false;
                end
                if verboseTog
                    disp(['            - ' rdPars.varnames{var} ' speed status ' speedStr.prev ' -> ' speedStr.curr])
                end
            end
        end % for var=stayers

        % test 2.5 -- check for high/low time derivative changes
        if algOptions.speedFussy && algOptions.speedFussyLevel == 2 % only force new regimes when fussy level is 2
            [isPres presIx] = ismember(epochPos,newTSlist(:,1));
            if isPres && newTSlist(presIx,3) ~= 0 % then small derivative test made a new epoch, so must force a new regime here
                newRegime = true;
	%             speedChange.flag = true;
	%             speedChange.num  = speedChange.num + 1;
	%             speedChange.var(speedChange.num)  = newTSlist(presIx,2);
                var = newTSlist(presIx,2);
                if newTSlist(presIx,3) == 1 % speeding up
	%                 speedChange.source(speedChange.num) = 0; % normal (or whatever)
	%                 speedChange.dest(speedChange.num)   = 1; % fast
                    if verboseTog
                        disp(['            - ' rdPars.varnames{var} ' speed status -> fast (high time derivative)'])
                    end
                else % slowing down
	%                 speedChange.source(speedChange.num)   = 1; % fast
	%                 speedChange.dest(speedChange.num)   =   0; % normal (or whatever)
                    if verboseTog
                        disp(['            - ' rdPars.varnames{var} ' speed status no longer fast (low time derivative)'])
                    end
                end
                speedChange.newreg = false; % don't set this to true here otherwise the whole reason for the new regime will not be recorded in the speed change
            end
        end
        
        
        % test 3 -- don't start a new regime when a var leaves if there are neither any dynamic variables nor qsPars
        %           (including non-near-constant slow variables that will become classified as qsPars at the end of the regime)
        if (regimeNum > 0 & ~isempty(leavers)) && (~isempty(regimes(regimeNum).dynVars) | ~isempty(regimes(regimeNum).qsPars) ...
                | ~isempty(setdiff(regimes(regimeNum).slowVars,regimes(regimeNum).constVars)))
            if algOptions.fastLeaveFussy % don't start new regime if only fast vars are leaving
                if ~isempty(setdiff(leavers,regimes(regimeNum).fastVars))
                    newRegime = true;
                    optStr = ' due to fussiness over fast leaving variables option being ON';
                end
            else % start new regime when any fast vars leave
                newRegime = true;
                optStr    = '';
            end
            if verboseTog && newRegime
                dispStr = '';
                for var=leavers
                    dispStr = [dispStr ' ' rdPars.varnames{var}];
                end
                disp(['            * new regime due to leaver(s)' dispStr optStr])
            end
        end
        if newRegime & speedChange.newreg & speedChange.flag
            speedChange.flag = false; % don't bother with any forced speed status changes for some occasions of "new regime"
            if verboseTog
                disp('            - cancelled speed changes due to new regime')
            end
        end

        % test 4 -- if a variable is staying almost constant, not because it's slow in the above sense, but because it's close to its instantaneous target,
        %   then make it 'effectively' slow --> `near-constant`
        newRegime_temp = false; % temp holder for updated `newRegime` flag. Want to know if previous tests set this throughout the for loop before updating it here.
        for varIx = [setdiff(rdPars.inputsIx{primaryFocVarIx},passiveVars),primaryFocVarIx]
            if ismember(varIx,union(dynVars.curr,dynCross.curr)) | varIx == primaryFocVarIx % var must be intrinsically dynamic and active or be the primary focused variable
                near_const = true; % initial value only
                if regimeNum > 0 & ~newRegime
                    [isp poslo_reg] = ismember(1,times <= regimes(regimeNum).timeInt(1));
                else
                    poslo_reg = poslo.curr;
                end
                poshi_reg = poshi.curr;
                if poshi_reg > poslo_reg
                    posIntervals = {[poslo_reg, poshi_reg]};
                else % regime for cyclic data is split over the cycle ends
                    [isp poshi_reg] = ismember(1,times <= regimes(regimeNum).timeInt(2));
                    posIntervals = {[poslo_reg,poshi_reg], [poslo.curr,poshi.curr]};
                end
                numIntervals = length(posIntervals);
                for intIx=1:numIntervals
                    posStep = floor((posIntervals{intIx}(2)-posIntervals{intIx}(1)) / numSamples);
                    if posStep == 0
                        numSamples_temp = 0;
                        near_const = false; % don't consider such a small regime
                    else
                        numSamples_temp = numSamples;
                    end
                    vLine_init = vars(posIntervals{intIx}(1),:);
                    for sample=1:numSamples_temp % if zero then will skip this
                        vLine_temp = vars(posIntervals{intIx}(1)+sample*posStep,:);
                        if abs(vLine_temp(varIx) - vLine_init(varIx)) > rdPars.varChangeThresh * abs(varBounds(varIx,2)-varBounds(varIx,1))
                            near_const = false;
                            break
                        end
                    end
                    if ~near_const % don't bother searching in next interval
                        break
                    end
                end
                if varIx == primaryFocVarIx && ~algOptions.speedFussy
                    if near_const & verboseTog
                        disp(['            - ' rdPars.varnames{varIx} ' varies by <= ' num2str(rdPars.varChangeThresh*100,'%.2f') '% over regime so far -> `near-constant`'])
                        disp( '              but this is the primary focused variable, so making no formal change to variable`s classification')
                    end
                else
                    if near_const
                        if verboseTog
                            disp(['            - ' rdPars.varnames{varIx} ' varies by <= ' num2str(rdPars.varChangeThresh*100,'%.2f') '% over regime so far -> `near-constant`'])
                        end
                        % DEVELOPMENT OPTION: Don't let speedFussy option affect near-constant variables if fussy level = 1
                        if algOptions.speedFussy && (algOptions.speedFussyLevel == 2 & ~const.prev(varIx))
                            newRegime_temp = true;
                            if verboseTog
                                disp(['            * new regime because ' rdPars.varnames{varIx} ' -> near-constant status, and wasn`t near-constant last epoch & speedFussy set at level 2'])
                            end
                        else
                            newRegime_temp = false;
                        end
                        if ( newRegime | newRegime_temp | (~newRegime & const.prev(varIx)) ) % && regimeNum > 0 ??
                            constVars = [constVars, varIx];
                            const.curr(varIx) = 1;
                        else
                            if ~(newRegime | newRegime_temp) && regimeNum > 0
                                regimes(regimeNum).constVars = setdiff(regimes(regimeNum).constVars,varIx);
                                if verboseTog
                                    disp(['            - ... but not new regime, so removing ' rdPars.varnames{varIx} ' from list'])
                                end
                            elseif verboseTog
                                disp(['            - ... but ' rdPars.varnames{varIx} ' wasn`t near-constant before, so not adding to list'])
                            end
                        end
                    else
                        if ~newRegime & const.prev(varIx) & regimeNum > 0 % if inside regime but for previous epochs it was near_const then cancel status
                            % DEVELOPMENT OPTION: see above ...
                            if algOptions.speedFussy && algOptions.speedFussyLevel == 2
                                newRegime_temp = true;
                                if verboseTog
                                    disp(['            - ' rdPars.varnames{varIx} ' no longer near-constant in this regime'])
                                    disp( '            * new regime because speedFussy option set to level 2')
                                end
                            else % keeping this regime & must remove from list
                                regimes(regimeNum).constVars = setdiff(regimes(regimeNum).constVars,varIx);
                                if verboseTog
                                    disp(['            - ' rdPars.varnames{varIx} ' no longer `near-constant`, removing from list'])
                                end
                            end
                        end
                    end
                end % if varIx == primaryFocVarIx
            end
        end
        newRegime = newRegime | newRegime_temp; % update this after possible change in speedFussy option's code in above test
    else
        newRegime = true;
        if verboseTog
            disp('            * new regime because it`s the first epoch of a non-cycle (or only one epoch exists)')
        end
        % test for high/low time derivatives
        if algOptions.speedFussy
            [isPres presIx] = ismember(epochPos,newTSlist(:,1));
            if isPres && newTSlist(presIx,3) ~= 0 % then small derivative test made a new epoch, so note speed change (already new regime!)
                dProceed = false;
                % if primary foc var is fast + not near-constant AND target var is slowing down, then change of speed status
                dVarRel = abs((vLine.curr(primaryFocVarIx) - vLine.next(primaryFocVarIx))/varMaxIntervals(varIx));
                if ismember(primaryFocVarIx,fastVars) & ~ismember(primaryFocVarIx,constVars)
                    if dVarRel >= dt*rdPars.derivThresh & newTSlist(presIx,3) == -1 % target var is slowing down
                        dProceed = true;
                        speedChange.source(speedChange.num)   = 1; % fast
                        speedChange.dest(speedChange.num)   =   0; % normal (or whatever)
                        if verboseTog
                            disp(['            - ' rdPars.varnames{speedChange.var(speedChange.num)} ' speed status no longer fast (low time derivative)'])
                        end
                    end
                else % primary foc var is not fast or is fast but is near-constant, AND target var is speeding up -> change of speed status
                    if dVarRel < dt*rdPars.derivThresh & newTSlist(presIx,3) == 1
                        dProceed = true;
                        speedChange.source(speedChange.num) = 0; % normal (or whatever)
                        speedChange.dest(speedChange.num)   = 1; % fast
                        if verboseTog
                            disp(['            - ' rdPars.varnames{speedChange.var(speedChange.num)} ' speed status -> fast (high time derivative)'])
                        end
                    end
                end
                
                if dProceed
                    speedChange.flag = true;
                    speedChange.num  = speedChange.num + 1;
                    speedChange.var(speedChange.num)  = newTSlist(presIx,2);
                    speedChange.newreg = false; % don't set this to true here otherwise the whole reason for the new regime will not be recorded in the speed change
                end
            end
        end
    end
    
    % After tests, update speeds.curr, slowVars and fastVars now that possible speed changes have occurred or been forced
    if speedChange.flag
        for j=1:speedChange.num
            speeds.curr(speedChange.var(j))=speedChange.dest(j); % update
        end
        [p i]     = find(speeds.curr==-1); % p unused
        slowVars  = intersect( union(dynVars.curr,dynCross.curr), i);
        [p i]     = find(speeds.curr== 1); % p unused
        fastVars  = intersect( union(dynVars.curr,dynCross.curr), i);
    end
    
    
    % Update regimes
    if newRegime
        % update `current` (soon to be `old`) regime's dynamic dimension,
        %     if the regime exists already (e.g. for pos == 1)
        % version "A" is the simple sum of the dynamic vars (inc.
        %     cross-multiplying vars)
        % version "B" takes into account fast vars that are adiabatically
        %     eliminated, and the slow vars that can be approximated with constants
        % version "C" takes into account fast vars that are adiabatically eliminated,
        %     the slow vars that can be approximated with constants,
        %     and the near-constant vars that do not change significantly over the epoch
        %     (except, for the latter two, the `initial` value needs to be known from the previous regime)
        % The `1 +` takes into account the equation for the primary focused variable
        if regimeNum > 0
            regimes(regimeNum).dimensionA = 1 + length(regimes(regimeNum).dynVars) + length(regimes(regimeNum).dynCross);
            regimes(regimeNum).dimensionB = 1 + length(setdiff(regimes(regimeNum).dynVars,union(regimes(regimeNum).fastVars,regimes(regimeNum).slowVars))) ...
                                              + length(setdiff(regimes(regimeNum).dynCross,union(regimes(regimeNum).fastVars,regimes(regimeNum).slowVars)));
            regimes(regimeNum).dimensionC = 1 + length(setdiff(regimes(regimeNum).dynVars,union(regimes(regimeNum).fastVars,union(regimes(regimeNum).slowVars,regimes(regimeNum).constVars)))) ...
                                              + length(setdiff(regimes(regimeNum).dynCross,union(regimes(regimeNum).fastVars,union(regimes(regimeNum).slowVars,regimes(regimeNum).constVars))));
            if isempty(regimes(regimeNum).qsPars)
                regimes(regimeNum).qsPars = setdiff(regimes(regimeNum).slowVars,regimes(regimeNum).constVars); % this may be empty, but qsPars was already empty so doesn't matter
            end
        end

        regimeNum = regimeNum + 1;
        doneEpochs = doneEpochs + 1; % this epoch has been incorporated into a regime
        regimes(regimeNum).epochs     = epochPos;
        regimes(regimeNum).dynVars    = dynVars.curr;
        regimes(regimeNum).dynCross   = dynCross.curr;
        regimes(regimeNum).nonDynVars = nonDynVars.curr;
        regimes(regimeNum).fastVars   = fastVars;
        regimes(regimeNum).slowVars   = slowVars;
        regimes(regimeNum).constVars  = constVars;
        regimes(regimeNum).vLines     = {vLine.curr};
        regimes(regimeNum).qsPars     = qsPars;
        if regimeNum > 1
            if regimes(regimeNum-1).timeInt(2) > tlo.curr
                if global_period == 0
                    fprintf('RegimeDat:  Error: global_period is zero when needed (new regime %i)\n',regimeNum)
                end
                if ~isCycle
                    disp('RegimeDat:  Error: Not in a cycle, but epoch times out of order')
                    return
                end
                tlo.curr = tlo.curr + global_period;
                thi.curr = thi.curr + global_period;
            end
        end
        regimes(regimeNum).timeInt = [tlo.curr, thi.curr];
    elseif regimeNum > 0
        doneEpochs = doneEpochs + 1; % this epoch has been incorporated into a regime
        if verboseTog & numEps > 1
            disp('            - still in same regime')
        end
        regimes(regimeNum).epochs       = [ regimes(regimeNum).epochs, epochPos ];
        regimes(regimeNum).dynVars      = union(regimes(regimeNum).dynVars, dynVars.curr); % extend to currently potentiated, previously active vars
        regimes(regimeNum).dynCross     = union(regimes(regimeNum).dynCross, dynCross.curr); % extend
        regimes(regimeNum).nonDynVars   = union(regimes(regimeNum).nonDynVars, nonDynVars.curr); % extend to previously active non-dyn vars
        regimes(regimeNum).vLines       = { regimes(regimeNum).vLines{:}, vLine.curr };
        regimes(regimeNum).constVars    = union(regimes(regimeNum).constVars, constVars);
        if speedChange.flag % currently only get here for a non-new-regime fast -> slow change, so update slow & fast vars
            regimes(regimeNum).slowVars = slowVars; % currently overwrite them all, but could use other speedChange fields to be specific
            regimes(regimeNum).fastVars = fastVars;
        end
        if regimes(regimeNum).timeInt(1) >= thi.curr
            if global_period == 0
                fprintf('RegimeDat:  Error: global_period is zero when needed (old regime %i)\n',regimeNum)
            end
            if ~isCycle
                disp('RegimeDat:  Error: Not in a cycle, but epoch times out of order')
                return
            end
            thi.curr = thi.curr + global_period;
        end
        regimes(regimeNum).timeInt = [regimes(regimeNum).timeInt(1), thi.curr]; % extend
        regimes(regimeNum).qsPars  = setdiff(union( regimes(regimeNum).qsPars, qsPars ), regimes(regimeNum).dynVars); % extend, but take out any that are in dynVars if that was extended
        % the following will keep being rechecked as the time interval grows to its full size, which will not create inconsistencies
        if ~isempty(regimes(regimeNum).slowVars) % check relative length of regime vs. slow time scale
            regDuration = regimes(regimeNum).timeInt(2) - regimes(regimeNum).timeInt(1);
            for var=regimes(regimeNum).slowVars
                sumG1_av = 0;
                for vL=1:length(regimes(regimeNum).vLines)
                    sumG1 = SumGamma1(DEqns_comp{rdPars.DEixMap(var)}{DE_GAMMA1i},regimes(regimeNum).vLines{vL}); % speed, so 1/(C*sumG1(1)) is a time scale
                    sumG1_av = sumG1_av + sumG1(1);
                end
                speed_av = (sumG1_av / DEqns_comp{rdPars.DEixMap(var)}{DE_CFACi}) / length(regimes(regimeNum).vLines);
                if regDuration*speed_av < rdPars.tScaleThresh | regDuration*speed_av > 1/rdPars.tScaleThresh % multiply sumG1(1) because it's a speed not a time scale
                    % then this `slow` variable will actually change O(1) over this `long` epoch, so better upgrade the variable to normal if algOptions.longFussy set
                    if algOptions.longFussy
                        regimes(regimeNum).slowVars = setdiff( regimes(regimeNum).slowVars, var ); % remove it from slow list
                        speedChange.flag   = true; % update this for next epoch around loop
                        speedChange.var    = var;
                        speedChange.newreg = false;
                        speedChange.source = -1; % was slow
                        speedChange.dest   = 0; % force to normal speed
                        if verboseTog
                            fprintf('            - %s was forced to change from slow -> normal speed because of long duration of regime\n',rdPars.varnames{var})
                        end
                    else
                        if verboseTog
                            fprintf('            - %s will change O(1) over long epoch, but not changing status because longFussy option not set\n',rdPars.varnames{var})
                        end
                    end
                    break % don't bother checking any more vLines
                end
            end
        end

    else % else this epoch belongs to an 'earlier' regime that will start at the end of the loop, so skip
        if verboseTog
            disp('            - skipping epoch until new regime found')
        end
    end
    
    if doneEpochs == numEps
        break % for loop -- finish search
    end
    
    const.prev = const.curr; % transfer this to the next epoch
end

% Finish updating regime information for final regime (update current regime's dynamic dimension, if regime exists)
if regimeNum > 0
    regimes(regimeNum).dimensionA = 1 + length(regimes(regimeNum).dynVars) + length(regimes(regimeNum).dynCross);
    regimes(regimeNum).dimensionB = 1 + length(setdiff(regimes(regimeNum).dynVars,union(regimes(regimeNum).fastVars,regimes(regimeNum).slowVars))) ...
                                      + length(setdiff(regimes(regimeNum).dynCross,union(regimes(regimeNum).fastVars,regimes(regimeNum).slowVars)));
    regimes(regimeNum).dimensionC = 1 + length(setdiff(regimes(regimeNum).dynVars,union(regimes(regimeNum).fastVars,union(regimes(regimeNum).slowVars,regimes(regimeNum).constVars)))) ...
                                      + length(setdiff(regimes(regimeNum).dynCross,union(regimes(regimeNum).fastVars,union(regimes(regimeNum).slowVars,regimes(regimeNum).constVars))));
    if isempty(regimes(regimeNum).qsPars)
        regimes(regimeNum).qsPars = setdiff(regimes(regimeNum).slowVars,regimes(regimeNum).constVars); % this may be empty, but qsPars was already empty so doesn't matter
    end
end

if numEps == 1
    if algOptions.speedFussy
        disp('RegimeDet:  Only one epoch passed, and speedFussy algorithm did not identify')
        disp('            further partitions of the orbit. Therefore, the epochs became a regime.')
    else
        disp('RegimeDet:  Only one epoch passed, and speedFussy option not set.')
        disp('            Therefore, the epoch became a regime.')
    end
end

if exist( 'regimes' ) ~= 1
    disp('RegimeDet:  No regimes created. Returning to StepNet ...')
    return
end

if isCycle
    startTime = regimes(1).timeInt(1); % updated from passed parameter since isCycle can cause actual first regime to begin later
else
    startTime = rdPars.startTime;
end


%%%%% Build return structure
regimeStruct.regimes    = regimes;
regimeStruct.varnames   = rdPars.varnames;
regimeStruct.numEps     = numEps;
regimeStruct.numRegs    = regimeNum;
regimeStruct.transSeq   = loopTS;
regimeStruct.isCycle    = isCycle;
regimeStruct.startTime  = startTime; % possibly updated because of isCycle
regimeStruct.period     = global_period; % zero if ~isCycle
regimeStruct.algOptions = algOptions;
regimeStruct.tScaleThresh = rdPars.tScaleThresh;
regimeStruct.focVarIx     = primaryFocVarIx;

result = {1, regimeStruct}; % 1 = success
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                            AUXILIARY FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function result = isAbnormalSpeed(inputDE,focDE,varDataLine,tScaleThresh)
% only accepts "compiled" DEs
global DE_GAMMA1i DE_CFACi

inputSg1 = SumGamma1(inputDE{DE_GAMMA1i},varDataLine); % input to 'foc' var's DE
focSg1 = SumGamma1(focDE{DE_GAMMA1i},varDataLine);
input_tscale = inputDE{DE_CFACi}/inputSg1(1);
foc_tscale = focDE{DE_CFACi}/focSg1(1);
tRat = abs(input_tscale / foc_tscale);
if tRat > tScaleThresh
    result = [-1, tRat, 0]; % 'input is slow'
elseif tRat < 1/tScaleThresh
    result = [1, 1/tRat, 1]; % 'input is fast'
else
    if tRat < 1
        result = [0, 1/tRat, 1]; % 'normal' (third param here indicates whether regular or reciprocal of ratio returned)
    else
        result = [0, tRat, 0]; % 'normal'
    end
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function qsfp = GetQsfpVal(thisDE,varDataLine,onlyActives)
% only accepts "compiled" thisDE
global DE_NAMEi DE_GAMMA1i DE_GAMMA2i DE_ACTSWi DE_ISTAUFFILEi DE_TAURECIPi DE_GAMVARNAMEi ...
    DE_GAMVARPOWi DE_ISTGTFFILEi DE_TARGETi DE_INTVARNAMEi DE_INTVARPOWi LARGEBOUND
% V0eqn in terms of conductances and currents (for example of voltage equations):
%  curly braces imply summation
%  define {all currents} = {conductances*revpots} + {direct currents}
%  so V0 = sum {all currents} / sum {conductances}
if nargin == 2
    onlyActives = false;
end

sumG1 = SumGamma1(thisDE{DE_GAMMA1i},varDataLine,onlyActives);
if sumG1(1) ~= 0
    qsfp = ( sumG1(2) + SumGamma2(thisDE{DE_GAMMA2i},varDataLine,onlyActives) ) / sumG1(1);
else
    disp('RegimeDet:  Warning. Division by zero in GetQsfpVal()')
    qsfp = LARGEBOUND;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sum = SumGamma1(Gamma1,varDataLine,onlyActives)
% only accepts "compiled" Gamma1
% gam1term order: actSw, filefuncflag, taurecip, var, power, filefuncflag, target [, intvar, power]
global DE_NAMEi DE_GAMMA1i DE_GAMMA2i DE_ACTSWi DE_ISTAUFFILEi DE_TAURECIPi DE_GAMVARNAMEi ...
    DE_GAMVARPOWi DE_ISTGTFFILEi DE_TARGETi DE_INTVARNAMEi DE_INTVARPOWi
if nargin == 2
    onlyActives = false;
end
sum1 = 0;
sum2 = 0;
for g1t=1:length(Gamma1)
    g1term = Gamma1{g1t};
    if ~onlyActives || ( onlyActives & g1term{DE_ACTSWi} )
        varix = g1term{DE_GAMVARNAMEi};
        if varix == 0
            varVal = 1;
        else
            varVal = varDataLine(varix);
        end
        tau_recipVal = GetTauRVal(g1term{DE_ISTAUFFILEi},g1term{DE_TAURECIPi},...
            varVal,g1term{DE_GAMVARPOWi});
        targetVal = GetTargVal(g1term{DE_ISTGTFFILEi},g1term{DE_TARGETi},varVal,varDataLine);
        if g1term{DE_INTVARNAMEi} > 0
            sum1 = sum1 + tau_recipVal * varDataLine(g1term{DE_INTVARNAMEi})^g1term{DE_INTVARPOWi};
            sum2 = sum2 + tau_recipVal * varDataLine(g1term{DE_INTVARNAMEi})^g1term{DE_INTVARPOWi} * targetVal;
        else
            sum1 = sum1 + tau_recipVal;
            sum2 = sum2 + tau_recipVal * targetVal;
        end
    end
end
sum = [sum1, sum2];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sum = SumGamma2(Gamma2,varDataLine,onlyActives)
% only accepts "compiled" Gamma2
% gam2term order: actSw, filefuncflag, taurecip, var, power
global DE_NAMEi DE_GAMMA1i DE_GAMMA2i DE_ACTSWi DE_ISTAUFFILEi DE_TAURECIPi DE_GAMVARNAMEi ...
    DE_GAMVARPOWi DE_ISTGTFFILEi DE_TARGETi DE_INTVARNAMEi DE_INTVARPOWi
if nargin == 2
    onlyActives = false;
end
sum = 0;
for g2t=1:length(Gamma2)
    g2term = Gamma2{g2t};
    if ~onlyActives || ( onlyActives & g2term{DE_ACTSWi} )
        tau_recipVal = GetTauRVal(g2term{DE_ISTAUFFILEi},g2term{DE_TAURECIPi},...
            varDataLine(g2term{DE_GAMVARNAMEi}),g2term{DE_GAMVARPOWi});
        sum = sum + tau_recipVal;
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = GetTargVal(sourceType,source,varVal,vLine)
switch sourceType
    case 0
		result = source;
    case 1
        result = eval( [ source '(' num2str(varVal) ')' ]);
    case 2
        result = vLine(source);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = GetTauRVal(sourceType,parOrFunc,varVal,power)
if sourceType == 1
    result = eval( [ parOrFunc '(' num2str(varVal) ')' ]);
else
    result = parOrFunc * varVal^power;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% useful as a dummy argument for KeyPressFcn callback, to prevent
% echo of key presses to command window
function DoNothing()
return
