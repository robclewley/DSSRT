function attRes = AttEst_regimes(regimeStruct, varData, seqTimes, DEqns, DEpars, t0_TS, Tperiod, attEstParams, dScaleThresh, ...
                         tScaleThresh, numExt, numInt, varnames, actsIxMap, inputsIx, varBounds, figHandle, verboseTog)
% Attractor Estimation procedure for a single local variable using reduced regimes
% Version 1.1, November 2005
% (c) Robert Clewley, Cornell University
%
% Notes about variables/parameters (including those for main auxiliary functions,
%   CompressEst(), NeighbHdEst() ):
% focusVarIx   defines a primary variable of focus (absolute index in set 'varnames')
% TSfocusSet   is the set of all focusVarIx-dependent variables that are in focus in the Transition
%               Sequence if this option is used (includes focusVarIx as first element), otherwise empty
% dVmax        is the largest acceptable range for focusVarIx's variable when searching
% Vinres       is the variable's inner-layer resolution when searching for the inverse Phi map domain
%               (when the use of the infinity values for V-dependent variables is switched off)
% dScaleThresh is the dominance scale threshold (1/epsilon) used to generate the functional network
% Vpert        defines the +/- distance from V(t0) used as an initial target domain
% tScaleThresh defines the scale threshold (1/gamma) for the fast time scale of fast V-dependent variables
% delta        is the time tolerance at chosen critical epoch/regime changeover,
%               but is currently unused (will replace Vpert)
% Vknown_t1    original V(t1) from known orbit used by neighbourhood estimator
% Vend         V(t2) targets used by neighbourhood estimator

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEVELOPMENT NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 1.1, November 2005
%   Fixed GetPsiVals function to support gap-junction coupling.
%
% Version 1.03, June 2005
%   Slight alteration to way in with PrepareFocVars avoids focused variable
%   in loop.
%
% Version 1.02, August 2004
%   Added noShooting and compareActs threshold parameters as user-passed parameter in attEstParams
%   Fixed apparent error in array reference code for arguments to CompareActs calls (see comments in code)
%   Fixed infSwitchOn not being reset to true for VendIx = 2 in first interval check of NeighbHdEst
%   Fixed compareActs to work for domain of validity calculation
%   Added spaces between every 10 dots printed for time intervals checked
%
% Version 1.00 (based on v.1.41 of AttEst_epochs.m)
%   Option to change output verbosity on the `escape` key interrupt now added
%   Will remove numExt, numInt, t/dScaleThresh, varnames, actsIxMap, inputsIx, varBounds,
%     and pass rdPars instead
%   t0_TS is now used to get correct offset against times associated with passed TS,
%     because regimes may not start at the first epoch
%   CompareActs() is completely different
%
% Pre-v.1.00
%  See notes for AttEst_epochs.m version 1.41 and before
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            CONSTANTS and INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global figHand maxVresFrac diagnosticOn diagnosticOn_def
GlobalConstants

%%%%% Internal parameters %%%%%
%% for initial uniformly-spaced tests in V interval at t1
maxVresFrac = 25;

%% for ignoring some changes in Actives set during perturbations
ignoreActs_flag = false; % THIS OPTION IS IN DEVELOPMENT (NOT CURRENTLY USED)

%% diagnostic information switch (for development purposes), used e.g. for tuning maxVresFrac
diagnosticOn_def = true; % but will only activate if verboseTog = true
screenWrap = 70; % in characters (for displaying dots per interval done with diagnosticOn_def but ~verboseTog)

success = false; % success of this AttEst function call
regimes = {}; % focused on the variable chosen
attRes = {success, regimes, 0, 0, 0}; % default result (error status)

caOpts.pt1 = attEstParams.posThresh1;
caOpts.pt2 = attEstParams.posThresh2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             PASSED-PARAMETER CHECKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('AttEst:  This is the `regimes` version of AttEst')

if nargin ~= 18
    beep
    disp('AttEst:  Internal error! Wrong number of arguments passed. Expected 18.')
    return
end
figHand = figHandle; % for use globally by subroutines
tseq = regimeStruct.transSeq;
TSname = tseq{TS_NAME};
TSlength = length(tseq{TS_TSEQ});

if Tperiod <= 0
    beep
    disp('AttEst:  Internal error - Tperiod must be greater than zero.')
    return
end

TSfocusSet = tseq{TS_FOCUS};
if ~isempty(TSfocusSet)
    focusVarIx = TSfocusSet(1);
    if length(TSfocusSet) > 1
        for tfix=2:length(TSfocusSet)
            if TSfocusSet(tfix) < 0 | TSfocusSet(tfix) > length(varnames)
                disp('AttEst:  Internal error! TSfocusSet entry out of range...')
                return
            end
        end
    end
else
    disp('AttEst:  Warning - focus set is empty. Attractor Estimation is')
    disp('          not possible for unfocused transition sequences')
    return
end
if focusVarIx < 0 | focusVarIx > length(varnames)
    beep
    disp('AttEst:  Internal error! Index for variable of focus is out of range.')
    return
end
if numExt + numInt ~= length(varnames)
    beep
    disp('AttEst:  Discrepancy between length of `varnames` paramter and')
    disp('           expected number of internal and external variables.')
    return
end
numTot = numExt + numInt;

%% delta is currently unused
% if delta <= 0 | delta >= Tperiod
%     beep
%     disp('AttEst:  delta must be greater than zero and less than Tperiod.')
%     return
% end

diagnosticOn = diagnosticOn_def & verboseTog;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   REGIME INITIALIZATION (including BOLSTERING)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('AttEst:  Initializing regime data ... ')

% Time step dt is according to assumed FIXED step size from integrated data
dt = (round(seqTimes(2)/TIMERES)-round(seqTimes(1)/TIMERES))*TIMERES; % accurate to TIMERES resolution

allActNames = {}; % names of all active vars appearing in focusVarIx's Acts set
if verboseTog
    fprintf('AttEst:  All Actives names for focused variable:  ')
end
for entry = 1:length( actsIxMap{focusVarIx} )
    newName = varnames{actsIxMap{focusVarIx}(entry)};
    if verboseTog
        fprintf('%s ',newName)
    end
    allActNames = {allActNames{:} newName};
end
if verboseTog
    disp(' ')
end

numRegimes = regimeStruct.numRegs; % does NOT include a repeated final entry in the case of cycles (cf. for transition sequences)
if numRegimes == 0
    disp(' ')
    disp('AttEst:  No regimes present. Cannot proceed.')
    return
end

% Build dummy final regime from first if cycle, otherwise from TS
if regimeStruct.isCycle
    regimeStruct.regimes(numRegimes+1) = regimeStruct.regimes(1);
    regimeStruct.regimes(numRegimes+1).timeInt = regimeStruct.regimes(numRegimes+1).timeInt + Tperiod;
else
    regimeStruct.regimes(numRegimes+1) = regimeStruct.regimes(1); % most of the entries are irrelevant for this dummy entry
    TS_dummy = regimeStruct.transSeq{TS_TSEQ}{TSlength}; % this is the dummy TS entry
    regimeStruct.regimes(numRegimes+1).timeInt = [TS_dummy{TSEQ_TIME}; TS_dummy{TSEQ_TIME}+1];
    regimeStruct.regimes(numRegimes+1).dynVars = union(TS_dummy{TSEQ_ACTS}{focusVarIx},...
        regimeStruct.regimes(numRegimes+1).qsPars); % doesn't really matter that they may not all be dynamic here
    regimeStruct.regimes(numRegimes+1).epochs  = TSlength;
end
numRegimes = numRegimes+1;

numVarsFocused = length(TSfocusSet); % not to be confused with numActsFocused, concerning actives
% Prepare derivative thresholds for deciding resolution of bolster regimes
varMaxIntervals = zeros(1,numVarsFocused);
for varIx = 1:numVarsFocused
    absVarIx = TSfocusSet(varIx);
    varMaxIntervals(varIx) = varBounds(absVarIx,2) - varBounds(absVarIx,1);
end

tEpBeginMap = zeros(1,TSlength);
for epIx=1:TSlength
    tEpBeginMap(epIx) = tseq{TS_TSEQ}{epIx}{TSEQ_TIME}; % relative times of the beginnings of each epoch (rel. to t0_TS)
end

% Diagnositc counters for bolstered regimes
numBlong = 0;
numBshort = 0;
numBmedium = 0;

numActs = length(tseq{TS_TSEQ}{1}{TSEQ_ACTS}); % arbitrarily look at first TS to see how many entries
numRegsTot = 0;
t0 = regimeStruct.regimes(1).timeInt(1);
for pos = 1:numRegimes % numRegimes is not total # of time intervals used in analysis, because we might bolster the set...
    if verboseTog
        fprintf('          Regime %2i',pos)
    end
    reg_temp = regimeStruct.regimes(pos);
    % Note: [tlo, thi] must be time intervals *relative* to t0. For absolute times, add t0
    tlo = reg_temp.timeInt(1) - t0;
    if pos < numRegimes
        thi = reg_temp.timeInt(2) - t0;
    else
        thi = regimeStruct.regimes(1).timeInt(2)+Tperiod;
        if thi < tlo
            beep
            disp(' ')
            disp(['           Specified period smaller than regime ' num2str(pos) '`s interval'])
            disp('            length in transition sequence. Cannot proceed.')
            return
        end
    end

    if verboseTog
        fprintf(' over rel. time [%3.2f, %3.2f], [%3.2f, %3.2f] absolute\n',tlo,thi,t0+tlo,t0+thi)
    end

    %%%%%% Find which actives are present for the whole regime (including bolster intervals)
    acts = union(reg_temp.dynVars,reg_temp.nonDynVars);
    actNames = {}; % initialize
    allActives = cell(1,numActs);
    for epIx = reg_temp.epochs
        for varIx=1:numActs
            allActives{varIx} = union(allActives{varIx}, tseq{TS_TSEQ}{epIx}{TSEQ_ACTS}{varIx});
        end
    end
    if verboseTog
        fprintf('            Actives for primary focused variable:  ')
    end
    for entry = 1:length( acts )
        newName = varnames{acts(entry)};
        if verboseTog
            fprintf('%s ',newName);
        end
        actNames = {actNames{:} newName};
    end
    if verboseTog
        fprintf('\n')
    end
    
    %%%%%% Add regular regime to regime list, but ...
    %   Flag the need for additional bolster regimes in this [tlo,thi],
    %       depending on maximum of derivatives of all focused variables!
    regTlo = tlo; % current lower regime time
    [isp regTloPos] = ismember(1,seqTimes <= regTlo+t0); % was FindPosInT(regTlo,t0,seqTimes);
    tloPos = regTloPos; % true at this initial point
    [isp thiPos] = ismember(1,seqTimes <= thi+t0); % was FindPosInT(thi,t0,seqTimes);
    origFlag = false;
    dummyVInterval = [0,0]; % this is a placeholder for the V interval of the attractor cross-section
    if pos == numRegimes % final regime is special: see above note
        numRegsTot = numRegsTot + 1;
        regimes = {regimes{:} {numRegsTot; [tlo, thi]; [tloPos, thiPos]; allActives; actNames; dummyVInterval; true; reg_temp.qsPars} };
        % force origFlag = true so that red vertical line appears at end of attractor estimate diagram.
        continue
    end
    addRegimes = true; % initial value (adding bolster intervals to the original regime interval)
    bolstered = false; % flag to indicate whether bolster intervals were added
    searchMarkPos = -1; % dummy value indicating "currently not in use"
    tsearch_old = regTlo; % initial value (start at tlo, and then whatever regTlo is set subsequently)
    tsearchPos_old = regTloPos; % initial value
    origIntPos = [tloPos, thiPos]; % original time interval for regime, in terms of positions
    varData_old = varData(tsearchPos_old,:);
    dVarRel = zeros(1,numVarsFocused);
    while addRegimes
        origFlag = false;
        tsearch_new    = tsearch_old + attEstParams.searchRes*dt;
        tsearchPos_new = tsearchPos_old + attEstParams.searchRes;
        if tsearchPos_new >= thiPos % then end of an original regime (not considered to be adding a bolster interval here)
            if ismember( regTloPos, origIntPos )
                origFlag = true; % part of an original (non-bolster) regime
            end
            if regTloPos < thiPos
                numRegsTot = numRegsTot + 1;
                regimes = {regimes{:} {numRegsTot; [regTlo, thi]; [regTloPos, thiPos]; allActives; actNames; ...
                            dummyVInterval; origFlag; reg_temp.qsPars} };
            else
                if verboseTog | diagnosticOn
                    disp('            Warning while bolstering: `consecutive` time positions were in fact not!')
                    disp('              (skipping)')
                end
            end
            addRegimes = false; % no more time left in the original regime interval [tlo, thi]
            continue
        end
        if tsearchPos_new == searchMarkPos % reached a lo-res marker set earlier
            if ismember( regTloPos, origIntPos )
                origFlag = true; % part of an original (non-bolster) regime
            end
            if regTloPos < tsearchPos_new
                numRegsTot = numRegsTot + 1;
                numBlong = numBlong + 1;
                regimes = {regimes{:} {numRegsTot; [regTlo, tsearch_new]; [regTloPos, tsearchPos_new]; allActives; ...
                            actNames; dummyVInterval; origFlag; reg_temp.qsPars}};
            else
                if verboseTog | diagnosticOn
                    disp('            Warning while bolstering: `consecutive` time positions were in fact not!')
                    disp('              (skipping)')
                end
            end
            bolstered = true;
            regTlo = tsearch_new;
            regTloPos = tsearchPos_new;
            searchMarkPos = -1;
        end
        varData_new = varData(tsearchPos_new,:);
        % approximate time derivatives
        for varIx = 1:numVarsFocused
            thisVar = TSfocusSet(varIx);
            dVarRel(varIx) = abs((varData_new(thisVar) - varData_old(thisVar))/varMaxIntervals(varIx));
        end
        % test size of derivatives
        if max(dVarRel) >=  attEstParams.searchRes*dt*attEstParams.derivThresh
            % allocate high resolution --> regime boundary here at tsearch position
            % close up the present regime if it wasn't a high res bolster regime
            nextRegTlo    = tsearch_old; % default value
            nextRegTloPos = tsearchPos_old; % default value
            if regTloPos ~= tsearchPos_old
                if ismember( regTloPos, origIntPos )
                    origFlag = true; % part of an original (non-bolster) regime
                end
                if regTloPos < tsearchPos_old
                    numRegsTot = numRegsTot + 1;
                    numBmedium = numBmedium + 1;
                    regimes = {regimes{:} {numRegsTot; [regTlo, tsearch_old]; [regTloPos, tsearchPos_old]; ...
                                allActives; actNames; dummyVInterval; origFlag; reg_temp.qsPars}};
                else
                    if verboseTog | diagnosticOn
                        disp('            Warning while bolstering: `consecutive` time positions were in fact not!')
                        disp('              (skipping)')
                    end
                    nextRegTlo    = regTlo;
                    nextRegTloPos = regTloPos;
                end
            end
            origFlag = false; % default value here (may be changed below)
            % add high res regime
            if ismember( nextRegTloPos, origIntPos )
                origFlag = true; % part of an original (non-bolster) regime
            end
            if nextRegTloPos < tsearchPos_new
                numRegsTot = numRegsTot + 1;
                numBshort = numBshort + 1;
                regimes = {regimes{:} {numRegsTot; [nextRegTlo, tsearch_new]; [nextRegTloPos, tsearchPos_new]; ...
                            allActives; actNames; dummyVInterval; origFlag; reg_temp.qsPars}};
            else
                if verboseTog | diagnosticOn
                    disp('            Warning while bolstering: `consecutive` time positions were in fact not!')
                    disp('              (skipping)')
                end
            end
            regTlo = tsearch_new; % for next iteration
            regTloPos = tsearchPos_new;
            searchMarkPos = -1; % reset this, since we had to add a high-res regime in the meantime
            bolstered = true;
        else
            % if it hasn't already been set, then set searchMarkPos ...
            % allocating low resolution, i.e. bolster regime end at tsearchPos_old + lowResMultiple * searchRes,
            % pending search through remainder of the [tlo,thi] interval
            if searchMarkPos == -1
                searchMarkPos = tsearchPos_old + attEstParams.lowResMultiple * attEstParams.searchRes;
            end
            % if search stepping overruns original thi boundary before this is reached,
            % then this will be ignored
        end
        varData_old    = varData_new;
        tsearch_old    = tsearch_new;
        tsearchPos_old = tsearchPos_new;
    end % while addRegimes
    if verboseTog & bolstered
            fprintf('            Bolster intervals were added (%i short, %i medium, %i long)\n', numBshort, numBmedium, numBlong)
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               REGIME ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numRegs = numRegsTot - 1;

if attEstParams.noShooting
    disp('AttEst:  Calculating domain of validity for regimes')
    disp('          (Not using pseudo-shooting estimation)')
else
    disp('AttEst:  Calculating attractor estimate for regimes')
    disp('          (Using pseudo-shooting estimation)')
end

fprintf('AttEst:  %i sample points will be analyzed.\n',numRegs)
disp('         This may take several minutes (press <escape> to cancel or switch output verbosity) ...')
if verboseTog
    fprintf(['          Search step size in ' varnames{focusVarIx} ' is %.4f\n'],attEstParams.Vinres)
end

actsFocused = cell(1,numVarsFocused); % e.g. for `E_demo`, actsFocused is the absolute indices corresponding to {'mu1', 'nu1', 'si', 'iu1', 'lu1'}
numActsFocused = zeros(1,numVarsFocused);
trivialCandActs = zeros(1,numVarsFocused); % this is also used in neighbourhood estimate
for focusSetIx=1:numVarsFocused
	absIx = TSfocusSet(focusSetIx);
    if absIx > numExt
        continue % candidate actives only exist for externals
    end
    varFocIps     = inputsIx{absIx};
    thisActsIxMap = actsIxMap{absIx};
    lenAF = 0; % length of actsFocused counter
	for ipIx = varFocIps
        if ismember(ipIx,thisActsIxMap) % then its a candidate active variable, so must include
            lenAF = lenAF + 1;
            actsFocused{focusSetIx} = [actsFocused{focusSetIx}, ipIx];
        end
	end
    numActsFocused(focusSetIx) = lenAF;
    trivialCandActs(focusSetIx) = ( lenAF <= 1 );
end


%% Get dependency information for focusVarIx variable 'V' for UNCOMPILED DEqns, and for others in TSfocusSet
% This is used only by NeighbHdEst()
focusVarNames = cell(1,numVarsFocused);
for focusSetIx=1:numVarsFocused
    VdepInfo{focusSetIx} = GetDepInfo(TSfocusSet(focusSetIx), varnames, inputsIx, DEqns);
    focusVarNames{focusSetIx} = varnames{TSfocusSet(focusSetIx)}; % used for compiling correct equations    
end

%% Work out which equation belongs to focusVarIx and others in TSfocusSet, by creating an index map
%% Also, "compile" names -> values where possible for DEqns
deqnIxMap = zeros(1,numTot); % (only for focused variables) - zero entries will denote no corresponding differential equation present
allDEixMap = zeros(1,numTot); % (for all variables)
numEqns = length(DEqns);
DEqns_comp = cell(1,numEqns);
for eqn=1:numEqns
    % create map entry
    [ispresent1 pos1] = ismember(DEqns{eqn}{DE_NAMEi},focusVarNames); % pos1 is relative to focusVarNames
    [ispresent2 pos2] = ismember(DEqns{eqn}{DE_NAMEi},varnames); % pos2 is an absolute index
    if ispresent1
        deqnIxMap(TSfocusSet(pos1)) = eqn;
    end
    if ispresent2
        allDEixMap(pos2) = eqn;
    end
    % compile this eqn
    DEqns_comp{eqn} = CompileGammas(DEqns{eqn}, DEpars, varnames);
end

%% Create map from [Gamma1, Gamma2] order of input variables to the set of absolute indices (the order of varnames)
absGamIxMap = cell(1,numVarsFocused);
for focusSetIx=1:numVarsFocused
    absGamIxMap_temp = [];
    thisDEix = deqnIxMap(TSfocusSet(focusSetIx));
    if thisDEix > 0
        for g1t=1:length(DEqns_comp{thisDEix}{DE_GAMMA1i})
            absGamIxMap_temp = [ absGamIxMap_temp DEqns_comp{thisDEix}{DE_GAMMA1i}{g1t}{DE_GAMVARNAMEi} ];
        end
        for g2t=1:length(DEqns_comp{thisDEix}{DE_GAMMA2i})
            absGamIxMap_temp = [ absGamIxMap_temp DEqns_comp{thisDEix}{DE_GAMMA2i}{g2t}{DE_GAMVARNAMEi} ];
        end
        absGamIxMap{focusSetIx} = absGamIxMap_temp;
    end
end

%% Get primary focused variable's equation
primaryEqnIx = deqnIxMap(focusVarIx);
if primaryEqnIx == 0
    disp('AttEst:  Internal error - primary focused variable does not have an associated differential equation')
    return
end
primaryDE = DEqns{primaryEqnIx}; % UNCOMPILED version
% Get derivative information for UNCOMPILED DEqns relative to primary focused var
derivInfo = GetDerivInfo(focusVarIx, varnames, primaryEqnIx, DEqns, DEpars); % call this function with UNCOMPILED DEqns only

% PhiMaps will hold coefficient and offset for affine linear Phi map of each regime *forward* in time
PhiMaps = zeros(numRegs,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%       MAIN LOOP THROUGH REGIMES       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

horizDotsDone = 0; % only used if diagnosticOn_def and ~verboseTog
if ~verboseTog & diagnosticOn_def
    fprintf('         ')
end

% Step through regimes backwards from the final one
for pos = 1:numRegs % can multiply by an integer to go around a cycle more times, if a cycle is being studied
    if regimeStruct.isCycle
        prevRegPos = mod(-pos-1,numRegs)+1;
        regPos     = mod(prevRegPos,numRegs)+1;
        nextRegPos = mod(regPos,numRegs)+1;
    else
        prevRegPos = numRegs-pos; % will go to zero! so be careful using it
        regPos     = prevRegPos + 1;
        nextRegPos = regPos + 1;
    end
    regime   = regimes{regPos}; % searches backwards!
    regAbsId = regime{1}; % absolute ID number of regime (according to order established above)
    nextRegAllActs = regimes{nextRegPos}{4}; % the one forward in time -- wraps around, assuming periodic orbit!
     %%% THE ABOVE LINE WILL MESS-UP FOR NON-CYCLIC ORBITS W/O ADDITION OF A FINAL DUMMY REGIME
    regAllActs     = regime{4}; % active sets for all variables
    regActNames    = regime{5}; % for screen display use (primary focused variable only)
    qsPars         = regime{8}; % quasi-static bifurcation parameters that are not active in regime
    tlo   = regime{2}(1);
    thi   = regime{2}(2);
    poslo = regime{3}(1);
    poshi = regime{3}(2);
    tData = seqTimes(poslo:poshi);
    vData = varData(poslo:poshi,:);

    % original actives sets for all vars in epochs at t1 and t2
    epLo  = max(find(tEpBeginMap <= tlo+(t0-t0_TS)+SMALLBOUND)); % need smallbound correction for accumulated rounding errors
    epHi  = max(find(tEpBeginMap <=  thi+(t0-t0_TS)+SMALLBOUND));
    if epLo == TSlength
        epLo = max(find(tEpBeginMap<=tlo+(t0-t0_TS)+SMALLBOUND-Tperiod));
    end
    if epHi == TSlength
        epHi = max(find(tEpBeginMap<=thi+(t0-t0_TS)+SMALLBOUND-Tperiod));
    end
    local_acts_t1 = tseq{TS_TSEQ}{ epLo }{TSEQ_ACTS};
    local_acts_t2 = tseq{TS_TSEQ}{ epHi }{TSEQ_ACTS};

    % DIAGNOSTICS OPTION
    % Use this section to add scheduled changes in verbosity, or to place breakpoints
%     if regAbsId < 12
%         verboseTog = true;
%         diagnosticOn = diagnosticOn_def;
%     end

    if ~verboseTog & diagnosticOn_def % then diagnosticOn == false but _def == true
        if horizDotsDone == screenWrap
            horizDotsDone = 0;
            fprintf('\n         ')
        end
        if mod(horizDotsDone,10) == 0 && horizDotsDone > 0
            fprintf(' ')
        end
        fprintf('.')
        horizDotsDone = horizDotsDone + 1;
    end
    % check if user is trying to interrupt AttEst or change output verbosity
    drawnow
    if get(figHandle,'CurrentCharacter') == 27
        disp(' ')
        ButtonName = questdlg('Quit or change verbosity?','AttEst interrupted','Quit','Change','Continue','Continue');
        switch ButtonName
            case 'Quit'
                disp('AttEst:  Attractor estimate operation cancelled by user')
                result{5} = -1; % cancel code
                return
            case 'Change'
                if ~verboseTog & diagnosticOn_def % then diagnosticOn == false but _def == true
                    fprintf('\n')
                end
                verboseTog   = ~verboseTog;
                diagnosticOn = diagnosticOn_def & verboseTog;
                disp('AttEst:  Output verbosity changed. Press any key over main window to continue ...')
            case 'Continue'
                % do nothing
                if ~verboseTog & diagnosticOn_def % then diagnosticOn == false but _def == true
                    fprintf('\n')
                end
                disp('AttEst:  Press any key over main window to continue ...')
        end
        %% Wait for user command over figure
        keyokay = false;
        while ~keyokay
            w = 0;
            while w~=1
                w = waitforbuttonpress;
            end
            keypress = get(figHandle,'CurrentCharacter');
            if ischar(keypress) & ~ismember(keypress, [27, 14, 15]) % ensures ESCAPE no longer current char
                keyokay = true;
            end
        end
        if ~verboseTog
            disp('         Continuing ...')
            if diagnosticOn_def
                fprintf('         ')
            end
        end
    end
    
    %%%%%%%%% INITIAL V TARGET INTERVAL TO SHOOT FOR
    %%% TEMPORARY %%%
    if regPos == numRegs % then this is last regime in sequence = first to analyze, so set the target here
        lastVval = varData(poshi,focusVarIx); % and use this value at end of function for error cf. V_Phi_mapped
        % change initial perturbation to be maximal if calculating DoV (override)
        if attEstParams.noShooting
            initVinterval = [varBounds(focusVarIx,1), varBounds(focusVarIx,2)];
        else
            initVinterval = [lastVval-initVpert,lastVval+initVpert]; % initial value only 
        end
        Vend = initVinterval;
        firstInterval_flag = true;
    else
        firstInterval_flag = false;
    end
    %%% TEMPORARY END

    %%% set Actives Switches according to current actives for focused variables
    DEqns_comp = SetActivesSwitch(DEqns_comp, TSfocusSet, deqnIxMap, numExt, actsIxMap, regAllActs);

    if verboseTog
        fprintf('\n')
        if diagnosticOn
            disp('**********************************************************************')
        end
        fprintf('AttEst:  Beginning analysis of sample interval %i\n',regAbsId);
        if pos > numRegs % replace this later with explicitly known cycle number
            disp( '          ... another time around a cycle')
        end
        fprintf('          Time interval is [%.2f, %.2f] relative, [%.2f, %.2f] absolute\n',tlo,thi,tlo+t0,thi+t0);
        fprintf('            %s is active\n',regActNames{:});
        if regPos == numRegs % then first regime analyzed
            boundedStr = ' (before bounds checking)';
        else
            boundedStr = '';
        end
        fprintf(['          ' varnames{focusVarIx} ' target at t2 is [%.3f, %.3f]' boundedStr '\n'],Vend(1),Vend(2));
        disp('          Starting compression estimate')
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Compression Estimate %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     compressResult = CompressEst(focusVarIx, varnames, primaryEqnIx, DEqns_comp, DEpars, derivInfo, dt, tlo, thi, ...
%         vData, tData, ignoreActs_flag, verboseTog);
    if attEstParams.noShooting
        compressResult = [1 0 0 0];
    else
        compressResult = CompressEst(focusVarIx, varnames, primaryEqnIx, DEqns_comp, DEpars, derivInfo, dt, tlo, thi, ...
                vData, tData, false, verboseTog);
    end
    % the values returned here are used by NeighbHdEst()

    PhiMaps(regAbsId,:) = [compressResult(1) compressResult(2)];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Neighbourhood Estimate at t1 %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isempty(compressResult) % then proceed (arbitrary check for error status from CompressEst)
        % Vend = V(t2) target interval (or at t1 if noShooting)
        Vt1 = vData(1,focusVarIx); % Vt1 = V(tlo) from known orbit
        Vt2 = vData(poshi-poslo+1,focusVarIx);
        if verboseTog & ~attEstParams.noShooting
            disp('          Compression estimate completed')
            disp(['          Estimating ' varnames{focusVarIx} ' attractor cross-section at absolute time ' num2str(t0+tlo,'%.3f')])
            fprintf('          %s(t1) = %.3f, %s(t2) = %.3f, coefficient = %.3f, offset = %.3f\n', varnames{focusVarIx}, Vt1, ...
                varnames{focusVarIx}, Vt2, compressResult(1), compressResult(2) )
        end
        if diagnosticOn & ~verboseTog
            disp(' ')
            disp(   '  *********************')
            fprintf('         Regime %i\n',regAbsId)
            disp(   '  *********************')
            disp(' ')
        end
        if regPos == 2
            pause(0.1)
        end
        neighbResult = NeighbHdEst(TSfocusSet, varnames, deqnIxMap, allDEixMap, DEqns_comp, DEpars, trivialCandActs, ...
            VdepInfo, varBounds, tScaleThresh, tlo, thi, vData, dt, compressResult, numExt, caOpts, ...
            Vend, attEstParams.Vinres, attEstParams.dVmax, dScaleThresh, regAllActs, nextRegAllActs, absGamIxMap, ...
            actsIxMap, local_acts_t1, local_acts_t2, qsPars, firstInterval_flag, ignoreActs_flag, verboseTog);
        if isempty(neighbResult{1})
            disp('AttEst:  Operation cancelled. Returning ...')
            return
        end
        Vinterval    = neighbResult{1};
        verboseTog   = neighbResult{3}; % possibly updated from user interrupt
        diagnosticOn = neighbResult{4}; % possibly updated from user interrupt
        if firstInterval_flag
            if sum(initVinterval == neighbResult{2}) ~= 2
                initVinterval = neighbResult{2}; % this had been trimmed to maintain Actives at t2
                if verboseTog
                    fprintf('          Initial target V interval at time %.3f trimmed down to maintain Actives\n',t0+thi)
                end
            end
        end
        if verboseTog
            if attEstParams.noShooting
                disp('          Domain of Validity cross-section completed')
                fprintf(['          ' varnames{focusVarIx} ' domain of validity interval = [ %.3f, %.3f ]\n'],Vinterval(1), Vinterval(2));
            else
                disp('          Attractor estimate cross-section completed')
                fprintf(['          ' varnames{focusVarIx} ' attractor estimate interval = [ %.3f, %.3f ]\n'],Vinterval(1), Vinterval(2));
            end
        end
        regimes{regPos}{6} = Vinterval;
        Vend = Vinterval; % for use with the next epoch
    else
        beep
        fprintf('          Compression estimate failed for regime %i\n',regAbsId)
        return
    end % if
end % for pos

% set this interval in dummy final regime (numRegsTot) to be the initial target interval (for plotting purposes)
regimes{numRegsTot}{6} = initVinterval;

% show overall contraction or expansion of orbits from t0 to t_period on the t0 cross section
PhiMapTotal = [1 0]; % intial value for coefficient, offset
for epIx = 1:numRegs % this time go forwards in time
    % Concatentation of Phi_2 to Phi_1: " a2( a1 v + b1 ) + b2 = ( a2 a1 ) v + ( a2 b1 + b2 ) "
    PhiMapTotal = PhiMaps(epIx,1) * PhiMapTotal + [0 PhiMaps(epIx,2)];
end

if ~verboseTog & diagnosticOn_def % then diagnosticOn == false but _def == true
    fprintf('\n')
end

V_Phi_mapped = PhiMapTotal(1) * Vinterval + PhiMapTotal(2);
if ~attEstParams.noShooting
	fprintf('AttEst:  Attractor cross-section [ %.6f, %.6f ] in variable %s at initial t0 = %.3f\n', Vinterval(1), Vinterval(2), varnames{focusVarIx}, t0)
	fprintf('                         maps to [ %.6f, %.6f ] at final t = %.3f (forwards in time)\n', V_Phi_mapped(1), V_Phi_mapped(2), t0+Tperiod)
	fprintf('                         by the `Phi_total` map %s |-> %.5f * %s + %.5f\n', varnames{focusVarIx}, PhiMapTotal(1), varnames{focusVarIx}, PhiMapTotal(2))
	fprintf('         Final %s(%.3f) = %.6f for comparison with concatenated linear map `Phi_total` for accumulated error.\n', varnames{focusVarIx}, t0+Tperiod, lastVval)
end

success = true; % will have already returned from function if failed
attRes = {success, regimes, numRegs, numRegsTot, dt};
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                 MAIN ESTIMATOR FUNCTIONS                   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% B calculation (in the tstep for loop, in CompressEst() below)
% e.g. (actsSw is actives switch) ... using the theory in FN_theory.pdf
%   B =   actsSw(act_miX) * mdot(V,m) * Psi_m / m
%       + actsSw(act_niX) * ndot(V,n) * Psi_n / n
%       + actsSw(act_siX) * sdot(V,s) * Psi_s / s
% But have to also add the zero-derivative terms which were ignored above
%     ... the formal definition of B is
%   B = sum { actSw(i) * d var_i / dt * Psi_i / var_i }
%     where i ranges over all inputs to the current equation in focus
% Note: differential equation data structure, per equation entry:
%   { subjIx, {actSw + gam1terms}, {actSw + gam2terms} }
%
function result = CompressEst(focusVarIx, varnames, eqnIx, DEqns, DEpars, derivInfo, dt, ...
              t1, t2, varData, tData, ignoreActs_flag, verboseTog)
% Function to calculate compression estimate (using reduced model as determined by DSSRT)
%  of Phi on original orbit, and judge error `x`
% (c) 2003 Robert Clewley
global DE_NAMEi DE_GAMMA1i DE_GAMMA2i DE_ACTSWi DE_ISTAUFFILEi DE_TAURECIPi DE_GAMVARNAMEi ...
    DE_GAMVARPOWi DE_ISTGTFFILEi DE_TARGETi DE_INTVARNAMEi DE_INTVARPOWi diagnosticOn

result = []; % initial, temporary value (error status)

%%% Initializations
Vt1 = varData(1,focusVarIx);
Tperiod = t2-t1;
Vt2 = varData(length(tData),focusVarIx);
primaryDE = DEqns{eqnIx}; % compiled version

%%% Calculate instantaneous target of focusVarIx at t1
V0t1 = GetQsfpVal(primaryDE,varData(1,:));

%%% Determine number of steps in time data
steps=round((t2-t1)/dt)+1;
lenVarData = length(varData(:,1)); % arbitrarily check 1st column
if steps > lenVarData
    if verboseTog | diagnosticOn
        fprintf('CompressEst:  Debugging info: steps > lenVarData by %i (temporarily fixed)\n',steps-lenVarData);
    end
    steps = lenVarData; % stops rounding errors when funny step sizes were used in integration data
end

%%% Sum of tau_recip's (e.g. sum of conductances) over whole period
Gamma1 = primaryDE{DE_GAMMA1i};
Gamma2 = primaryDE{DE_GAMMA2i};
numGam1Terms = derivInfo{1};
numGam2Terms = derivInfo{2};

%%% Integral as finite sum over each partial period
integral_sum = 0;
sumData      = zeros(numGam1Terms,steps);
% Build sumData array
for g1t=1:numGam1Terms
    if Gamma1{g1t}{DE_ACTSWi} | ignoreActs_flag % then term is active, so it contributes (or we're ignoring actives set)
        for yi=1:steps
            if Gamma1{g1t}{DE_INTVARNAMEi} > 0
                sumData(g1t,yi) = GetTauRVal(Gamma1{g1t}{DE_ISTAUFFILEi}, Gamma1{g1t}{DE_TAURECIPi}, varData(yi, Gamma1{g1t}{DE_GAMVARNAMEi}), Gamma1{g1t}{DE_GAMVARPOWi}) ...
                    * varData(yi, Gamma1{g1t}{DE_INTVARNAMEi})^Gamma1{g1t}{DE_INTVARPOWi};
            else
                sumData(g1t,yi) = GetTauRVal(Gamma1{g1t}{DE_ISTAUFFILEi}, Gamma1{g1t}{DE_TAURECIPi}, varData(yi, Gamma1{g1t}{DE_GAMVARNAMEi}), Gamma1{g1t}{DE_GAMVARPOWi});
            end
        end
    end % else make no entry (remains zero)
end

partialSum = 0; % running partial sum
for tstep=steps:-1:1 % go from end backwards (more efficient re-use of partial sums)
    % Calculate partial sum
    partialSum = partialSum + sum(sumData(:,tstep));
    Ghat = partialSum * dt;
    varDataLine = varData(tstep,:);
    PsiVals = GetPsiVals(varDataLine,primaryDE,false); % get signed Psi values (i.e. par unsigned = false)
    derivVals = GetDerivVals(varDataLine,DEqns,derivInfo);
    B = 0;
    for derivIx=1:numGam1Terms
        g1term = Gamma1{derivIx};
        varVal = varDataLine(g1term{DE_GAMVARNAMEi});
        B = B + g1term{DE_ACTSWi} * derivVals(derivIx) * PsiVals(derivIx) / varVal;
    end
    for derivIx=1:numGam2Terms
        g2term = Gamma2{derivIx};
        varVal = varDataLine(g2term{DE_GAMVARNAMEi});
        B = B + g2term{DE_ACTSWi} * derivVals(derivIx+numGam1Terms) * PsiVals(derivIx+numGam1Terms) / varVal;
    end
    % Add to the integral that measures compression over the regime
    integral_sum = integral_sum + B * exp(-Ghat) * dt;
end

% At end of the tstep for loop, Ghat_t1_t2 is the final Ghat calculated.
Ghat_t1_t2 = Ghat;

% Contraction estimate -> contraction variable value estimate at t2
Xt2_est = (Vt1-V0t1)*exp(-Ghat_t1_t2) - integral_sum;

% Calculate instantaneous target of focusVarIx at t2
V0t2 = GetQsfpVal(primaryDE,varData(steps,:)); % steps = t2 position

% Corrected V(t2) value
Vt2_est = V0t2 + Xt2_est;
x = Vt2 - Vt2_est;

% offset (second return argument) incorporates regime correction x
result = [exp(-Ghat_t1_t2), V0t2 - V0t1*exp(-Ghat_t1_t2) - integral_sum + x, V0t1, V0t2];
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function result = NeighbHdEst(focusSet, varnames, deqnIxMap, allDEixMap, DEqns, DEpars, trivialCandActs, ...
            VdepInfo, varBounds, tScaleThresh, t1, t2, varData, dt, compressResult, numExt, ...
            caOpts, Vend, Vinres, dVmax, dScaleThresh, acts_t1, acts_t2, absGamIxMap, ...
            absActsIxMap, local_acts_t1, local_acts_t2, qsPars, firstInterval_flag, ignoreActs_flag, verboseTog)
% Function to calculate estimate of inverse Phi map's domain (using reduced regime model determined by DSSRT)
%  over one regime. (Compare with numerically obtained neighbourhood by regular shooting.)
% (c) 2003-2004, Robert Clewley
%   May 2004: Updated to work with reduced dynamical regimes
%   Jan 2004: Added check for changes in all focusSet variables' active sets as part of attractor estimate search
%   Dec 2003: Added general term structure for differential equations ... to compute attractor estimate for non-voltage variables
global DE_NAMEi DE_GAMMA1i DE_GAMMA2i DE_ACTSWi DE_ISTAUFFILEi DE_TAURECIPi DE_GAMVARNAMEi ...
    DE_GAMVARPOWi DE_ISTGTFFILEi DE_TARGETi DE_INTVARNAMEi DE_INTVARPOWi LARGEBOUND SMALLBOUND
global figHand maxVresFrac diagnosticOn diagnosticOn_def

result = {[],[], 0, 0}; % default value (error status)

%%%%% Initializations
primaryFocIx = focusSet(1);
numFocused   = length(focusSet);
lenVarData   = length(varData(:,1)); % arbitrarily check 1st column
expGhat      = compressResult(1);
offset       = compressResult(2);
% Thus, Phi map over the regime is Vt2 = Vt1 * expGhat + offset
Vknown_t1    = varData(1,primaryFocIx);
Vknown_t2    = varData(lenVarData,primaryFocIx);

% Precision check for compression estimate (can avoid needless searching)
if expGhat <= SMALLBOUND % DSSRT precision level
    % then assume variable is extremely flat here (and extremely compressed), and Vt1 is best approximated by 'offset', without need for attractor mapping
    if verboseTog | diagnosticOn
        disp('          No calculations done - compression ratio too high and orbit is very compressed in this regime')
    end
    if firstInterval_flag % make sure we don't return an out-of-bounds target Vend (this code copied from below)
        % note changes to exclusion of lower and upper bound for tests 1 and 4 (compared to code's copy source)
        bdTest = [Vend(1) < varBounds(primaryFocIx,1), Vend(1) >= varBounds(primaryFocIx,2), ...
                Vend(2) <= varBounds(primaryFocIx,1), Vend(2) > varBounds(primaryFocIx,2)];
        if sum(bdTest) ~= 0
            if verboseTog | diagnosticOn
                disp('          Initial target V interval Vend outside variable bounds in first regime analyzed ...')
            end
            if bdTest(1)
                Vend(1) = varBounds(primaryFocIx,1);
            end
            if bdTest(2) % irregular out of bounds
                disp('             Warning: Initial target V interval low entry above high bound!')
                Vend(1) = varBounds(primaryFocIx,1);
            end
            if bdTest(3) % irregular out of bounds
                disp('             Warning: Initial target V interval high entry below low bound!')
                Vend(2) = varBounds(primaryFocIx,2);
            end
            if bdTest(4)
                Vend(1) = varBounds(primaryFocIx,1);
            end
            if verboseTog | diagnosticOn
                disp('          -> Vend adjusted')
            end
        end
    end % else leave Vend alone        
    result = { [max(varBounds(primaryFocIx,1),Vknown_t1-SMALLBOUND), min(varBounds(primaryFocIx,2),Vknown_t1+SMALLBOUND)], Vend, verboseTog, diagnosticOn };
    return
end

% Check legitimacy of initial Vend using binary search, if this is first regime of sequence
if firstInterval_flag % check passed V interval `Vend` at t2, and trim down if necessary
	bdTest = [Vend(1) <= varBounds(primaryFocIx,1), Vend(1) >= varBounds(primaryFocIx,2), ...
            Vend(2) <= varBounds(primaryFocIx,1), Vend(2) >= varBounds(primaryFocIx,2)];
    if sum(bdTest) ~= 0
        if verboseTog | diagnosticOn
            disp('          Initial target V interval Vend outside variable bounds in first regime analyzed ...')
        end
        if bdTest(1)
            Vend(1) = varBounds(primaryFocIx,1);
        end
        if bdTest(2) % irregular out of bounds
            disp('             Warning: Initial target V interval low entry above high bound!')
            Vend(1) = varBounds(primaryFocIx,1);
        end
        if bdTest(3) % irregular out of bounds
            disp('             Warning: Initial target V interval low entry above high bound!')
            Vend(2) = varBounds(primaryFocIx,2);
        end
        if bdTest(4)
            Vend(2) = varBounds(primaryFocIx,2);
        end
        if verboseTog | diagnosticOn
            disp('          -> Vend adjusted')
        end
	end

    fail = [ true, true ]; % assume "fail" state initially, for convenience
    VendOrig = Vend;
    VendEpoint = Vend; % starting end points of the V interval
    VendCpoint = [ varData(lenVarData,primaryFocIx), varData(lenVarData,primaryFocIx) ]; % centre point position of V interval
    % while loop assumes not all focused vars have trivial CandActs, otherwise will end up in an infinite while loop
    for VendIx = 1:2
        infSwitchOn = true; % until search refines to the inner layer
        while fail(VendIx)
            success = ones(1,numFocused); % reset for this test Vend interval
            focusedInfo = PrepareFocVars( tScaleThresh, varData(lenVarData,:), primaryFocIx, Vend(VendIx), focusSet, ...
                DEqns, allDEixMap, infSwitchOn );
%             fastFocVars = focusedInfo.fastFocVars;
            for focusSetIx = 1:numFocused
                if trivialCandActs(focusSetIx) % then nothing can change in the Acts set for this variable
                    continue
                end
                absFocusedIx = focusSet(focusSetIx);
                if absFocusedIx > numExt % only do this loop for external variables
                    continue
                end
                focEqnIx = deqnIxMap(absFocusedIx);
                if focEqnIx == 0
                    continue % then no equation associated, and nothing to check. so skip ahead
                end
                abs_acts_t2 = absActsIxMap{absFocusedIx}( acts_t2{absFocusedIx} );
                %% Prepare estimated varDataLine and find fast V-dependent variables at start of regime
                pestInfo = PrepareEstimates( focusSetIx, primaryFocIx, focusedInfo.varData, VdepInfo, ...
                    tScaleThresh, focusedInfo.focVarSpeeds(focusSetIx), focEqnIx, DEqns, allDEixMap, infSwitchOn );

                %% NEW to April 04: Check that grouping of timescales hasn't changed
%                 focusedInfo_update = PrepareFocVars( tScaleThresh, pestInfo.varDataEst, primaryFocIx, ...
%                     pestInfo.varDataEst(primaryFocIx), focusSet, DEqns, allDEixMap, infSwitchOn );
% %                 if ~isempty(setxor( intersect(focusedInfo.fastFocVars,abs_acts_t2) ,intersect(focusedInfo_update.fastFocVars,abs_acts_t2)))
%                 if ~isempty(setxor( focusedInfo.fastFocVars, focusedInfo_update.fastFocVars ))
%                     success(focusSetIx) = 0;
%                     if diagnosticOn
%                         fprintf('             Grouping of timescales changed for focusSetIx %i\n',focusSetIx)
%                     end
%                     break % for loop
%                 end

                %% Check that set of V-dependent (and fastest) actives doesn't change at t2 due to perturbation
                %  ( order of Psi data (and thus E_b) is the same as the Gam1, Gam2 order )
                E_b_t2    = ScaleThreshSort( GetPsiVals( pestInfo.varDataEst, DEqns{focEqnIx}, true ), dScaleThresh );
                absE_b_t2 = absGamIxMap{focusSetIx}(E_b_t2);
                % don't pass qsPars to this call because we're working at t2!
                match = CompareActs(absE_b_t2, abs_acts_t2, VdepInfo{focusSetIx}, pestInfo.fastestRelVars, [], ...
                    absActsIxMap{absFocusedIx}(local_acts_t2{absFocusedIx}), caOpts);
                % WAS:
%                 match = CompareActs(absE_b_t2, abs_acts_t2, VdepInfo{focusSetIx}, pestInfo.fastestRelVars, [], ...
%                     absActsIxMap{focusSetIx}(local_acts_t2{absFocusedIx}), caOpts);
%            but focusSetIx as arg to absActsIxMap returned too small a set to reference local_acts_t2{}
                if ~match(1) % mismatch, exit `for` loop
                    success(focusSetIx) = 0;
                    break % for loop
                end
            end % for focusSetIx

            if sum(success) == numFocused % then all passed for this Vend
                if Vend(VendIx) == VendOrig(VendIx) % then search is over, passed first time (because initial VendEpoint = Vend)
                    fail(VendIx) = false;
                else % passed, but haven't necessarily narrowed interval small enough
                    VendCpoint(VendIx) = Vend(VendIx);
                    if abs( VendCpoint(VendIx) - VendEpoint(VendIx) ) <= SMALLBOUND % then we're close enough
                        fail(VendIx) = false; % next line will be the while loop, so no need to use 'continue'
                    else % change Vend(VendIx)
                        Vend(VendIx) = ( Vend(VendIx) + VendEpoint(VendIx) ) / 2;
                    end
                end
            else % else one of the focusSetIx failed
                if infSwitchOn && abs(Vend(VendIx) - varData(lenVarData,primaryFocIx)) < Vinres*1.1 % need *1.1 because VendCpoint could be set one Vinres from Vknown_t1, and so this will never quite be true (would if <=) after asymptoting through repeated halving in the failure case -- eventually causing VendCpoint == VendEpoint error in failure case
                    infSwitchOn = false; % so close to orbit that using _infinity values for V-dep vars may be inaccurate
                            % ... and might prevent even the original orbit value from passing the above tests (hence an infinite loop!)
                    VendEpoint(VendIx) = Vend(VendIx);
                    Vend(VendIx) = ( VendCpoint(VendIx) + Vend(VendIx) ) / 2;
                    if VendEpoint(VendIx) == VendCpoint(VendIx)
                        disp('             In first regime analyzed (checks at t2):')
                        disp('             ERROR: VendE == VendC after failing in binary search, while close to Vinres')
                        return
                    end
                    %% DO NOTHING NOW, BECAUSE WE MUST FIND A PASSING POINT NOT ON-ORBIT BEFORE CONTINUING
                else % we'll continue - so change Vend(VendIx)
                    VendEpoint(VendIx) = Vend(VendIx);
                    Vend(VendIx) = ( VendCpoint(VendIx) + Vend(VendIx) ) / 2;
                    if VendEpoint(VendIx) == VendCpoint(VendIx)
                        disp('             In first regime analyzed (checks at t2):')
                        disp('             ERROR: VendE == VendC after failing in binary search')
                        return
                    end
                end
            end % if sum(success)
        end % while fail(VendIx)
    end % for VendIx
end % if firstInterval_flag

firstVend = Vend; % a record of what Vend is, in case it was changed for the first regime of an attractor estimate above
changedVend = [0,0]; % diagnostic

% Adjust Vend (V(t2) target interval) so that it maps back to a valid interval at t1
bounds_t1 = [ max( varBounds(primaryFocIx,1), Vknown_t1-dVmax ), min( varBounds(primaryFocIx,2), Vknown_t1+dVmax ) ];
if expGhat == 1 & offset == 0
    Vend_at_t1 = bounds_t1;
    Vend       = Vend_at_t1;
else
    Vend_at_t1 = (Vend - offset) / expGhat; % map interval back to t1 using Phi inverse
end
if Vend_at_t1(1) > bounds_t1(2) | Vend_at_t1(2) < bounds_t1(1)
    disp(   '          Fatal error! Something messed up at start of regime -- Vend mapped to t1 no longer contains known V(t1)')
    fprintf('            V(t1) = %.4f, Vend_at_t1 = (%.4f, %.4f)\n', Vknown_t1, Vend_at_t1(1), Vend_at_t1(2))
    beep
    return
end
if Vend_at_t1(1) < bounds_t1(1)
    Vend(1) = bounds_t1(1)*expGhat + offset; % still guaranteed to be OK for Actives at t2
    Vend_at_t1(1) = bounds_t1(1);
    if diagnosticOn
        disp('          Target V interval (lo) at t2 trimmed to fit into maximum search interval at t1 (specified by dVmax)')
    end
    changedVend(1) = changedVend(1)+1;
end
if Vend_at_t1(2) > bounds_t1(2)
    Vend(2) = bounds_t1(2)*expGhat + offset; % still guaranteed to be OK for Actives at t2
    Vend_at_t1(2) = bounds_t1(2);
    if diagnosticOn
        disp('          Target V interval (hi) at t2 trimmed to fit into maximum search interval at t1 (specified by dVmax)')
    end
    changedVend(2) = changedVend(2)+1;
end

%%%%% Search for inverse-Phi map's domain
%
% Check points inside inner layer around original orbit at t1.
% If the one-Vinres points fail, we do not need to search in that direction.
% Otherwise, search outwards from 'centre' (i.e. Vknown_t1) in medium-sized steps until we fail
%  The points of failure in each direction will be the maximum search interval for the binary search
%  in the small intervals now known to contain the first point of failure
%
%  stage 1: check inner layer points, stage 2: check points lower than Vknown_t1,
%  stage 3: check points higher than Vknown_t1, stage 4: end (dummy value to save on logic statements in while loop)
if diagnosticOn
    disp('          Starting pre-binary search stages:')
end
maxVres     = max(Vinres, (varBounds(primaryFocIx,2)-varBounds(primaryFocIx,1)) / maxVresFrac);
Vresolution = [ maxVres maxVres ];
numChecks   = floor(abs(Vend_at_t1 - Vknown_t1) ./ Vresolution);
for i=1:2
    if numChecks(i) > 30 % stops numChecks reaching huge values
        if diagnosticOn
            disp('          numChecks capped to 30')
        end
        Vresolution(i) = abs(Vend_at_t1(i) - Vknown_t1) / 30; % still fill the appropriate interval with checks, but fewer of them
        numChecks(i)   = 30;
    end
end
initialCheckVals = [ -Vinres, Vinres ]; % initial value
numInitChecks    = 2; % initial value

% Setup initial checked values for stage 1 (inner layer)
stageOneResult    = [1,1]; % for initial inner layer checks
for i=1:2
	if numChecks(i) == 0 % then no checks at Vresolution steps for that stage
        % but do the inner layer checks fit into Vend_at_t1?
        if abs(Vend_at_t1(i) - Vknown_t1) < Vinres % if not, then do tighter checks
            initialCheckVals(i) = Vend_at_t1(i) - Vknown_t1; % get relative vals to Vknown_t1
        end
        % and are the inner layer checks inside the variable's bounds?
        if initialCheckVals(i) + Vknown_t1 < bounds_t1(1)
            if abs(Vknown_t1 - bounds_t1(i)) <= SMALLBOUND
                stageOneResult(i) = 0;
                numInitChecks = numInitChecks - 1;
            else
                initialCheckVals(i) = bounds_t1(1) - Vknown_t1;
            end
        end
        if initialCheckVals(i) + Vknown_t1 > bounds_t1(2)
            if abs(Vknown_t1 - bounds_t1(i)) <= SMALLBOUND
                stageOneResult(i) = 0;
                numInitChecks = numInitChecks - 1;
            else
                initialCheckVals(i) = bounds_t1(2) - Vknown_t1;
            end
        end
	end
end
if sum(stageOneResult) ~= 0
    stagePositions = [1, 0, 0, 0]; % initial value - for positive entries, when initCheckIx reaches these, then we've started a different stage
else
    stagePositions = [0, 0, 0, 0];
end
otherInitVals    = {[],[]};

% Setup initial checked values for stages 2 and 3 (outer layer), if necessary
for stage=2:3
	if numChecks(stage-1) >= 1 % then there are checks to do in that stage (if only one then that will be done by stage 1)
        % stage 2 is decreasing direction, etc. so sign() function is a trick to get direction
        otherInitVals{stage-1}  = sign(stage-2.5) * [1:numChecks(stage-1)]*Vresolution(stage-1);
        initialCheckVals = [ initialCheckVals, otherInitVals{stage-1} ];
        stagePositions(stage) = numInitChecks + 1;
        numInitChecks = numInitChecks + numChecks(stage-1);
	end
end
stagePositions(4) = numInitChecks+1;
firstFailVals     = [stagePositions(2), stagePositions(3)]; % initialize first low and high fail values (first stage irrelevant)
initialCheckVals  = Vknown_t1 + initialCheckVals;

% Prepare relevant focus set for rest of checks
allAbsActs = cell(1,numFocused);
doFocusSet = [];
for focusSetIx = 1:numFocused
    if trivialCandActs(focusSetIx) % then nothing can change in the Acts set for this variable
        continue
    end
    absFocusedIx = focusSet(focusSetIx);
    focEqnIx = deqnIxMap(absFocusedIx);
    if focEqnIx == 0
        continue % then no equation associated, and nothing to check. so skip ahead
    end
    doFocusSet = [ doFocusSet focusSetIx ];
    if absFocusedIx > numExt % then an internal variable so no actives set needed
        allAbsActs{focusSetIx} = [];
    else
        allAbsActs{focusSetIx} = absActsIxMap{absFocusedIx}( acts_t1{absFocusedIx} );
    end
end

% Do the initial checks
num_steps = 0; % diagnostic, counting number of searches in V. Value is displayed on success in binary search loop
initCheckIx = 1;
stage = 1; % reset this (it was used above)
while initCheckIx <= numInitChecks
    num_steps = num_steps + 1; % diagnostic
    [ispres stage] = ismember(1, (initCheckIx >= stagePositions) .* (stagePositions>0) ); % update 'stage' value
    if ~ispres
        disp('       Internal error. stage position indicator messed up')
        beep
        return
    end
    if stage > 1 % can only possibly skip these stages
        infSwitchOn = true; % use _infinity values outside Vinres 'inner layer'
        if ~stageOneResult(stage-1) % then pass on this stage
            nextStage = stage+1;
            while nextStage <= 4
                if stagePositions(nextStage) == 0 % that stage is not being done
                    nextStage = nextStage + 1;
                else
                    initCheckIx = stagePositions(nextStage); % move on to next stage (if there's one left to do)
                    break % while nextStage loop
                end
            end
            if diagnosticOn
                fprintf('            Skipped stage %i because corresponding stage one value failed\n',stage)
            end
            continue % to next while initCheckIx for next stage
        end
    else % stage = 1 -- inner layer check
        infSwitchOn = false; % don't use _infinity values for V-dependent variables this close to original orbit
    end
    focusedInfo = PrepareFocVars( tScaleThresh, varData(1,:), primaryFocIx, initialCheckVals(initCheckIx), ...
        focusSet, DEqns, allDEixMap, infSwitchOn );
%     if focusedInfo.varData(3) > - 50
%         pause(0.1)
%     end
	for focusSetIx = doFocusSet
        absFocusedIx = focusSet(focusSetIx);
        if absFocusedIx > numExt % then it's an internal variable so don't do this loop
            continue
        end
        focEqnIx = deqnIxMap(absFocusedIx);
        abs_acts_t1 = allAbsActs{focusSetIx};
        %% Prepare estimated varDataLine and find fast V-dependent variables at start of regime
        pestInfo = PrepareEstimates( focusSetIx, primaryFocIx, focusedInfo.varData, VdepInfo, tScaleThresh, ...
            focusedInfo.focVarSpeeds(focusSetIx), focEqnIx, DEqns, allDEixMap, infSwitchOn );

        %% NEW to April 04: Check that grouping of timescales hasn't changed
%         focusedInfo_update = PrepareFocVars( tScaleThresh, pestInfo.varDataEst, primaryFocIx, ...
%             pestInfo.varDataEst(primaryFocIx), focusSet, DEqns, allDEixMap, infSwitchOn );
% %         if ~isempty(setxor( intersect(focusedInfo.fastFocVars,abs_acts_t1),intersect(focusedInfo_update.fastFocVars,abs_acts_t1)))
%         if ~isempty(setxor( focusedInfo.fastFocVars, focusedInfo_update.fastFocVars ))
%             match(1) = false;
%             match(2) = 100;
%         else
            match(1) = true;
            match(2) = 0;
%         end

        %% Check that set of V-dependent (and fastest) actives doesn't change at t1 due to perturbation
        % order of Psi data (and thus E_b) is the same as the Gam1, Gam2 order
        E_b_t1       = ScaleThreshSort( GetPsiVals( pestInfo.varDataEst, DEqns{focEqnIx}, true ), dScaleThresh );
        absE_b_t1    = absGamIxMap{focusSetIx}(E_b_t1);
        if match(1)
            if ignoreActs_flag
                currFast = intersect(pestInfo.fastestRelVars,abs_acts_t1);
                if ~isempty(currFast)
                    match = CompareActs(absE_b_t1, abs_acts_t1, VdepInfo{focusSetIx}, setdiff(pestInfo.fastestRelVars, currFast), ...
                        qsPars, absActsIxMap{absFocusedIx}(local_acts_t1{absFocusedIx}), caOpts); % changed arg of absActsIxMap from focusSetIx to absFocusedIx
                else
                    match = CompareActs(absE_b_t1, abs_acts_t1, VdepInfo{focusSetIx}, pestInfo.fastestRelVars, ...
                        qsPars, absActsIxMap{absFocusedIx}(local_acts_t1{absFocusedIx}), caOpts); % changed arg of absActsIxMap from focusSetIx to absFocusedIx
                end
            else
                match = CompareActs(absE_b_t1, abs_acts_t1, VdepInfo{focusSetIx}, pestInfo.fastestRelVars, ...
                    qsPars, absActsIxMap{absFocusedIx}(local_acts_t1{absFocusedIx}), caOpts); % changed arg of absActsIxMap from focusSetIx to absFocusedIx
            end
        end

        if ~match(1) % mismatch (nothing to do if pass!)
            if stage == 1 % for the inner layer checks (initCheckIx is either 1 or 2)
                stageOneResult(initCheckIx) = 0;
                if (verboseTog | diagnosticOn) & (~stageOneResult(1) | ~stageOneResult(2))
                    if initCheckIx == 2
                        dirStr = 'high';
                    else
                        dirStr = 'low';
                    end
                    disp(['            Conditions failed for inner-layer check at ' num2str(initialCheckVals(initCheckIx),'%.4f') ' in ' dirStr ' direction'])
                    fprintf('            -> actives set at t1 violated - reason %i\n', match(2))
                    if match(2)==100 & diagnosticOn
                        disp('             (grouping of timescales changed)')
                    end
                    DisplayAbsActs(absE_b_t1, 15, 'E_b_t1', varnames)
                    DisplayAbsActs(abs_acts_t1, 15, 'acts_t1', varnames)
                end
            else % stage >= 2
                firstFailVals(stage-1) = initCheckIx; % mark first failing value's index, in whichever direction
                nextStage = stage+1;
                while nextStage <= 4
                    if stagePositions(nextStage) == 0 % that stage is not being done
                        nextStage = nextStage + 1;
                    else
                        initCheckIx = stagePositions(nextStage); % move on to next stage (if there's one left to do)
                        break % while loop
                    end
                end
                break % for focusSetIx loop, bypassing stage-passing update
            end
        end % if ~match(1)
    end % for focusSetIx
    initCheckIx = 1 + initCheckIx;
    if stage >= 2 % then all focusSetIx's passed match for previous initCheckIx at this stage
        % find out if we did the last initCheckIx of this stage
        nextStage = stage+1;
        while nextStage <= 4
            if stagePositions(nextStage) == 0 % that stage is not being done
                nextStage = nextStage + 1;
            else
                if initCheckIx == stagePositions(nextStage) % then all focusSetIx's passed match for stage 2 or 3
                    firstFailVals(stage-1) = initCheckIx; % indicates didn't fail by setting outside of stage's index range
                    break % while loop
                else % we didn't do the last initCheckIx of this stage, so carry on
                    nextStage = nextStage + 1;
                end
            end
        end
    end % else carry on with stage 1
end % while initCheckIx
if diagnosticOn
    disp('          Pre-binary search stages completed')
end

%%% Determine whether binary searches are necessary, and initialize search intervals
% If we get this far, we know that at least one of the two one-Vinres points passed
% Set Vend to the forward Phi map of the last passing V values, and update binary search interval
if expGhat == 1 & offset == 0
    VendCpoint = [ Vknown_t1, Vknown_t1 ];
    VendEpoint = Vend; % will be an interval at t1 (hasn't been changed since defined)
else
    VendCpoint = [ Vknown_t2, Vknown_t2 ]; % default centre point position of V interval
    VendEpoint = Vend; % starting end points of the V interval
end
if diagnosticOn
    VEt1 = (VendEpoint - offset)/expGhat;
    VCt1 = (VendCpoint - offset)/expGhat;
    disp('          Determining whether binary searches are necessary:')
    fprintf('             Vinres     = %.4f\n',Vinres)
    fprintf('             Vknown_t1  = %.4f\n',Vknown_t1)
    fprintf('             Vknown_t2  = %.4f\n',varData(lenVarData,primaryFocIx))
    fprintf('             Vend(lo)   = %.4f, Vend(hi) = %.4f\n',Vend(1),Vend(2))
    fprintf('             exp(-Ghat) = %.4f\n',expGhat)
    fprintf('             VendE (t1) = (%.4f, %.4f)\n', VEt1(1), VEt1(2))
    fprintf('             VendC (t1) = (%.4f, %.4f)\n', VCt1(1), VCt1(2))
end
failBinSearch = [ 1, 1 ]; % for convenience, assume "fail" state of search initially, i.e. that must continue search
% recall that, by definition of Phi map,   varData(lenVarData,primaryFocIx) == Vknown_t1*expGhat + offset
for stage=2:3
    if stagePositions(stage) > 0
        if firstFailVals(stage-1) == stagePositions(stage) % only the inner layer point might have passed
            if stageOneResult(stage-1) % the inner layer point passed
                % sign of Vinres = sign(stage-2.5) ... i.e. stage 2 is decreasing direction, etc.
                VendCpoint(stage-1) = initialCheckVals(stage-1)*expGhat + offset; % set centre point to last passed value
                VendEpoint(stage-1) = (initialCheckVals(stage-1) + sign(stage-2.5)*Vresolution(stage-1))*expGhat + offset; % a known failing point
                Vend(stage-1) = ( VendEpoint(stage-1) + VendCpoint(stage-1) ) / 2; % test this val next, in binary search
                changedVend(stage-1) = changedVend(stage-1)+1;
                if diagnosticOn
                    fprintf('          Stage %i attempted and only inner layer point at %.4f passed\n',stage,initialCheckVals(stage-1))
                end
            else % none passed in this direction from Vknown_t1 - but we HAVE to find a closer passing point
                % leave centre point on orbit
                VendEpoint(stage-1) = (initialCheckVals(stage-1) + sign(stage-2.5)*Vinres)*expGhat + offset; % set end point to a known failing point
                Vend(stage-1) = ( VendEpoint(stage-1) + VendCpoint(stage-1) ) / 2; % test this val next, in binary search
                changedVend(stage-1) = changedVend(stage-1)+1;
                if diagnosticOn
                    fprintf('          Stage %i attempted and not even inner layer point at %.4f passed\n',stage,initialCheckVals(stage-1))
                end
            end
        elseif firstFailVals(stage-1) > stagePositions(stage) % then two cases still need to be distinguished
            nextStage = stage+1;
            while nextStage <= 4
                if stagePositions(nextStage) == 0 % that stage was not done
                    nextStage = nextStage + 1;
                else
                    if firstFailVals(stage-1) < stagePositions(nextStage) % then actually failed during stage
                        if diagnosticOn
                            fprintf('          Failed during stage %i\n',stage)
                        end
                        % so map last passing value to VendCpoint at t2 for binary search
                        VendCpoint(stage-1) = initialCheckVals(firstFailVals(stage-1)-1)*expGhat + offset;
                        % end point should be the first failing value
                        VendEpoint(stage-1) = initialCheckVals(firstFailVals(stage-1))*expGhat + offset;
                        Vend(stage-1)       = ( VendEpoint(stage-1) + VendCpoint(stage-1) ) / 2; % test this val next, in binary search
                        changedVend(stage-1) = changedVend(stage-1)+1;
                    else % firstFailVals(stage-1) == stagePositions(nextStage) only indicates that whole stage passed
                        if diagnosticOn
                            fprintf('          Whole stage %i passed\n',stage)
                        end
                        % Vend(stage-1) doesn't need to be altered (assumes it wasn't tested due to perfectly spaced test intervals reaching that far)
                        VendCpoint(stage-1) = initialCheckVals(firstFailVals(stage-1)-1)*expGhat + offset; % last passing value
                        VendEpoint(stage-1) = Vend(stage-1); % end of searchable interval (unknown pass/fail)
                    end
                    break % while loop
                end
            end % while
        end
    else % stage was not even attempted ... only the inner layer point might have passed
        if stageOneResult(stage-1) % the inner layer point passed
            % sign of Vinres = sign(stage-2.5) ... i.e. stage 2 is decreasing direction, etc.
            if abs(Vknown_t1 - Vend_at_t1(stage-1)) <= SMALLBOUND % a very small difference - so don't do binary search
                failBinSearch(stage-1)=0;
                if diagnosticOn
                    fprintf('          Stage %i not attempted and only inner layer point at %.4f passed\n',stage,initialCheckVals(stage-1))
                    disp(   '           -> target interval mapped to t1 is close enough to V(t1) - not attempting binary search')
                end
            else
                % leave VendCpoint on orbit: setting to initialCheckVals(...) could make it = Vend_at_t1(...) if Vend_at_t1 was already in inner layer
                VendEpoint(stage-1) = Vend(stage-1);
                % leave Vend(stage-1) at VendEpoint
                if diagnosticOn
                    fprintf('          Stage %i not attempted and only inner layer point at %.4f passed\n',stage,initialCheckVals(stage-1))
                end
            end
        else % none passed in this direction from Vknown_t1 - but we HAVE to find a closer passing point
            % leave centre point on orbit
            VendEpoint(stage-1) = initialCheckVals(stage-1)*expGhat + offset; % set end point to a known failing point
            Vend(stage-1) = ( VendEpoint(stage-1) + VendCpoint(stage-1) ) / 2; % test this val next, in binary search
            changedVend(stage-1) = changedVend(stage-1)+1;
            if diagnosticOn
                fprintf('          Stage %i not attempted: inner layer point at %.4f also failed\n',stage,initialCheckVals(stage-1))
            end
        end
    end
end

%%% Now do binary search, one in each direction away from known Vt2, but checking Actives at t1
if diagnosticOn
    disp(' ')
    disp('          Beginning binary search:')
end
if expGhat == 1 & offset == 0
    % update it from known first failing values for stages 2 and 3
    Vend_at_t1 = VendEpoint;
%     Vend = Vend_at_t1; % do all tests on Vend at t1
    VendCpoint = [Vknown_t1 Vknown_t1];
else % carry on as before
    Vend_at_t1 = (Vend - offset) / expGhat; % map interval back to t1 using Phi inverse
end
if Vend_at_t1(1) > bounds_t1(2) | Vend_at_t1(2) < bounds_t1(1) | Vknown_t1 < Vend_at_t1(1) | Vknown_t1 > Vend_at_t1(2)
    disp(   '          Fatal error! Something messed up -- Vend mapped to t1 no longer contains known V(t1)')
    fprintf('            V(t1) = %.4f, Vend_at_t1 = (%.4f, %.4f)\n', Vknown_t1, Vend_at_t1(1), Vend_at_t1(2))
    beep
    return
end
prev_Vend = [ -LARGEBOUND, LARGEBOUND ];
if diagnosticOn
    if VendEpoint(1) == VendCpoint(1)
        disp('            VendEpoint(1) == VendCpoint(1) at beginning of binary search (not necessarily an error)')
    end
    if VendEpoint(2) == VendCpoint(2)
        disp('            VendEpoint(2) == VendCpoint(2) at beginning of binary search (not necessarily an error)')
    end
    fprintf( '            ( number of searches before binary search = %i )\n',num_steps)
end
% while loop assumes not all focused vars have trivial CandActs, otherwise will end up in an infinite while loop
% VendIx == 1 => stage 2 (decreasing from Vknown_t1), VendIx == 2 => stage 3 (increasing from Vknown_t1)
num_steps_before = num_steps; % diagnostic
numBinSteps = [0 0];
dirStr = {'low','high'};
for VendIx = 1:2
    if diagnosticOn
        disp(['          Starting ' dirStr{VendIx} ' direction ...'])
    end
    lastFailInfo = cell(1,numFocused); % for most recent failure during the binary search
    while failBinSearch(VendIx)
        if abs(prev_Vend(VendIx) - Vend(VendIx)) <= SMALLBOUND % then haven't made significant change in search, so quit it
            failBinSearch(VendIx) = 0;
            if diagnosticOn
                disp('           No significant change made between successive steps - stopping binary search in this direction')
                fprintf('             after %i steps ( %i total )\n', numBinSteps(VendIx), num_steps)
                disp(['            Done in ' dirStr{VendIx} ' direction'])
            end
            continue % while
        end
        prev_Vend(VendIx) = Vend(VendIx);
        num_steps = num_steps + 1; % diagnostic
        numBinSteps(VendIx) = numBinSteps(VendIx) + 1; % diagnostic
        if diagnosticOn
            drawnow
            if get(figHand,'CurrentCharacter') == 27
                disp(' ')
                ButtonName = questdlg('Quit or change verbosity?','AttEst interrupted','Quit','Change','Continue','Continue');
                switch ButtonName
                    case 'Quit'
                        disp('AttEst:  Attractor estimate operation cancelled by user')
                        result{5} = -1; % cancel code
                        return
                    case 'Change'
                        verboseTog   = ~verboseTog;
                        diagnosticOn = diagnosticOn_def & verboseTog;
                        disp('AttEst:  Output verbosity changed. Press any key over main window to continue ...')
                    case 'Continue'
                        % do nothing
                        disp('AttEst:  Press any key over main window to continue ...')
                end
                %%%% Wait for user command over figure
                keyokay = false;
                while ~keyokay
                    w = 0;
                    while w~=1
                        w = waitforbuttonpress;
                    end
                    keypress = get(figHand,'CurrentCharacter');
                    if ischar(keypress) & keypress ~= 27 % ensures ESCAPE no longer current char
                        keyokay = true;
                    end
                end
                if ~verboseTog
                    disp('         Continuing ...')
                    if diagnosticOn_def
                        fprintf('         ')
                    end
                end
            end
            if numBinSteps(VendIx) > 30
                beep
                disp('         Too many binary search steps (30) - something`s probably going wrong')
                disp('          (press ESCAPE to break out)')
                pause(5)
            end
        end
        success = ones(1,numFocused);
        Vend_at_t1(VendIx) = (Vend(VendIx) - offset) / expGhat; % map interval back to t1 using Phi inverse
        % now restrict it to variable's bounds
        if Vend_at_t1(1) > bounds_t1(2) | Vend_at_t1(2) < bounds_t1(1)
            disp(  ['          Fatal error! Something messed up in ' dirStr{VendIx} ' direction'])
            disp(   '            Vend mapped to t1 no longer contains known V(t1)')
            fprintf('            V(t1) = %.4f, Vend_at_t1 = (%.4f, %.4f)\n', Vknown_t1, Vend_at_t1(1), Vend_at_t1(2))
            beep
            return
        end
        if Vend_at_t1(1) < bounds_t1(1)
            Vend(1) = bounds_t1(1)*expGhat + offset; % still guaranteed to be OK for Actives at t2
            Vend_at_t1(1) = bounds_t1(1);
            if verboseTog | diagnosticOn
                disp('          Target V interval (lo) at t2 trimmed to fit into maximum search interval at t1 (specified by dVmax)')
                fprintf('            V(t1) = %.4f, Vend_at_t1 = (%.4f, %.4f)\n', Vknown_t1, Vend_at_t1(1), Vend_at_t1(2))
            end
        end
        if Vend_at_t1(2) > bounds_t1(2)
            Vend(2) = bounds_t1(2)*expGhat + offset; % still guaranteed to be OK for Actives at t2
            Vend_at_t1(2) = bounds_t1(2);
            if verboseTog | diagnosticOn
                disp('          Target V interval (hi) at t2 trimmed to fit into maximum search interval at t1 (specified by dVmax)')
                fprintf('            V(t1) = %.4f, Vend_at_t1 = (%.4f, %.4f)\n', Vknown_t1, Vend_at_t1(1), Vend_at_t1(2))
            end
        end
        if abs(Vend_at_t1(VendIx)-Vknown_t1) <= Vinres % we're in the inner layer
            infSwitchOn = false;
        else
            infSwitchOn = true;
        end
        focusedInfo = PrepareFocVars( tScaleThresh, varData(1,:), primaryFocIx, Vend_at_t1(VendIx), focusSet, ...
            DEqns, allDEixMap, infSwitchOn );
        for focusSetIx = doFocusSet
            absFocusedIx = focusSet(focusSetIx);
            if absFocusedIx > numExt % then it's an internal variable so don't do this loop
                continue
            end
            focEqnIx = deqnIxMap(absFocusedIx);
            abs_acts_t1 = allAbsActs{focusSetIx};
            pestInfo  = PrepareEstimates( focusSetIx, primaryFocIx, focusedInfo.varData, VdepInfo, tScaleThresh, ...
                focusedInfo.focVarSpeeds(focusSetIx), focEqnIx, DEqns, allDEixMap, infSwitchOn );

            %% NEW to April 04: Check that grouping of timescales hasn't changed
%             focusedInfo_update = PrepareFocVars( tScaleThresh, pestInfo.varDataEst, primaryFocIx, ...
%                 pestInfo.varDataEst(primaryFocIx), focusSet, DEqns, allDEixMap, infSwitchOn );
% %             if ~isempty(setxor( intersect(focusedInfo.fastFocVars,abs_acts_t1), intersect(focusedInfo_update.fastFocVars,abs_acts_t1)))
%             if ~isempty(setxor( focusedInfo.fastFocVars, focusedInfo_update.fastFocVars ))
%                 match(1) = false;
%                 match(2) = 100;
%             else
                match(1) = true;
                match(2) = 0;
%             end

            %% Check that set of V-dependent (and fastest) actives doesn't change at t1 due to perturbation
            % ( order of Psi data (and thus E_b) is the same as the Gam1, Gam2 order )
            E_b_t1    = ScaleThreshSort( GetPsiVals( pestInfo.varDataEst, DEqns{focEqnIx}, true ) , dScaleThresh );
            absE_b_t1 = absGamIxMap{focusSetIx}(E_b_t1);

            if match(1)
                if ignoreActs_flag
                    currFast = intersect(pestInfo.fastestRelVars,abs_acts_t1);
                    if ~isempty(currFast)
                        match = CompareActs(absE_b_t1, abs_acts_t1, VdepInfo{focusSetIx}, setdiff(pestInfo.fastestRelVars, currFast), ...
                            qsPars, absActsIxMap{absFocusedIx}(local_acts_t1{absFocusedIx}), caOpts); % changed arg of absActsIxMap from focusSetIx to absFocusedIx
                    else
                        match = CompareActs(absE_b_t1, abs_acts_t1, VdepInfo{focusSetIx}, pestInfo.fastestRelVars, ...
                            qsPars, absActsIxMap{absFocusedIx}(local_acts_t1{absFocusedIx}), caOpts); % changed arg of absActsIxMap from focusSetIx to absFocusedIx
                    end
                else
                    match = CompareActs(absE_b_t1, abs_acts_t1, VdepInfo{focusSetIx}, pestInfo.fastestRelVars, ...
                        qsPars, absActsIxMap{absFocusedIx}(local_acts_t1{absFocusedIx}), caOpts); % changed arg of absActsIxMap from focusSetIx to absFocusedIx
                end
            end

            if ~match(1) % then failed
                lastFailInfo{focusSetIx} = {absE_b_t1, abs_acts_t1, match(2)};
                success(focusSetIx) = 0;
                break % for focusSetIx
            end
        end % for focusSetIx

        if sum(success) == numFocused % then all passed for this Vend
            % but haven't necessarily narrowed interval to resolution of Vinres, in which case change Vend(VendIx)
            VendCpoint(VendIx) = Vend(VendIx);
            VC_at_t1(VendIx) = ( VendCpoint(VendIx) - offset ) / expGhat;
            VE_at_t1(VendIx) = ( VendEpoint(VendIx) - offset ) / expGhat;
            if abs( VC_at_t1(VendIx) - VE_at_t1(VendIx) ) <= SMALLBOUND % then we're close enough
                failBinSearch(VendIx) = false;
                if diagnosticOn
                    fprintf('            Close enough: total number of searches = %i\n',num_steps)
                end
                if verboseTog
                    if ~isempty(lastFailInfo{focusSetIx})
                        fprintf('            Actives set at t1 violated - reason %i\n',lastFailInfo{focusSetIx}{3})
                        if lastFailInfo{focusSetIx}{3}==100 & diagnosticOn
                            disp('             (grouping of timescales changed)')
                        end
                        DisplayAbsActs(lastFailInfo{focusSetIx}{1}, 15, 'E_b_t1', varnames)
                        DisplayAbsActs(lastFailInfo{focusSetIx}{2}, 15, 'acts_t1', varnames)
                    else
                        disp(   '            Reached precision limit')
                    end
                    disp(['            Done in ' dirStr{VendIx} ' direction'])
                end
                continue % while loop
            else
                Vend(VendIx) = ( Vend(VendIx) + VendEpoint(VendIx) ) / 2;
                changedVend(VendIx) = changedVend(VendIx)+1;
            end
        else % failed for some one of the focusSetIx
            if infSwitchOn && abs(Vend_at_t1(VendIx) - Vknown_t1) < Vinres*1.1 % need *1.1 because VendCpoint could be set one Vinres from Vknown_t1, and so this will never quite be true (would if <=) after asymptoting through repeated halving in the failure case -- eventually causing VendCpoint == VendEpoint error in failure case
                infSwitchOn = false; % so close to orbit that using _infinity values for V-dep vars may be inaccurate
                        % ... and might prevent even the original orbit value from passing the above tests (hence an infinite loop!)
                VendEpoint(VendIx) = Vend(VendIx);
                Vend(VendIx) = ( VendCpoint(VendIx) + Vend(VendIx) ) / 2;
                changedVend(VendIx) = changedVend(VendIx)+1;
                if VendEpoint(VendIx) == VendCpoint(VendIx)
                    fprintf('             ERROR: VendE == VendC after failing in binary search, in inner layer, after %i binary search steps\n',num_steps-num_steps_before)
                    return
                end
                % DO NOTHING BECAUSE WE MUST FIND A PASSING POINT NOT ON-ORBIT BEFORE CONTINUING
            else % we'll continue - so change Vend(VendIx)
                VendEpoint(VendIx) = Vend(VendIx);
                Vend(VendIx) = ( VendCpoint(VendIx) + Vend(VendIx) ) / 2;
                changedVend(VendIx) = changedVend(VendIx)+1; % diagnostic
                if VendEpoint(VendIx) == VendCpoint(VendIx)
%                     fprintf('             ERROR: VendE == VendC after failing in binary search after %i binary search steps\n',num_steps-num_steps_before)
%                     return
                    continue % now trap in the <= SMALLBOUND check at beginning of loop
                end
            end
        end % if sum(success)        
    end % while
end % for VendIx

% compose final result to return
if diagnosticOn
    disp('          Finished this regime:')
    fprintf('           Attractor size = %.6f, changed Vend (%i, %i) # of times\n', Vend_at_t1(2) - Vend_at_t1(1), changedVend(1), changedVend(2))
end
result = {Vend_at_t1, firstVend, verboseTog, diagnosticOn};
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                     AUXILIARY FUNCTIONS                    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% copy of this function also exists in StepNet.m
function DEqns = SetActivesSwitch(DEqns, focusSet, deqnIxMap, numExt, actsIxMap, allActs)
% do this only on COMPILED Gamma sets, for storing in either compiled or uncompiled equations
global DE_NAMEi DE_GAMMA1i DE_GAMMA2i DE_ACTSWi DE_ISTAUFFILEi DE_TAURECIPi DE_GAMVARNAMEi ...
    DE_GAMVARPOWi DE_ISTGTFFILEi DE_TARGETi DE_INTVARNAMEi DE_INTVARPOWi

for focusSetIx = 1:length(focusSet)
    thisAbsIx = focusSet(focusSetIx);
    thisDEix = deqnIxMap(thisAbsIx);
    if thisDEix > 0 && thisAbsIx <= numExt % internal variables (even those with DE's) don't have candidate actives
        thisAbsActs = actsIxMap{thisAbsIx}(allActs{thisAbsIx});
        thisGamma1 = DEqns{thisDEix}{DE_GAMMA1i};
        thisGamma2 = DEqns{thisDEix}{DE_GAMMA2i};
        lenG1 = length(thisGamma1);
        lenG2 = length(thisGamma2);
        for g1t=1:lenG1
            DEqns{thisDEix}{DE_GAMMA1i}{g1t}{DE_ACTSWi} = ismember( thisGamma1{g1t}{DE_GAMVARNAMEi}, thisAbsActs );
        end
        for g2t=1:lenG2
            DEqns{thisDEix}{DE_GAMMA2i}{g2t}{DE_ACTSWi} = ismember( thisGamma2{g2t}{DE_GAMVARNAMEi}, thisAbsActs );
        end
    end % else nothing to do
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function result = GetDepInfo(focusVarIx, varnames, inputsIx, DEqns)
%%% Identify which of V's input variables are themselves directly V-dependent (ignores their speeds)
% Method:
%  Find those inputs to V in dEqn, and select those that have V-dependent tau_recip's (controlling their time-constant)
%  Use these to know which entries in varData to change to correspond to the
%   estimated V values during perturbation
global DE_NAMEi DE_GAMMA1i DE_GAMMA2i DE_ACTSWi DE_ISTAUFFILEi DE_TAURECIPi DE_GAMVARNAMEi ...
    DE_GAMVARPOWi DE_ISTGTFFILEi DE_TARGETi DE_INTVARNAMEi DE_INTVARPOWi

VdepList = []; % List of V-dependent variable indices (absolute indices in `varnames`)
Vname = varnames{focusVarIx};
Vinputs = inputsIx{focusVarIx};
for eqnIx=1:length(DEqns)
    dEqn = DEqns{eqnIx};
    eqnSubjName = dEqn{DE_NAMEi};
    [p eqnSubjNameIx] = ismember( eqnSubjName, varnames ); % p unused, known to be true
    if ismember( eqnSubjNameIx, Vinputs )
        % go through Gamma1 and Gamma2 sets to check tau_recip arguments & possible gating variables of the terms
        %  for presence of V (often this is just for one term, in Gamma1 -- e.g. for gating variables
        %   m, n, h in Hodgkin-Huxley)
        Gamma1 = dEqn{DE_GAMMA1i};
        Gamma2 = dEqn{DE_GAMMA2i};
        for g1t=1:length(Gamma1)
            g1term = Gamma1{g1t};
            if strcmp( g1term{DE_GAMVARNAMEi}, Vname ) % then source Deqn has direct dependence on focusVarIx
                [p depVarIx] = ismember( eqnSubjName, varnames ); % get position in varData column, p unused
                VdepList = [ VdepList, depVarIx ];
            end
        end
        for g2t=1:length(Gamma2)
            g2term = Gamma2{g2t};
            if strcmp( g2term{DE_GAMVARNAMEi}, Vname ) % then source Deqn has direct dependence on focusVarIx
                [p depVarIx] = ismember( eqnSubjName, varnames ); % get position in varData column, p unused
                VdepList = [ VdepList, depVarIx ];
            end
        end
    end
end

result = VdepList; % absolute indices
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function thisDE_new = CompileGammas(thisDE_old, DEpars, varnames)
% gam1term order: actSw, filefuncflag, taurecip, var, power, filefuncflag, target [, intvar, power]
% gam2term order: actSw, filefuncflag, taurecip, var, power
global DE_NAMEi DE_GAMMA1i DE_GAMMA2i DE_ACTSWi DE_ISTAUFFILEi DE_TAURECIPi DE_GAMVARNAMEi ...
    DE_GAMVARPOWi DE_ISTGTFFILEi DE_TARGETi DE_INTVARNAMEi DE_INTVARPOWi

Gamma1 = thisDE_old{DE_GAMMA1i}; % to become the index- and par- ready versions
Gamma2 = thisDE_old{DE_GAMMA2i};

if length(Gamma1) + length(Gamma2) == 0
    beep
    disp(' CompileGammas:  FATAL ERROR - no terms in focused equation!')
    return
end

for g1t=1:length(Gamma1)
    g1term = Gamma1{g1t};
    if ~g1term{DE_ISTAUFFILEi}
        result = LookupDEpar( g1term{DE_TAURECIPi}, DEpars);
        if result{2}
            Gamma1{g1t}{DE_TAURECIPi} = result{1};
        else
            disp(' CompileGammas:  Bad DEpar lookup')
            return
        end
    end
    [p ix] = ismember(g1term{DE_GAMVARNAMEi}, varnames);
    Gamma1{g1t}{DE_GAMVARNAMEi} = ix; % replace name with ix in varnames
    if g1term{DE_ISTGTFFILEi} == 0
        result = LookupDEpar( g1term{DE_TARGETi}, DEpars );
        if result{2}
            Gamma1{g1t}{DE_TARGETi} = result{1};
        else
            fprintf('CompileGammas:  Bad DEpar lookup for %s\n',g1term{DE_TARGETi})
            return
        end
    elseif g1term{DE_ISTGTFFILEi} == 2
        [p result] = ismember(Gamma1{g1t}{DE_TARGETi}, varnames);
        if p
            Gamma1{g1t}{DE_TARGETi} = result;
        else
            fprintf('CompileGammas:  Bad variable lookup for %s\n',g1term{DE_TARGETi})
            return
        end
    end
    [p ix] = ismember(g1term{DE_INTVARNAMEi}, varnames);
    Gamma1{g1t}{DE_INTVARNAMEi} = ix; % replace name with ix in varnames (o/w 0 if none present)
end

for g2t=1:length(Gamma2)
    g2term = Gamma2{g2t};
    if ~g2term{DE_ISTAUFFILEi}
        result = LookupDEpar( g2term{DE_TAURECIPi}, DEpars);
        if result{2}
            Gamma2{g2t}{DE_TAURECIPi} = result{1};
        else
            disp(' CompileGammas:  Bad DEpar lookup')
            return
        end
    end
    [p ix] = ismember(g2term{DE_GAMVARNAMEi}, varnames);
    Gamma2{g2t}{DE_GAMVARNAMEi} = ix; % replace name with ix in varnames (o/w 0 if none present)
end

thisDE_new = thisDE_old;
thisDE_new{DE_GAMMA1i} = Gamma1;
thisDE_new{DE_GAMMA2i} = Gamma2;
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function result = PrepareFocVars( focFastScaleThresh, varDataLine, Vix, Vval, focusSet, DEqns, allDEixMap, infSwitchOn )
% This function updates varDataLine with estimates for all fast primaryVar-dependent variables (relative to primary var)
global DE_GAMMA1i DE_GAMMA2i DE_CFACi

varDataLine_temp = varDataLine;
varDataLine_temp(Vix) = Vval; % the primary focused variable perturbed
varDataLine_est = varDataLine_temp; % initial value
fastFocVars  = [];
focVarSpeeds = zeros(1,length(focusSet));
if infSwitchOn
	sumG1prim = SumGamma1(DEqns{allDEixMap(Vix)}{DE_GAMMA1i},varDataLine_temp);
    focVarSpeeds(1) = abs(sumG1prim(1)*DEqns{allDEixMap(Vix)}{DE_CFACi}); % primary var speed
% 	for focIx=2:length(focusSet) % don't include first entry (which is "self", ie. primary var)
	for focIx=1:length(focusSet)
        focVarIx = focusSet(focIx);
        % new if statement follows
        if focVarIx ~= Vix
            focDEix  = allDEixMap(focVarIx);
            if focDEix > 0
				sumG1foc = SumGamma1(DEqns{focDEix}{DE_GAMMA1i},varDataLine_temp);
                focVarSpeeds(focIx) = abs(sumG1foc(1)*DEqns{focDEix}{DE_CFACi}); % foc var speed at original orbit except w/ primary var perturbed
				if focVarSpeeds(focIx) >= focVarSpeeds(1) * focFastScaleThresh
		            fastFocVars = [ fastFocVars focVarIx ];
                    varDataLine_est(focVarIx) = GetQsfpVal(DEqns{focDEix},varDataLine_temp); % for fast focused vars, use their asymptotic value
				end
            end
        end
	end
end

result.varData = varDataLine_est;
result.fastFocVars = fastFocVars;
result.focVarSpeeds = focVarSpeeds;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If either V-dependent variable in a cross-multiplying external x internal term is in the
%  fast set then consider the corresponding external variable the fast variable ...
% In which case we ask: is the cross-multiplying term in the Diff. Eq. in the fast set?
% If so, this is its associated external variable in the same term
% But the assoc. ex. var. might not be V-dependent! so we have a new "fast variable"
%  version of VdepList for use by CompareActs that can include non V-dependent external vars
%  that have fast V-dependent internal vars associated with them
% and replace this entry with that external variable index in absFastVars
% don't bother updating fastVars because we won't use that again
%
function result = PrepareEstimates(focusVarIx, primaryFocusVarIx, varDataLine, VdepInfo, ...
        fastScaleThresh, focVarSpeed, focEqnIx, DEqns, allDEixMap, infSwitchOn)
% This function forms the estimated varData line by getting infinity-values for V-dependent targets
% and finds the fastest V-dependent variables for use by CompareActs
global DE_NAMEi DE_GAMMA1i DE_GAMMA2i DE_CFACi DE_ACTSWi DE_ISTAUFFILEi DE_TAURECIPi DE_GAMVARNAMEi ...
    DE_GAMVARPOWi DE_ISTGTFFILEi DE_TARGETi DE_INTVARNAMEi DE_INTVARPOWi

Gamma1 = DEqns{focEqnIx}{DE_GAMMA1i};
Gamma2 = DEqns{focEqnIx}{DE_GAMMA2i};
varDataLine_est  = varDataLine;
varDataLine_temp = varDataLine; % keep a copy of original while it gets updated

% Form estimated varDataLine using all V-dep vars if they have faster timescales than that of primary var
% (depIx is an absolute index)
VdepList    = VdepInfo{focusVarIx};
varSpeeds   = zeros(1,length(VdepList));
absFastVars = [];
for dv = 1:length(VdepList)
    depIx         = VdepList(dv); % assume that there is a DE associated with this focVar-dependent variable
    sumG1dep      = SumGamma1(DEqns{allDEixMap(depIx)}{DE_GAMMA1i},varDataLine_est);
    depVarSpeed   = abs(sumG1dep(1)*DEqns{allDEixMap(depIx)}{DE_CFACi}); % this is using any var_infinity values established by PrepareFocVars, so is eval'd at perturbed point in *all* variables now
    varSpeeds(dv) = depVarSpeed;
    if infSwitchOn && depVarSpeed >= focVarSpeed * fastScaleThresh
        absFastVars = [absFastVars, depIx];
    end
%     if infSwitchOn && depVarSpeed >= focVarSpeed * fastScaleThresh
%         varDataLine_est(depIx) = GetQsfpVal(DEqns{allDEixMap(depIx)},varDataLine_temp); % estimated value for that variable
%     end % else leave varDataLine alone
end

if infSwitchOn && length(VdepList) > 0
    % DEVELOPMENT OPTION: DO NOT SCALE-THRESH THESE (SEE COMMENTED VERSION ABOVE)
    fastestVars = ScaleThreshSort(varSpeeds, fastScaleThresh);
    % only make absFastestVars from those varSpeeds that are sufficiently faster than focVarSpeed
    absFastestVars = intersect(VdepList(fastestVars),absFastVars); % array for absolute variable indices of fastestVars relative to focused variable
	for g1t = 1:length(Gamma1) % will have to do this for Gamma2 too later on
        [ pres1 pix1 ] = ismember( Gamma1{g1t}{DE_INTVARNAMEi}, absFastVars );
        [ pres2 pix2 ] = ismember( Gamma1{g1t}{DE_INTVARNAMEi}, absFastestVars );
        if pres1
            absFastVars(pix1) = Gamma1{g1t}{DE_GAMVARNAMEi};
        end
        if pres2
            absFastestVars(pix2) = Gamma1{g1t}{DE_GAMVARNAMEi};
        end
	end
	% Remove any repetitions in absFastVars from any cross-mult. internals that have been substituted
%     [aUnique aixs] = unique( absFastVars ); % aUnique is sorted, so no use to us (??)
%     absFastVars = absFastVars( setdiff([1:length(absFastVars)], setdiff([1:length(absFastVars)],aixs) ) );
    absFastVars    = unique( absFastVars ); % order doesn't matter now!
    absFastestVars = unique( absFastestVars );
else
    absFastestVars = [];
end

% Only the fastest vars get subs'd with their Qsfp values, if we're not in the 'inner layer'
if infSwitchOn
	for depIx = VdepList
        if ismember(depIx,absFastestVars)
            varDataLine_est(depIx) = GetQsfpVal(DEqns{allDEixMap(depIx)},varDataLine_temp); % estimated value for that variable
        end
	end
end

result.varDataEst     = varDataLine_est;
result.fastRelVars    = absFastVars; % "relative" to focused variable
result.fastestRelVars = absFastestVars;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    disp('AttEst:  Warning. Division by zero in GetQsfpVal()')
    qsfp = LARGEBOUND;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    varix = g1term{DE_GAMVARNAMEi};
    if varix == 0
        varix = g1term{DE_TARGETi};
        varVal = 1;
    else
        varVal = varDataLine(varix);
    end
    if ~onlyActives || ( onlyActives & g1term{DE_ACTSWi} )
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        tau_recipVal = GetTauRVal(g2term{DE_ISTAUFFILEi},g2term{DE_TAURECIPi},varDataLine(g2term{DE_GAMVARNAMEi}),g2term{DE_GAMVARPOWi});
        sum = sum + tau_recipVal;
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Psi calculation:
% curly braces {} indicate implied summing over set
% for Gamma1 inputs:
%   Psi_i = gi*pi*si^pi* ( Ei*{gj * sj^pj} - {gj * sj^pj * Ej} - {currents} ) / {conductances}^2
%                                 gj_sum         G1j_sum            sumG2            sumG1
%    (where j ~= i are other conductances)
%    for h type (internal) inactivation variable multiplying m, Ei -> Ei*h in above
%  for Gamma2 inputs:
%   Psi_i = gi*pi*si^(pi-1) / {conductances}
%
% A copy of this function is also present inside FuncNet.m
%
function Psi_data = GetPsiVals(varDataLine, thisDE, unsigned)
% calculates Psi influence strength values for each term (Gamma1 and Gamma2) of dEqn differential equation
global DE_NAMEi DE_GAMMA1i DE_GAMMA2i DE_ACTSWi DE_ISTAUFFILEi DE_TAURECIPi DE_GAMVARNAMEi ...
    DE_GAMVARPOWi DE_ISTGTFFILEi DE_TARGETi DE_INTVARNAMEi DE_INTVARPOWi

Gamma1 = thisDE{DE_GAMMA1i}; % index- and par-ready version ('compiled')
numGamma1Terms = length(Gamma1);
Gamma2 = thisDE{DE_GAMMA2i}; % index- and par-ready version ('compiled')
numGamma2Terms = length(Gamma2);
lenPsiData = numGamma1Terms + numGamma2Terms; % number of terms in differential equation
Psi_data = zeros(1,lenPsiData);

sumG1temp = SumGamma1(Gamma1,varDataLine);
sumG1 = sumG1temp(1);
sumG2 = SumGamma2(Gamma2,varDataLine);

if numGamma1Terms == 1
    g1term = Gamma1{1};
    varix = g1term{DE_GAMVARNAMEi};
    if varix == 0
        varix = g1term{DE_TARGETi};
        varVal = 1;
    else
        varVal = varDataLine(varix);
    end
    tau_recipVal = GetTauRVal(g1term{DE_ISTAUFFILEi},g1term{DE_TAURECIPi},varVal,g1term{DE_GAMVARPOWi});
    intIx = g1term{DE_INTVARNAMEi};
    if intIx > 0
        intVarValP = varDataLine(intIx)^g1term{DE_INTVARPOWi};
    else
        intVarValP = 1;
    end
    Psi_data(1) = g1term{DE_GAMVARPOWi}*tau_recipVal*intVarValP * ( - sumG2 ) / (sumG1*sumG1);
else
    varVals = zeros(1,numGamma1Terms);
    intVarValsP = zeros(1,numGamma1Terms); % power already included, where applicable
    for i=1:numGamma1Terms
        varix = Gamma1{i}{DE_GAMVARNAMEi};
        if varix == 0
            varix = Gamma1{i}{DE_TARGETi};
            varVals(i) = 1;
        else
            varVals(i) = varDataLine(varix);
        end
        intIx = Gamma1{i}{DE_INTVARNAMEi};
        if intIx > 0
            intVarValsP(i) = varDataLine(intIx)^Gamma1{i}{DE_INTVARPOWi};
        else
            intVarValsP(i) = 1;
        end
    end
    for pix=1:numGamma1Terms
        g1term = Gamma1{pix};
        tau_recipVal = GetTauRVal(g1term{DE_ISTAUFFILEi},g1term{DE_TAURECIPi},varVals(pix),g1term{DE_GAMVARPOWi});
        if g1term{DE_ISTGTFFILEi} == 2
            % target refers to another external variable, and 'source'
            % index will be in the domain of varDataLine, not VarVals
            targetVal = GetTargVal(g1term{DE_ISTGTFFILEi},g1term{DE_TARGETi},varVals(pix),varDataLine);
        else
            targetVal = GetTargVal(g1term{DE_ISTGTFFILEi},g1term{DE_TARGETi},varVals(pix),varVals);
        end
        gj_sum = sumG1 - tau_recipVal * intVarValsP(pix);
        G1j_sum = sumG1temp(2) - tau_recipVal * intVarValsP(pix) * targetVal;
        Psi_data(pix) = g1term{DE_GAMVARPOWi}*tau_recipVal*intVarValsP(pix) * ( targetVal * gj_sum - G1j_sum - sumG2 ) / (sumG1*sumG1);
    end
end

for pix=1:numGamma2Terms
    g2term = Gamma2{pix};
    % note the '-1' in the power to tau_recipVal
    tau_recipVal = GetTauRVal(g2term{DE_ISTAUFFILEi},g2term{DE_TAURECIPi},varDataLine(g2term{DE_GAMVARNAMEi}),g2term{DE_GAMVARPOWi}-1);
    % the following value for Psi is only valid if the Gamma2 term is NOT a file function
    Psi_data(numGamma1Terms+pix) = g2term{DE_GAMVARPOWi}*tau_recipVal / sumG1;
end

if unsigned
    Psi_data = abs( Psi_data );
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function derivInfo = GetDerivInfo(focusVarIx,varnames,thisDEix,DEqns,DEpars)
% processes *UNCOMPILED* DEqns only
global DE_NAMEi DE_GAMMA1i DE_GAMMA2i DE_ACTSWi DE_ISTAUFFILEi DE_TAURECIPi DE_GAMVARNAMEi ...
    DE_GAMVARPOWi DE_ISTGTFFILEi DE_TARGETi DE_INTVARNAMEi DE_INTVARPOWi
dEqn = DEqns{thisDEix};
Gamma1 = dEqn{DE_GAMMA1i};
Gamma2 = dEqn{DE_GAMMA2i};
numGamma1Terms = length(Gamma1);
numGamma2Terms = length(Gamma2);
numDerivs = numGamma1Terms + numGamma2Terms;

% only calculate derivs for those inputs to focusVarIx which have their own DEs -- all others are zero
inputListG1 = {};
inputListG2 = {};
for ip=1:numGamma1Terms
    inputListG1 = {inputListG1{:}, Gamma1{ip}{DE_GAMVARNAMEi}};
end

for ip=1:numGamma2Terms
    inputListG2 = {inputListG2{:}, Gamma2{ip}{DE_GAMVARNAMEi}};
end

derivList = cell(1,numDerivs);
for eqnIx=1:length(DEqns)
    name = DEqns{eqnIx}{DE_NAMEi};
    [ p1 ix1 ] = ismember( name, inputListG1 );
    if p1
        [v vix] = ismember(name,varnames);
        derivList{ix1} = {vix, eqnIx};
    else
        [ p2 ix2 ] = ismember( name, inputListG2 );
        if p2
            [v vix] = ismember(name,varnames);
            derivList{ix2} = {vix, eqnIx};
        end
    end
end

derivInfo = {numGamma1Terms, numGamma2Terms, numDerivs, derivList};
return
% first two returned vals are just freebie information from this function
% (observe that only 3rd and 4th entries are used in GetDerivVals() function)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function derivs = GetDerivVals(varDataLine,DEqns,derivInfo)
global DE_NAMEi DE_GAMMA1i DE_GAMMA2i DE_CFACi DE_ACTSWi DE_ISTAUFFILEi DE_TAURECIPi DE_GAMVARNAMEi ...
    DE_GAMVARPOWi DE_ISTGTFFILEi DE_TARGETi DE_INTVARNAMEi DE_INTVARPOWi
% derivList = derivInfo{4};
% d var_i / dt = ( Sum_j{ gj*sj^pj*(Ej-var_i) } + Sum_k{ gk*sk^pk } ) / C
derivs = zeros(1,derivInfo{3});
for derivIx=1:derivInfo{3}
    if ~isempty(derivInfo{4}{derivIx})
        Cfac = DEqns{derivInfo{4}{derivIx}{2}}{DE_CFACi};
        totG1temp = SumGamma1( DEqns{derivInfo{4}{derivIx}{2}}{DE_GAMMA1i}, varDataLine );
        derivs(derivIx) = (totG1temp(2) - totG1temp(1) * varDataLine(derivInfo{4}{derivIx}{1}) + ...
            SumGamma2( DEqns{derivInfo{4}{derivIx}{2}}{DE_GAMMA2i}, varDataLine ) ) / Cfac;
    end % else already 0
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function result = ScaleThreshSort(dataVals, thresh)
% This function is used to scale-threshold sort Psi values, and timescales of V-dependent variables,
%  in NeighbHdEst()
global LARGEBOUND

numVals = length(dataVals);
[sort_vals_rev ix_vals_rev] = sort(dataVals);
sort_vals = zeros(1,numVals);
ix_vals = zeros(1,numVals);
for n = 1:numVals
    sort_vals(n) = sort_vals_rev(numVals-n+1);
    ix_vals(n) = ix_vals_rev(numVals-n+1); % map to tell us what got moved where after sorting!
end

ratios_vals = zeros(numVals-1,1);
for i = numVals-1 : -1 : 1
    if sort_vals(i+1) == 0
        ratios_vals(i) = LARGEBOUND; % should be a big enough representation of infinity!
    else
        ratios_vals(i) = sort_vals(1)/sort_vals(i+1);
    end
end
result = [ix_vals(1)]; % first entry is definitely most dominant!
for i = 1:numVals-1
    if ratios_vals(i) <= thresh % then within scale tolerance of most dominant variable
        result = [result ix_vals(i+1)]; % so include it in the result set
    else
        break % for
    end
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% CompareActs()
% This function needs to check the variables which have File Function tau_reciprocal
%  and target values for the terms in their differential equation. This will not
%  include internal variables that cross-multiply tau_reciprocal values (such as a sodium
%  inactivation gating variable 'h') since they won't appear in the Actives reference set.
% We only check for V-dependent variables (because of short-term independence of instantaneous target on other variables)
%  in perturbed Actives set, comparing to the reference set of Actives (for the unperturbed dynamics)
% Also, need a general rule for whether to check one or more of these together, depending
%  on which ones were active in the previous regime. So far, we're picking out the fastest
%  V-dep variables to focus on the most...
%
function result = CompareActs(actsPert, regimeActs, VdepIxList, VdepFastIxList, qsPars, local_absActs, caOpts)
% NOTE - order of entries in these actives sets is crucial!
% As of May 04 actsRef is relabelled regimeActs -- true actsRef is now local_absActs
result = [false, 0]; % initial and default value (this result indicates *mismatch* between actsPert and actsRef)

% Note: VdepIxList necessarily includes VdepFastIxList
actsP = [];
actsR = [];

% This is too harsh!
% if ~isempty(intersect(actsPert,qsPars)) % no qsPars are allowed to join the actives set
%     result(2) = 1;
%     return % reason 1
% end

% don't count qsPars that are in regimeActs (possible if RegimeDet said they are slow)
qsPars_ref = setdiff(qsPars,regimeActs);

lenApert = length(actsPert);
lenAref  = length(local_absActs); % use local_absActs as the reference set

qsParsPres = intersect(actsPert,qsPars_ref);
for var=qsParsPres
    [pres pos] = ismember(var, actsPert); % pres always true here!
    if ismember(var,VdepFastIxList)
        if pos <= max(lenApert-2,1)
            result(2) = 1;
        end
    else
        if lenApert < 3
            if ~isempty(setdiff(actsPert,local_absActs)) % if pert Acts is this small (2) then just make sure both vars are in set
                result(2) = 1;
                return
            end
        else
            if pos <= max(min(lenApert-1,3),1)
                result(2) = 1; % mismatch
                return
            end
        end
    end
end

%% Form reduced Actives sets, with only V-dependent variables included (in original order)
% Don't do this using `intersect` so that we preserve the order of the entries
% for aix = 1:lenApert
%     if ismember(actsPert(aix),VdepFastIxList)
%         actsP = [actsP actsPert(aix)]; % keep most active at front of list
%     end
% end
% for aix = 1:lenAref
%     if ismember(local_absActs(aix),VdepFastIxList)
%         actsR = [actsR local_absActs(aix)]; % keep most active at front of list
%     end
% end

% In this version, the order of these doesn't matter
actsP = intersect(actsPert,VdepFastIxList);
actsR = intersect(local_absActs,VdepFastIxList); % note, not actsRef

% if a fast var lost from actsPert then mismatch
% if ~isempty(intersect(setdiff(actsR,actsP),VdepFastIxList))
%     result(2) = 2; % mismatch
%     return
% end

% if a fast var joins actsPert (that's not in the regime's Acts) then it's a mismatch
%  if position difference >= attEstParams.posThresh1 (default 2) or length of local acts -1 if smaller
actsR = union(actsR,intersect(regimeActs,VdepFastIxList)); % add fast regime vars not present locally to local set
discrep = setxor(actsP,actsR);
if ~isempty(discrep)
    for var=discrep
        [isp1 pos1] = ismember(var,local_absActs);
        if isp1 == 0 % then isp2 will or will not be true, depending on whether var is in not regimeActs or it is (respectively)
            if ismember(var, regimeActs) % then var is present in regimeActs, but inactive in local epoch, so effectively pos > length(local_absActs)
                pos1 = lenAref + 1;
            else
                pos1 = lenAref + 100; % put it far down!
            end
        end
        [isp2 pos2] = ismember(var,actsPert);
        if isp2 == 0 % then isp1 was true or isp1 was false but var is in regimeActs
            pos2 = lenApert + 100; % put it far down!
        end
        posDiff = abs(pos1-pos2);
        if posDiff >= min(caOpts.pt1,length(local_absActs)-1)
            result(2) = 3; % mismatch
            return
        end
    end
else
    if ~isempty(actsR)
        for var = actsR % fast var
            [isp1 pos1] = ismember(var,local_absActs);
            [isp2 pos2] = ismember(var,actsPert);
            posDiff = abs(pos1-pos2);
            if posDiff >= min(caOpts.pt2,length(local_absActs)-1)
                result(2) = 4; % mismatch
                return
            end
        end
    end
end

result(1) = true; % is a match if got this far!
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = LookupDEpar(parName,DEpars)
result = {0,false}; % default value (error status)
for parIx=1:length(DEpars)
    if strcmp(parName,DEpars{parIx}{1})
        result = {DEpars{parIx}{2},true};
        break
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = GetTauRVal(sourceType,parOrFunc,varVal,power)
if sourceType == 1
    result = eval( [ parOrFunc '(' num2str(varVal) ')' ]);
else
    result = parOrFunc * varVal^power;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DisplayAbsActs(acts, numSpaces, actsStr, varnames)
for sp=1:numSpaces
    fprintf(' ')
end
fprintf([actsStr ' = '])
for ai = 1:length(acts)
    fprintf('%s',varnames{acts(ai)})
    if ai ~= length(acts)
        fprintf(', ')
    end
end
fprintf('\n')
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% useful as a dummy argument for KeyPressFcn callback, to prevent
% echo of key presses to command window
function DoNothing()
return
