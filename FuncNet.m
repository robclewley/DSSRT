function result = FuncNet(networkObjectDir, CFGfilename, thresh, varnames, DEqns, DEpars, ...
    start_time, stop_time, hobj, figHandle, varBounds, varCache, timeCache, CFGinfo, updateNetState)
% Post-processing tool for dimension-reduction analysis of ODEs
% via derivation of a Functional Network. Intended to be called
% from StepNet.m only.
%
% Usage: FuncNet( networkObjectDir, CFGfilename, scaleThreshold,
%          varnames, DEqns, DEpars, startTime, stopTime, netObjects,
%          figHandle, varBounds, varCache, timeCache, CFGinfo,
%          updateNetState )
% returns a cell array containing lots of funky things that are
% too complicated to explain in this help listing. The final
% argument is optional. Read the documentation in the
% Documentation folder for more information.
%
% Version 3.0 November 2005
% (c) 2002-2004 Robert Clewley, Center for Biodynamics, Boston University
% (c) 2005 Robert Clewley, Department of Mathematics, Cornell University

%%%%%%%%%%%%%%%%%%%%%%%% DEVELOPMENT NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 3.0 November 2005
%   Added support for gap junctions in ODE system.
%   Removed useTermsFlag option and made it permanent -- FuncNet now
%    produces data for both Psi influence strength values and term
%    sizes.
%   Cleaned up unused definitions and other minor issues.
%
% Version 2.07 July 2004
%   Added useTermsFlag development (possibly permanent) option to calculate
%    dominance using term sizes, rather than influence strengths.
%
% Version 2.06 May 2004
%   Tidied up syntax and presentation. Renamed `FNS`, `FNO` etc. for new
%    `DSSRT` naming scheme.
%   Added check for 0 entries in absGamIxMapActs and Pots
%
% Version 2.05 February 2004
%   Bug fixes to calculation of Psi values for potentials.
%   Bug fix for storage of info in `domInfo` when num_acts == 1.
%   Some vectorization and other optimization of loops.
%   DEqns argument is now expected to have been compiled already.
%
% Version 2.0 January 2004
%   No longer supports read from file. That data should only be passed via
%    the data cache. All actives and potentials are calculated directly from
%    the variable data.
%
% Version 1.35 January 2004
%   Added 'drawnow' figure refresh commands in main loop
%   Added user-interrupt 'escape' key command
%   varCache is still not being used. Expected to be useful for scale-
%    threshold searches in future versions.
%
% Version 1.3 December 2003
%   Changed order of internal and external variables
%
% Version 1.2 November 2003
%   Added parameter to main function call: networkObjectDir
%   Uses proper MATLAB global variables
%   inputsIx has a range of all system variables, hence offsets to it
%    for external variables' inputs require + numInt
%   References to domGraphIx removed, thus returned value 'result' has
%    different entries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: DEpars argument is currently unused
result = {}; % initial (and default) return value

switch nargin
    case 13
        CFGinfo = {};
        CFGinfoPassed = false;
        updateNetState = true;
    case 14
        if ~isempty(CFGinfo)
            CFGinfoPassed = true;
        else
            CFGinfoPassed = false;
        end
        updateNetState = true;
    case 15
        if ~isempty(CFGinfo)
            CFGinfoPassed = true;
        else
            CFGinfoPassed = false;
        end
	otherwise
        disp('FuncNet:  Internal error! Incorrect number of arguments passed.')
        disp('           Expected 13, 14, or 15.')
        return
end

%%%% Load global constants
%% GLOBAL CONSTANTS (common to all DSSRT system files)
% June 2005

global DE_NAMEi DE_GAMMA1i DE_GAMMA2i DE_CFACi
global DE_ACTSWi DE_ISTAUFFILEi DE_TAURECIPi DE_GAMVARNAMEi
global DE_GAMVARPOWi DE_ISTGTFFILEi DE_TARGETi DE_INTVARNAMEi DE_INTVARPOWi
global OBJ_STATE OBJ_HANDLE OBJ_TYPE OBJ_DEPVAR OBJ_ACTSMAP OBJ_POTSMAP OBJ_STATELT OBJ_LABEL
global OBJ_OCOORD OBJ_LCOORD OBJ_OBJVAR NODE LINK LINKOFFSTATE
global ESEQ_TIMES ESEQ_ACTIVES_BAR ESEQ_ACTIVES ESEQ_POTENTS_BAR ESEQ_POTENTS
global TS_TSEQ TS_NAME TS_VALX TS_VALY TS_ERRS TS_PADS TS_FOCUS TS_PASS TS_ENDP
global TSEQ_ACTS TSEQ_TIME TSEQ_POSN
global LARGEBOUND SMALLBOUND TIMERES ST_ERR NOT_AP POTENT ACTIVE NUMPARS

% differential equation cell array indices
DE_NAMEi   = 1;
DE_GAMMA1i = 2;
DE_GAMMA2i = 3;
DE_CFACi   = 4;

DE_ACTSWi      = 1;
DE_ISTAUFFILEi = 2;
DE_TAURECIPi   = 3;
DE_GAMVARNAMEi = 4;
DE_GAMVARPOWi  = 5;
DE_ISTGTFFILEi = 6;
DE_TARGETi     = 7;
DE_INTVARNAMEi = 8;
DE_INTVARPOWi  = 9;

% object indices
OBJ_STATE   = 1;
OBJ_HANDLE  = 2;
OBJ_TYPE    = 3;
OBJ_DEPVAR  = 4;
OBJ_ACTSMAP = 5;
OBJ_POTSMAP = 6;
OBJ_STATELT = 7;
OBJ_LABEL   = 8;
OBJ_OCOORD  = 9;
OBJ_LCOORD  =10;
OBJ_OBJVAR  =11; % for link objects only

LINKOFFSTATE = 100; % arbitrary value!

% object types
NODE = 1;
LINK = 2;

% object state maps
ST_ERR  =-1;
NOT_AP  = 0;
POTENT  = 1;
ACTIVE  = 2;

% curly E set sequence indices (actives and potentials sets)
ESEQ_TIMES       = 1;
ESEQ_ACTIVES_BAR = 2;
ESEQ_ACTIVES     = 3;
ESEQ_POTENTS_BAR = 4;
ESEQ_POTENTS     = 5;

% transition sequence indices
TS_TSEQ  = 1;
TS_NAME  = 2;
TS_VALX  = 3;
TS_VALY  = 4;
TS_ERRS  = 5;
TS_PADS  = 6;
TS_FOCUS = 7;
TS_PASS  = 8;
TS_ENDP  = 9;

% inside a transition sequence
TSEQ_ACTS = 1;
TSEQ_TIME = 2;
TSEQ_POSN = 3;

% number of parameters in parameter .par files
NUMPARS = 17;

% DSSRT 'infinity' (e.g. used in Psi ratios)
LARGEBOUND = 1e6; % e.g. ratio to use when divisor for scale ratios is zero, or effective upper limit on variable size

% time resolution
TIMERES = 1e-6;

% maximum precision for DSSRT
SMALLBOUND = 1e-4;

%%%% END of common global constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%% Essential bounds checking
if stop_time <= start_time
    disp('FuncNet:  Internal error. Stop time argument passed was smaller than start time.')
    return
end

if thresh <= 1
    disp('FuncNet:  Internal error. Argument scale threshold must be > 1')
    return
end

%%%% CFG setup
if ~CFGinfoPassed
	setup = ReadConfigFile(CFGfilename,networkObjectDir);
	
	if ~isempty(setup)
        fprintf('FuncNet:  CFG file `%s` read successfully\n',CFGfilename)
	else
        fprintf('FuncNet:  Error reading CFG file `%s`\n',CFGfilename)
        return
	end
else
    setup = CFGinfo;
end

numInt = setup{1}; % total number of internal variables in the system
numExt = setup{2}; % total number of external variables in the system
numTot = numInt + numExt; % total number of variables in the system

% intVarNames = setup{3};   % unused
extVarNames = setup{4};
allVarNames = setup{5}; % the concatentation [ extVarNames, intVarNames ] except in terms of cell arrays

inputsIx         = setup{6};   % input indices (absolute) per system variable
candActsIxMapAbs = setup{7};   % map from relative indices of observable's actives -> absolute ix
candPotsIxMapAbs = setup{8};   % map from relative indices of observable's potentials -> absolute ix

ixExtSetAbs = [1:numExt];        % set of absolute indices for internal variables
ixIntSetAbs = [numExt+1:numTot]; % set of absolute indices for external variables


% get number of internal, external variables, and valid potentials connected to each of the
% *external* variables in allVarNames (-> network connectivity)
internals_per_ev  = zeros(1,numExt); % initialisations
externals_per_ev  = zeros(1,numExt);
actives_per_ev    = zeros(1,numExt);
potentials_per_ev = zeros(1,numExt);
for ev = 1:numExt
    for i=1:length(inputsIx{ev})
        if ismember(inputsIx{ev}(i),ixIntSetAbs)
			internals_per_ev(ev)  = internals_per_ev(ev) + 1;
        elseif ismember(inputsIx{ev}(i),ixExtSetAbs)
			externals_per_ev(ev)  = externals_per_ev(ev) + 1;
        end
    end
    actives_per_ev(ev)    = length(candActsIxMapAbs{ev});
    potentials_per_ev(ev) = length(candPotsIxMapAbs{ev});
end
% total_internals = sum(internals_per_ev);
% total_externals = sum(externals_per_ev);
% total_actives = sum(actives_per_ev);
% total_potents = sum(potentials_per_ev);

% map of candidate potentials' indices to
% relative positions in 'pot'
candPotsIxMapRel = cell(1,numExt);
for ev = 1:numExt % for every external var
    for i=1:length(candPotsIxMapAbs{ev})
        [isIn ix] = ismember(candPotsIxMapAbs{ev}(i), candActsIxMapAbs{ev});
        if isIn
            candPotsIxMapRel{ev} = [candPotsIxMapRel{ev} ix];
        end
    end
end

% for each external variable, we need a
% map from candidate active / pot. *internal* vars -> index for additional node state maps (if applicable)
statemapIxMaps = cell(1,numExt);
for ev = 1:numExt
	statemapIxMapA = zeros(1,numInt);
	statemapIxMapP = zeros(1,numInt);
    nextIxPosA = 1; % which state map in entry of hobj (one for every internal candidate active / pot. var)
    nextIxPosP = 1;
	for intIx = numExt+1:numTot % search internals
        if ismember(intIx, candActsIxMapAbs{ev}) % then there'll be a additional state map
            statemapIxMapA(intIx-numExt) = nextIxPosA;
            nextIxPosA = nextIxPosA + 1;
        end
        if ismember(intIx, candPotsIxMapAbs{ev}) % then there'll be a additional state map
            statemapIxMapP(intIx-numExt) = nextIxPosP;
            nextIxPosP = nextIxPosP + 1;
        end
    end
    statemapIxMaps{ev} = {statemapIxMapA, statemapIxMapP};
end

% initialize size of passed network objects
numObjects = length(hobj);

% Build map from source variable position (1:numTot) to the ext_var
% destination (1:numExt) -- map called `whosevars`
% e.g. whosevars_inv = {[1 3 5 10 11], [2 4 6 9 12], [7], [8], [], []}
% e.g. whosevars = [1 2 1 2 1 2 3 4 2 1 1 2] ... zero indicates no-one
% We won't care about internal variables in inputsIx, so we have to miss out
% those first entries in order to be left with whosevars_inv
whosevars_inv = cell(1,numExt);
for ev=1:numExt
    whosevars_inv{ev} = inputsIx{ev}; % inverse of whosevars map, for external variables only
end

if isempty(setup{9}) % then setup came from ReadConfigFile(), which didn't return a ninth entry for whosevars
	whosevars = zeros(1,numTot); % corresponds to entries in varnames
	for sv = 1:numTot
%         offset=0;
		for ev = 1:numExt
            if intersect(sv,whosevars_inv{ev})
                whosevars(sv)=ev;
            end
		end
	end
else
    whosevars = setup{9};
end

% find initial position in data given by the start_time argument if using the data file
if start_time >= 0
    if start_time < timeCache(1) % then specified start time was smaller than that included in data file
        disp('FuncNet:  Specified start time was smaller than time of first entry')
        disp('            in data file. Using first entry instead.')
        start_time = timeCache(1);
        start_pos = 1;
        start_time_found = true;  % currently not used, but consistent with else statement
    else
        [start_time_found start_pos] = ismember(1,timeCache <= start_time);
        if start_time_found
            start_time = timeCache(start_pos);
        else
            disp('FuncNet:  Start time too large: not found in variable data')
            return
        end
    end
else
    disp('FuncNet:  Start time must be a non-negative real number!')
    return
end

[stop_time_found stop_pos] = ismember(1,timeCache <= stop_time);
if ~stop_time_found
    disp('FuncNet:  Stop time too large: not found in variable data')
    return
else
    stop_time = timeCache(stop_pos); % just to be absolutely accurate about the time;
end

fprintf('FuncNet:  Processing variable data from times %.2f -> %.2f\n',start_time,stop_time)
disp(   '           This may take several minutes (press <escape> to cancel) ...')
% varCacheLen = length(varCache); % number of steps present in the cache


%%%% Work out actives and potentials sets
numEqns   = length(DEqns);
deqnIxMap = zeros(1,numExt); % zero entries will denote no corresponding differential equation present
for varIx = 1:numExt
	for eqnix=1:numEqns % create map entries
        if strcmp(DEqns{eqnix}{DE_NAMEi},varnames{varIx})
            deqnIxMap(varIx) = eqnix;
            break % inner for loop: done!
        end
	end
end

% initialize curly E sets -- THESE HOLD ABSOLUTE INDICES i.e. 'i' in inputsIx{numInt+i}
Eseq          = cell(1,stop_pos-start_pos+1); % holds the temporal sequence of the curly E sets in a cell array
% all_E_b_old   = cell(1,stop_pos-start_pos+1); % holder for previous curly E_b sets
% all_E_old     = cell(1,stop_pos-start_pos+1); % etc.
% all_E_pot_old = cell(1,stop_pos-start_pos+1);
Eseq_term          = cell(1,stop_pos-start_pos+1); % holds the temporal sequence of the curly E sets in a cell array
% all_E_b_old_term   = cell(1,stop_pos-start_pos+1); % holder for previous curly E_b sets
% all_E_old_term     = cell(1,stop_pos-start_pos+1); % etc.
% all_E_pot_old_term = cell(1,stop_pos-start_pos+1);

% initialize cell array to hold all diagram states, and dominance graphs, ratios data
diag_states = cell(1,stop_pos-start_pos+1);
domGraph    = zeros(numExt,stop_pos-start_pos+1,3);
domInfo     = cell(numExt,stop_pos-start_pos+1);
diag_states_term = cell(1,stop_pos-start_pos+1);
domGraph_term    = zeros(numExt,stop_pos-start_pos+1,3);
domInfo_term     = cell(numExt,stop_pos-start_pos+1);

%% Create map from [Gamma1, Gamma2] order of input variables to the set of relative indices of
% candidate actives / potentials to each ext_var
relGamIxMapActs = cell(1,numExt);
relGamIxMapPots = cell(1,numExt);
absGamIxMap = cell(1,numExt); % build this first
for evIx=1:numExt
    absGamIxMap_temp = [];
    thisDEix = deqnIxMap(evIx);
    if thisDEix > 0
        thisDE = DEqns{thisDEix};
        for g1t=1:length(thisDE{DE_GAMMA1i})
            g1 = thisDE{DE_GAMMA1i}{g1t};
            if g1{DE_ISTGTFFILEi} == 2
                absGamIxMap_temp = [ absGamIxMap_temp g1{DE_TARGETi} ];
            else
                absGamIxMap_temp = [ absGamIxMap_temp g1{DE_GAMVARNAMEi} ];
            end
        end
        for g2t=1:length(thisDE{DE_GAMMA2i})
            absGamIxMap_temp = [ absGamIxMap_temp thisDE{DE_GAMMA2i}{g2t}{DE_GAMVARNAMEi} ];
        end
        absGamIxMap{evIx} = absGamIxMap_temp;
    end
end

% abs_ix = candActsIxMapAbs{ev}( rel_ix ) etc.
for ev=1:numExt
    thisAbsGamMap = absGamIxMap{ev};
    for i = 1:length(thisAbsGamMap)
        [ispres aix] = ismember( thisAbsGamMap(i), candActsIxMapAbs{ev} );
        if ispres
            relGamIxMapActs{ev}(aix) = i;
        end
        [ispres pix] = ismember( thisAbsGamMap(i), candPotsIxMapAbs{ev} );
        if ispres
            relGamIxMapPots{ev}(pix) = i;
        end
    end
    if ~isempty(find(relGamIxMapActs{ev}==0, 1))
        fprintf('FuncNet:  Error in differential equation setup for external variable %s\n          -- terms in `relGamIxMapActs` are zero.\n',varnames{ev})
        return
    end
    if ~isempty(find(relGamIxMapPots{ev}==0, 1))
        fprintf('FuncNet:  Error in differential equation setup for external variable %s\n          -- terms in `relGamIxMapPots` are zero.\n',varnames{ev})
        return
    end
end
% Note that relGamIxMapActs and relGamIxMapPots may be identical if candidate actives and potentials are always the same.
% Having two separate maps is really just for backwards compatibility right now.

interruptCounter = 0;
drawnow

for tpos = start_pos:stop_pos

    interruptCounter = interruptCounter + 1; % check interrupt status every 100 steps
    if interruptCounter == 100
        interruptCounter = 0;
        % check if user is trying to interrupt
        drawnow
        if get(figHandle,'CurrentCharacter') == 27 % escape key
            beep
            disp('FuncNet:  Functional network computation cancelled by user')
            result = {-1}; % cancel code
            return
        end
    end

    var_row = varCache(tpos,:);

    % initialise this time step's E sets
    all_E_b     = cell(1,numExt); % actives, including internals
    all_E       = cell(1,numExt); %      ... excluding internals
    all_E_pot_b = cell(1,numExt); % potentials, including internals
    all_E_pot   = cell(1,numExt); %      ... excluding internals
    all_E_b_term     = cell(1,numExt); % actives, including internals
    all_E_term       = cell(1,numExt); %      ... excluding internals
    all_E_pot_b_term = cell(1,numExt); % potentials, including internals
    all_E_pot_term   = cell(1,numExt); %      ... excluding internals


    for ext_var = 1:numExt
        num_ints = internals_per_ev(ext_var);
%         num_exts = externals_per_ev(ext_var);  % unused
%         num_ipts = num_ints + num_exts; % unused
        num_acts = actives_per_ev(ext_var);
        num_pots = potentials_per_ev(ext_var);

        %%%% BEGIN { DETERMINE ACTIVES }
        if num_acts > 0
            if num_acts == 1
                allIps = inputsIx{ext_var}; % WAS extIps = inputsIx{ext_var}( inputsIx{ext_var} <= numExt );
                if length(allIps) > 1
                    disp('FuncNet:  Internal error: extIps not of length 1!')
                    return
                end
                [isAct ixa] = ismember( allIps, candActsIxMapAbs{ext_var} ); % ixa is a relative index
                % don't need act, sort_act etc., but in case they are
                % later, they are defined here for consistency
                act = 1; % define this to be 1
                sort_act = 1;
                ixa_term = ixa;
                act_term = 1;
                sort_act_term = 1;
            else
                % is there a differential equation for this ext_var?
                if deqnIxMap(ext_var) > 0
                    [act, act_term] = GetDominanceVals(var_row,DEqns{deqnIxMap(ext_var)},ext_var,true); % unsigned
                    % act is not yet in the order that this code expects: currently in Gamma1, Gamma2 order ...
                    % so we map it to the correct range
                    act = act( relGamIxMapActs{ext_var} ); % will use act_gam in potentials calculation
                    act_term = act_term( relGamIxMapActs{ext_var} );
                else
                    act = ones(num_acts); % dummy values
                    act_term = ones(num_acts);
                end
                [sort_act_rev, ixa_rev] = sort(act);
                [sort_act_term_rev, ixa_term_rev] = sort(act_term);
                % reverse order to get descending order
                sort_act = sort_act_rev(num_acts:-1:1);
                sort_act_term = sort_act_term_rev(num_acts:-1:1);
                ixa = ixa_rev(num_acts:-1:1); % map to tell us what got moved where after sorting
                ixa_term = ixa_term_rev(num_acts:-1:1);
            end
        
            % initialize curly E sets for ext_var
    		E_b = [ixa(1)]; % first entry is undoubtedly most dominant
			if candActsIxMapAbs{ext_var}(ixa(1)) <= numExt % then is an external var
                E = [ixa(1)];
			else
                E = [];
            end
            E_b_term = [ixa_term(1)];
            if candActsIxMapAbs{ext_var}(ixa_term(1)) <= numExt % then is an external var
                E_term = [ixa_term(1)];
			else
                E_term = [];
            end
        else
            E_b = [];
            E   = [];
            E_b_term = [];
            E_term = [];
        end
        E_pot_b = [];
        E_pot   = [];
        E_pot_b_term = [];
        E_pot_term = [];

        if num_acts > 1
            % build ordered set of ratios for actives (all vs. rank #1)
            % (i.e. convert ordered influence strengths into ratio form)
			ratios_act = zeros(num_acts-1,1);
            ratios_act_term = zeros(num_acts-1,1);
			for i = num_acts-1 : -1 : 1
                if sort_act(i+1) == 0
                    ratios_act(i) = LARGEBOUND; % should be a big enough representation of infinity!
                else
                    ratios_act(i) = sort_act(1)/sort_act(i+1);
                end
                if sort_act_term(i+1) == 0
                    ratios_act_term(i) = LARGEBOUND; % should be a big enough representation of infinity!
                else
                    ratios_act_term(i) = sort_act_term(1)/sort_act_term(i+1);
                end
            end

            % add element to domGraph
            domGraph(ext_var,tpos-start_pos+1,:) = [timeCache(tpos), ratios_act(1), candActsIxMapAbs{ext_var}(ixa(1))];
            domGraph_term(ext_var,tpos-start_pos+1,:) = [timeCache(tpos), ratios_act_term(1), candActsIxMapAbs{ext_var}(ixa_term(1))];

            for i = 1:num_acts-1
                if ratios_act(i) <= thresh % then within scale tolerance of most dominant variable
                    E_b = [E_b ixa(i+1)]; % so include it in the \bar{curly E} set
                else
                    break % for
                end
            end % for

            for i = 1:num_acts-1
                if ratios_act_term(i) <= thresh % then within scale tolerance of most dominant variable
                    E_b_term = [E_b_term ixa_term(i+1)]; % so include it in the \bar{curly E} set
                else
                    break % for
                end
            end % for

            sort_act_dom = sort_act(1); % use this for domInfo record later (after pots re-use `sort_act` variable)
            sort_act_term_dom = sort_act_term(1);
            ixa_dom = ixa; % use this for domInfo record later (after pots re-use `ixa` variable)
            ixa_term_dom = ixa_term;

            E = E_b( candActsIxMapAbs{ext_var}(E_b) <= numExt ); % pick out only the external variables
            E_term = E_b_term( candActsIxMapAbs{ext_var}(E_b_term) <= numExt );
            
        elseif num_acts == 1
            ratios_act = [];
            sort_act_dom = sort_act;
            ixa_dom = ixa;
            if candActsIxMapAbs{ext_var}(ixa(1)) <= numExt % then is an external variable
                E = [ixa(1)];
            end
            E_b = [ixa(1)];
            ratios_act_term = [];
            sort_act_term_dom = sort_act_term;
            ixa_term_dom = ixa_term;
            if candActsIxMapAbs{ext_var}(ixa_term(1)) <= numExt % then is an external variable
                E_term = [ixa_term(1)];
            end
            E_b_term = [ixa_term(1)];
        elseif num_acts == 0
            ratios_act = [];
            sort_act_dom = [];
            ixa_dom = [];
            E = [];
            E_b = [];
            ratios_act_term = [];
            sort_act_term_dom = [];
            ixa_term_dom = [];
            E_term = [];
            E_b_term = [];
        end % if num_acts > 1
        %%%% END { DETERMINE ACTIVES }

        %%%% BEGIN { DETERMINE POTENTIALS }
        if (num_pots > 1 || (num_pots == 1 && num_ints > 0)) && deqnIxMap(ext_var) > 0
            % either > 1 potentials, or there are internal variables involved if there's only one potential
            % compare potential interactions at time t against scale of dominant
            % set E_b (i.e. include internal variables)
            cand_act_ix = candPotsIxMapRel{ext_var}([1:num_pots]); % candPotsIxMapRel{ext_var} is not *always* the identity map

            % to compare like for like, both cand_act_ix and E_b must be in terms of absolute indices
            E_bAbs = candActsIxMapAbs{ext_var}(E_b);
            E_b_termAbs = candActsIxMapAbs{ext_var}(E_b_term);

            % find set of indices for potentials to sort...
            % ...don't include any whose active values are already in the dominant set \bar{curly E}
            pot_ix = cand_act_ix(ismember( candActsIxMapAbs{ext_var}(cand_act_ix), E_bAbs ) == 0); % relative indices in candActsIxMap
            pot_ix_term = cand_act_ix(ismember( candActsIxMapAbs{ext_var}(cand_act_ix), E_b_termAbs ) == 0);
            % Note that, even though we use candActsIxMapAbs and not the Pots, the cand_act_ix set ensures we won't pick out
            %  any invalid potentials from the set of actives. Since cand_act_ix is defined to be indices relative to the set
            %  of possible actives, then if we use the Pots instead we will get index out of bounds.

            thisDE = DEqns{deqnIxMap(ext_var)};
            pot_varVals = GetPsiPotVals(var_row, thisDE, varBounds, candActsIxMapAbs{ext_var}(pot_ix));
            pot_term_varVals = GetTermPotVals(var_row, thisDE, varBounds, ext_var, candActsIxMapAbs{ext_var}(pot_ix_term));

            lpi = length(pot_ix);
            test_act_dom = cell(1, lpi); % store these for domInfo variable, later
            for n = 1 : lpi % find out what actives set would look like with each pot in there individually
                curr_potIx = pot_ix(n);
                % find rank #1 active influence strength with this potential's max Psi val and test this potential against it.
                test_vars = var_row;
                test_vars( candActsIxMapAbs{ext_var}(curr_potIx) ) = pot_varVals( relGamIxMapActs{ext_var}(curr_potIx) );
                [test_act, ignored] = GetDominanceVals(test_vars,DEqns{deqnIxMap(ext_var)},ext_var,true); % unsigned
                % test_act is not yet in the order that this code expects: currently in Gamma1, Gamma2 order ...
                % so we map it to the correct range
                test_act = test_act( relGamIxMapActs{ext_var} );
                test_act_dom{n} = test_act;
                % find actives for test set
                [sort_act_rev, ixa_rev] = sort(test_act);
                if ixa_rev(num_acts) == curr_potIx % then this potential is most dominant when activated, so we're done
                    E_pot_b = [E_pot_b curr_potIx];
                    if candActsIxMapAbs{ext_var}(curr_potIx) <= numExt % then is an external var
                        % here we use ActsIxMap and not Pots because pot_ix is defined relative to set of possible actives, even
                        % though it will only pick out valid potentials as specified by user
                        E_pot = [E_pot curr_potIx];
                    end
                else
                    % else continue to check ratios...
                    % reverse order to get descending order
                    sort_act = sort_act_rev(num_acts:-1:1);
                    ixa = ixa_rev(num_acts:-1:1); % map to tell us what got moved where after sorting
                    [ ispres ix_pot ] = ismember(curr_potIx, ixa); % ispres guaranteed true
                    if sort_act(1)/sort_act(ix_pot) <= thresh % then within scale tolerance set vs. active rank #1
                        E_pot_b = [E_pot_b curr_potIx];
                        if candActsIxMapAbs{ext_var}(curr_potIx) <= numExt % then is an external var
                            % here we use ActsIxMap and not Pots because pot_ix is defined relative to set of possible actives, even
                            % though it will only pick out valid potentials as specified by user
                            E_pot = [E_pot curr_potIx];
                        end
                    end
                end
            end % for n (curr_potIx)
            lpi_t = length(pot_ix_term);
            test_act_term_dom = cell(1, lpi_t); % store these for domInfo variable, later
            for n = 1 : lpi_t % find out what actives set would look like with each pot in there individually
                curr_potIx_term = pot_ix_term(n);
                % find rank #1 active influence strength with this
                % potential's max Psi val and test this potential against it.
                test_vars_term = var_row;
                test_vars_term( candActsIxMapAbs{ext_var}(curr_potIx_term) ) = pot_term_varVals( relGamIxMapActs{ext_var}(curr_potIx_term) );
                [ignored, test_act_term] = GetDominanceVals(test_vars_term,DEqns{deqnIxMap(ext_var)},ext_var,true); % unsigned
                % test_act is not yet in the order that this code expects: currently in Gamma1, Gamma2 order ...
                % so we map it to the correct range
                test_act_term = test_act_term( relGamIxMapActs{ext_var} );
                test_act_term_dom{n} = test_act_term;
                % find actives for test set
                [sort_act_term_rev, ixa_term_rev] = sort(test_act_term);
                if ixa_term_rev(num_acts) == curr_potIx_term % then this potential is most dominant when activated, so we're done
                    E_pot_b_term = [E_pot_b_term curr_potIx_term];
                    if candActsIxMapAbs{ext_var}(curr_potIx_term) <= numExt % then is an external var
                        % here we use ActsIxMap and not Pots because pot_ix is defined relative to set of possible actives, even
                        % though it will only pick out valid potentials as specified by user
                        E_pot_term = [E_pot_term curr_potIx_term];
                    end
                else
                    % else continue to check ratios...
                    % reverse order to get descending order
                    sort_act_term = sort_act_term_rev(num_acts:-1:1);
                    ixa_term = ixa_term_rev(num_acts:-1:1); % map to tell us what got moved where after sorting
                    [ ispres ix_pot ] = ismember(curr_potIx, ixa_term); % ispres guaranteed true
                    if sort_act_term(1)/sort_act_term(ix_pot) <= thresh % then within scale tolerance set vs. active rank #1
                        E_pot_b_term = [E_pot_b_term curr_potIx_term];
                        if candActsIxMapAbs{ext_var}(curr_potIx_term) <= numExt % then is an external var
                            % here we use ActsIxMap and not Pots because pot_ix is defined relative to set of possible actives, even
                            % though it will only pick out valid potentials as specified by user
                            E_pot_term = [E_pot_term curr_potIx_term];
                        end
                    end
                end
            end % for n (curr_potIx_term)
        else % else no potentials need considering, since either
            % num_pots == 0, or ==1 but only 1 variable acts on ext_var's
            % eqn, and it's an external variable!
            % do nothing, since E_pot == [] already
            test_act_dom = {};
            pot_ix = [];
            test_act_term_dom = {};
            pot_ix_term = [];
        end % if num_pots
        %%%% END { DETERMINE POTENTIALS }

        % update current set of E sets over all ext vars
		all_E_b{ext_var}     = E_b;
		all_E{ext_var}       = E;
		all_E_pot{ext_var}   = E_pot;
		all_E_pot_b{ext_var} = E_pot_b;
        domInfo{ext_var,tpos-start_pos+1} = { E_b, ratios_act, ixa_dom, sort_act_dom, test_act_dom, pot_ix };
        all_E_b_term{ext_var}     = E_b_term;
		all_E_term{ext_var}       = E_term;
		all_E_pot_term{ext_var}   = E_pot_term;
		all_E_pot_b_term{ext_var} = E_pot_b_term;
        domInfo_term{ext_var,tpos-start_pos+1} = { E_b_term, ratios_act_term, ixa_term_dom, sort_act_term_dom, test_act_term_dom, pot_ix_term };
    end % for ext_var

    %%%% Update sequence of curly sets
    Eseq{tpos-start_pos+1} = {timeCache(tpos), all_E_b, all_E, all_E_pot_b, all_E_pot};
    Eseq_term{tpos-start_pos+1} = {timeCache(tpos), all_E_b_term, all_E_term, all_E_pot_b_term, all_E_pot_term};

    %%%% BEGIN { DETERMINE CURRENT NETWORK STATE }
    if updateNetState
        % E_b is in order of observables in list extVarNames
        % ix = candActsIxMapAbs{ext_var}(tempEb(v)) is an absolute index of
        %   a variable 'v' that's active for ext_var
		% To update all objects, only have to go through the observables (external vars), and follow up
		%   their dependent states, updating the objects associated with those on the way.

        % currentActs = all_E_b = Eseq{tpos}{2}
        % currentPots = all_E_pot_b = Eseq{tpos}{4}
        
        % Update sequence of diagram states
        diag_states{tpos-start_pos+1} = GetDiagState(numExt, allVarNames, extVarNames, all_E_b, all_E_pot_b, ...
                    statemapIxMaps, ixIntSetAbs, candActsIxMapAbs, hobj, numObjects);
        diag_states_term{tpos-start_pos+1} = GetDiagState(numExt, allVarNames, extVarNames, all_E_b_term, all_E_pot_b_term, ...
                    statemapIxMaps, ixIntSetAbs, candActsIxMapAbs, hobj, numObjects);

    end % if updateNetState
    %%%% END { DETERMINE CURRENT FUNCTIONAL NETWORK STATE }
    
end % while

tstep = timeCache(2)-timeCache(1); % arbitrary to check 2nd and 1st positions, since numerical time-step assumed fixed
fprintf('FuncNet:  Completed generation of functional network\n')
fprintf('            at a time step of %.4f\n', tstep);

result = { diag_states, diag_states_term, Eseq, Eseq_term, numTot, whosevars, candActsIxMapAbs, candPotsIxMapAbs, ...
        allVarNames, stop_pos-start_pos+1, tstep, domGraph, domGraph_term, domInfo, domInfo_term, timeCache(start_pos:stop_pos), ...
        varCache(start_pos:stop_pos,:) };
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                Secondary functions only used by FuncNet
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function diag_state = GetDiagState(numExt, allVarNames, extVarNames, currentActs, currentPots, ...
            statemapIxMaps, ixIntSetAbs, candActsIxMapAbs, hobj, numObjects)
global OBJ_ACTSMAP OBJ_POTSMAP OBJ_STATE LINKOFFSTATE

diag_state = zeros(1, numObjects);
for ext_var = 1:numExt
    nodeName  = extVarNames{ext_var}; % each ext_var should have a node associated with it
    currEb    = currentActs{ext_var};
    currEp    = currentPots{ext_var};
    lenCurrEb = length(currEb);
    lenCurrEp = length(currEp);
    if lenCurrEb > 0 % DO ACTIVES FIRST
        nodeIx = FindObjIx(hobj, nodeName, 'Node'); % this is the index for this node's state in diag_state
           % ... as well as being the object handle in hobj
        if nodeIx == 0
            disp('FuncNet:  Node name not found among network objects. Probable error in configuration')
            return
        else
            for var = 1:lenCurrEb
                stateInt = diag_state(nodeIx); % may be updated each time round 'var' loop
                varIx = candActsIxMapAbs{ext_var}(currEb(var)); % absolute index of input variable var in allVarNames
                if ismember( varIx, ixIntSetAbs ) % then updating node (internal) state
                    % map old state to new state for node object
                    smIxMapA = statemapIxMaps{ext_var}{1}; % get index from varIx -> which map
                    mapAct = hobj{nodeIx}{OBJ_ACTSMAP}{smIxMapA(varIx-numExt)};
                    stateInt  = DoStateMap(stateInt, mapAct);
                    if stateInt >= 0
                        diag_state(nodeIx) = stateInt;
                    else
                        if stateInt == -2
                            fprintf('FuncNet:  Error updating diagram state at time position %i for observable %s',tpos,nodeName)
                            disp(   '           -- active`s state not in node`s internal StateMap')
                            return
                        elseif stateInt == -1
                            disp('FuncNet:  Network state mapping error during determination of actives for node')
                            return
                        else
                            disp('FuncNet:  Error in value returned by state mapping of actives for node')
                            return
                        end
                    end
                else % is in ixExtSetAbs, and updating a link state
                    % map old state to new state for link object
                    varName  = allVarNames{varIx};
                    varObjIx = FindObjIx(hobj, varName, 'Link', nodeName);
                    if varObjIx == 0
                        fprintf('FuncNet:  Failed to find variable name %s in object handles \n',varName)
                        return
                    end
                    if hobj{varObjIx}{OBJ_STATE} ~= LINKOFFSTATE % then update link state
                        stateExt = diag_state(varObjIx);
                        mapAct = hobj{varObjIx}{OBJ_ACTSMAP}{1}; % links only have one map
                        stateExt  = DoStateMap(stateExt, mapAct);
                        if stateExt >= 0
                            diag_state(varObjIx) = stateExt;
                        else
                            if stateExt == -2
                                fprintf('FuncNet:  Error updating diagram state at time position %i for observable %s\n',tpos,nodeName)
                                disp(   '           -- active`s state not in external link`s StateMap')
                                return
                            elseif stateExt == -1
                                disp('FuncNet:  Network state mapping error during determination of actives for link')
                                return
                            else
                                disp('FuncNet:  Error in value returned by state mapping of actives for link')
                                return
                            end
                        end
                    end % else leave the un-displayed state alone
                end
            end
        end
    end % if lenCurrEb > 0

    if lenCurrEp > 0 % DO POTENTIALS NEXT
        nodeIx = FindObjIx(hobj, nodeName, 'Node'); % this is the index for this node's state in diag_state
           % ... as well as being the object handle in hobj
        if nodeIx == 0
            disp('FuncNet:  Node name not found amongst network objects. Probable error in configuration')
            return
        else
            for var = 1:lenCurrEp
                stateInt = diag_state(nodeIx); % may be updated each time round 'var' loop
                varIx = candActsIxMapAbs{ext_var}(currEp(var)); % absolute index of input variable var in allVarNames
                if ismember( varIx, ixIntSetAbs ) % then updating node (internal) state
                    % map old state to new state for node object
                    smIxMapP = statemapIxMaps{ext_var}{2}; % get index from varIx -> which map
                    mapPot = hobj{nodeIx}{OBJ_POTSMAP}{smIxMapP(varIx-numExt)};
                    stateInt  = DoStateMap(stateInt, mapPot);
                    if stateInt >= 0
                        diag_state(nodeIx) = stateInt;
                    else
                        if stateInt == -2
                            fprintf('FuncNet:  Error updating diagram state at time position %i for observable %s',tpos,nodeName)
                            disp(   '           -- potential`s state not in node`s internal StateMap')
                            return
                        elseif stateInt == -1
                            disp('FuncNet:  Network state mapping error during determination of potentials for node')
                            return
                        else % ??
                            disp('FuncNet:  Error in value returned by state mapping of potentials for node')
                            return
                        end
                    end
                else % is in ixExtSetAbs, and updating a link state
                    % map old state to new state for link object
                    varName  = allVarNames{varIx};
                    varObjIx = FindObjIx(hobj, varName, 'Link', nodeName);
                    if varObjIx == 0
                        fprintf('FuncNet:  Failed to find variable name %s in object handles \n',varName)
                        return
                    end
                    if hobj{varObjIx}{OBJ_STATE} ~= LINKOFFSTATE % then update link state
                        stateExt = diag_state(varObjIx);
                        mapPot = hobj{varObjIx}{OBJ_POTSMAP}{1}; % links only have one map
                        stateExt  = DoStateMap(stateExt, mapPot);
                        if stateExt >= 0
                            diag_state(varObjIx) = stateExt;
                        else
                            if stateExt == -2
                                fprintf('FuncNet:  Error updating diagram state at time position %i for observable %s\n',tpos,nodeName)
                                disp(   '           -- potential`s state not in external link`s StateMap')
                                return
                            elseif stateExt == -1
                                disp('FuncNet:  Network state mapping error during determination of potentials for link')
                                return
                            else
                                disp('FuncNet:  Error in value returned by state mapping of potentials for link')
                                return
                            end
                        end
                    end % else leave the un-displayed state alone
                end
            end
        end
    end % if lenCurrEp > 0     
end % for ext_var


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = DoStateMap(state, stateMap)
result = -2; % initial and default value (signifies state not found in map)
lenMap = length(stateMap);

for entry = 1:lenMap
    if state == stateMap(entry,1)
        result = stateMap(entry,2);
    end
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function result = FindObjIx(objects, name, typeCf, objVarName)
% type is either 1 (node), or 2 (link)
% typeStr and typeCf are strings, either "Link" or "Node"
% name and objVarName is a string, referring to a dependent variable name
global OBJ_STATE OBJ_HANDLE OBJ_TYPE OBJ_DEPVAR OBJ_ACTSMAP OBJ_POTSMAP OBJ_STATELT
global OBJ_LABEL OBJ_OCOORD OBJ_LCOORD OBJ_OBJVAR LINKOFFSTATE

result = 0; % initial (and default) value
len = length(objects);

if nargin == 3
    objVarName = '';
end

for obj = 1:len
    check1 = false;
	check2 = false;
    depVarStr = objects{obj}{OBJ_DEPVAR};
    type      = objects{obj}{OBJ_TYPE};
    if type == 1
        typeStr = 'Node';
    else
        typeStr = 'Link';
        objVarCf = objects{obj}{OBJ_OBJVAR};
    end
    lenTypSt = length(typeStr);
    lenDepVS = length(depVarStr);
    if length(name) == lenDepVS
        if sum( name ~= depVarStr ) == 0 % then strings equal
            if type == 2 % then have to also check objVarName for links before choosing this entry
                if strcmp( objVarName, objVarCf ) % then strings equal
                    check1 = true;
                end
            else
                check1 = true;
            end
        end
    end
    if lenTypSt == 4 % then has correct length for a type string
        if sum( typeCf ~= typeStr ) == 0 % then strings equal
            check2 = true;
        end
    else
        disp('  FindObjIx: Passed invalid type string in arguments!')
    end
    if check1 && check2
        result = obj;
    end
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% A copy of this function is also present inside AttEst_*.m and StepNet.m
%
function [Psi_data, Term_data] = GetDominanceVals(varDataLine, thisDE, ext_var, unsigned)
% calculates Psi influence strength values for each term (Gamma1 and Gamma2) of dEqn differential equation
global DE_NAMEi DE_GAMMA1i DE_GAMMA2i DE_ACTSWi DE_ISTAUFFILEi DE_TAURECIPi DE_GAMVARNAMEi ...
    DE_GAMVARPOWi DE_ISTGTFFILEi DE_TARGETi DE_INTVARNAMEi DE_INTVARPOWi

Gamma1 = thisDE{DE_GAMMA1i}; % index- and par-ready version ('compiled')
numGamma1Terms = length(Gamma1);
Gamma2 = thisDE{DE_GAMMA2i}; % index- and par-ready version ('compiled')
numGamma2Terms = length(Gamma2);
lenDomData = numGamma1Terms + numGamma2Terms; % number of terms in differential equation
Psi_data = zeros(1,lenDomData);
Term_data = zeros(1,lenDomData);

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
    targetVal = GetTargVal(g1term{DE_ISTGTFFILEi},g1term{DE_TARGETi},varVal,varDataLine);
    intIx = g1term{DE_INTVARNAMEi};
    if intIx > 0
        Psi_data(1) = g1term{DE_GAMVARPOWi}*tau_recipVal* varDataLine(intIx)^g1term{DE_INTVARPOWi} * ( - sumG2 ) / (sumG1*sumG1);
        Term_data(1) = tau_recipVal * varDataLine(g1term{DE_INTVARNAMEi})^g1term{DE_INTVARPOWi} * (targetVal - varDataLine(ext_var));
    else
        Term_data(1) = tau_recipVal * (targetVal - varDataLine(ext_var));
        Psi_data(1) = g1term{DE_GAMVARPOWi}*tau_recipVal * ( - sumG2 ) / (sumG1*sumG1);
    end
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
        % When ISTGTFFILE != 2 then the final argument is not used anyway!
%         if g1term{DE_ISTGTFFILEi} == 2
%             % target refers to another external variable, and 'source'
%             % index will be in the domain of varDataLine, not VarVals
        targetVal = GetTargVal(g1term{DE_ISTGTFFILEi},g1term{DE_TARGETi},varVals(pix),varDataLine);
%         else
%             targetVal = GetTargVal(g1term{DE_ISTGTFFILEi},g1term{DE_TARGETi},varVals(pix),varVals);
%         end
        gj_sum = sumG1 - tau_recipVal * intVarValsP(pix);
        G1j_sum = sumG1temp(2) - tau_recipVal * intVarValsP(pix) * targetVal;
        Psi_data(pix) = g1term{DE_GAMVARPOWi}*tau_recipVal*intVarValsP(pix) * ( targetVal * gj_sum - G1j_sum - sumG2 ) / (sumG1*sumG1);
        Term_data(pix) = tau_recipVal * intVarValsP(pix) * (targetVal - varDataLine(ext_var));
    end
end

for pix=1:numGamma2Terms
    g2term = Gamma2{pix};
    % note the '-1' in the power to tau_recipVal for Psi value
    tau_recipVal = GetTauRVal(g2term{DE_ISTAUFFILEi},g2term{DE_TAURECIPi},varDataLine(g2term{DE_GAMVARNAMEi}),g2term{DE_GAMVARPOWi}-1);
    % the following value for Psi is only valid if the Gamma2 term is NOT a file function
    Psi_data(numGamma1Terms+pix) = g2term{DE_GAMVARPOWi}*tau_recipVal / sumG1;
    % effectively restore the power of the variable in the expression for
    % the original term:
    Term_data(numGamma1Terms+pix) = tau_recipVal*varDataLine(g2term{DE_GAMVARNAMEi});
end

if unsigned
    Psi_data = abs( Psi_data );
    Term_data = abs( Term_data );
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function var_data = GetPsiPotVals(varDataLine_orig, thisDE, varBds, pot_ixAbs)
% calculates variable values for maximum Psi influence strength values (potentials) for candidate potential terms (Gamma1 and Gamma2) of dEqn differential equation
global DE_NAMEi DE_GAMMA1i DE_GAMMA2i DE_ACTSWi DE_ISTAUFFILEi DE_TAURECIPi DE_GAMVARNAMEi ...
    DE_GAMVARPOWi DE_ISTGTFFILEi DE_TARGETi DE_INTVARNAMEi DE_INTVARPOWi LARGEBOUND

Gamma1 = thisDE{DE_GAMMA1i}; % index- and par-ready version ('compiled')
Gamma2 = thisDE{DE_GAMMA2i}; % index- and par-ready version ('compiled')
numGamma1Terms = length(Gamma1);
numGamma2Terms = length(Gamma2);
var_data   = zeros(1,numGamma1Terms + numGamma2Terms); % different order to varDataLine_orig
sumG1temp  = SumGamma1(Gamma1,varDataLine_orig);
% sumG1_orig = sumG1temp(1);

if numGamma1Terms == 1
    var_data(1) = 0; % will make Psi -> infinity
else
    intVarValsP = zeros(1,numGamma1Terms); % power already included, where applicable
    for i=1:numGamma1Terms
        intIx = Gamma1{i}{DE_INTVARNAMEi};
        if intIx > 0
            intVarValsP(i) = varDataLine_orig(intIx)^Gamma1{i}{DE_INTVARPOWi};
        else
            intVarValsP(i) = 1;
        end
    end
    for pix=1:numGamma1Terms
        extIx = Gamma1{pix}{DE_GAMVARNAMEi};
        if extIx == 0
            extIx = Gamma1{pix}{DE_TARGETi};
        end
        if ~ismember(extIx, pot_ixAbs) % then this var is not a candidate potential, so ignore
            continue
        end
        varDataLine   = varDataLine_orig; % reset for each term
        g1term        = Gamma1{pix};
        gi            = GetTauRVal(g1term{DE_ISTAUFFILEi},g1term{DE_TAURECIPi},1,1);
        varDataLine(extIx) = 0; % only to get gj_sum right in a single sum
        gj_sum_temp   = SumGamma1(Gamma1,varDataLine);
        gj_sum        = gj_sum_temp(1);
        % greatest potential value is either at (the unique) turning point or interval end-points. quickest to maximize results from each.
        varMaxer      = [ varBds(extIx,1), varBds(extIx,2), max( varBds(extIx,1), min( varBds(extIx,2), ( gj_sum / gi )^(1/g1term{DE_GAMVARPOWi}) ) ) ];
        varMaxed      = [0 0 0];
        for i=1:3 % in varMaxer
            tau_recipVal  = GetTauRVal(g1term{DE_ISTAUFFILEi},g1term{DE_TAURECIPi},varMaxer(i),g1term{DE_GAMVARPOWi});
            targetVal     = GetTargVal(g1term{DE_ISTGTFFILEi},g1term{DE_TARGETi},varMaxer(i),varDataLine);
            varDataLine(extIx) = varMaxer(i);
            sumG1temp = SumGamma1(Gamma1,varDataLine);
            sumG1     = sumG1temp(1);
            sumG2     = SumGamma2(Gamma2,varDataLine);
            G1j_sum   = sumG1temp(2) - tau_recipVal * intVarValsP(pix) * targetVal;
            varMaxed(i) = tau_recipVal * ( targetVal * gj_sum - G1j_sum - sumG2 ) / (sumG1*sumG1);
        end
        [maxPsi maxIx] = max(abs(varMaxed));
        % was max(abs(varMaxed_term * g1term{DE_GAMVARPOWi} * intVarValsP(pix))) ... but
        % note that because the same positive power and intVarValsP multiplies each of these terms to be maxed, these can be dropped
        var_data(pix) = varMaxer(maxIx);
    end
end

for pix=1:numGamma2Terms % !!! Update this for internal values when they are introduced
%     g2term = Gamma2{pix};
    extIx =  Gamma2{pix}{DE_GAMVARNAMEi};
    if ~ismember(extIx, pot_ixAbs) % then this var is not a candidate potential, so ignore
        continue
    end
    var_data(numGamma1Terms+pix) = varBds(extIx,2);
end
% Note: non-file function Gamma 2 set psi formulae are monotonic in their variable, so are maximized at var val's max value.
% Therefore, the value for Psi calculated in main loop is only valid if the Gamma2 term is NOT a file function
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function var_data = GetTermPotVals(varDataLine_orig, thisDE, varBds, ext_var, pot_ixAbs)
% calculates variable values for maximum term size values (potentials) for candidate potential terms (Gamma1 and Gamma2) of dEqn differential equation
global DE_NAMEi DE_GAMMA1i DE_GAMMA2i DE_ACTSWi DE_ISTAUFFILEi DE_TAURECIPi DE_GAMVARNAMEi ...
    DE_GAMVARPOWi DE_ISTGTFFILEi DE_TARGETi DE_INTVARNAMEi DE_INTVARPOWi LARGEBOUND

Gamma1 = thisDE{DE_GAMMA1i}; % index- and par-ready version ('compiled')
Gamma2 = thisDE{DE_GAMMA2i}; % index- and par-ready version ('compiled')
numGamma1Terms = length(Gamma1);
numGamma2Terms = length(Gamma2);
var_data   = zeros(1,numGamma1Terms + numGamma2Terms); % different order to varDataLine_orig
sumG1temp  = SumGamma1(Gamma1,varDataLine_orig);
% sumG1_orig = sumG1temp(1);

if numGamma1Terms == 1
    var_data(1) = 0; % will make Psi -> infinity
else
    intVarValsP = zeros(1,numGamma1Terms); % power already included, where applicable
    for i=1:numGamma1Terms
        intIx = Gamma1{i}{DE_INTVARNAMEi};
        if intIx > 0
            intVarValsP(i) = varDataLine_orig(intIx)^Gamma1{i}{DE_INTVARPOWi};
        else
            intVarValsP(i) = 1;
        end
    end
    for pix=1:numGamma1Terms
        extIx = Gamma1{pix}{DE_GAMVARNAMEi};
        if extIx == 0
            extIx = Gamma1{pix}{DE_TARGETi};
        end
        if ~ismember(extIx, pot_ixAbs) % then this var is not a candidate potential, so ignore
            continue
        end
        varDataLine   = varDataLine_orig; % reset for each term
        g1term        = Gamma1{pix};
        % greatest potential value is at an interval end-point.
        varMaxer      = [varBds(extIx,1), varBds(extIx,2)];
        varMaxed      = [0 0];
        for i=1:2 % in varMaxer
            tau_recipVal  = GetTauRVal(g1term{DE_ISTAUFFILEi},g1term{DE_TAURECIPi},varMaxer(i),g1term{DE_GAMVARPOWi});
            targetVal     = GetTargVal(g1term{DE_ISTGTFFILEi},g1term{DE_TARGETi},varMaxer(i),varDataLine);
            varDataLine(extIx) = varMaxer(i);
            varMaxed(i) = tau_recipVal * intVarValsP(pix) * ( targetVal - varDataLine(ext_var) );
        end
        [maxTerm maxIx] = max(abs(varMaxed));
        var_data(pix) = varMaxer(maxIx);
    end
end

for pix=1:numGamma2Terms % !!! Update this for internal values when they are introduced
%     g2term = Gamma2{pix};
    extIx =  Gamma2{pix}{DE_GAMVARNAMEi};
    if ~ismember(extIx, pot_ixAbs) % then this var is not a candidate potential, so ignore
        continue
    end
    var_data(numGamma1Terms+pix) = varBds(extIx,2);
end
% Note: non-file function Gamma 2 set terms are monotonic in their variable, so are maximized at var val's max value.
% Therefore, the value for calculated in main loop above is only valid if the Gamma2 term is NOT a file function
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sum = SumGamma1(Gamma1,varDataLine)
% only accepts "compiled" Gamma1
% gam1term order: actSw, filefuncflag, taurecip, var, power, filefuncflag, target [, intvar, power]
global DE_NAMEi DE_GAMMA1i DE_GAMMA2i DE_ACTSWi DE_ISTAUFFILEi DE_TAURECIPi DE_GAMVARNAMEi ...
    DE_GAMVARPOWi DE_ISTGTFFILEi DE_TARGETi DE_INTVARNAMEi DE_INTVARPOWi

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
    tau_recipVal = GetTauRVal(g1term{DE_ISTAUFFILEi},g1term{DE_TAURECIPi},varVal,g1term{DE_GAMVARPOWi});
    intIx = g1term{DE_INTVARNAMEi};
    if intIx > 0
        intVarValP = varDataLine(intIx)^g1term{DE_INTVARPOWi};
    else
        intVarValP = 1;
    end
    targetVal = GetTargVal(g1term{DE_ISTGTFFILEi},g1term{DE_TARGETi},varVal,varDataLine);
    tempVal = tau_recipVal * intVarValP;
    sum1 = sum1 + tempVal;
    sum2 = sum2 + tempVal * targetVal;
end
sum = [sum1, sum2];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sum = SumGamma2(Gamma2,varDataLine)
% only accepts "compiled" Gamma2
% gam2term order: actSw, filefuncflag, taurecip, var, power
global DE_NAMEi DE_GAMMA1i DE_GAMMA2i DE_ACTSWi DE_ISTAUFFILEi DE_TAURECIPi DE_GAMVARNAMEi ...
    DE_GAMVARPOWi DE_ISTGTFFILEi DE_TARGETi DE_INTVARNAMEi DE_INTVARPOWi

sum = 0;
for g2t=1:length(Gamma2)
    g2term = Gamma2{g2t};
    varVal = varDataLine(g2term{DE_GAMVARNAMEi});
    tau_recipVal = GetTauRVal(g2term{DE_ISTAUFFILEi},g2term{DE_TAURECIPi},varVal,g2term{DE_GAMVARPOWi});
    sum = sum + tau_recipVal;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = LookupDEpar(parName,DEpars)
result = {0,false}; % default value
for parIx=1:length(DEpars)
    if strcmp(parName,DEpars{parIx}{1})
        result = {DEpars{parIx}{2},true};
        break
    end
    % should never get here!
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% useful as a dummy argument for KeyPressFcn callback, to prevent
% echo of key presses to command window
function DoNothing()
return
