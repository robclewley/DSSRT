function setup = ReadConfigFile(filename,networkObjDir)
% Function to interpret Functional Network setup `.cfg` file
%
% Version 1.11, (c) Robert Clewley, Center for Biodynamics, Boston University, 2004
%   and Department of Mathematics, Cornell University, 2005.

%%%%%%% Version history
% v1.11 November 2005.
%  Fixed bug in Gamma1 target = external variable.
%
% v1.10 June 2005.
%  Added optional CFAC spec, for each differential equation, to allow
%    capacitance to change from unity. If omitted, CFAC = 1.
%  Added support for Gamma1 target to be an external variable, to support
%    electrical synapses.
%
% v1.00 May 2004.
%  Created.

%%%% MUST HAVE THE FOLLOWING SPECIFICATIONS PRESENT, IN ANY ORDER (EXCEPT NODE* SPECS), IN THE .cfg FILE
%
% VARSINT <varname> [ <varname> ... ]
% VARSEXT <varname> [ <varname> ... ]
% INPUTS <dest_var> [ <varname> ... ]
% BOUNDS <varname> <lo_val> <hi_val>
% UNITBOUNDS <varname> [ <varname> ... ]
% DEQNS <eqn_subject_name> [<eqn_subject_name> ...]
% GAM1TERM <eqn_subject_name> <taureciprocal> <extvar> <power> <target> [<intvar> <power>]
% GAM2TERM <eqn_subject_name> <taureciprocal> <extvar> <power>
% CFAC <eqn_subject_name> <value>
% DEPAR <parname> <value>
% NODE <obslabel> <x-centre> <y-centre> <radius>
% NODELABEL <obslabel> <name> <x-textpos> <y-textpos> <textsize>
% NODEADDSTATES <obslabel> <numstates> <style> [<style> ... ]
% NODESTATEMAP_ACTS <obslabel> [<statemapfunction> ... ]
% NODESTATEMAP_POTS <obslabel> [<statemapfunction> ... ]
% LINK <source_obslabel> <dest_obslabel> <x-tail> <y-tail> <x-head> <y-head> [ <visibility-switch> ]
% VBAR <varlabel> <x-pos> <y-botpos> <y-height> <variablelo_val> <variablehi_val> <log-scale-switch>

setup = {}; % initial and temporary value (error exit status)
errorFlag = 0;

if nargin ~= 2
    errorFlag = true;
end

if isempty(filename) | exist(filename,'file') ~= 2
    errorFlag = true;
else
    if ~strcmp('.cfg',filename(length(filename)-3:length(filename)))
        errorFlag = true;
    end
end

if exist( networkObjDir, 'dir' ) ~= 7
    errorFlag = true;
end

if errorFlag % errorFlag only used here
    beep
    disp('ReadConfigFile:  Need a valid configuration filename with extension `.cfg`')
    disp('                   and a valid network object directory')
    return
end

fid = fopen(filename,'r');

%%% INITIALIZATIONS AND GLOBAL VARIABLES
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
global LARGEBOUND SMALLBOUND TIMERES LINKOFFSTATE ST_ERR NOT_AP POTENT ACTIVE NUMPARS

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



global VARSINTi VARSEXTi INPUTSi CANDACTSi CANDPOTSi NODEi NODELABELi NODEADDSTATESi
global NODESTATEMAP_ACTSi NODESTATEMAP_POTSi LINKi VBARi DEQNSi GAM1TERMi GAM2TERMi DEPARi
global BOUNDSi UNITBOUNDSi CFACi

% indices for setup_done array
VARSINTi       = 1;
VARSEXTi       = 2;
INPUTSi        = 3;
NODEi          = 4;
NODELABELi     = 5;
NODEADDSTATESi = 6;
NODESTATEMAP_ACTSi = 7;
NODESTATEMAP_POTSi = 8;
LINKi          = 9;
VBARi          = 10;
DEQNSi         = 11;
GAM1TERMi      = 12;
GAM2TERMi      = 13;
DEPARi         = 14;
BOUNDSi        = 15;
UNITBOUNDSi    = 16;
CFACi          = 17;

numCommands = 17;

commands = {'VARSINT', 'VARSEXT', 'INPUTS', 'NODE', 'NODELABEL', ...
       'NODEADDSTATES', 'NODESTATEMAP_ACTS', 'NODESTATEMAP_POTS', ...
       'LINK', 'VBAR', 'DEQNS', 'GAM1TERM', 'GAM2TERM', 'DEPAR', ...
       'BOUNDS', 'UNITBOUNDS', 'CFAC'};

%%%% Other initializations
rawArgs = cell(1,numCommands);
comArgs = cell(1,numCommands);
% The specifications must each be preceded by `command_str` at the beginning of a line
command_str = ''; % was '@FN ' in earlier versions, but turned out not to be useful
ignoreOccSet = [12,13,14,17]; % don't need at least one of these commands

% setup_done array holds counts for how many times each specification has been found
setup_done = zeros(1,numCommands);
diagObjects = {};
vBarObjects = {};

%%%% READ FILE, GET INFO ON DSSRT COMMANDS CONTAINED IN IT
% initialize rawArgs, comArgs
for c = 1:numCommands
    rawArgs{c} = {};
    comArgs{c} = {};
end

file_done = false;
lines = 0;

while ~file_done
    fline = fgetl(fid);
    if ~ischar(fline) % then EOF
        file_done = true;
        break
    end
    
    isDSSRTcom = IsCfgCommand(fline, commands, command_str);
    ctest = isDSSRTcom{1};
    args  = isDSSRTcom{2};
    [isCom comIx] = ismember(1,ctest); % presumes '1' can only appear once in ctest!
    if isCom % then a specification was found
        setup_done(comIx) = setup_done(comIx) + 1; % value n indicates we found specification with index comIx n times
        rawArgs{comIx} = {rawArgs{comIx}{:} args};
    end
end

fclose(fid);

%%%% COMPILE COMMANDS FOUND
compile = true; % initial value only
maxNumOccs = 0; % initial value for maximum number of specification occurrences (to get necessary array size for comLens)
for c=1:numCommands
    numOccs = setup_done(c);
    if numOccs == 0 % then no occurrences of a certain specification were found
        if ~ismember(c,ignoreOccSet) % do not worry about occurrences of DEpars, Gam1Terms, Gam2Terms, Cfac
            compile = false;
        end
    end
    if numOccs > maxNumOccs
        maxNumOccs = numOccs;
    end
end

if compile % then all relevant specifications found at least once ...

    comLens = zeros(numCommands,maxNumOccs);

    for c = 1:numCommands
        for i=1:setup_done(c)
            processed = ProcessCfgArgs(rawArgs{c}{i});
            comLens(c,i) = processed{1};
            comArgs{c} = { comArgs{c}{:} processed{2} };
        end
    end

    if setup_done(VARSINTi) == 1
        numInt = comLens(VARSINTi,1);
    else
        disp('ReadConfigFile:  More than one occurrence of VARSINT found')
        return
    end
    if setup_done(VARSEXTi) == 1
        numExt = comLens(VARSEXTi,1);
    else
        disp('ReadConfigFile:  More than one occurrence of VARSEXT found')
        return
    end
    numTot = numExt+numInt;

    % check that there's no overlap between external and internal variables
    extVarNames = comArgs{VARSEXTi}{1};
    intVarNames = comArgs{VARSINTi}{1};
    if isempty( intersect( extVarNames, intVarNames ) )
        allVarNames = { extVarNames{:}, intVarNames{:} };
    else
        disp('ReadConfigFile:  There is overlap between internal and external variable names!')
        return
    end

    % there should be numTot number of INPUTS specs
    if setup_done(INPUTSi) ~= numTot
        disp(   'ReadConfigFile:  Incorrect number of INPUTS specifications.')
        fprintf('                 Found %i/%i.\n', setup_done(INPUTSi), numTot)
        return
    end
    
    % Bounds setup
    if setup_done(BOUNDSi) + setup_done(UNITBOUNDSi) == 0
        disp('ReadConfigFile:  Bounds must be specified for every variable')
        return
    end
    
    % Differential Equations setup
    if setup_done(DEQNSi) ~= 1
        disp('ReadConfigFile:  There must be exactly one DEQNS specification')
        return
    end
    numDEpars = setup_done(DEPARi);
    numGam1Terms = setup_done(GAM1TERMi);
    numGam2Terms = setup_done(GAM2TERMi);
    numCfacs = setup_done(CFACi);

    % Diagram setup
    numNodes = setup_done(NODEi);
    numLinks = setup_done(LINKi);
    numObjects = numNodes + numLinks;
    
    % There should be numNodes # of the following specifications:
    % NODELABEL, NODEADDSTATES, NODESTATEMAP_ACTS, NODESTATEMAP_POTS
    if setup_done(NODELABELi) ~= numNodes | setup_done(NODEADDSTATESi) ~= numNodes | setup_done(NODESTATEMAP_ACTSi) ~= numNodes | setup_done(NODESTATEMAP_POTSi) ~= numNodes
        disp(   'ReadConfigFile:  Incorrect number of NODELABEL, NODEADDSTATES, NODESTATEMAP_ACTS,')
        disp(   '                   NODESTATEMAP_POTS specifications.')
        fprintf('                 %i are required. You specified %i, %i, %i, and %i, for the above specs, respectively.\n', ...
                 numNodes, setup_done(NODELABELi), setup_done(NODEADDSTATESi), setup_done(NODESTATEMAP_ACTSi),  setup_done(NODESTATEMAP_POTSi) )
        return
    end
    
    if numNodes ~= numExt
        disp(   'ReadConfigFile:  Incorrect number of NODE specifications')
        fprintf('                   Expected %i\n',numExt)
        return
    end

    numVbars = setup_done(VBARi);
    if numVbars ~= numTot
        fprintf('ReadConfigFile: Incorrect number of VBAR specifications. Expected %i\n',numTot)
        return
    end

    inputsIx = cell(1,numTot); % max # possible
    candActsIx = cell(1,numExt); % max # possible
    candPotsIx = cell(1,numExt); % max # possible
    totNumObsInputs = 0; % number of connections between observable variables
    for i = 1:numTot
        % INPUTS ... can be internal or external variables
        numArgs = comLens(INPUTSi,i);
        % notice in Fetch... call below, 3rd arg is allVarNames, not extVarNames, because
        %   Inputs are specified for all external *and* internal variables
        result = FetchArgsSubjVars(comArgs{INPUTSi}{i}, numArgs, allVarNames, allVarNames, 'an', 'all', 'INPUT', true);
        if ~isempty(result)
            subjectIx = result{1};
            argsLen = length(result{2});
            if argsLen > 0
                seenIx = zeros(1,numTot);
                for argnum = 1:argsLen
                    if ismember(result{2}(argnum),[1:numTot])
                        if seenIx(result{2}(argnum)) == 0
                            seenIx(result{2}(argnum))=1;
                        else
                            disp('ReadConfigFile:  Fatal Error. Arguments to INPUTS repeated')
                            return
                        end
                    end
                    if ismember(subjectIx,[1:numExt]) & ismember(result{2}(argnum),[1:numExt])
                        % then it's an external variable -> external variable link
                        % (need to know this total to compare with # of declared LINKs)
                        totNumObsInputs = totNumObsInputs + 1;
                    end
                end
            end
            inputsIx{subjectIx} = result{2};
        else
            disp('ReadConfigFile:  Fatal Error (INPUTS). Returning.')
            return
        end
    end

    
    %%%% BOUNDS ... for all variables
    numBdComs = setup_done(BOUNDSi);
    numUBdComs = setup_done(UNITBOUNDSi);
    seenVarList = [];
    varBounds = zeros(numTot,2);
    % unit bounds
    for i=1:numUBdComs
        subjNameList = comArgs{UNITBOUNDSi}{i};
        listLen = comLens(UNITBOUNDSi,i); % = length(subjNameList);
        for subjIx=1:listLen
            thisSubjName = subjNameList{subjIx};
            thisVarIx = VarnameIxMap(thisSubjName,allVarNames);
            if thisVarIx > 0
                if ~ismember(thisVarIx, seenVarList)
                    seenVarList = [seenVarList, thisVarIx];
                    varBounds(thisVarIx,1) = 0;
                    varBounds(thisVarIx,2) = 1;
                else
                    fprintf('ReadConfigFile:  Variable names can only be declared in one bounds specification (%s)\n',thisSubjName)
                    return
                end
            else
                fprintf('ReadConfigFile:  Variable name %s in UNITBOUNDS specification not known\n',thisSubjName)
                return
            end
        end
    end
    
    % other bounds
    for i=1:numBdComs
        argsList = comArgs{BOUNDSi}{i};
        if comLens(BOUNDSi,i) ~= 3 % length(argsList)
            disp('ReadConfigFile:  Incorrect number of parameters to BOUNDS command. Expected 3.')
            return
        end
        subjName = argsList{1};
        thisVarIx = VarnameIxMap(subjName,allVarNames);
        if thisVarIx > 0
            if ~ismember(thisVarIx, seenVarList)
                loBdStr = argsList{2};
                hiBdStr = argsList{3};
                if isNum( loBdStr ) & isNum( hiBdStr )
                    loBdVal = str2num(loBdStr);
                    hiBdVal = str2num(hiBdStr);
                    if loBdVal == hiBdVal
                        fprintf('ReadConfigFile:  Upper and lower bounds cannot be equal for variable %s\n',subjName)
                        return
                    end
                    seenVarList = [seenVarList, thisVarIx];
                    varBounds(thisVarIx,1) = min(loBdVal,hiBdVal); % ignores order given!
                    varBounds(thisVarIx,2) = max(loBdVal,hiBdVal);
                else
                    fprintf('ReadConfigFile: Second and third parameters to BOUNDS specification for variable %s\n',subjName)
                    disp(   '                  must be numerical values')
                    return
                end
            else
                fprintf('ReadConfigFile:  Variable names can only be declared in one bounds specification (%s)\n',subjName)
                return
            end
        else
            fprintf('ReadConfigFile:  Variable name %s in BOUNDS specification not known\n',subjName)
            return
        end
    end
    
    if length(seenVarList) ~= numTot
        disp('ReadConfigFile:  Not all variable names have bounds specified.')
        return
    end


    %%%% DIFFERENTIAL EQUATIONS SETUP
    % Differential Equations
    numDEqns = comLens(DEQNSi,1);
    DEqns = cell(1,numDEqns);
    for i=1:numDEqns
        eqnSubjName = comArgs{DEQNSi}{1}{i};
        eqnSubjIx = VarnameIxMap(eqnSubjName,allVarNames);
        if eqnSubjIx == 0
            disp('ReadConfigFile:  Arguments to DEQNS must be internal or external variable names')
            return
        end
        DEqns{i} = {eqnSubjName, {}, {}};
    end
    
    % DE pars
    DEpars = {};
    seenParList = {};
    if numDEpars > 0
        for i=1:numDEpars
            numArgs = comLens(DEPARi,i);
            if numArgs ~=2
                disp('ReadConfigFile:  Only two arguments expected for DEPAR specification')
                return
            end
            parName = comArgs{DEPARi}{i}{1};
            if ismember(parName,seenParList)
                disp(['ReadConfigFile:  Parameter ' parName ' repeated. Expected unique occurrence.'])
                return
            else
                if ~ismember( parName, allVarNames )
                    seenParList = {seenParList{:}, parName};
                else
                    disp(['ReadConfigFile:  Parameter name ' parName ' clashes with variable names.'])
                    return
                end
            end
            parVal  = comArgs{DEPARi}{i}{2};
            if isNum(parVal)
                DEpars  = { DEpars{:}, {parName, str2num(parVal)} };
            else
                fprintf('ReadConfigFile:  Occurrence %i of DEPAR specification had non-numeric value\n',i)
                return
            end
        end
    end
    
    % Gam1Terms
    Gam1Terms = {};
    %seenInputList = {};
    if numGam1Terms > 0
        for i=1:numGam1Terms
            numArgs = comLens(GAM1TERMi,i);
            if numArgs ~= 5 & numArgs ~= 7
                fprintf('ReadConfigFile:  GAM1TERM #%i\n',i)
                disp('                 Wrong number of arguments passed to GAM1TERM specification. Expected 5 or 7.')
                return
            end
            eqnSubjName = comArgs{GAM1TERMi}{i}{1}; % equation subject name
            isEqnSubjIx = false;
            eqnSubjIx = 0;
            for eqnIx=1:numDEqns
                if strcmp(eqnSubjName,DEqns{eqnIx}{1})
                    isEqnSubjIx = true;
                    thisEqnIx = eqnIx;
                    [p eqnSubjIx] = ismember(eqnSubjName,allVarNames); % p not used
                    break
                end
            end
            if ~isEqnSubjIx % then not a var name associated with an equation
                fprintf('ReadConfigFile:  GAM1TERM #%i\n',i)
                disp('                 Subject of GAM1TERM (arg 1) must be the name of a system variable')
                disp('                   associated with a differential equation.')
                return
            end
            % should check that each input variable appears only once in Gam1Terms (and also once out of both Gamma sets)
	
            arg2name = comArgs{GAM1TERMi}{i}{2};
            arg2isDEpar = isDEpar(arg2name,DEpars);
            arg2isFileFunc = exist([ networkObjDir '/' arg2name '.m'],'file')==2;
            if arg2isDEpar & arg2isFileFunc
                fprintf('ReadConfigFile:  GAM1TERM #%i\n',i)
                disp('                 2nd argument to GAM1TERM cannot be both a DE parameter and file-specified function.')
                return
            end
            if ~arg2isDEpar & ~arg2isFileFunc
                fprintf('ReadConfigFile:  GAM1TERM #%i\n',i)
                disp('                 2nd argument to GAM1TERM must be either a specified DE parameter or a file-specified function.')
                return
            end
	
            arg3name = comArgs{GAM1TERMi}{i}{3};
            arg3isDEpar = isDEpar(arg3name,DEpars);
            [proceed1 objIx] = ismember(arg3name, allVarNames);
            if proceed1
                proceed2 = ismember(objIx,inputsIx{eqnSubjIx});
            else
                % other alternative is for arg3 to be a DE parameter
                proceed2 = arg3isDEpar;
            end
            if ~proceed2
                fprintf('ReadConfigFile:  GAM1TERM #%i\n',i)
                disp('ReadConfigFile:  3rd argument to GAM1TERM must be the name of a')
                disp('               variable input to this subject variable, or a DE parameter')
                return
            end
	
            arg4 = comArgs{GAM1TERMi}{i}{4};
            power = 0;
            if isNum( arg4 )
                arg4val = str2num(arg4);
                if arg4val >= 0 & arg4val == round(arg4val) % then non-negative integer
                    power = arg4val;
                else
                    fprintf('ReadConfigFile:  GAM1TERM #%i\n',i)
                    disp('                 4th argument to GAM1TERM must be a non-negative integer.')
                    return
                end
            end
	
            arg5name = comArgs{GAM1TERMi}{i}{5};
            arg5isDEpar = isDEpar(arg5name,DEpars);
            arg5isFileFunc = exist([ networkObjDir '/' arg5name '.m'],'file')==2;
            arg5isExtVar = ismember(arg5name, extVarNames); % e.g. for "electrical coupling"
            numtrue = sum([arg5isDEpar, arg5isFileFunc, arg5isExtVar]);
            if numtrue > 1
                fprintf('ReadConfigFile:  GAM1TERM #%i\n',i)
                disp('                 5th argument to GAM1TERM cannot be both a DE parameter / variable and file-specified function.')
                return
            end
            if numtrue == 0
                fprintf('ReadConfigFile:  GAM1TERM #%i\n',i)
                disp('                 5th argument to GAM1TERM must be either a specified DE parameter, external variable, or a file-specified function.')
                return
            end
            if arg5isDEpar % DE parameter
                if arg3isDEpar
                    fprintf('ReadConfigFile:   GAM1TERM #%i\n', i)
                    disp('                 5th argument cannot be a DE parameter is 3rd argument is a DE parameter')
                    return
                end
                arg5sourceType = 0;
            elseif arg5isFileFunc % file function
                arg5sourceType = 1;
            else % external variable
                if strcmp(eqnSubjName,arg5name)
                    fprintf('ReadConfigFile:  GAM1TERM #%i\n',i)
                    disp('                 5th argument to GAM1TERM must not be the name of the subject equation`s dependent variable.')
                    return                    
                else
                    [proceed1 objIx] = ismember(arg5name, allVarNames);
                    proceed2 = ismember(objIx,inputsIx{eqnSubjIx});
                    if ~proceed2
                        fprintf('ReadConfigFile:  GAM1TERM #%i\n',i)
                        disp('ReadConfigFile:  5th argument to GAM1TERM is not the name of a')
                        disp('               variable input to this subject variable.')
                        return
                    end
                    arg5sourceType = 2;
                end
            end
	
            intVarName = '';
            intVarPower = 0;
            if numArgs == 7
                arg6name = comArgs{GAM1TERMi}{i}{6};
                [proceed1 objIx] = ismember(arg6name, intVarNames);
                proceed2 = false;
                if proceed1 % then arg6 is an internal var name
                    if ismember(objIx+numExt,inputsIx{eqnSubjIx})
                        proceed2 = true;
                    end
                end
                if ~proceed2
                    fprintf('ReadConfigFile:  GAM1TERM #%i\n',i)
                    disp('                 6th argument to GAM1TERM must be the name of an internal')
                    disp('                   variable present in inputs for this subject variable.')
                    return
                else
                    intVarName = arg6name;
                    arg7 = comArgs{GAM1TERMi}{i}{7};
                    if isNum( arg7 )
                        arg7val = str2num(arg7);
                        if arg7val > 0 & arg7val == round(arg7val) % then natural number
                            intVarPower = arg7val;
                        else
                            fprintf('ReadConfigFile:  GAM1TERM #%i\n',i)
                            disp('                 7th argument to GAM1TERM must be a natural number.')
                            return
                        end
                    end
                end
            end
            % gam1term order:
            % actSw, filefuncflag, taurecip, var, power, filefuncflag, target [, intvar, power]
            DEqns{thisEqnIx}{DE_GAMMA1i} = { DEqns{thisEqnIx}{DE_GAMMA1i}{:}, {0, arg2isFileFunc, arg2name, arg3name, power, arg5sourceType, arg5name, intVarName, intVarPower } };
        end % for i
    end % if

    % Gam2Terms
    Gam2Terms = {};
    %seenInputList = {};
    if numGam2Terms > 0
        for i=1:numGam2Terms
            numArgs = comLens(GAM2TERMi,i);
            if numArgs ~= 4
                fprintf('ReadConfigFile:  GAM2TERM #%i\n',i)
                disp('                 Wrong number of arguments passed to GAM2TERM specification. Expected 4.')
                return
            end
            eqnSubjName = comArgs{GAM2TERMi}{i}{1}; % equation subject name
            isEqnSubjIx = false;
            eqnSubjIx = 0;
            for eqnIx=1:numDEqns
                if strcmp(eqnSubjName,DEqns{eqnIx}{1})
                    isEqnSubjIx = true;
                    thisEqnIx = eqnIx;
                    [p eqnSubjIx] = ismember(eqnSubjName,allVarNames); % p not used
                    break
                end
            end
            if ~isEqnSubjIx % then not a var name associated with an equation
                fprintf('ReadConfigFile:  GAM2TERM #%i\n',i)
                disp('                 Subject of GAM2TERM (arg 1) must be the name of a system variable')
                disp('                   associated with a differential equation.')
                return
            end
            % should check that each input variable appears only once in Gam2Terms (and also once out of both Gamma sets)
	
            arg2name = comArgs{GAM2TERMi}{i}{2};
            arg2isDEpar = isDEpar(arg2name,DEpars);
            arg2isFileFunc = exist([ networkObjDir '/' arg2name '.m'],'file')==2;
            if arg2isDEpar & arg2isFileFunc
                fprintf('ReadConfigFile:  GAM2TERM #%i\n',i)
                disp('                 2nd argument to GAM2TERM cannot be both a DE parameter and file-specified function.')
                return
            end
            if ~arg2isDEpar & ~arg2isFileFunc
                fprintf('ReadConfigFile:  GAM2TERM #%i\n',i)
                disp('                 2nd argument to GAM2TERM must be either a specified DE parameter or a file-specified function.')
                return
            end
	
            arg3name = comArgs{GAM2TERMi}{i}{3};
            [proceed1 objIx] = ismember(arg3name, allVarNames);
            proceed2 = false;
            if ismember(objIx,inputsIx{eqnSubjIx})
                proceed2 = true;
            end
            if ~proceed2
                fprintf('ReadConfigFile:  GAM2TERM #%i\n',i)
                disp('                 3rd argument to GAM2TERM must be the name of a')
                disp('                   variable input to this subject variable.')
                return
            end
	
            arg4 = comArgs{GAM2TERMi}{i}{4};
            power = 0;
            if isNum( arg4 )
                arg4val = str2num(arg4);
                if arg4val >= 0 & arg4val == round(arg4val) % then non-negative integer
                    power = arg4val;
                else
                    fprintf('ReadConfigFile:  GAM2TERM #%i\n',i)
                    disp('                 4th argument to GAM2TERM must be a non-negative integer.')
                    return
                end
            end
            % gam2term order:
            % actSw, filefuncflag, taurecip, var, power
            DEqns{thisEqnIx}{DE_GAMMA2i} = { DEqns{thisEqnIx}{DE_GAMMA2i}{:}, {0, arg2isFileFunc, arg2name, arg3name, power } };
        end % for i
    end % if

    % Cfac (optional, defaults to 1 for each DE)
    seenEqnList = [];
    cfacs = ones(1,numDEqns); % defaults
    if numCfacs > 0
        for i=1:numCfacs
            numArgs = comLens(CFACi,i);
            if numArgs ~= 2
                fprintf('ReadConfigFile:  CFAC #%i\n',i)
                disp('                 Wrong number of arguments passed to CFAC specification. Expected 2.')
                return
            end
            eqnSubjName = comArgs{CFACi}{i}{1}; % equation subject name
            isEqnSubjIx = false; % initial value
            eqnSubjIx = 0;
            for eqnIx=1:numDEqns
                if strcmp(eqnSubjName,DEqns{eqnIx}{1})
                    isEqnSubjIx = true;
                    thisEqnIx = eqnIx;
                    if ismember(eqnIx, seenEqnList)
                        fprintf('ReadConfigFile:  CFAC #%i\n',i)
                        disp(['                 Duplicated CFAC command for Diff Eq in variable ' eqnSubjName])
                        return
                    else
                        seenEqnList = [ seenEqnList, eqnIx ];
                    end
                    [p eqnSubjIx] = ismember(eqnSubjName,allVarNames); % p not used
                    break
                end
            end
            if ~isEqnSubjIx % then not a var name associated with an equation
                fprintf('ReadConfigFile:  CFAC #%i\n',i)
                disp('                 Subject of CFAC (arg 1) must be the name of a system variable')
                disp('                   associated with a differential equation.')
                return
            end
            % should check that each input variable appears only once
	
            arg2name = comArgs{CFACi}{i}{2};
            lookupResult = LookupDEpar(arg2name,DEpars);
            arg2isDEpar = lookupResult{2};
            if ~arg2isDEpar
                fprintf('ReadConfigFile:  CFAC #%i\n',i)
                disp('                 2nd argument to CFAC must be a DE parameter.')
                return
            end
            if lookupResult{1} ~= 0
                cfacs(thisEqnIx) = lookupResult{1}; % 'compile' this value right now
            else
                fprintf('ReadConfigFile:  CFAC #%i\n',i)
                disp('                 2nd argument: DE parameter value must not be 0.')
                return
            end
        end
    end
    % update DEqns
    for eqn = 1:numDEqns
        DEqns{eqn}{DE_CFACi} = cfacs(eqn);
    end
    
    
    % For each Eqn, there should be no more than INPUTS{eqn} number of GAM1TERMs + GAM2TERMs
    % Also, for fewer than that number, those present in INPUTS that are not in a Gamma term
    % may only be internal variables. Note that for electrical synapses the
    % INPUT may appear as a target value.
    for eqn=1:numDEqns
        eqnSubjName = DEqns{eqn}{DE_NAMEi};
        [present eqnSubjIx] = ismember(eqnSubjName, allVarNames);
        subjIps = inputsIx{eqnSubjIx};
        expectedNum = length(subjIps);
        Gam1 = DEqns{eqn}{DE_GAMMA1i};
        Gam2 = DEqns{eqn}{DE_GAMMA2i};
        numGam1Terms = length(Gam1);
        numGam2Terms = length(Gam2);
        numGamTot = numGam1Terms + numGam2Terms;
        if expectedNum ~= numGamTot
            if expectedNum < numGamTot
                disp(   'ReadConfigFile:  Incorrect number of GAM1TERM or GAM2TERM specifications')
                disp(  ['                   for equation ' eqnSubjName '.'])
                fprintf('                 Expected no more than %i specs. Found %i.\n',expectedNum,numGamTot);
                return
            else % >  and so there are additional checks to be made
                for ipIx=1:expectedNum
                    ipNameFound = false;
                    ipName = allVarNames{subjIps(ipIx)};
                    for g1t=1:numGam1Terms
                        g1term = Gam1{g1t};
                        break_out = false;
                        if ismember( ipName, g1term{DE_GAMVARNAMEi} )
                            ipNameFound = true;
                            break_out = true;
                        end
                        if ismember( ipName, g1term{DE_TARGETi} )
                            if break_out
                                disp('ReadConfigFile:  Cannot match input name to both gating variable and target value')
                                return
                            else
                                ipNameFound = true;
                                break_out = true;
                            end
                        end
                        if break_out
                            break % for g1t
                        end
                    end
                    if ~ipNameFound
                        for g2t=1:numGam2Terms
                            g2term = Gam2{g2t};
                            if ismember( ipName, g2term{DE_GAMVARNAMEi} )
                                ipNameFound = true;
                                break % for g2t
                            end
                        end 
                    end
                    if ~ipNameFound % then it's in input list but not Gamma sets, so it had better be an internal variable
                        if ~ismember( ipName, intVarNames )
                            disp(['ReadConfigFile:  Inputs to equation ' eqnSubjName ' without GAM1TERM or'])
                            disp( '                   GAM2TERM specs must be internal variables to the subject')
                            disp( '                   (or, in the case of electrical synapses, the target value)')
                            return
                        end
                    end
                end % for ipIx
            end
        end
    end

    deqnIxMap = zeros(1,numTot); % zero entries will denote no corresponding differential equation present
    for varIx = 1:numTot
        for eqnix=1:numDEqns % create map entries
            if strcmp(DEqns{eqnix}{DE_NAMEi},allVarNames{varIx})
                deqnIxMap(varIx) = eqnix;
                break % inner for loop: done!
            end
        end
    end

    %%% Candidate actives / potentials
    % prepare list of cross-multiplying internal variables to avoid for candidate actives
    crossMultIntList = [];
    for eqnix=1:numDEqns
        Gamma1 = DEqns{eqnix}{DE_GAMMA1i}; % will have to do this for Gamma2 later on (when it's allowed to have cross-mult terms)
        for g1t = 1:length(Gamma1)
            internalName = Gamma1{g1t}{DE_INTVARNAMEi};
            [ pres internalIx ] = ismember( internalName, allVarNames ); % will catch empty internalName
            if pres && ~ismember(internalIx, crossMultIntList)
                crossMultIntList = [ crossMultIntList internalIx ];
            end
        end
	end
    % compile candidate actives for the observables having associated differential equations
    for i=1:numExt
        if deqnIxMap(i) > 0 % then a DE exists for this external variable, so it has candidate actives
            candActsIx{i} = setdiff( inputsIx{i}, crossMultIntList ); % but mask any cross-multiplying internal variables
            potsList = [];
            for cAix = candActsIx{i}
%                 if cAix <= numExt && deqnIxMap(cAix) > 0 % then is an observable that can vary in time - therefore a valid potential if not active
%                     % for now, don't care about potentials for internal variables (even if they are dummy ones that vary in time but don't have a DE)
				if deqnIxMap(cAix) > 0 % then can vary in time - therefore a valid potential if not active
                    potsList = [ potsList cAix ];
                end
            end
            candPotsIx{i} = potsList;
        end
    end

	%%%% DIAGRAM SETUP
	diagObjects = cell(numObjects,1);

    % Nodes...
    seenIx = zeros(1,numExt); % numNodes == numExt already verified
    for i = 1:numNodes
        % NODE
        numArgs = comLens(NODEi,i);
        args = comArgs{NODEi}{i}; %  numArgs, extVarNames, allVarNames, 'an external', 'all', 'INPUT', true);
        subjectName = args{1};
        [p subjIx] = ismember( subjectName, extVarNames );
        if ~p
            disp(['ReadConfigFile:  Invalid external variable name ' subjectName ' given in NODE specification'])
            return
        else
            if seenIx(subjIx) == 1 % then seen before
                disp(['ReadConfigFile:  More than one occurrence of ' subjectName ' for NODE specification'])
                return
            else
                seenIx(subjIx) = 1;
            end
        end
        if numArgs ~= 4
            disp( 'ReadConfigFile:  Wrong number of arguments to NODE specification')
            disp(['                   for subject variable ' subjectName])
            return
        end
        if ismember(subjectName, extVarNames)
            subjectIx = VarnameIxMap(subjectName, extVarNames);
            if isNum(args{2}) & isNum(args{3}) & isNum(args{4})
                objCoords = [str2num(args{2}), str2num(args{3}), str2num(args{4})]; % circle x, y, size
            else
                disp( 'ReadConfigFile: NODE arguments 2-4 contain non-numeric characters')
                disp(['                  for subject variable ' subjectName])
                return
            end
        else
            disp( 'ReadConfigFile:  First argument to NODE must be an external variable name')
            disp(['                   for subject variable ' subjectName])
            return
        end
        
        % NODELABEL
        numArgs = comLens(NODELABELi,i);
        args = comArgs{NODELABELi}{i};
        subjectName = args{1};
        [p pIx] = ismember( subjectName, extVarNames );
        if pIx ~= subjIx
            fprintf('ReadConfigFile:  NODELABEL spec %i`s subject expected to be %s\n',i,allVarNames{subjIx})
            return
        end
        if numArgs ~= 5
            disp( 'ReadConfigFile:  Wrong number of arguments to NODELABEL specification')
            disp(['                   for subject variable ' subjectName])
            return
        end
        if ismember(subjectName, extVarNames)
            subjectIx = VarnameIxMap(subjectName, extVarNames);
            labelStr = args{2};
            if isNum(args{3}) & isNum(args{4}) & isNum(args{5})
                labelCoords = [str2num(args{3}), str2num(args{4}), str2num(args{5})]; % circle x, y, size
            else
                disp( 'ReadConfigFile: NODELABEL arguments 3-5 contain non-numeric characters')
                disp(['                  for subject variable ' subjectName])
                return
            end
        else
            disp( 'ReadConfigFile:  First argument to NODELABEL must be an external variable name')
            disp(['                   for subject variable ' subjectName])
            return
        end
        
        % NODEADDSTATES
        numArgs = comLens(NODEADDSTATESi,i);
        args = comArgs{NODEADDSTATESi}{i};
        subjectName = args{1};
        [p pIx] = ismember( subjectName, extVarNames );
        if pIx ~= subjIx
            fprintf('ReadConfigFile:  NODEADDSTATES spec %i`s subject expected to be %s\n',i,allVarNames{subjIx})
            return
        end
        if numArgs < 2
            disp( 'ReadConfigFile:  Wrong number of arguments to NODEADDSTATES specification')
            disp(['                   for subject variable ' subjectName])
            return
        end
        if ismember(subjectName, extVarNames)
            subjectIx = VarnameIxMap(subjectName, extVarNames);
            if isNum(args{2})
                numAddedStates = str2num(args{2});
            else
                disp( 'ReadConfigFile:  NODEADDSTATES argument 2 was non-numeric')
                disp(['                   for subject variable ' subjectName])
                return
            end
            if numAddedStates ~= numArgs-2
                disp( 'ReadConfigFile:  Wrong number of arguments to NODEADDSTATES specification')
                disp(['                   for subject variable ' subjectName])
                return
            end
            stateLineTypes = cell( numArgs-2, 1);
            for statenum = 1:numArgs-2
                stateLineTypes{statenum} = args{2+statenum};
            end
        else
            disp( 'ReadConfigFile:  First argument to NODEADDSTATES must be an external variable name')
            disp(['                   for subject variable ' subjectName])
        end
        
        % NODESTATEMAP_ACTS
        numArgs = comLens(NODESTATEMAP_ACTSi,i);
        args = comArgs{NODESTATEMAP_ACTSi}{i};
        subjectName = args{1};
        [p pIx] = ismember( subjectName, extVarNames );
        if pIx ~= subjIx
            fprintf('ReadConfigFile:  NODESTATEMAP_ACTS spec %i`s subject expected to be %s\n',i,allVarNames{subjIx})
            return
        end
        if ismember(subjectName, extVarNames)
            subjectIx = VarnameIxMap(subjectName, extVarNames);
            tempActs = candActsIx{subjectIx};
            numIntActs = 0;
            for a = 1:length(tempActs)
                if tempActs(a) > numExt % then is an internal variable
                    numIntActs = numIntActs + 1;
                end
            end
            if numArgs ~= numIntActs + 1
                disp( 'ReadConfigFile:  Wrong number of arguments to NODESTATEMAP_ACTS specification')
                disp(['                   for subject variable ' subjectName])
                return
            end
            stateMapActs = {};
            if numIntActs == 0
                stateMapActs = { emptyStateMap };
            else
                for p = 1:numIntActs
                    stateMapActs = {stateMapActs{:}, eval(args{1+p})};
                end
            end
        else
            disp( 'ReadConfigFile:  First argument to NODESTATEMAP_ACTS must be an external variable name')
            disp(['                   for subject variable ' subjectName])
        end
        
        % NODESTATEMAP_POTS
        numArgs = comLens(NODESTATEMAP_POTSi,i);
        args = comArgs{NODESTATEMAP_POTSi}{i};
        subjectName = args{1};
        [p pIx] = ismember( subjectName, extVarNames );
        if pIx ~= subjIx
            fprintf('ReadConfigFile:  NODESTATEMAP_POTS spec %i`s subject expected to be %s\n',i,allVarNames{subjIx})
            return
        end
        if ismember(subjectName, extVarNames)
            subjectIx = VarnameIxMap(subjectName, extVarNames);
            tempPots = candPotsIx{subjectIx};
            numIntPots = 0;
            for a = 1:length(tempPots)
                if tempPots(a) > numExt % then is an internal variable
                    numIntPots = numIntPots + 1;
                end
            end
            if numArgs ~= numIntPots + 1
                disp( 'ReadConfigFile:  Wrong number of arguments to NODESTATEMAP_POTS specification')
                disp(['                   for subject variable ' subjectName])
                return
            end
            stateMapPots = {};
            if numIntPots == 0
                stateMapPots = { emptyStateMap };
            else
                for p = 1:numIntPots
                    stateMapPots = {stateMapPots{:}, eval(args{1+p})};
                end
            end
        else
            disp( 'ReadConfigFile:  First argument to NODESTATEMAP_POTS must be an external variable name')
            disp(['                   for subject variable ' subjectName])
        end
        
        diagObjects{i} = InitNetworkObj('Node', subjectName, objCoords, stateMapActs, stateMapPots, ...
                                        stateLineTypes, labelStr, labelCoords );
    end

    % Links ...
    if totNumObsInputs ~= numLinks
        disp('ReadConfigFile:  Incorrect number of LINKS specified. Require one for every interaction')
        disp('                   specified in an INPUTS command')
        return
    end
    for i = 1:numLinks
        % LINK
        numArgs = comLens(LINKi,i);
        args = comArgs{LINKi}{i};
        if ~( numArgs == 6 | numArgs == 7 )
            disp( 'ReadConfigFile:  Wrong number of arguments to LINK specification')
            disp(['                   number ' num2str(i) ])
            return
        end
        if numArgs == 7
            argTemp = args{7};
            if isNum(argTemp)
                if ~(str2num(argTemp) == 0 | str2num(argTemp) == 1)
                    disp( 'ReadConfigFile:  7th argument to LINK specification must be 0 or 1')
                    disp(['                   (LINK specification # ' num2str(i) ')' ])
                    return
                end
            else
                disp( 'ReadConfigFile:  7th argument to LINK specification must be 0 or 1')
                disp(['                   (LINK specification # ' num2str(i) ')' ])
                return
            end
        end
        subjectName = args{1};
        [p subjectIx] = ismember( subjectName, extVarNames );
        if p
            objectName = args{2};
            if ismember(objectName, extVarNames)
                objectIx = VarnameIxMap(objectName, extVarNames);
            else
                disp( 'ReadConfigFile:  Second argument to LINK must be a valid external variable name')
                disp(['                   (LINK specification # ' num2str(i) ')' ])
                return
            end
        else
            disp( 'ReadConfigFile:  First argument to LINK must be an external variable name')
            disp(['                   (LINK specification # ' num2str(i) ')' ])
            return
        end
        if isNum(args{3}) & isNum(args{4}) & isNum(args{5}) & isNum(args{6})
            linkCoords = [str2num(args{3}), str2num(args{4}), str2num(args{5}), str2num(args{6})];
        else
            disp( 'ReadConfigFile:  LINK coords contained non-numerics')
            disp(['                   (LINK specification # ' num2str(i) ')' ])
            return
        end
        % argument order to InitNetworkObj: typeStr, dependentVarStr, objCoords, stateMapA, stateMapP, ...
        %                                   stateLineTypes, labelStr, labelCoords, linkObjVarStr
        if numArgs == 6
            diagObjects{i+numNodes} = InitNetworkObj('Link', subjectName, linkCoords, {}, {}, {}, '', [], objectName, false);
        else % numArgs == 7
            diagObjects{i+numNodes} = InitNetworkObj('Link', subjectName, linkCoords, {}, {}, {}, '', [], objectName, str2num(args{7}));
        end
    end

    % VBAR SETUP
    vBarObjects = cell(numVbars,1); % numVbars == numTot has been verified already
    seenIx = zeros(1,numVbars);
    for i = 1:numVbars
        args = comArgs{VBARi}{i};
        numArgs = comLens(VBARi,i);
        if numArgs ~= 7
            disp( 'ReadConfigFile:  Wrong number of arguments to VBAR specification')
            disp(['                   (VBAR specification # ' num2str(i) ')' ])
            return
        end
        [p posIx] = ismember( args{1}, allVarNames );
        if p % then valid variable name
            if seenIx(posIx) == 0
                seenIx(posIx) = 1;
                if isNum(args{2}) & isNum(args{3}) & isNum(args{4}) & isNum(args{5}) & isNum(args{6}) & isNum(args{7})
                    vBarObjects{posIx} = InitVBar(args{1}, str2num(args{2}), str2num(args{3}), str2num(args{4}), ...
                        str2num(args{5}), str2num(args{6}), str2num(args{7}));
                else
                    disp( 'ReadConfigFile:  VBAR parameters 2-7 contained non-numerics')
                    disp(['                   (VBAR specification # ' num2str(i) ')' ])
                    return
                end
            else
                fprintf('ReadConfigFile:  Repeated VBAR specification for variable %s\n',args{1})
                return
            end
        else
            fprintf('ReadConfigFile:  Invalid variable name %s given to VBAR specification\n',args{1})
            return
        end
    end

else

    fprintf('ReadConfigFile:  Fatal Error. DSSRT specifications missing in file %s, or syntax messed up:\n',filename)
    fprintf('                 Incomplete specifications: ')
    for c = 1:numCommands
        if setup_done(c) == 0
            fprintf('%s ', commands{c})
        end
    end
    fprintf(['\n                 Commands must start "' command_str '" on a new line in the CFG file \n'])
    return

end

% all indices in returned values named "*Ix" are ABSOLUTE with respect to
% the list allVarNames
setup = {numInt; numExt; intVarNames; extVarNames; allVarNames; inputsIx;
        candActsIx; candPotsIx; diagObjects; vBarObjects; DEqns; DEpars; varBounds};
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function result = IsCfgCommand(line, commands, commandStr)
numCommands = length(commands);
args = '';
comStrLen = length(commandStr);
% This routine DOES NOT ignore leading whitespace: specification must start with
% command_str on a new line
len = length(line);
test = zeros(1, numCommands);
for c = 1:numCommands
    clens(c) = length(commands{c});
end

result = {test, args}; % default initial value signifies 'not a specification'
if len > comStrLen
    if comStrLen > 0
        first = line(1:comStrLen);
    else
        first = '';
    end
    rest = line(comStrLen+1:len);
    lenrest = len - comStrLen;

    %if sum( commandStr ~= first ) > 0 % then strings not equal
    if ~strcmp( commandStr, first )
        return
    else % line is a potential specification for us
        % test each specification to see if line corresponds to it
        for c = 1:numCommands
            if lenrest >= clens(c)
                if lenrest == clens(c) % then specification has been given no arguments
                    test(c) = (  sum( rest(1:clens(c)) ~= commands{c} ) == 0  ); % 1 if specification recognised
                else % have to check there's a space after the specification before its arguments
                    test(c) = (  sum( rest(1:clens(c)+1) ~= [commands{c} ' '] ) == 0  ); % 1 if specification recognised
                end
                if test(c) & clens(c) < lenrest
                    args = rest(clens(c)+2:lenrest);
                end
            end
        end
    end
else
    return
end

if sum(test) > 1
    % should only be one successful match in the list!!!
    disp('IsCFGcommand:  fatal error!')
    return
end
result = {test, args};
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function result = ProcessCfgArgs(argStrIn)
% ignores leading whitespace in string
argStrLen = length(argStrIn)+1;
argStr = [argStrIn ' ']; % ensure whitespace at end to properly update args in loop below
result = {};
arg = []; % placeholder for each argument's formation
args = {}; % final argument cell array
nargs = 0;
newSpace = false;

for pos = 1:argStrLen
    if ~isspace(argStr(pos))
        if ~newSpace % then have just started a new argument
            newSpace = true;
            nargs = nargs + 1;
        end
        if ischar(argStr(pos))
            arg = [arg argStr(pos)];
        else
            fprintf('ProcessCfgArgs:  Non-character found in argument %i\n',nargs)
            return
        end
    else
        if newSpace % then have finished with an argument, so update list
            newSpace = false;
            args = {args{:}, arg};
            arg = [];
        end
    end
end
result = {nargs, args};
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function index = VarnameIxMap(name, varnames)
% assumes varnames is a cell array, containing no repetitions!
[present index] = ismember(name,varnames);
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function result = FetchArgsSubjVars(Args, numArgs, VarNames1, VarNames2, VarTitle1, VarTitle2, comName, sorted);
% for commands of the form COMMAND <subject> <var> [<var> ... ]
result = {};
badArgs = false;
if numArgs > 0
    subjectName = Args{1};
    if ismember(subjectName, VarNames1)
        subjectIx   = VarnameIxMap(subjectName, VarNames1);
        if numArgs > 1
            unsortedIxs = zeros(1,numArgs-1);
            for a = 1:numArgs-1 % remaining arguments to the specification
                candidateName = Args{1+a};
                if ismember(candidateName, VarNames2)
                    unsortedIxs(a) = VarnameIxMap(candidateName, VarNames2);
                else
                    fprintf('ReadConfigFile:  Argument to %s with subject %s not in list of %s variables\n', comName, subjectName,VarTitle2)
                    badArgs = true;
                    break % for
                end
            end
            if ~badArgs
                if sorted
                    result = {subjectIx, sort(unsortedIxs)};
                else
                    result = {subjectIx, unsortedIxs};
                end
            end
        else
            result = {subjectIx, []};
        end
    else
        fprintf('ReadConfigFile:  First argument to %s must be %s variable name\n',comName,VarTitle1)
    end
else
    fprintf('ReadConfigFile:  At least one argument needed for %s specification\n',comName)
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Function currently defunct
function result = FetchArgsAllVars(Args, numArgs, VarNames, VarTitle, comName, sorted);
% for specifications of the form DOMGRAPH <var> [<var> ... ]
result = {};
badArgs = false;

if numArgs ~= 0
    unsortedIxs = zeros(1,numArgs-1);
    for a = 1:numArgs % remaining arguments to the specification
        candidateName = Args{a};
        if ismember(candidateName, VarNames)
            unsortedIxs(a) = VarnameIxMap(candidateName, VarNames);
        else
            fprintf('ReadConfigFile:  Argument to %s not in list of %s variables\n', comName,VarTitle)
            badArgs = true;
            break % for
        end
    end
    if ~badArgs
        if sorted
            result = sort(unsortedIxs);
        else
            result = unsortedIxs;
        end
    end
else
    fprintf('ReadConfigFile:  At least one argument needed for %s specification\n',comName)
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function objdata = InitNetworkObj(typeStr, dependentVarStr, objCoords, stateMapA, stateMapP, ...
                                   stateLineTypes, labelStr, labelCoords, linkObjVarStr, linkOff)
global NODE LINK OBJ_STATE OBJ_HANDLE OBJ_TYPE OBJ_DEPVAR OBJ_ACTSMAP OBJ_POTSMAP OBJ_STATELT OBJ_LABEL
global OBJ_OCOORD OBJ_LCOORD OBJ_OBJVAR LINKOFFSTATE

OK             = 0;
NO_DEP_VAR     = 1;
NO_RECOG_TYPE  = 2;
WRONG_NUM_ARGS = 3;
MAPS_NEEDED    = 4;
errStatus = OK;

objdata = cell(1,11);
lenTypeStr = length(typeStr);
if lenTypeStr == 4
    if sum( typeStr ~= 'Node' ) == 0 % then strings equivalent
        type = NODE;
    elseif sum( typeStr ~= 'Link' ) == 0
        type = LINK;
    else
        type = 0;
    end
else
    errStatus = NO_RECOG_TYPE;
end

if nargin <= 3 & ~errStatus
    if type == LINK
        stateMapA = linkStateMapActs;
        stateMapP = linkStateMapPots;
        labelStr = '';
        labelCoords = [0 0 0];
    else
        errStatus = MAPS_NEEDED;
    end
elseif nargin <= 4
    if type == LINK
        stateMapP = linkStateMapPots;
        labelStr = '';
        labelCoords = [0 0 0];
    else
        errStatus = MAPS_NEEDED;
    end
end

if nargin > 8 & ~errStatus & type == LINK
    if isempty(stateMapA)
        stateMapA = linkStateMapActs;
    end
    if isempty(stateMapP)
        stateMapP = linkStateMapPots;
    end
    if isempty(labelCoords)
        labelCoords = [0 0 0];
    end
    if nargin == 9
        linkOff = false;
    end % else linkOff was set in the argument list
end

if ~errStatus
    if nargin == 2
        errStatus = NO_DEP_VAR;
    elseif nargin == 1 | nargin > 10
        errStatus = WRONG_NUM_ARGS;
    end
    if nargin > 8 & type == NODE
        disp('InitNetworkObj:  Must not pass link object variable in arguments when object type is NODE.')
        errstatus = WRONG_NUM_ARGS;
    end
    if nargin == 8
        if type == LINK
            disp('InitNetworkObj:  Not enough arguments for object type is LINK.')
            errstatus = WRONG_NUM_ARGS;
        else
            linkObjVarStr = '';
            linkOff = false;
        end
    end
    if nargin == 9 & type == LINK
        linkOff = false;
    end
end

switch errStatus
    case NO_DEP_VAR
        disp('InitNetworkObj:  Nodes must be supplied with name of dependent variable')
    case NO_RECOG_TYPE
        disp('InitNetworkObj:  Supplied object type not recognised -- must be "Link" or "Node"')
    case WRONG_NUM_ARGS
        disp('InitNetworkObj:  Wrong number of arguments passed to InitNetworkObj')
    case MAPS_NEEDED
        disp('InitNetworkObj:  Explicit StateMaps must be passed to InitNetworkObj for Link objects')
    case OK
        if type == LINK
            if ~linkOff % normal -- display the link!
        		objdata{OBJ_STATE}   = 0; % initial current state for the displayed graphics object associated with this object
            else
                objdata{OBJ_STATE}   = LINKOFFSTATE; % represents off state
            end
        else
            objdata{OBJ_STATE}   = 0; % initial current state for the displayed graphics object associated with this object
        end
        objdata{OBJ_HANDLE}  = 0; % initial (dummy) object handle for the circle or arrow graphic
        objdata{OBJ_TYPE}    = type;
        objdata{OBJ_DEPVAR}  = dependentVarStr;
        objdata{OBJ_ACTSMAP} = stateMapA;
        objdata{OBJ_POTSMAP} = stateMapP;
        objdata{OBJ_STATELT} = stateLineTypes;
        objdata{OBJ_LABEL}   = labelStr;
        objdata{OBJ_OCOORD}  = objCoords;
        objdata{OBJ_LCOORD}  = labelCoords;
        objdata{OBJ_OBJVAR}  = linkObjVarStr;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function vbdata = InitVBar(name, xb, yb, yh, vlo, vhi, can_log);

vbdata = cell(10,1);

switch nargin
    case 7
        % do nothing, okay!
    case 6
        can_log = false;
    otherwise
        disp('InitVBar:  Wrong number of arguments passed')
        return
end

if can_log & (vlo ~= 0 | vhi ~= 1)
    disp('InitVBar:  Only gating variables belonging to (0,1) can have log scale')
    disp('            option turned on. Setting option to off.')
    can_log = false;
end

vbdata{1} = 0;
vbdata{2} = 0;
vbdata{3} = xb; % bar x value
vbdata{4} = yb; % bar y value (of bottom)
vbdata{5} = vlo; % lowest value shown
vbdata{6} = vhi-vlo; % greatest value shown as relative offset from lowest
vbdata{7} = yh; % bar height (y coord)
vbdata{8} = name; % string!
vbdata{9} = can_log; % allowed to be converted to log scale (for 0..1 gating vars)
   % this val can be 2 for scaling the top end of the range rather than the
   % bottom
vbdata{10} = 0; % initial (dummy) object handle for the dot used to signify 'on log scaling'

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function result = isNum(data)
% parameter `data` is a string
result = false;
pointCount = 0;
num = [char(48:57),'.','-'];
for i=1:length(data)
    if data(i) == '.'
        pointCount = pointCount + 1;
        if pointCount > 1
            return % false
        end
    end
    if ~ismember(data(i),num)
        return % false
    end
end
result = true; % can only get here if all chars were numeric
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = isDEpar(name, DEpars)
result = false;
for i=1:length(DEpars)
    if strcmp(name,DEpars{i}{1})
        result = true;
        break
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
-1.000000e+000
