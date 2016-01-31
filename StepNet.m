function StepNet(StepNetOptions)
% StepNet is the Graphical User Interface for DSSRT: integrating the Functional
%  Network, Transition Sequence, Regime Determination, and Attractor Estimation
%  utilities. These are post-processing tools for dimension-reduction analysis
%  of numerically integrated ODEs, using the dominant-scale method.
%
% Version 2.0, November 2005
% (c) Robert Clewley, Center for Biodynamics, Boston University, and
%                     Department of Mathematics, Cornell University.
%
% For help using this program see `DSSRT_documentation.rtf` and the theory
%   paper (in PDF format), both in the Documentation folder.
% Some bugs and general inconveniences are still to be expected in this developmental
%   version, so *please* send error reports and comments to rhc28@cornell.edu
%
% Geometric drawing and text functions by F. Auger, March 1999, with minor additions
%   and modifications by R. Clewley. These also include contributions by
%   Fred M. Staudaher.  Contact: f.auger@ieee.org, freds@packet.net

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEVELOPMENT NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 2.0, November 2005
%   Added option to plot sequence of phase planes after regime determination
%   Minor bug fixes (esp. 'speed fussy' strict option was not accessible)
%   Updates to support calculation of dominance strengths for systems with
%     electrical synapses.
%   Added switching of functional network diagram between that for
%     influence strengths and that for term sizes, along with associated
%     changes in functionality that use the dominant scale "E sequence"
%     information, etc.
%   Added option to restrict forward and backward "identical diagram state"
%     searches to a focused variable.
%   Added option to restrict comparison of transition sequences to a
%     focused variable.
%
% Version 1.84, September 2004
%   For consistency with theory papers, showing dScaleThresh and tScaleThresh in reciprocal form in read-outs
%   Added dialog boxes for user selection of DoV or Att Est (noShooting param) for domain of validity /
%     attractor estimate (command 'A'), and internal CompareActs() function's position thresholds
%     (although the latter are currently left commented-out as a development option)
%
% Version 1.83, July 2004
%   Added scale separation, and error calculation for reduced models to command '`'
%   Minor adjustments to messages, etc.
%   Removed automatic inclusion of all passive inputs in calculation of local model's target value,
%     in command '~'.
%   Added development option to compute functional network using term sizes (e.g. currents) rather than
%     Psi values. See global option useTermsFlag.
%   Known bug: When functional net re-calculated, most recently calculated trans. seq. still thinks it's
%     "fresh" and command 'X' suggests its data is still being displayed in the main window, which it is not.
%
% Version 1.82, May 2004
%   Made time scale / speed plot take absolute value of SumGamma1(1).
%   Added relative end position to transition sequence cell array.
%   Added extra parameters for RegimeDet, and updated DisplayRegimes accordingly.
%   Added a user confirmation stage before recomputing a functional network.
%   Added argument to DSSRT call to pass to StepNet, for use by parameters file (overriding `defaults`)
%   Removed ignorePassive option from command 'C' (create TS) -- seems unnecessary with RegimeDet now.
%   Added error line info when parameter file is bad.
%   Fixed minor error when deleting TS's in the presence of a freshTS that's retained
%     (updates its TS ref #)
%   Added COMPILE_ONLY option in DSSRTuseroptions and added underscores to all par names.
%
% Version 1.80, April 2004
%   System now known as DSSRT: Dominant-Scale System Reduction Tool. Other renamings include
%      `.fno` -> `.cfg` and `.fnp` -> `.par`
%   `sigRes` parameter has been replaced by `tScaleThresh`, and `scaleThresh` is now `dScaleThresh`.
%   Command '3' now changes time-scale threshold for AttEst and regimes (also saved in .par).
%   Command '`' now also shows time-scale information in console window.
%   Command '!' plots timescales or speeds over a marked interval in a separate figure.
%   Command '~' now plots, over a marked interval:
%     1) A focused variable & its asymptotic 'target' value (original and reduced for epochs)
%     2) All influence strength and additive input term values for that focused variable.
%   'About' message moved to command '*'.
%   Fixed a minor bug concerning attractor estimates when epoch markers are exactly at the ends of the
%     marked time interval.
%   Added general colour code help info for functional network diagram (e.g. for internal states),
%     specific to any ODE system that requires it, by placing `DSSRT_colourcodehelp.txt` in ODE folder.
%   New format for passed parameters to StepNet: StepNetOptions (a structure).
%   Command 'I' (change .cfg file) has been removed.
%   Command 'I' will redisplay reduced-regime information and interactively plot phase diagrams.
%   Command 'R' now finds reduced dynamical regimes. 'A' is now for AttEst, and '7' is now for toggling
%     `ignore pots/acts order` in transition sequences, etc.
%   Transition sequences can now be created that ignore passive variables (i.e. those without
%     associated differential equations), and addition TS_PASS entry created. Focused passives are now
%     screened out in command 'C'.
%   Removed pointless leading spaces from all console output :-)
%   GetNextFilename() now used by `save FN window as figure` command '$'. EndStripStr() updated (and
%     backwards compatible), and EndRemainStr() function added.
%   Updated FindIdentState(), GenTransSeq(), CompareTransSeq() to use verboseTog.
%
% Version 1.75, March 2004
%   Command 'G' now displays information about potential dominance values.
%   Command '`' added to display variable information at time t.
%   Included title in AttFig based on main window name.
%   Added trap for division by zero in all m_infinity, etc. functions for Hodgkin-Huxley activation eqns
%
% Version 1.74, February 2004
%   Minor changes.
%
% Version 1.73, January 2004
%   Changed focused transition sequences to include all obsIx-dependent variables too (see AttEst.m)
%   It is now obligatory to use a focused transition sequence before calling AttEst.
%   A new 'bounds' specification has been added to the .fno file. This is currently used to
%    estimate derivatives for bolster epochs, but will also be used in the automatic generation of
%    `potential` influence strength values in later versions.
%   A new command to refresh the main StepNet window has been added ('0', the character for zero),
%    in case of occasional Matlab figure mess-ups (although this doesn't yet fix the most serious ones).
%   Attractor estimates can be re-displayed using '&'.
%   Updated to work with new FuncNet.m v 2.0
%
% Version 1.72, December 2003
%   Changed order of internals and externals in .fno file, so that XPP auxiliaries can be used to add
%    non-differential equation generated variables at the end of the internals.
%   Added view of freshest transition sequence in main window, stepped through using SPACEBAR.
%   Command 'I' has become defunct, and is buggy anyway (see below). It will be removed in the future.
%   Bug: if .fno and variable data .dat files are renamed, reloading _FN.dat files referring to
%    the old names means that changing .dat and .fno files using 'I' and 'O' commands will not be
%    possible. 'O' will fail but StepNet recovers, 'I' will crash the program.
%   Removed "@FN " from the beginning of all .fno specifications
%   For AttEst, added focus option to transition sequences to allow less cluttered, more relevant
%     set of epochs for focused variable. This currently focuses on only one variable.
%
% Version 1.6, November 2003
%   Changed order of calling this function. StepNetInit.m has been replaced by FNSinit.m is now run
%    from FNS.m script, as is StepNet. This keeps the user-modified code parts more distinct from
%    internal routines.
%   Incorporated AttFig.m, GenTransSeq.m, CompareTransSeq.m, and ViewTransSeq.m into this file
%   Changed FuncNet.m so that it no longer returns domGraphIx entry
%   inputsIx map now has range of all system variables
%   Double buffering of figure animation added, for smoother graphics
%   AttEst now works for "non-voltage" variables with arbitrary number of inputs to the equation of study
%    (although it doesn't yet work correctly for such variables)
%   VarnameIxMap() made efficient
%   Selection of variables for variable of focus and variable to display are now done in terms of names
%   No bug fixes yet
%   All references to sigRes and sigResSet (for scale threshold continuation) are currently defunct
%   Compatability issue: older versions of MATLAB do not support the predefined boolean constants
%   Bug: Command 'P' says .fno file already set and will not update, but FNS seems to reload it anyway!
%
% Version 1.5, October 2003
%   Added automatic local attractor estimation for a single (voltage) variable -- assumes 7 currents
%     (this is implemented in AttEst.m)
%   Added save of observable of focus (obsIx) and variable to view (viewVarIx) to be saved with _FN file
%   No bug fixes yet
%
% Version 1.4, May 2003
%   Fixed a bug when an observable (node) has no outputs to another node
%   `Potential` arrow plots (in green) are commented out on lines #2850, 3543
%   Bug: when 'P' command changes param file, time bar is not correctly restored
%    and subsequent movement in time is not displayed with the vertical marker (nor are the markers)
%
% Version 1.3
%   Came and went. Very quickly.
%
% Version 1.2
%    Relative time display added.
%    Transition sequences now store the relative time of the sequence positions
%     (which are shown when using ViewTransSeq utility)
%
% Version 1.1
%    This version incorporates fully general network specification in the FNO file
%
% Version 1.0, January 2003
%    domGraphIx is defunct!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CONSTANTS (global index names, etc.)
GlobalConstants

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                       INITIAL SETUP                          %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin ~= 1
    disp('StepNet:  Incorrect number of arguments passed (1 expected). Quitting ...')
    return
end

if isempty(StepNetOptions.rootPath) || isempty(StepNetOptions.networkObjectDir)
    disp('StepNet:  Directory information for network objects improperly passed')
    return
end
networkObjectDir = StepNetOptions.networkObjectDir; % more convenient name
DSSRT_ROOTPATH = StepNetOptions.rootPath; % more convenient name
cd(DSSRT_ROOTPATH)
addpath(DSSRT_ROOTPATH)
addpath([DSSRT_ROOTPATH networkObjectDir])

% set default verbosity level
verboseTog = StepNetOptions.startVerbose;
if verboseTog, verbT_str = 'ON';  else  verbT_str = 'OFF';  end
disp(['StepNet:  Verbose messaging is ' verbT_str ]);

% temporary initializations
CFGfilename  = '';
CFGfilenameDisp = '';
varDataFilename = '';
varDataFilenameDisp = '';
CFGfilenameSet = false;
varDataFilenameSet = false;
availStTime = 0;
availSpTime = 0;
attEstParsSet = false;

% Attractor estimate parameter defaults (for neural equations)
attEstDefaults.delta  = 0.02; % currently unused
attEstDefaults.dVmax  = 145;
attEstDefaults.Vinres = 0.1; % step size for nhd. est. search
attEstDefaults.Vpert  = 25; % initial pert for chosen first epoch
attEstDefaults.searchRes = 2; % search resolution (also lowest bolster epoch interval size (highest resolution))
attEstDefaults.lowResMultiple = 30; % this value * searchRes
attEstDefaults.derivThresh = 0.3; % threshold for determining resolution of bolster epochs

% Use default parameters file, if present
if ~isempty(StepNetOptions.parFile)
    parFile = StepNetOptions.parFile;
else
    parFile = 'defaults';
end
if exist([ networkObjectDir '/' parFile '.par' ], 'file') == 2 % then file of params exists
    fprintf('StepNet:  Parameters file `%s.par` found. Loading ...\n',parFile)
    pars = ReadParamsFile( [ networkObjectDir '/' parFile '.par' ] );
    if ~isempty(pars)
        if verboseTog
            disp('StepNet:  Parameters loaded successfully')
        end
        CFGfilename        = [ networkObjectDir '/' pars{1} ];
        CFGfilenameDisp    = pars{1};
        CFGfilenameSet     = true;
        nOut               = pars{6};
        nOutSet            = true;
        varDataFilename    = [ networkObjectDir '/' pars{2} ];
        resultVAT = GetVarsAndTimes([DSSRT_ROOTPATH networkObjectDir],varDataFilename,CFGfilename,CFGfilenameSet,verboseTog);
        if isempty(resultVAT)
            disp('StepNet:  Fatal error. GetVarsAndTimes returned nothing')
            return
        end
        availStTime = resultVAT{1};
        availSpTime = resultVAT{2};
        numTot      = resultVAT{3};
        varnames    = resultVAT{4};
        varsAll     = resultVAT{5};
        seqTimesAll = resultVAT{6};
        vars        = varsAll(1:nOut:length(seqTimesAll),:);
        seqTimes    = seqTimesAll(1:nOut:length(seqTimesAll));
        numInt      = resultVAT{7};
        numExt      = resultVAT{8};
        actsIxMap   = resultVAT{9};
        potsIxMap   = resultVAT{10};
        inputsIx    = resultVAT{11};
        hobj        = resultVAT{12};
        hvb         = resultVAT{13};
        DEqns       = resultVAT{14};
        DEpars      = resultVAT{15};
        varBounds   = resultVAT{16};
        if isempty(actsIxMap) % arbitrarily use this one from CFG file to test CFG file read success
            beep
            disp('StepNet:  Cannot continue. CFG file error') % error message already will have been displayed
            return
        end
        if availSpTime == 0 % then varDataFilename messed up
            beep
            disp('StepNet:  Internal error. Reading available stop time messed up from varDataFilename')
            varDataFilename = ''; % reset these
            availStTime = 0;
            startTimeSet = false;
            stopTimeSet  = false;
        else
            varDataFilenameDisp   = pars{2};
            varDataFilenameSet    = true;
            startTime          = pars{3};
            stopTime           = pars{4};
            if ~ ( startTime >= availStTime && startTime < min(stopTime,availSpTime) )
                startTime = availStTime;
            end
            startTimeSet       = true;
            if ~ ( stopTime <= availSpTime && stopTime > max(startTime, availStTime) )
                stopTime = availSpTime;
            end
            stopTimeSet        = true;
        end
        dScaleThresh        = pars{5};
        dScaleThreshSet     = true;
        ignoreActsPotsOrder = pars{7};
        guessPeriod         = pars{8};
        marginPC            = pars{9};
        tScaleThresh        = pars{10};
        tScaleThreshSet     = true;
        attEstParams        = pars{11};
        attEstParsSet       = true;
        parsInitSet         = true;
    else
        disp('StepNet:  Default parameters file messed up')
        startTimeSet    = false;
        stopTimeSet     = false;
        dScaleThreshSet = false;
        nOutSet         = false;
        tScaleThreshSet = false;
        parsInitSet     = false;
    end
else
    parsInitSet = false;
end

% Quit if no parameters file specified
if ~parsInitSet
    disp('StepNet:  Missing default parameters file')
    disp('StepNet:  Please specify a parameters file ... [NOT YET IMPLEMENTED]')
    disp('          ... quitting')
    return
end

% Other default initializations
varViewIx = 1;
TSmatchTog = 2;
obsIx    = 1; % default to first variable of whatever system
obsIxSet = true;
attractorEst.exist = false;
attractorEst.number = 0;
attractorEst.attEstList = {};
varChangeThresh_def = 0.10; % = 10%
derivThresh_def = 0.1;
PSIS = 0;
TERMS = 1;
diagSw = PSIS;  % which functional net to display

% Compile Differential equations
numEqns = length(DEqns);
% (produce index map)
allDEixMap = zeros(1,numTot); % (for all variables)                        
DEqns_comp = cell(1,numEqns);
for eqn=1:numEqns
    % create map entry
    [ispresent posfound] = ismember(DEqns{eqn}{DE_NAMEi},varnames); % posfound is an absolute index
    if ispresent
        allDEixMap(posfound) = eqn;
    end
    % compile this eqn
    DEresult = CompileGammas(DEqns{eqn}, DEpars, varnames);
    if isempty(DEresult)
        disp('Fatal error')
        return
    else
        DEqns_comp{eqn} = DEresult;
    end
end
DEqnsCompiled = true;
if verboseTog
    disp('StepNet:  Compiled DEs successfully')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%            Network object setup            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure setup
figTitle1 = 'StepNet >> ';
figTitleName = 'No Functional Net';
DSSRTfigHandle = figure('NumberTitle','off','Name',[ figTitle1 figTitleName],...
            'Position',[10, 40, 550, 700], 'MenuBar','none', 'KeyPressFcn', 'DoNothing');
viewFigHandle = -1; % default initial value for ViewTransSeq window
diagSw_handle = -1;

% Double buffering setup for smoother graphics
set(DSSRTfigHandle,'RendererMode','manual')
set(DSSRTfigHandle,'Renderer','painters')
set(DSSRTfigHandle,'DoubleBuffer','on')

if parsInitSet
	numobjects = length(hobj); % number of network objects (observables)
	if numobjects == 0
        disp('StepNet:   Network graphics failed to initialize. No objects to display!')
        return
	end
	
	%% Get windows set up and display default network state
	initFNresult = InitFuncNet(DSSRTfigHandle, false, hobj);
	if isempty(initFNresult)
        disp('StepNet:  Internal error: empty `initFNresult`. Cannot proceed ...')
        disp('           Perhaps check that a directory for network object methods is present and correct')
        return
	end
	
	handVVaxes = initFNresult{1}; % Variable viewer subplot
    handTSaxes = initFNresult{2}; % Transition Sequence subplot
	handTBaxes = initFNresult{3}; % Time box subplot
	handFNaxes = initFNresult{4}; % Functional net subplot
	tBarHand   = initFNresult{5}{1};
	tTickHand  = initFNresult{5}{2};
	
	numVBs = length(hvb);

	%%%% Set up initial stuff on screen...
	% display text for relevant objects, using second argument (initialize) == true
	axes(handFNaxes) % set to current axes
	for obj = 1:numobjects
        DrawNet(hobj{obj}, true); % discard returned (null) handle for initialization
        hobj{obj}{OBJ_HANDLE} = DrawNet( hobj{obj} ); % draw inactive state
	end

	% display initialised vertical bars
	for vb = 1:numVBs
        hvbc = hvb{vb};
        xb = hvbc{3};
        yb = hvbc{4};
        yh = hvbc{7};
        hvb{vb}{1} = plot([xb xb],[yb yb+yh],'b-');
        hvb{vb}{2} = plot([xb-0.01 xb+0.01], [yb yb],'k-','LineWidth',2);
	end
	
	% display first variable in VarViewBox, if CFGfilenameSet
	% (i.e. if data was loaded correctly)
	if CFGfilenameSet
        PlotVBox(handVVaxes,seqTimes,vars(:,varViewIx),varViewIx,varnames)
        axes(handFNaxes)
	end
    
%     test_inits = ones(1,numobjects); % initial values of test results for a state change
	diag_state_old = zeros(1, numobjects);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%        SETUP MAIN LOOP FOR USER INTERACTION          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% log scale variables toggle
log_scale = false; % default initial value
if log_scale,  log_str = 'ON';  else   log_str = 'OFF';  end
lmag = 15;
toggle_ls = true; % shows whether log_scale has just been updated (for graphics changes)
fprintf(['StepNet:  Log variable scaling ' log_str ' (magnification set to %i)\n'],lmag);

% whether to care or not about order of actives in a curly E set when
% shooting for periodic states
if ignoreActsPotsOrder,  iAO_str = 'ON';  else  iAO_str = 'OFF';  end
disp(['StepNet:  Ignore order of actives/potentials in network state is ' iAO_str ]);

switch TSmatchTog
    case 0
        TSm_str = 'PADDING and SWAPPING both OFF';
    case 1
        TSm_str = 'PADDING ON, SWAPPING OFF';
    case 2
        TSm_str = 'PADDING ON, SWAPPING ON';
end

% startup info message
if verboseTog
    disp('StepNet:  Put the pointer over the figure window in order to issue')
    disp('            keyboard commands and press H, * or _ for help')
end

% other initial settings, flags
keypress  = 'a'; % dummy initial keypress (impossible because real keypresses are first capitalized)
resetView = true;
existNet  = false;
clearFNet = false;
freshTS   = {false, 0, 0, 0, 0};
epochTimesSet = false;
epochTimes = [];
markLset   = false;
markRset   = false;
markLeft   = 0; % initial value out of range of legal position in diagram state position
markRight  = 0; % initial value out of range of legal position in diagram state position
markLchanged  = false;
markRchanged  = false;
existTransSeq = false;
TransSeqs     = {}; % initial setting only
refreshDiagState = false;
forceNetRedraw = false;
existFNpars   = CFGfilenameSet & varDataFilenameSet & nOutSet & startTimeSet & stopTimeSet & dScaleThreshSet;
existRegime   = false;

% FNdataFname    = '';
FNdataFilenameDisp = '';
FNdataFnameSet = false;
% TSdataFname    = '';
TSdataFilenameDisp = '';
TSdataFnameSet = false;

pos = 1; % initial value for time position
doQuit = false; % initial value
errStatus = false; % initial value
while ~doQuit && ~errStatus % give me two reasons not to quit right now
    switch keypress
        case 'H' % command help
            % currently-used keyboard command characters (AND all the alphabet is used up):
            % 012345678_.,<>+-=;:"[]{}\|/?`~!@#$%^&*
            % other keys used:   ESCAPE KEY, SPACE BAR
            % keys still available: BACKSPACE, ARROW KEYS, ENTER KEY, 9()'
            % (Codes: BSPC = 8, ENTER = 13) ... TAB is not recognized!
            fprintf('\n')
            disp('************************************************************************************')
			disp('StepNet:         Single key commands (active only over StepNet figure window)')
            disp('************************************************************************************')
			disp('Frame changes:   N(ext +1)   B(ack -1)  , (+10)  . (-10)')
            disp('                 < (go to start)  > (go to end)')
			disp('General:         ^ (toggle log-scaled vars)      $ (save screenshot)  Q(uit)  H(elp)')
            disp('                 _ (func. net colour-code help)  * (about DSSRT)')
            disp('                 T (go to time)                  Y (go to position)')
            disp('                 + (write parameter file)        = (show parameters specified)')
            disp('                 P (read parameter file)         ? (toggle diagnostics verbosity)')
            disp('                 0 [zero] (refresh main window)  ESC(ape from current computation)')
            disp('Functional net:  1 (change start time)           2 (change stop time)')
            disp('                 3 (change scale thresholds)     4 (change resampling-step)')
            disp('                 5 (change observable of focus)  6 (change variable to view)')
            disp('                 U (compute functional net)      @ (list current influence scales)')
            disp('                 W (write func. net to file)     O (specify raw data file name)')
            disp('                 E (read func. net from file)    - (switch between Psis and terms)')
            disp('                 % (show graph of dominances)    ` (list current variable info)')
            disp('                 ~ (show graph of variable targets and influence strengths / terms)')
            disp('                 ! (show graph of variable timescales or speeds)')
            disp('File system:     { (list var data files)         } (list CFG & PAR files)')
            disp('                 [ (list FN data files)          ] (list TS data files)')
            disp('Position marks:  J (set left mark)               K (set right mark)')
            disp('                 ; (jump to left mark)           : (jump to right mark)')
            disp('                 M (show marks set)              / (clear marks)')
            disp('Transition seqs: D (find next identical state of actives, backwards)')
            disp('                 F (find next identical state of actives, forwards)')
            disp('                 C(alculate transition sequence from left mark -> right mark)')
            disp('                 V(iew transition sequences)     X (list transition sequences)')
            disp('                 S(ave transition sequences)     L(oad transition sequences)')
            disp('                 SPACEBAR (cycle through epochs in `fresh` transition sequence)')
            disp('                 \ (delete 1 or all transition sequences in memory)')
            disp('                 " (rename and/or change x or y data values of a trans. sequence)')
            disp('                 Z (compare & pad transition sequences in memory)')
            disp('                 # (toggle: padding/swap algorithm level)')
            disp('                 7 (toggle: ignore order for actives/potentials set)')
            disp('                 8 (toggle: ignore potentials for transition sequence)')
            disp('                 | (plot transition sequence y data versus x data)')
            disp('Regimes:         R (find dynamical regimes)      I (show regime info w/ diagrams)')
            disp('Attractor / DoV: A (calculate local attractor / DoV for one variable)')
            disp('                 & (redisplay previously calculated local attractor / DoV)')
            disp('************************************************************************************')
            fprintf('\n')
        case '_' % colour code help
            if existNet
                colFilename = [ networkObjectDir '/' 'DSSRT_colourcodehelp.txt'];
                if exist(colFilename,'file') == 2
                    fid = fopen(colFilename,'r');
                    fileLines = {};
                    eofile = false;
                    if fid ~= -1
                        while ~eofile
                            fileLine = fgetl(fid);
                            if fileLine == -1
                                eofile = true;
                            else
                                fileLines = {fileLines{:}, fileLine};
                            end
                        end
                        proceed = true;
                    else
                        disp('StepNet:  Error opening `DSSRT_colourcodehelp.txt` in ODE system folder')
                        proceed = false;
                    end
                    if isempty(fileLines)
                        disp('StepNet:  Error with `DSSRT_colourcodehelp.txt`: no text contents!')
                        proceed = false;
                    end
                    if proceed
                        fprintf('\n')
                        disp('************************************************************************************')
                        disp('StepNet:  Help for colour codes in functional network figure')
                        for i = 1:length(fileLines)
                            fprintf([' ' fileLines{i} '\n'])
                        end
                        disp('************************************************************************************')
                        fprintf('\n')
                    end
                    fclose(fid)
                else
                    disp('StepNet:  No colour code help available')
                    if verboseTog
                        disp('           (no file `DSSRT_colourcodehelp.txt` in ODE system folder)')
                    end
                end
            else
                disp('StepNet:  No functional network calculated!')
            end
        case '*' % about FuncNet
            fprintf('\n')
            disp(' ********************************** About DSSRT ***************************************')
            help StepNet
            disp(' ************************************************************************************')
            fprintf('\n')
        case '0' % [zero] refresh main StepNet window (in case of axis mess-up)
			set(handVVaxes,'Position', [0.1 0.885 0.85 0.07]);
            axis(handVVaxes,'on')            
            set(handTSaxes,'Position', [0.1 0.84 0.85 0.01]);
            axis(handTSaxes,'off')
            set(handTBaxes,'Position', [0.1 0.8 0.85 0.04]);
			set(handFNaxes,'Position', [0.05 0.04 0.9 0.71] );
			set(handFNaxes,'PlotBoxAspectRatioMode','manual')
			set(handFNaxes,'XLimMode','manual')
			set(handFNaxes,'YLimMode','manual')
            axis(handFNaxes,'off')
            axes(handFNaxes)
            refreshDiagState = true;
            if verboseTog
                disp('StepNet:  Refreshed main StepNet window')
            end
        case 'N' % next step
            if existNet
                pos = pos + 1;
                refreshDiagState = true;
            else
                if verboseTog
                    disp('StepNet:  No functional network calculated!')
                end
            end
        case 'B' % previous step
            if existNet
                pos = pos - 1;
                refreshDiagState = true;
            else
                if verboseTog
                    disp('StepNet:  No functional network calculated!')
                end
            end
        case '.' % -10 steps
            if existNet
                pos = pos + 10;
                refreshDiagState = true;
            else
                if verboseTog
                    disp('StepNet:  No functional network calculated!')
                end
            end
        case ',' % +10 steps
            if existNet
                pos = pos - 10;
                refreshDiagState = true;
            else
                if verboseTog
                    disp('StepNet:  No functional network calculated!')
                end
            end
        case '}' % list .cfg and .par system setup files in local directory
            disp('StepNet:  CFG system setup files available to load using `I` command:')
            dir( [ networkObjectDir '/*.cfg' ] )
            disp('StepNet:  PAR parameter files available to load using `P` command:')
            dir( [ networkObjectDir '/*.par' ] )
        case '{' % list .dat (potential) variable data files in local directory
            disp('StepNet:  Variable data files available to load using `O` command:')
            dir( [ networkObjectDir '/*.dat' ] )
        case '[' % list _FN.mat functional net data files in local directory
            disp('StepNet:  Functional net data files available to load using `E` command:')
            dir( [ networkObjectDir '/*_FN.mat' ] )
        case ']' % list _TS.mat transition sequence files in local directory
            disp('StepNet:  Transition sequence data files available to load using `L` command:')
            dir( [ networkObjectDir '/*_TS.mat' ] )
        case '$' % save figure snapshot to disc
            if existNet
                date_string = strrep(datestr(now,1),'-','');
                root_string = strcat(networkObjectDir,'/Figures/',EndStripStr(CFGfilenameDisp,'.',true),'_', ...
                     EndStripStr(varDataFilenameDisp,'.',true));
                figFilename = GetNextFilename([root_string '-' date_string],'fig');
                saveas(DSSRTfigHandle, figFilename)
                if verboseTog
                    fprintf('StepNet:  Figure saved for step %i as `%s`\n',pos,RootStripStr(figFilename))
                end
            else
                disp('StepNet:  No functional network calculated!')
            end
        case '?' % toggle verbose messaging
            verboseTog = not(verboseTog);
            if verboseTog, verbT_str = 'ON';  else  verbT_str = 'OFF';  end
            disp(['StepNet:  Verbose messaging is ' verbT_str ]);
        case '^' % toggle log scaling for appropriate variables
            log_scale = not(log_scale);
            toggle_ls = true;
            if log_scale
                log_str = ['ON with magnification ' num2str(lmag)];
            else
                log_str = 'OFF';
            end
            fprintf(['StepNet:  Log variable scaling ' log_str '\n'],log_scale);
            if verboseTog
                disp('            (the dot marks the end of the slider that has greater magnification)')
            end
            refreshDiagState = true; % use this to force update of diagram to show dots for log-scaled variables
        case 'T' % go to time
            if existNet
                result = DialogBox('New time in sequence',seqTimes(pos),verboseTog);
                if ~isempty(result)
                    resultVal = str2num(result);
                    pos = FindPosFromTime(resultVal,seqTimes,numSteps);
                    fprintf('           New time set to %.3f\n',seqTimes(pos))
                    refreshDiagState = true;
                else
                    if verboseTog
                        disp('           Cancelled')
                    end
                end
            else
                disp('StepNet:  No functional network calculated!')
            end
        case 'Y' % go to position
            if existNet
                result = DialogBox('New position in sequence',pos,verboseTog);
                if ~isempty(result)
                    resultVal = str2num(result);
                    if round(resultVal) == resultVal
                        if resultVal > 0 && resultVal <= numSteps
                            pos = resultVal;
                            refreshDiagState = true;
                        else
                            fprintf('StepNet:  Position out of range [1,%i]\n',numSteps)
                        end
                    else
                        disp('StepNet:  Position must be an integer')
                    end
                else
                    if verboseTog
                        disp('           Cancelled')
                    end
                end
            else
                if verboseTog
                    disp('StepNet:  No functional network calculated!')
                end
            end
        case 'J' % set mark for left hand of sequence
            if existNet
                if ~markRset
                    markLeft = pos;
                    markLset = true;
                    markLchanged = true;
                    refreshDiagState = true;
                    if verboseTog
                        fprintf('StepNet:  Left mark set at position %i\n',markLeft)
                    end
                else
                    if pos >= markRight
                        fprintf('StepNet:  Left mark must be to the left of the right mark at position %i\n',markRight)
                    else
                        markLeft = pos;
                        markLset = true;
                        markLchanged = true;
                        refreshDiagState = true;
                        if verboseTog
                            fprintf('StepNet:  Left mark set at position %i\n',markLeft)
                        end
                    end
                end
            else
                if verboseTog
                    disp('StepNet:  No functional network calculated!')
                end
            end
        case 'K' % set mark for right
            if existNet
                if ~markLset
                    markRight = pos;
                    markRset = true;
                    markRchanged = true;
                    refreshDiagState = true;
                    if verboseTog
                        fprintf('StepNet:  Right mark set at position %i\n',markRight)
                    end
                else
                    if pos <= markLeft
                        fprintf('StepNet:  Right mark must be to the right of the left mark at position %i\n',markLeft)
                    else
                        markRight = pos;
                        markRset = true;
                        markRchanged = true;
                        refreshDiagState = true;
                        if verboseTog
                            fprintf('StepNet:  Right mark set at position %i\n',markRight)
                        end
                    end
                end
            else
                if verboseTog
                    disp('StepNet:  No functional network calculated!')
                end
            end
        case 'M' % show marks set
            if existNet
                if markLset && ~markRset
                    fprintf('StepNet:  Left mark is set at position %i\n' ,markLeft);
                end
                if markRset && ~markLset
                    fprintf('StepNet:  Right mark is set at position %i\n' ,markRight);
                end
                if markLset && markRset
                    fprintf('StepNet:  Left mark is set at position  %i\n' ,markLeft);
                    fprintf('          Right mark is set at position %i\n' ,markRight);
                end
                if ~markLset && ~markRset
                    disp('StepNet:  No marks set')
                end
            else
                if verboseTog
                    disp('StepNet:  No functional network calculated')
                end
            end
        case '/' % clear marks
            if existNet
                if markLset || markRset
                    if verboseTog
                        disp('StepNet:  Clear marks? (see dialog box)')
                    end
                    ButtonName = questdlg('Are you sure?','Clear marks?','Yes','No','No');
                    switch ButtonName
                        case 'Yes'
                            if markLset
                                markLchanged = true;
                            end
                            if markRset
                                markRchanged = true;
                            end
                            markLset = false;
                            markRset = false;
                            markLeft = 0;
                            markRight = 0;
                            refreshDiagState = true;
                            disp('StepNet:  Cleared marks')
                        otherwise
                            if verboseTog
                                disp('           Cancelled')
                            end
                    end
                end
            else
                if verboseTog
                    disp('StepNet:  No functional network, so no marks exist')
                end
            end
        case '&' % redisplay local attractor basin / domain of validity estimate
            if existNet
                if attractorEst.exist
                    disp('StepNet:  Choose an attractor estimate / domain of validity to re-display (see dialog box)')
                    disp('           Data sets in memory:')
                    for ae = 1 : attractorEst.number
                        attEstData = attractorEst.attEstList{ae};
                        fprintf( '             #%2i - focused on variable %s, from time %.4f\n', ae, varnames{attEstData.obsIx}, attEstData.t0)
                        fprintf( '                   over time period of length %.3f, with %i sample epochs\n', attEstData.Tperiod, attEstData.numEpochs)
                        if attEstData.noShooting && attEstData.version == 1
                            disp('                  - domain of validity data')
                        else
                            disp('                  - attractor estimate data')
                        end
                        if attEstData.version == 1 % then used regimes version (NEW to DSSRT v.1.12)
                            disp('                   (used regimes version of AttEst.m)')
                        else % == 0
                            disp('                   (used epochs version of AttEst.m)')
                        end
                    end
                    result = DialogBox(['data set # (from 1 -> ' num2str(attractorEst.number) ')'],attractorEst.number,verboseTog);
                    if ~isempty(result)
                        attIx = str2num(result);
                        if attIx > 0 && attIx <= attractorEst.number % then proceed
                            attEstData = attractorEst.attEstList{attIx};
                            fignum = AttFig(attEstData); % note, this currently over-writes any old fignum handle...
                            fprintf('StepNet:  Re-displaying data set #%i in graph',attIx)
                            if verboseTog
                                disp(['            for variable ' varnames{attEstData.obsIx} ' starting at time ' num2str(attEstData.t0,'%.3f')]);
                            else
                                fprintf('\n')
                            end
                            
                        end
                    end
                else
                    disp('StepNet:  No attractor estimate / domain of validity data sets in memory')
                end
            else
                if verboseTog
                    disp('StepNet:  No functional network calculated')
                end
            end
        case 'A' % calculate local attractor basin estimate / domain of validity
            if existNet
                if ~tScaleThreshSet
                    disp('StepNet:  Time scale threshold must be set before running AttEst')
                    proceed = false;
                else
                    proceed = true;
                end
                if proceed
                    ButtonName = questdlg('Use regimes (default) or epochs version?','Regimes or epochs version','Regimes','Epochs','Regimes');
                    switch ButtonName
                        case 'Regimes'
                            if existRegime
                                proceed = true;
                                if verboseTog
                                    disp('StepNet:  Using the regimes version ...')
                                end
                                attest_ver = 1;
                            else
                                proceed = false;
                                disp('StepNet:  No reduced regimes available - calculate them first')
                            end
                        case 'Epochs'
                            if existTransSeq
                                proceed = true;
                                if verboseTog
                                    disp('StepNet:  Using the epochs version ...')
                                end
                                attest_ver = 0;
                            else
                                proceed = false;
                                disp('StepNet:  No transition sequence available - calculate epochs first')
                            end
                        otherwise
                            proceed = false;
                    end
                end
                if proceed && attest_ver == 0
                    %%%%%%%%%%%%%
                    %% not sure it's necessary to have the marks set at the period -- these positions should just come from the TS
                    %%%%%%%%%%%%%%
                    if markLset && markRset
                        disp('StepNet:  Choose a transition sequence containing no padding')
                        disp('           on which to base attractor calculation (see dialog box)')
                        numTSeqs = length(TransSeqs);
                        proceed = false; % initial default value
                        result = DialogBox(['a trans seq # to base attractor calculation (from 1 -> ' num2str(numTSeqs) ')'],numTSeqs,verboseTog);
                        if ~isempty(result)
                            TSIx = str2num(result);
                            if TSIx > 0 && TSIx <= numTSeqs % then proceed
                                if verboseTog
                                    fprintf('StepNet:  Using transition sequence #%i\n',TSIx)
                                end
                                if ~isempty(TransSeqs{TSIx}{TS_FOCUS})
                                    proceed = true;
                                else
                                    disp('StepNet:  Transition sequence selected is not focused')
                                    disp('           (Use the `focus` option when creating transition sequence)')
                                end
                            else
                                disp('           Invalid sequence number entered')
                            end
                        end
                    else
                        proceed = false;
                        disp('StepNet:  Use left and right markers to demark time limits')
                    end
                end
                if proceed
                    if attest_ver == 0
                        obsIx = TransSeqs{TSIx}{TS_FOCUS}(1);
                    else
                        obsIx = regimeStruct.focVarIx;
                    end
                    focVarName = varnames{obsIx};
                    
                    % defaults and dVmax
                    proceed = false; % reset
                    % check obsIx bounds first, and create defaults
                    obsIxInt = varBounds(obsIx,2)-varBounds(obsIx,1);
                    
                    if attEstParsSet && attEstParams.dVmax >= obsIxInt % e.g. if loaded from a parameters file for a different obsIx
                        % (this wasn't checked there because the bound interval wasn't available)
                        attEstParsSet = false;
                        if verboseTog
                            disp('StepNet:  Previous dVmax value was outside focused variable`s bounds.')
                            disp('          It has been reset to the size of the bounds, and the')
                            disp('            other parameters adjusted according to defaults')
                        end
                    end
                    
                    if ~attEstParsSet
                        % re-adjust defaults to current system if not already determined these params
                        attEstDefaults.dVmax  = 0.99 * obsIxInt;
                        attEstDefaults.Vinres = 0.1/150 * obsIxInt;
                        attEstDefaults.Vpert  = 40/150 * obsIxInt;
                        attEstDefaults.posThresh1 = 2;
                        attEstDefaults.posThresh2 = 2;
                        attEstParams = attEstDefaults;
                    end

                    % determine noShooting (DoV or Att Est) if using regimes version
                    if attest_ver == 1
                        if verboseTog
                            disp('StepNet:  Calculate domain of validity or attractor estimate? (see dialog box)')
                        end
                        ButtonName = questdlg('Calculate Domain of Validity or Attractor Estimate?','DoV or AttEst ...','DoV','AttEst','DoV');
                        switch ButtonName
                            case 'DoV'
                                attEstParams.noShooting = true;
                                proceed = true;
                            case 'AttEst'
                                attEstParams.noShooting = false;
                                proceed = true;
                            otherwise
                                proceed = false;
                        end
                    else
                        attEstParams.noShooting = false;
                        proceed = true;
                    end
                end

                if proceed
                    if attest_ver == 0 || (attest_ver == 1 && ~attEstParams.noShooting) % then Att. Est. calc needs these params
                        proceed = false;
                        result = DialogBox(['maximum perturbation in ' focVarName],num2str(attEstParams.dVmax),verboseTog);
                        if ~isempty(result)
                            if isNum(result)
                                dVmaxNew = str2num(result);
                                if dVmaxNew > 0 && dVmaxNew < obsIxInt
                                    attEstParams.dVmax = dVmaxNew;
                                    proceed = true;
                                    if verboseTog
                                        disp(['StepNet:  Using max perturbation ' result])
                                    end
                                else
                                    if verboseTog
                                        disp('          Value must be positive and less than focused variable`s bounds')
                                    else
                                        disp('StepNet:  Value must be positive and less than focused variable`s bounds')
                                    end
                                end
                            end
                        end

                        if proceed % Vpert will become delta instead
                            proceed = false; % reset
                            result = DialogBox(['initial perturbation in ' focVarName],num2str(attEstParams.Vpert),verboseTog);
                            if ~isempty(result)
                                if isNum(result)
                                    VpertNew = str2num(result);
                                    if VpertNew > 0 && VpertNew < attEstParams.dVmax
                                        attEstParams.Vpert = VpertNew;
                                        proceed = true;
                                        if verboseTog
                                            disp(['StepNet:  Using initial perturbation ' result])
                                        end
                                    else
                                        if verboseTog
                                            disp('          Value must be less than dVmax and positive')
                                        else
                                            disp('StepNet:  Value must be less than dVmax and positive')
                                        end
                                    end
                                end
                            end
                        end
                    else
                        attEstParams.dVmax = obsIxInt;
                        attEstParams.Vpert = 1; % this value will be ignored in AttEst.m
                    end
                end

                if proceed % Vinres
                    proceed = false; % reset
                    result = DialogBox(['inner-layer resolution in ' focVarName],num2str(attEstParams.Vinres),verboseTog);
                    if ~isempty(result)
                        if isNum(result)
                            VinresNew = str2num(result);
                            if VinresNew > 0 && VinresNew < attEstParams.dVmax
                                attEstParams.Vinres = VinresNew;
                                proceed = true;
                                if verboseTog
                                    disp(['StepNet:  Using perturbation step ' result])
                                end
                            else
                                if verboseTog
                                    disp('          Value must be less than dVmax and positive')
                                else
                                    disp('StepNet:  Value must be less than dVmax and positive')
                                end
                            end
                        end
                    end
                end
                
                if proceed % searchRes
                    proceed = false; % reset
                    result = DialogBox('search resolution (in time steps)',num2str(attEstParams.searchRes),verboseTog);
                    if ~isempty(result)
                        if isNum(result)
                            sResNew = str2num(result);
                            if sResNew > 0 && sResNew < 500 && (round(sResNew) == sResNew)
                                attEstParams.searchRes = sResNew;
                                proceed = true;
                                if verboseTog
                                    disp(['StepNet:  Using search resolution of ' result])
                                end
                            else
                                if verboseTog
                                    disp('          Value must be a positive integer < 500')
                                else
                                    disp('StepNet:  Value must be a positive integer < 500')
                                end
                            end
                        end
                    end
                end
                
                if proceed % lowResMultiple
                    proceed = false; % reset
                    result = DialogBox('low resolution value (multiple of search resolution)',num2str(attEstParams.lowResMultiple),verboseTog);
                    if ~isempty(result)
                        if isNum(result)
                            lResNew = str2num(result);
                            if lResNew > 0 && lResNew < 100 && (round(lResNew) == lResNew)
                                attEstParams.lowResMultiple = lResNew;
                                proceed = true;
                                if verboseTog
                                    disp(['StepNet:  Using low resolution multiple of ' result])
                                end
                            else
                                if verboseTog
                                    disp('          Value must be a positive integer < 100')
                                else
                                    disp('StepNet:  Value must be a positive integer < 100')
                                end
                            end
                        end
                    end
                end
                
                if proceed % derivThresh
                    proceed = false; % reset
                    result = DialogBox('derivative threshold',num2str(attEstParams.derivThresh),verboseTog);
                    if ~isempty(result)
                        if isNum(result)
                            dThNew = str2num(result);
                            if dThNew > 0 && dThNew < 100
                                attEstParams.derivThresh = dThNew;
                                proceed = true;
                                if verboseTog
                                    disp(['StepNet:  Using derivative threshold of ' result])
                                end
                            else
                                if verboseTog
                                    disp('          Value must be a positive real < 100')
                                else
                                    disp('StepNet:  Value must be a positive real < 100')
                                end
                            end
                        end
                    end
                end

%                 if proceed % posThresh1
%                     proceed = false; % reset
%                     result = DialogBox('position threshold 1 (positive int)',num2str(attEstParams.posThresh1),verboseTog);
%                     if ~isempty(result)
%                         if isNum(result)
%                             pt1 = str2num(result);
%                             if pt1 > 0 && pt1 < 10 && (round(pt1) == pt1)
%                                 attEstParams.posThresh1 = pt1;
%                                 proceed = true;
%                                 if verboseTog
%                                     disp(['StepNet:  Using position threshold 1 = ' result])
%                                 end
%                             else
%                                 if verboseTog
%                                     disp('          Value must be a positive integer < 10')
%                                 else
%                                     disp('StepNet:  Value must be a positive integer < 10')
%                                 end
%                             end
%                         end
%                     end
%                 end
% 
%                 if proceed % posThresh2
%                     proceed = false; % reset
%                     result = DialogBox('position threshold 2 (positive int)',num2str(attEstParams.posThresh2),verboseTog);
%                     if ~isempty(result)
%                         if isNum(result)
%                             pt2 = str2num(result);
%                             if pt2 > 0 && pt2 < 10 && (round(pt2) == pt2)
%                                 attEstParams.posThresh2 = pt2;
%                                 proceed = true;
%                                 if verboseTog
%                                     disp(['StepNet:  Using position threshold 2 = ' result])
%                                 end
%                             else
%                                 if verboseTog
%                                     disp('          Value must be a positive integer < 10')
%                                 else
%                                     disp('StepNet:  Value must be a positive integer < 10')
%                                 end
%                             end
%                         end
%                     end
%                 end
                
                if proceed
                    proceed = false; % reset
                    attEstParsSet = true;
                    if verboseTog
                        disp('StepNet:  Proceed to calculate attractor estimate? (see dialog box)')
                    end
                    ButtonName = questdlg('Proceed to calculation?','Parameters set...','Yes','No','Yes');
                    switch ButtonName
                        case 'Yes'
                            proceed = true;
                    end
                end

                if proceed && DEqnsCompiled % the DEqnsCompiled is an internal assertion -- it should always be true
                    refresh(DSSRTfigHandle)
                    drawnow
                    if attest_ver == 0 || (attest_ver == 1 && ~regimeStruct.isCycle)
                        lenTS    = length(TransSeqs{TSIx}{TS_TSEQ});
                        poslo    = TransSeqs{TSIx}{TS_TSEQ}{1}{TSEQ_POSN};
                        poshi    = TransSeqs{TSIx}{TS_TSEQ}{lenTS}{TSEQ_POSN};
                        t0_TS    = seqTimes(poslo);
                        t0       = t0_TS;
                        Tperiod  = seqTimes(poshi) - seqTimes(poslo);
                        posPd    = poshi - poslo; % length of `period` (but not a cycle) in time steps
                        attVars  = vars(poslo:poshi,:);
                        attTimes = seqTimes(poslo:poshi);
                    else % attest_ver = 1 and isCycle
                        tlo      = regimeStruct.regimes(1).timeInt(1);
                        thi      = regimeStruct.regimes(regimeStruct.numRegs).timeInt(2);
                        Tperiod  = thi - tlo;
                        t0_TS    = seqTimes(regimeStruct.transSeq{TS_TSEQ}{1}{TSEQ_POSN});
                        t0       = tlo;
                        [isp poslo] = ismember(1, seqTimes <= tlo);
                        [isp poshi] = ismember(1, seqTimes <= thi);
                        max_t_found = seqTimes(poshi);
                        if abs(max_t_found-thi) > SMALLBOUND % then thi is out of bounds for seqTimes!
                            attVars_temp  = vars(poslo:poshi,:);
                            attTimes_temp = seqTimes(poslo:poshi);
                            if t0-Tperiod > 0 % fill in rest of variable data from previous cycle
                                poshi_temp = poslo-1;
                                [isp poslo_temp] = ismember(1, seqTimes <= max_t_found-Tperiod);
                                attVars  = [attVars_temp; vars(poslo_temp:poshi_temp,:)];
                                attTimes = [attTimes_temp; seqTimes(poslo_temp:poshi_temp)+Tperiod];
                                posPd = length(attTimes);
                            else
                                disp('StepNet: Error - not enough time data available to run AttEst!')
                                return
                            end
                        else
                            attVars  = vars(poslo:poshi,:);
                            attTimes = seqTimes(poslo:poshi);
                            posPd = poshi - poslo; % length of period in time steps
                        end
                    end

                    % NB -- only pass UNcompiled DEqns to AttEst, because it currently needs to work with them like that...

                    % DIAGNOSTIC OPTION
                    % t_start = cputime;
                    if attest_ver == 0
                        attRes = AttEst_epochs(TransSeqs{TSIx}, attVars, attTimes, DEqns, DEpars, t0_TS, Tperiod, attEstParams, ...
                            dScaleThresh, tScaleThresh, numExt, numInt, varnames, actsIxMap, inputsIx, varBounds, ...
                            DSSRTfigHandle, verboseTog);
                    else % == 1
                        attRes = AttEst_regimes(regimeStruct, attVars, attTimes, DEqns, DEpars, t0_TS, Tperiod, attEstParams, ...
                            dScaleThresh, tScaleThresh, numExt, numInt, varnames, actsIxMap, inputsIx, varBounds, ...
                            DSSRTfigHandle, verboseTog);
                    end
                    % t_total = cputime - t_start
                    % END DIAGNOSTIC
                    
                    if attRes{1} % success
                        disp(' ')
                        disp('StepNet:  Calculation successful')
                        beep % possibly make this only happen if t_total > 30 sec
                        epochs = attRes{2}; % this is also `regimes` for the regimes version!
                        numEpochs = attRes{3};
                        numEpochsTot = attRes{4}; % usually just 1+numEpochs, but leaving the option open for them to be different in the future
                        dt = attRes{5};
                        attEstData.dt = dt;
                        attEstData.numEpochs = numEpochsTot;
                        attEstData.Tperiod = Tperiod;
                        attPlotData = cell(1,numEpochsTot);
                        for ep = 1:numEpochsTot
                            epoch = epochs{ep};
                            epochLabel = num2str(epoch{1});
                            graphToffset = 0;
                            % correct for out-of-bounds start/end positions if epoch boundary on end of marked time interval
                            if poslo + epoch{3}(1) <= 0
                                epochVlp = 1;
                                if poslo + epoch{3}(2) > length(seqTimes)
                                    epochVrp = length(seqTimes);
                                end
                            else
                                if poslo + epoch{3}(1) > length(seqTimes)
                                    if attest_ver == 1 && regimeStruct.isCycle
                                        epochVlp = poslo + epoch{3}(1) - posPd;
                                        epochVrp = poslo + epoch{3}(2) - posPd;
                                        graphToffset = Tperiod;
                                    else
                                        epochVlp = length(seqTimes);
                                        epochVrp = epochVlp;
                                    end
                                else
                                    if poslo + epoch{3}(2) > length(seqTimes)
                                        if attest_ver == 1 && regimeStruct.isCycle
                                            epochVlp = poslo + epoch{3}(1) - posPd; % reset
                                            epochVrp = poslo + epoch{3}(2) - posPd;
                                            graphToffset = Tperiod;
                                        else
                                            epochVrp = epochVlp;
                                        end
                                    else
                                        epochVlp = poslo + epoch{3}(1);
                                        epochVrp = poslo + epoch{3}(2);
                                    end                                    
                                end
                            end
                            if epochVlp > epochVrp
                                disp('StepNet:  Internal error: epochVlp > epochVrp!')
                                beep
                                return
                            end
                            epochVdataPosns = epochVlp:epochVrp;
                            % WAS: epochTlength = epoch{2}(2) - epoch{2}(1); % may be inaccurate w.r.t. small position changes if epoch boundary at end of T interval
                            epochTdata = seqTimes(epochVdataPosns) - t0 + graphToffset; % offset to t0 = 0
                            epochTlength = epochTdata(length(epochTdata)) - epochTdata(1);
                            epochVdata = vars(epochVdataPosns,obsIx); % will be generalized later!
                            epochActs = epoch{4}; % probably remain unused
                            epochActNames = epoch{5};
                            epochAttInterval = epoch{6};
                            epochOrigFlag = epoch{7};
                            attPlotData{ep} = {epochLabel, epochTlength, epochTdata, epochVdata, epochActNames, ...
                                    epochAttInterval, epochOrigFlag };
                        end % for
                        if attEstParams.noShooting
                            attEstData.nameStr = ['Domain of Validity in variable ' varnames{obsIx} ' starting at time ' num2str(t0,'%.3f')];
                        else
                            attEstData.nameStr = ['Attractor Estimate in variable ' varnames{obsIx} ' starting at time ' num2str(t0,'%.3f')];
                        end
                        attEstData.plotData = attPlotData;
                        attEstData.obsIx = obsIx;
                        attEstData.t0 = t0;
                        attEstData.figTitle = strrep(figTitleName,'_','\_');
                        attEstData.version  = attest_ver; %% NEW to DSSRT v1.12
                        fignum = AttFig(attEstData);
                        if attEstParams.noShooting
                            disp('StepNet:  Displaying domain of validity data in graph')
                        else
                            disp('StepNet:  Displaying attractor estimate data in graph')
                        end
                        disp(['            for variable ' varnames{obsIx} ' starting at time ' num2str(t0,'%.3f')]);
                        attractorEst.exist = true;
                        attractorEst.number = attractorEst.number + 1;
                        attractorEst.attEstList{attractorEst.number} = attEstData;
                        if verboseTog
                            disp('           You can re-display this graph later by pressing `&`')
                        end
                    elseif attRes{5} == -1
                        % do nothing, operation cancelled by user, message already put out by AttEst
                    else
                        disp(' ')
                        disp('StepNet:  Calculation not completed')
                    end % if AttEst succeed
                elseif ~DEqnsCompiled
                        beep
                        disp('StepNet:  Internal error: differential equations should have been compiled already')
                        disp('            Cannot run AttEst')
                end % final `if proceed`
                if ~proceed % catch the cancel case
                    if verboseTog
                        disp('           Calculation cancelled')
                    end
                end
            else
                if verboseTog
                    disp('StepNet:  No functional network calculated!')
                end
            end
        case ' ' % step forward through current transition sequence if "fresh"
            if existNet
                if epochTimesSet
                    foundEp = false;
                    absTimeOff = freshTS{4};
                    for ep=1:length(epochTimes)
                        epTabs = epochTimes(ep) + absTimeOff;
                        if epTabs > seqTimes(pos) % find next epoch start time rel. to current position
                            newTime = epTabs;
                            foundEp = true;
                            break
                        end
                    end
                    if ~foundEp
                        newTime = epochTimes(1) + absTimeOff; % go back to start of cycle
                    end
                    pos = FindPosFromTime(newTime,seqTimes,numSteps);
                    refreshDiagState = true;
                else
                    if verboseTog
                        disp('StepNet:  No fresh transition sequence available')
                    end
                end
            else
                if verboseTog
                    disp('StepNet:  No functional network calculated!')
                end
            end
        case 'D' % find next identical state of actives backwards
            if existNet
                if pos > 1
                    if verboseTog
                        fprintf('StepNet:  External variable in focus is %s ...\n',varnames{obsIx})
                        disp(   '           Restrict search to this variable? (see dialog box)')
                    end
                    ButtonName = questdlg(['Restrict search to focused variable ' varnames{obsIx} '?'], ...
                        'Find identical state...', 'Tightly', 'Loosely', 'None', 'Loosely');
                    switch ButtonName
                        case 'Tightly'
                            rtype = 2;
                            proceed = true;
                        case 'Loosely'
                            rtype = 1;
                            proceed = true;
                        case 'None'
                            proceed = true;
                            rtype = 0;
                        otherwise
                            proceed = false;
                    end
                    if proceed
                        if rtype == 1
                            FSfocusId = [obsIx];
                            % then add all external obsIx-dependents to
                            % list
                            for dix=1:numExt
                                if ismember(FSfocusId(1),inputsIx{dix})
                                    FSfocusId = [FSfocusId, dix];
                                end
                            end
                        elseif rtype == 2
                            FSfocusId = [obsIx];
                        else
                            FSfocusId = [];
                        end
                        if verboseTog
                            fprintf(    'StepNet:  Searching for identical state of actives, backwards in time from pos %i\n',pos)
                            disp(     [ '             with "ignore order of actives" ' iAO_str ])
                            disp(     [ '             with restriction level to focused variable: `' ButtonName '`' ])
                        end
                        result = str2num(DialogBox('Guess Period',guessPeriod,verboseTog));
                        if result > 0 && result < seqTimes(numSteps)
                            guessPeriod = result;
                            result = str2num(DialogBox('margin %',marginPC,verboseTog));
                            if result > 0 && result < 1
                                marginPC = result;
                                refreshDiagState = true;
                                newPos = FindIdentState(Eseq,pos,ESEQ_ACTIVES_BAR,-1,ignoreActsPotsOrder,guessPeriod,...
                                    marginPC,FSfocusId,verboseTog); % initial condition is the \bar{curly E} set for the current time
                                if newPos ~= -1
                                    % don't update period searching backwards unless make the search find the *last* identical
                                    %  state going backwards that is contiguous with the first one found
%                                     guessPeriod = round( abs(newPos - pos)*timeStep );
                                    pos = newPos;
                                    fprintf('            Found identical state at pos %i\n',newPos)
%                                     fprintf('            Guess period updated to %.3f\n',guessPeriod)
                                else
                                    disp(   '            No identical state found backwards')
                                end
                            else
                                disp('StepNet:  Value not a valid percentage in (0,1)')
                            end
                        else
                            fprintf('StepNet:  Guess period out of time range (0,%.3f)\n',seqTimes(numSteps))
                        end
                    else
                        if verboseTog
                            disp('          Cancelled')
                        end
                    end
                else
                    if verboseTog
                        disp('StepNet:  Try a position > 1 for search')
                    end
                end
            else
                if verboseTog
                    disp('StepNet:  No functional network calculated!')
                end
            end
        case 'F' % find next identical state of actives forwards
            if existNet
                if pos < numSteps
                    if verboseTog
                        fprintf('StepNet:  External variable in focus is %s ...\n',varnames{obsIx})
                        disp(   '           Restrict search to this variable? (see dialog box)')
                    end
                    ButtonName = questdlg(['Restrict search to focused variable ' varnames{obsIx} '?'], ...
                        'Find identical state...', 'Tightly', 'Loosely', 'None', 'Loosely');
                    switch ButtonName
                        case 'Tightly'
                            rtype = 2;
                            proceed = true;
                        case 'Loosely'
                            rtype = 1;
                            proceed = true;
                        case 'None'
                            proceed = true;
                            rtype = 0;
                        otherwise
                            proceed = false;
                    end
                    if proceed
                        if rtype == 1
                            FSfocusId = [obsIx];
                            % then add all external obsIx-dependents to list
                            for dix=1:numExt
                                if ismember(FSfocusId(1),inputsIx{dix})
                                    FSfocusId = [FSfocusId, dix];
                                end
                            end
                        elseif rtype == 2
                            FSfocusId = [obsIx];
                        else
                            FSfocusId = [];
                        end
                        if verboseTog
                            fprintf(    'StepNet:  Searching for identical state of actives, forwards in time from pos %i\n',pos)
                            disp(     [ '             with "ignore order of actives" ' iAO_str ])
                            disp(     [ '             with restriction level to focused variable: `' ButtonName '`' ])
                        end
                        result = str2num(DialogBox('Guess Period',guessPeriod,verboseTog));
                        if result > 0 && result < seqTimes(numSteps)
                            guessPeriod = result;
                            result = str2num(DialogBox('margin percentage  (decimal in (0,1))',marginPC,verboseTog));
                            if result > 0 && result < 1
                                marginPC = result;
                                newPos = FindIdentState(Eseq,pos,ESEQ_ACTIVES_BAR,1,ignoreActsPotsOrder,guessPeriod,...
                                    marginPC,FSfocusId,verboseTog); % initial condition is the \bar{curly E} set for the current time
                                if newPos ~= -1
                                    guessPeriod = round( abs(newPos - pos)*timeStep );
                                    pos = newPos;
                                    refreshDiagState = true;
                                    fprintf('            Found identical state at pos %i\n',newPos)
                                    fprintf('            Guess period updated to %.3f\n',guessPeriod)
                                else
                                    disp(   '            No identical state found forwards')
                                end
                            else
                                disp('StepNet:  Value not a valid percentage in (0,1)')
                            end
                        else
                            fprintf('StepNet:  Guess period out of time range (0,%.3f)\n',seqTimes(numSteps))
                        end
                    else
                        if verboseTog
                            disp('          Cancelled')
                        end
                    end
                else
                    if verboseTog
                        fprintf('StepNet:  Try a position < %i for search\n',numSteps)
                    end
                end
            else
                if verboseTog
                    disp('StepNet:  No functional network calculated!')
                end
            end
        case '#' % change padding/swap algorithm level
            switch TSmatchTog
                case 2
                    TSmatchTog = 0;
                    TSm_str = 'PADDING and SWAPPING both OFF';
                case 0
                    TSmatchTog = 1;
                    TSm_str = 'PADDING ON, SWAPPING OFF';
                case 1
                    TSmatchTog = 2;
                    TSm_str = 'PADDING ON, SWAPPING ON';
			end
            disp(['StepNet:  Trans Seq matching algorithm levels: ' TSm_str ])
        case '7' % toggle switch for ignoring order of actives/potentials in diag state
            ignoreActsPotsOrder = not(ignoreActsPotsOrder);
            if ignoreActsPotsOrder,  iAO_str = 'ON';  else  iAO_str = 'OFF';  end
            disp(['StepNet:  Ignore order of actives/potentials in network state is ' iAO_str ]);
        case '8'
            disp('StepNet:  Command not yet implemented. Currently only ACTIVE variables are considered in transition sequences')
        case 'Z' % compare transition sequences stored in memory
            if existTransSeq
                numTSeqs = length(TransSeqs);
                longestTransSeqIx = 0;
                longestTransSeqLe = 0;
                for ts = 1:numTSeqs
                    lenTS = length(TransSeqs{ts}{TS_TSEQ});
                    if lenTS > longestTransSeqLe % then found a longer sequence
                        longestTransSeqIx = ts;
                        longestTransSeqLe = lenTS;
                    end
                end
                fprintf('StepNet:  Longest transition seq is #%i with %i steps\n',longestTransSeqIx,longestTransSeqLe)
%                 compBaseLen = longestTransSeqLe; % base length for comparing all the sequences
%                 if verboseTog
%                     disp('StepNet:  Currently can only compare 2 sequences of identical length')
%                 end
                numTSeqs = length(TransSeqs);
                result = DialogBox(['First trans seq # to compare (from 1 -> ' num2str(numTSeqs) ')'],1,verboseTog);
                if ~isempty(result)
                    TSoneIx = str2num(result);
                    if TSoneIx > 0 && TSoneIx <= numTSeqs % then okay
                        TSoneDone = true;
                    else
                        TSoneDone = false;
                        disp('StepNet:  Invalid sequence number entered')
                    end
                else
                    TSoneDone = false;
                end
                if TSoneDone
                    result = DialogBox(['Second trans seq # to compare (from 1 -> ' num2str(numTSeqs) ')'],numTSeqs,verboseTog);
                    if ~isempty(result)
                        TStwoIx = str2num(result);
                        if TStwoIx > 0 && TStwoIx <= numTSeqs && TStwoIx ~= TSoneIx % then okay
                            TStwoDone = true;
                        else
                            TStwoDone = false;
                            disp('StepNet:  Invalid sequence number entered (must also be different from first)')
                        end
                    else
                        TStwoDone = false;
                    end
                end
                if TSoneDone && TStwoDone
                    if verboseTog
                        fprintf('StepNet:  External variable in focus is %s ...\n',varnames{obsIx})
                        disp(   '           Restrict transition sequence to this variable? (see dialog box)')
                    end
                    ButtonName = questdlg(['Restrict transition sequence to focused variable ' varnames{obsIx} '?'], ...
                        'Generate transition sequence...', 'Tightly', 'Loosely', 'None', 'Loosely');
                    switch ButtonName
                        case 'Tightly'
                            TSfocusId = [obsIx];
                        case 'Loosely'
                            TSfocusId = [obsIx];
                            % then add all obsIx-dependents (including internals) to list
                            for dix=1:numTot % allows self-referential variables (i.e. those that directly input to themselves!)
                                if ismember(TSfocusId(1),inputsIx{dix})
                                    TSfocusId = [TSfocusId, dix];
                                end
                            end
                        case 'None'
                            TSfocusId = [];
                        otherwise
                            TSfocusId = [];
                    end
                    TSone = TransSeqs{TSoneIx};
                    TStwo = TransSeqs{TStwoIx};
                    if TSone{TS_FOCUS} == TStwo{TS_FOCUS}
                        % Could WARN if sequence lengths are more than, e.g., 5 different, that calc time
                        % may be GRRRRREAT due to padding algorithm, if padTimes > 2
                        lenTS1 = length(TSone{TS_TSEQ});
                        lenTS2 = length(TStwo{TS_TSEQ});
                        if lenTS1 < lenTS2 % this stuff only used if padding/swapping formed a new sequence
                            TSnewName = [TSone{TS_NAME} ' padded, vs. TS#' num2str(TStwoIx) ': # of errors in value Y'];
                            TSvalX = TSone{TS_VALX};
                        else
                            TSnewName = [TStwo{TS_NAME} ' padded, vs. TS#' num2str(TSoneIx) ': # of errors in value Y'];
                            TSvalX = TStwo{TS_VALX};
                        end
                        matchAlg = [ (TSmatchTog==1 | TSmatchTog==2), 2, TSmatchTog==2 ];
                        result = CompareTransSeq(TSone{TS_TSEQ}, TStwo{TS_TSEQ}, ignoreActsPotsOrder,...
                                             true, matchAlg, TSfocusId, verboseTog);
                        if ~isempty(result)
                            numErrs = result{1}(1);
                            if ~isempty(result{2}) % then padding/swapping formed a new sequence
%                                 numTSeqs = length(TransSeqs);
                                TransSeqs = { TransSeqs{:} {result{2}, TSnewName, TSvalX, numErrs, ...
                                            result{1}(2:numErrs+1), result{3}, TSone{TS_FOCUS}, [], 0 } };
                                disp('StepNet:  New transition sequence created from the originally SHORTER, now padded, sequence')
                            else % update compared sequence with errors, error positions (assumed to be 2nd sequence)
                                % we'll only be here if the trans seqs. were already the same length (no padding occurred)
                                % so we know that TSnewName will be based on that of TStwo, as we want...
                                % first, if 2nd already has error / padding info, ask if overwrite or create new seq.
                                if ~isempty(TStwo{TS_ERRS})
                                    ButtonName = questdlg('Overwrite existing padding / error info in second sequence (no to create new trans. seq.)?');
                                    switch ButtonName
                                        case 'Yes'
                                            TransSeqs{TStwoIx}{TS_PADS} = {};
                                            TransSeqs{TStwoIx}{TS_VALY} = numErrs;
                                            TransSeqs{TStwoIx}{TS_ERRS} = result{1}(2:numErrs+1);
                                            TransSeqs{TStwoIx}{TS_NAME} = [ EndStripStr(TStwo{TS_NAME},' compared to TS#') ' compared to TS#' num2str(TSoneIx) ];
                                        otherwise
                                            TransSeqs = { TransSeqs{:} {TStwo{TS_TSEQ}, EndStripStr(TStwo{TS_NAME},' compared to TS#'),...
                                                            TSvalX, numErrs, result{1}(2:numErrs+1), result{3}, TStwo{TS_FOCUS}, [], 0 } };
                                            fprintf('StepNet:  New transition sequence created from the second sequence, %s\n',TStwo{TS_NAME})
                                    end
                                else % can go ahead and write over 2nd sequence info because it was already blank
                                    TransSeqs{TStwoIx}{TS_VALY} = numErrs;
                                    TransSeqs{TStwoIx}{TS_ERRS} = result{1}(2:numErrs+1);
                                    TransSeqs{TStwoIx}{TS_NAME} = [ EndStripStr(TransSeqs{TStwoIx}{TS_NAME},' compared to TS#') ' compared to TS#' num2str(TSoneIx) ];
                                end
                            end
                            fprintf('StepNet:  Comparison results between transition seqs #%i:%s and #%i:%s\n', ...
                                    TSoneIx,TransSeqs{TSoneIx}{TS_NAME},TStwoIx,TransSeqs{TStwoIx}{TS_NAME});
                            if numErrs > 0
                                fprintf('            There were %i position errors, at positions:\n',numErrs);
                                fprintf('             ')
                                for errnum = 1:numErrs
                                    fprintf('%4i ',result{1}(errnum+1)); % error positions in the remainder of array result{1}
                                end
                                fprintf('\n')
                            else
                                disp(   '            *** Sequences are identical ***')
                            end
                        else
                            disp('StepNet:  Internal error in comparing sequences')
                        end
                    else
                        disp('StepNet:  Mismatch between focus of the chosen transition sequences')
                        disp('           Compare cancelled')
                    end
                else
                    disp('           Compare cancelled')
                end % else return control
            else
                if verboseTog
                    disp('StepNet:  No transition sequences to compare')
                end
            end
        case 'C' % compute transition sequence
            if existNet
                if markLset && markRset
                    TSfocusId = [];
                    if obsIxSet
                        if isempty(inputsIx{obsIx}) || allDEixMap(obsIx) == 0
                            disp('StepNet:  For transition sequences, must focus on a non-passive external')
                            disp('            variable with at least one input')
                            proceed = false;
                        else
                            proceed = true;
                        end
                        if proceed
                            if verboseTog
                                fprintf('StepNet:  External variable in focus is %s ...\n',varnames{obsIx})
                                disp(   '           Restrict transition sequence to this variable? (see dialog box)')
                            end
                            ButtonName = questdlg(['Restrict transition sequence to focused variable ' varnames{obsIx} '?'], ...
                                'Generate transition sequence...', 'Tightly', 'Loosely', 'None', 'Loosely');
                            switch ButtonName
                                case 'Tightly'
                                    rtype = 2;
                                    proceed = true;
                                case 'Loosely'
                                    rtype = 1;
                                    proceed = true;
                                case 'None'
                                    proceed = true;
                                    rtype = 0;
                                otherwise
                                    proceed = false;
                            end
                        end
                        % NOT SURE THIS IS NECESSARY NOW THAT THERE IS RegimeDet.m
%                         if proceed
%                             ButtonName = questdlg(['Ignore passive variables?'],'Generate transition sequence...','Yes','No','No');
%                             switch ButtonName
%                                 case 'Yes'
%                                     ignorePassive = true;
%                                     ignorePassiveStr = '';
%                                     proceed = true;
%                                 case 'No'
%                                     ignorePassive = false;
%                                     ignorePassiveStr = 'not ';
%                                     proceed = true;
%                                 otherwise
%                                     proceed = false;
%                             end
%                         end
                        ignorePassive = false;
                        ignorePassiveStr = 'not ';
                    else
                        disp('           No focused variable set to restrict')
                    end
                    if proceed && DEqnsCompiled
                        if rtype == 1
                            TSfocusId = [obsIx];
                            % then add all obsIx-dependents (including internals) to list
                            for dix=1:numTot % allows self-referential variables (i.e. those that directly input to themselves!)
                                if ismember(TSfocusId(1),inputsIx{dix})
                                    TSfocusId = [TSfocusId, dix];
                                end
                            end
                        elseif rtype == 2
                            TSfocusId = [obsIx];
                        end
                        if ignorePassive
                            passiveVars = find(allDEixMap==0); % in range 0..numTot
                            passiveVars = passiveVars(passiveVars<=numExt); % we will de-select only external passive vars
                        else
                            passiveVars = [];
                        end
                        if verboseTog
                            fprintf('StepNet:  Generating transition sequence from position %i to %i\n',markLeft,markRight)
                            fprintf('           %signoring passive variables\n',ignorePassiveStr)
                        end
                        drawnow
%                         % TEMP
%                         for i = 1:length(TSfocusId)
%                             fprintf('TSfocus var = %s\n', varnames{TSfocusId(i)})
%                         end
%                         for i = 1:length(passiveVars)
%                             fprintf('passive var = %s\n', varnames{passiveVars(i)})
%                         end
%                         disp('Candidate ix map for PDs_V:')
%                         for i =1:numTot
%                             if strcmp(varnames{i}, 'PDs_V')
%                                 for j = 1:length(candActsIxMapAbs{i})
%                                     fprintf(' %s\n', varnames{candActsIxMapAbs{i}(j)})
%                                 end
%                             end
%                         end
                        result = GenTransSeq(Eseq,markLeft,markRight,false,ignoreActsPotsOrder,TSfocusId, ...
                                             candActsIxMapAbs,passiveVars,verboseTog);
                        if ~isempty(result)
                            existTransSeq = true;
                            numTSeqs = length(TransSeqs);
                            if diagSw == PSIS
                                ts_suffix = 'Psi';
                            else
                                ts_suffix = 'term';
                            end
                            defTSname = ['TS_' ts_suffix '_' num2str(numTSeqs+1)];
                            TSname = DialogBox('Name for transition sequence',defTSname,verboseTog);
                            if isempty(TSname)
                                if verboseTog
                                    disp('StepNet:  Using default name for transition sequence')
                                end
                                TSname = defTSname;
                            end
                            if length(TSfocusId) > 0
                                TSname = [TSname ' - focus on ' varnames{TSfocusId(1)}];
                            end
                            TSvarValDef = 0;
                            TSvarValStr = DialogBox('Dependent variable value X associated with sequence',TSvarValDef,verboseTog);
                            if isempty(TSvarValStr)
                                TSvarValX = TSvarValDef;
                                if verboseTog
                                    disp('StepNet:  Using default variable X value')
                                end
                            else
                                TSvarValX = str2num(TSvarValStr);
                            end
                            TSvarValDef = 0;
                            TSvarValStr = DialogBox('Dependent variable value Y associated with sequence',TSvarValDef,verboseTog);
                            if isempty(TSvarValStr)
                                TSvarValY = TSvarValDef;
                                if verboseTog
                                    disp('StepNet:  Using default variable Y value')
                                end
                            else
                                TSvarValY = str2num(TSvarValStr);
                            end
                            TSerrpos  = []; % for new sequences this is empty -- used when comparing
                            TSpadding = {}; % for new sequences this is empty -- used when comparing
                            TransSeqs = { TransSeqs{:}, {result, TSname, TSvarValX, TSvarValY, TSerrpos, TSpadding, ...
                                        TSfocusId, passiveVars, markRight-markLeft} };
                            if freshTS{3} > 0 & ishandle(freshTS{3})
                                delete(freshTS{3})
                            end
                            freshTS = {true, length(TransSeqs), 0, seqTimes(markLeft), TSfocusId}; % length(TS) is = the TS id #
                            epochTimesSet = true;
                            if verboseTog
                                disp('StepNet:  Transition sequence created')
                            end
                        else
                            disp('StepNet:  Error -- no transition sequence returned')
                        end
                    else
                        if verboseTog
                            disp('StepNet:  Transition sequence generation cancelled')
                        end
                    end % if proceed
                else
                    disp('StepNet:  Left and right marks are not both set')
                end
            else
                disp('StepNet:  No functional network present from which to calculate transition sequence')
            end
        case 'R' % find dynamical regimes from transition sequences and variable data, timescales, etc.
            proceed = false; % initial value
            if existNet
                if existTransSeq & freshTS{3} > 0 & ishandle(freshTS{3}) % && markLset && markRset
                    if isempty(TransSeqs{freshTS{2}}{TS_PASS}) && isempty(TransSeqs{freshTS{2}}{TS_PADS})
                        proceed = true;
                    else
                        disp('StepNet:  Cannot find regimes if passive variables ignored in fresh transition sequence,')
                        disp('            or transition sequence has been padded to compare with another.')
                        disp('            Recompute transition sequence with `ignore` option off')
                        proceed = false;
                    end
%                     thisTS = TransSeqs{freshTS{2}}{TS_TSEQ}; % get out the actual transition seq. of curly E sets
%                     inLeft = freshTS{4} - seqTimes(markLeft);
%                     inRight = seqTimes(markRight) - freshTS{4} - thisTS{length(thisTS)}{TSEQ_TIME};
%                     if inLeft>=0 && inRight>=0 % then marks set inside or at the freshTS limits
%                         proceed = true;
%                     end
                else
                    disp('StepNet:  No fresh transition sequence available')
                end
            else
                disp('StepNet:  No functional network present')
            end
            if proceed
                ButtonName = questdlg('Is the transition sequence for a limit cycle?','Limit cycle?','Yes','No','Cancel','Yes');
                switch ButtonName
                    case 'Yes'
                        isCycle = true;
                        isCycStr = 'for a cycle';
                    case 'No'
                        isCycle = false;
                        isCycStr = '';
                    otherwise
                        proceed = false;
                end
            end
            if proceed
                ButtonName = questdlg('Be fussy about speed changes?','Speed change fussiness','Yes, not strict','Yes, strict','No','No');
                switch ButtonName
                    case 'Yes, not strict'
                        algOptions.speedFussy = true;
                        algOptions.speedFussyLevel = 1;
                        spFussyStr = '';
                    case 'Yes, strict'
                        algOptions.speedFussy = true;
                        algOptions.speedFussyLevel = 2;
                        spFussyStr = '';
                    case 'No'
                        algOptions.speedFussy = false;
                        algOptions.speedFussyLevel = 0;
                        spFussyStr = 'not';
                    otherwise
                        proceed = false;
                end
            end
            if proceed
                ButtonName = questdlg('Be fussy about fast variables leaving?','Fast leavers fussiness','Yes','No','No');
                switch ButtonName
                    case 'Yes'
                        algOptions.fastLeaveFussy = true;
                        flFussyStr = '';
                    case 'No'
                        algOptions.fastLeaveFussy = false;
                        flFussyStr = 'not ';
                    otherwise
                        proceed = false;
                end
            end
            if proceed
                ButtonName = questdlg('Be fussy about long epochs with slow variables?','Long epoch fussiness','Yes','No','No');
                switch ButtonName
                    case 'Yes'
                        algOptions.longFussy = true;
                        lgFussyStr = '';
                    case 'No'
                        algOptions.longFussy = false;
                        lgFussyStr = 'not ';
                    otherwise
                        proceed = false;
                end
            end
            if proceed % small change threshold
                proceed = false; % reset
                result = DialogBox('small change threshold',num2str(varChangeThresh_def),verboseTog);
                if ~isempty(result)
                    if isNum(result)
                        vctNew = str2num(result);
                        if vctNew > 0 && vctNew < 100
                            rdPars.varChangeThresh = vctNew;
                            proceed = true;
                            if verboseTog
                                disp(['StepNet:  Using small change threshold of ' result])
                            end
                        else
                            if verboseTog
                                disp('          Value must be a positive real < 1')
                            else
                                disp('StepNet:  Value must be a positive real < 1')
                            end
                        end
                    end
                end
            end
            if proceed && algOptions.speedFussy % derivThresh
                proceed = false; % reset
                result = DialogBox('derivative threshold',num2str(derivThresh_def),verboseTog);
                if ~isempty(result)
                    if isNum(result)
                        dThNew = str2num(result);
                        if dThNew > 0 && dThNew < 100
                            rdPars.derivThresh = dThNew;
                            proceed = true;
                            if verboseTog
                                disp(['StepNet:  Using derivative threshold of ' result])
                            end
                        else
                            if verboseTog
                                disp('          Value must be a positive real < 100')
                            else
                                disp('StepNet:  Value must be a positive real < 100')
                            end
                        end
                    end
                end
            end
            
            if proceed
                if verboseTog
                    disp(['StepNet:  Calculating reduced dynamical regimes ' isCycStr ' ' spFussyStr ' being fussy'])
                    disp(['            about speed changes, ' flFussyStr 'being fussy about fast leaving variables,'])
                    disp(['            and ' lgFussyStr 'being fussy about long epochs.'])
                    disp( '            This may take a few moments (press <escape> to cancel) ...')
                end
                rdPars.TS           = TransSeqs{freshTS{2}};
                rdPars.Eseq         = Eseq;
                rdPars.startTime    = freshTS{4};
                rdPars.focusSet     = freshTS{5};
                rdPars.tScaleThresh = tScaleThresh;
                rdPars.varnames     = varnames;
                rdPars.DEixMap      = allDEixMap;
                rdPars.caIxMap      = candActsIxMapAbs;
                rdPars.inputsIx     = inputsIx;
                rdPars.numExt       = numExt;
                rdPars.numInt       = numInt;
                result = RegimeDet(rdPars, seqTimes, vars, DEqns_comp, DEpars, varBounds, isCycle, algOptions, ...
                                    verboseTog, false, DSSRTfigHandle);
                if result{1} % success
                    existRegime  = true; % currently only one can be held in memory at a time
                    regimeStruct = result{2};
                    DisplayRegimes(regimeStruct,verboseTog)
                    % INFO on returned structure's fields
                    % regimes, numEps, numRegs, varnames, transSeq, isCycle, period
                    % startTime  --> not necessarily the same as the value passed to RegimeDet
                    % tScaleThresh, algOptions, focVarIx
                    % phasePlotInfo --> only if phase planes plotted successfully
                    ButtonName = questdlg('Show phase planes for regimes?','Regime phase planes','Yes','No','Yes');
                    switch ButtonName
                        case 'Yes'
                            regimeStruct.phasePlotInfo = PlotPhaseDiagrams(regimeStruct, seqTimes, vars, varBounds, DEqns_comp, ...
                                allDEixMap, numExt, networkObjectDir, DSSRTfigHandle, verboseTog);
                        otherwise
                            regimeStruct.phasePlotInfo = [];
                    end
                else
                    if verboseTog
                        disp('StepNet:  Regimes not calculated')
                    end
                end
            else
                if verboseTog
                    disp('StepNet:  Cancelled')
                end
            end
        case 'S' % save all transition sequence(s) to file
            if existTransSeq
                disp('StepNet:  Save all transition sequence(s) to file ...')
                if TSdataFnameSet
                    result = DialogBox('WRITE: Trans. Seq. data file name',EndStripStr(TSdataFilenameDisp,'_TS.',true,7),verboseTog);
                else
                    result = DialogBox('WRITE: Trans. Seq. data file name',[],verboseTog);
                end
                if ~isempty(result)
                    result = EndStripStr(result,'_TS.',true,7);
                    TSwrite = true; % initial value
                    if exist([ networkObjectDir '/' result '_TS.mat'],'file') == 2
                        ButtonName = questdlg('File exists. Overwrite?');
                        switch ButtonName
                            case 'Yes'
                                TSwrite = true;
                            otherwise
                                TSwrite = false;
                        end
                    end
                    if TSwrite
                        TSdataFilenameDisp = [result '_TS.mat'];
                        TSdataFname        = [ networkObjectDir '/' TSdataFilenameDisp ];
                        TSdataFnameSet     = true;
                        refresh(DSSRTfigHandle)
                        dummyEntry = {};
                        TransSeqData = { TransSeqs, dummyEntry }; % this gives us a different varname
                           % to check the existence of when reloading to
                           % test success of file read, since we now delete
                           % this variable
                        save(TSdataFname, 'TransSeqData');
                        clear TransSeqData
                        disp('           Transition sequences data file written')
                    else
                        disp('           Transition sequences data file not written')
                    end
                else
                    if verboseTog
                        disp('           Cancelled')
                    end
                end
            else
                if verboseTog
                    disp('StepNet:  No transition sequences in memory')
                end
            end
        case 'L' % load transition sequence(s) from file
            disp('StepNet:  Load transition sequence(s) from file ...')
            readTSokay = false;
            TSappend = false;
            if TSdataFnameSet
                result = DialogBox('READ: Trans Seq data file name',EndStripStr(TSdataFilenameDisp,'_TS.',true,7),verboseTog);
            else
                result = DialogBox('READ: Trans Seq data file name','demo',verboseTog); % TEMPORARY
            end
            if ~isempty(result)
                result = EndStripStr(result,'_TS.',true,7);
                if exist([ networkObjectDir '/' result '_TS.mat'],'file') == 2
                    TSdataFilenameDisp = [result '_TS.mat'];
                    TSdataFname    = [ networkObjectDir '/' TSdataFilenameDisp ];
                    TSdataFnameSet = true;
                    if existTransSeq
                        ButtonName = questdlg('Overwrite or append sequences?','Overwrite or append transition sequences',...
                                              'Overwrite','Append','Append');
                        switch ButtonName
                            case 'Overwrite'
                                TSappend = false;
                            otherwise
                                TSappend = true;
                        end
                    end
                    refresh(DSSRTfigHandle)
                    load(TSdataFname);
                    if exist('TransSeqData','var')~=1 % then file didn't have correct variable name inside
                        fprintf('StepNet:  Unreadable TS data file %s\n',TSdataFilenameDisp)
                        readTSokay = false;
                    else
                        disp('           TS data file read successfully')
                        readTSokay = true;
                    end
                else
                    fprintf('StepNet:  Transition sequence file name %s_TS.mat does not exist\n',result)
                end
            else
                if verboseTog
                    disp('           Cancelled')
                end
            end
            if readTSokay
                if TSappend
                    lenTS = length(TransSeqs);
                    disp('StepNet:  WARNING. Currently, StepNet does not prevent clashing trans. seq. names')
                    if verboseTog
                        disp('           Check manually for this and use `"` command to change names.')
                        fprintf('           The appended sequences have TS#`s from %i and up\n',lenTS+1)
                    end
                    TransSeqs = {TransSeqs{:} TransSeqData{1}{:}};
                else
                    if existTransSeq && verboseTog
                        disp('StepNet:  Previous sequences overwritten')
                    end
                    TransSeqs = TransSeqData{1};
                end
                existTransSeq = true;
            end
        case 'X' % list transition sequence(s) in memory
            if existTransSeq
                numTseqs = length(TransSeqs);
                fprintf('StepNet:  Transition sequences in memory 1 -> %i ...\n',numTseqs)
                for ts = 1:numTseqs
                    currTS = TransSeqs{ts};
                    fprintf(    '        *  Trans. Seq. #%i: `%s`\n', ts, currTS{TS_NAME});
                    if freshTS{3} > 0
                        if freshTS{2} == ts
                            disp('              FRESH ( these events are displayed in the main window )')
                        end
                    end
                    fprintf(    '              Variable X value: %3.4f\n', currTS{TS_VALX});
                    fprintf(    '              Variable Y value: %3i\n', currTS{TS_VALY});
                    fprintf(    '              Length (# of events): %3i\n',length(currTS{TS_TSEQ}));
                    if ~isempty( currTS{TS_PASS} )
                        fprintf('              Passive variables ignored')
                    end
                    if ~isempty( currTS{TS_ERRS} )
                        fprintf('              Error positions vs. original compared sequence:\n')
                        fprintf('                ')
                        for errnum = 1:currTS{TS_VALY}
                            fprintf('%4i ',currTS{TS_ERRS}(errnum)); % error positions in the remainder of array result{1}
                        end
                        fprintf('\n')
                    end
                    if ~isempty( currTS{TS_PADS} )
                        fprintf('              Padded at positions:\n')
                        fprintf('                ')
                        for padPos = 1:length(currTS{TS_PADS}{1})
                            if currTS{TS_PADS}{1}(padPos) > 0 % unused paddings entered as 0
                                fprintf('%3i x%3i,',currTS{TS_PADS}{1}(padPos),currTS{TS_PADS}{2}(padPos));
                            end
                        end
                        fprintf('\n')
                    end
                    TSfocusId = currTS{TS_FOCUS};
                    if length(TSfocusId) > 0
                        fprintf(    '              Focused on variable %s\n',varnames{TSfocusId(1)})
                    end
                end
            else
                disp('StepNet:  No transition sequences to list')
            end
        case '\' % delete transition sequence(s) from memory
            if existTransSeq
                if verboseTog
                    disp('StepNet:  Delete ALL transition sequences? (see dialog box)')
                    beep
                end
                ButtonName = questdlg('Delete ALL trans. seqs.? (use `No` to delete individuals)','Delete transition sequences?','Yes','No','Cancel','Cancel');
                switch ButtonName
                    case 'Yes'
                        TransSeqs = {};
                        existTransSeq = false;
                        if freshTS{3} > 0
                            freshTS{1} = false; % retain previous TS info for possible future use
                            if ishandle(freshTS{3})
                                delete(freshTS{3})
                            end
                            freshTS{3} = 0; % reset TSbox points plot handle
                            freshTS{4} = 0;
                            epochTimesSet = false;
                        end
                        disp('StepNet:  All transition sequences deleted')
                    case 'No'
                        numTSeqs = length(TransSeqs);
                        result = DialogBox(['a trans seq # to delete (from 1 -> ' num2str(numTSeqs) ')'],numTSeqs,verboseTog);
                        if ~isempty(result)
                            TSdeleteIx = str2num(result);
                            if TSdeleteIx > 0 && TSdeleteIx <= numTSeqs % then proceed to delete
                                newTS = {};
                                if TSdeleteIx > 1
                                    for ts = 1:TSdeleteIx-1
                                        newTS = { newTS{:} TransSeqs{ts} };
                                    end
                                end
                                if TSdeleteIx < numTSeqs
                                    for ts = TSdeleteIx+1:numTSeqs
                                        newTS = { newTS{:} TransSeqs{ts} };
                                    end
                                end
                                TransSeqs = newTS;

                                if freshTS{2} == TSdeleteIx % then have deleted TS associated with displayed TS in TSbox
                                    if freshTS{3} > 0
                                        freshTS{1} = false; % retain previous TS info for possible future use
                                        if ishandle(freshTS{3})
                                            delete(freshTS{3})
                                        end
                                        freshTS{3} = 0; % reset TSbox points plot handle
                                        freshTS{4} = 0;
                                        epochTimesSet = false;
                                    end
                                elseif freshTS{2} > TSdeleteIx
                                    freshTS{2} = freshTS{2} - 1; % position got shifted by 1
                                end
                                if verboseTog
                                    fprintf('StepNet:  Deleted transition sequence #%i\n',TSdeleteIx)
                                end
                            else
                                disp('           Invalid sequence number entered')
                            end
                        else
                                disp('           Delete cancelled')
                        end
                    otherwise
                        disp('           Delete cancelled')
                end % switch
            else
                disp('StepNet:  No transition sequences to clear')
            end
        case '"' % rename transition sequence and var values
            if existTransSeq
                numTSeqs = length(TransSeqs);
                if verboseTog
                    disp('StepNet:  Rename transition sequence and/or change value data')
                    disp('            (Remember, you can use this command to reset padding/error comparison entries,')
                    disp('            and the extended name, later on by setting Y=0)')
                end
                result = DialogBox(['a trans seq # to change associated info (from 1 -> ' num2str(numTSeqs) ')'],numTSeqs,verboseTog);
                if ~isempty(result)
                    TSrenameIx = str2num(result);
                    if TSrenameIx > 0 && TSrenameIx <= numTSeqs % then proceed to rename etc.
                        % new name
                        newTSname = DialogBox(['a new name'],TransSeqs{TSrenameIx}{TS_NAME},verboseTog);
                        if ~isempty(newTSname)
                            TransSeqs{TSrenameIx}{TS_NAME} = newTSname;
                            if verboseTog
                                disp('           New name set')
                            end
                        else
                            if verboseTog
                                disp('           Name change cancelled')
                            end
                        end
                        % new X val
                        newTSXvalStr = DialogBox(['a new X value'],TransSeqs{TSrenameIx}{TS_VALX},verboseTog);
                        if ~isempty(newTSXvalStr)
                            newTSXval = str2num(newTSXvalStr);
                            if newTSXval ~= TransSeqs{TSrenameIx}{TS_VALX}
                                TransSeqs{TSrenameIx}{TS_VALX} = newTSXval;
                                if verboseTog
                                    disp('           New X val set')
                                end
                            else
                                if verboseTog
                                    disp('           X val unchanged')
                                end
                            end
                        else
                            if verboseTog
                                disp('           X-val change cancelled')
                            end
                        end
                        % new Y val
                        newTSYvalStr = DialogBox(['a new Y value'],TransSeqs{TSrenameIx}{TS_VALY},verboseTog);
                        if ~isempty(newTSYvalStr)
                            newTSYval = str2num(newTSYvalStr);
                            if newTSYval ~= TransSeqs{TSrenameIx}{TS_VALY}
                                TransSeqs{TSrenameIx}{TS_VALY} = newTSYval;
                                if verboseTog
                                    disp('           New Y val set')
                                end
                                if newTSYval == 0 && ( ~isempty(TransSeqs{TSrenameIx}{TS_ERRS}) || ~isempty(TransSeqs{TSrenameIx}{TS_PADS}) )
                                    ButtonName = questdlg('Confirm clear old padding/error comparison data?');
                                    switch ButtonName
                                        case 'Yes'
                                            TransSeqs{TSrenameIx}{TS_ERRS} = [];
                                            TransSeqs{TSrenameIx}{TS_PADS} = {};
                                            TransSeqs{TSrenameIx}{TS_NAME} = EndStripStr( TransSeqs{TSrenameIx}{TS_NAME} ,' compare to TS#');
                                            disp('StepNet:  Padding / error comparison data cleared, and name reset')
                                        otherwise
                                            % donkey
                                    end
                                end
                            else
                                if verboseTog
                                    disp('           Y val unchanged')
                                end
                            end
                        else
                            if verboseTog
                                disp('           Y-val change cancelled')
                            end
                        end
                    else
                        disp('           Invalid sequence number entered')
                    end
                else
                    if verboseTog
                        disp('           Cancelled')
                    end
                end
            else
                disp('StepNet:  No transition sequences to rename')
            end
        case '|' % plot transition sequences` Y data vs. X data
            if existTransSeq
                numTSeqs = length(TransSeqs);
                result = DialogBox(['a space-separated list of trans seq #`s to plot (#`s are 1 -> ' num2str(numTSeqs) ')'],numTSeqs,verboseTog);
                plotOK = false; % initial value
                if ~isempty(result)
                    if ~ismember(',',result) & ~ismember(';',result) % common separators that might have been used by accident!
                        TSplotIxList = sscanf(result, '%i');
                        if min(TSplotIxList) > 0 && max(TSplotIxList) <= numTSeqs % then proceed to plotting
                            plotOK = true;
                        else
                            disp('StepNet:  Invalid sequence number(s) entered')
                        end
                    else
                        disp('StepNet:  List should contain only spaces as separators')
                    end
                end
                if plotOK
                    numTS = length(TSplotIxList);
                    Xdata = zeros(1,numTS);
                    Ydata = zeros(1,numTS);
                    TSlistStr = '';
                    for ts = 1:numTS
                        TSix = TSplotIxList(ts);
                        Xdata(ts) = TransSeqs{TSix}{TS_VALX};
                        Ydata(ts) = TransSeqs{TSix}{TS_VALY};
                        TSlistStr = [TSlistStr num2str(TSix)];
                    end
                    plotTSfigHandle = figure('NumberTitle','off','Name', ...
                        [ 'Transition Sequences: Y vs X data -- ' figTitleName],...
                        'Position',[300, 100, 500, 270]);
                    plot( Xdata, Ydata, 'ko' );
                    title(['Transition Sequences' TSlistStr ': Y vs X data'])
                    figure(DSSRTfigHandle) % restore figure handle to main window
                else
                    disp('StepNet:  Plot cancelled')
                end
            else
                disp('StepNet:  No transition sequences from which to plot data')
            end
        case 'V' % view transition sequence(s) in memory graphically
            if existTransSeq && CFGfilenameSet && exist('varnames','var') && exist('actsIxMap') % etc. etc.
                longestTransSeqIx = 0;
                longestTransSeqLe = 0;
                for ts = 1:length(TransSeqs)
                    lenTS = length(TransSeqs{ts}{TS_TSEQ});
                    if lenTS > longestTransSeqLe % then found a longer sequence
                        longestTransSeqIx = ts;
                        longestTransSeqLe = lenTS;
                    end
                end
                fprintf('StepNet:  Info -- longest trans seq is #%i with %i steps\n',longestTransSeqIx,longestTransSeqLe)
                numTSeqs = length(TransSeqs);
                result = DialogBox(['a space-separated list of EQUAL-LENGTH trans seq #`s to view (#`s are 1 -> ' num2str(numTSeqs) ')'],numTSeqs,verboseTog);
                viewOK = false; % initial value
                if ~isempty(result)
                    if ~ismember(',',result) & ~ismember(';',result) % common separators that might have been used by accident!
                        TSviewIxList = sscanf(result, '%i');
                        if min(TSviewIxList) > 0 && max(TSviewIxList) <= numTSeqs % then proceed
                            firstLen = length(TransSeqs{TSviewIxList(1)}{TS_TSEQ});
                            testTSLens = zeros(1,length(TSviewIxList));
                            for i=2:length(TSviewIxList)
                                testTSLens(i) = ( firstLen ~= length(TransSeqs{TSviewIxList(i)}{TS_TSEQ}) );
                            end
                            if sum( testTSLens ) == 0 % then all sequences are of equal length
                                viewOK = true;
                            else
                                disp('StepNet:  Sequences were not all of same length. View cancelled')
                                if verboseTog
                                    disp('           Use `Z` command to compare sequences and pad the shorter')
                                    disp('            so that they are of the same length')
                                end
                            end
                        else
                            disp('StepNet:  Invalid sequence number(s) entered')
                        end
                    else
                        disp('StepNet:  List should contain only spaces as separators')
                    end
                end
                if viewOK % then proceed to viewing
                    numTS = length(TSviewIxList);
                    TSlistStr = '';
                    TSlistSeq = cell(1,numTS);
                    for ts = 1:numTS
                        TSix = TSviewIxList(ts);
                        TSlistSeq{ts} = TransSeqs{TSix};
                        TSlistStr = [TSlistStr [' ' num2str(TSix)] ];
                    end
                    if verboseTog
                        disp('StepNet:  Passing control over to ViewTransSeq...')
                    end
                    if ~ishandle(viewFigHandle)
                        viewFigHandle = figure('NumberTitle','off','Name', ...
                            [ 'Transition Sequences Viewer: -- ' figTitleName],...
                            'Position',[80, 100, 750, 500], 'MenuBar', 'none', 'KeyPressFcn', 'DoNothing');
                    else
                        set(viewFigHandle,'Name',[ 'Transition Sequences Viewer: -- ' figTitleName])
                        figure(viewFigHandle)
                    end
                    title(['Viewing Transition Sequences' TSlistStr])
                    ViewTransSeq(viewFigHandle, TSlistSeq, varnames, ignoreActsPotsOrder, ...
                        numInt, numExt, actsIxMap, potsIxMap) % later, add arg in case potentials are cared about
                    figure(DSSRTfigHandle) % restore figure handle to main window
                    if viewFigHandle >0 & ishandle(viewFigHandle)
                        delete(viewFigHandle)
					end
                    if verboseTog
                        disp('StepNet:  Returned from ViewTransSeq')
                    end
                else
                    disp('StepNet:  Transition sequence(s) view cancelled')
                end
            else
                if ~existTransSeq
                    disp('StepNet:  No transition sequences to view')
                elseif ~exist('varnames','var') || ~exist('actsIxMap','var') % etc.
                    disp('StepNet:  CFG file must first be specified in order to setup system,')
                    disp('            before running `View`')
                end
            end
        case '1' % change start time for functional net calculation
            if varDataFilenameSet
                if stopTimeSet
                    startTimeLimStr = [ 'stop time ' num2str(stopTime,'%3.2f')];
                else
                    startTimeLimStr = num2str(availSpTime,'%3.2f');
                end
                if startTimeSet
                    result = DialogBox(['New start time  ( < ' startTimeLimStr ' )'],startTime,verboseTog);
                else
                    result = DialogBox(['New start time  ( < ' startTimeLimStr ' )'],[],verboseTog);
                end
                if ~isempty(result)
                    resultVal = str2num(result);
                    if resultVal >= availStTime && resultVal < availSpTime
                        if ~stopTimeSet
                            startTime = resultVal;
                            startTimeSet = true;
                            if verboseTog
                                fprintf('           Start time set to %.2f\n',startTime);
                            end
                        else
                            if resultVal < stopTime
                                startTime = resultVal;
                                startTimeSet = true;
                                if verboseTog
                                    fprintf('           Start time set to %.2f\n',startTime);
                                end
                            else
                                disp('StepNet:  Invalid start time -- must be smaller than current stop time')
                            end
                        end
                    else
                        disp('StepNet:  Invalid start time -- check available start/stop limits using `=`')
                    end
                    existFNpars = CFGfilenameSet & varDataFilenameSet & nOutSet & startTimeSet & stopTimeSet & dScaleThreshSet;
                else
                    if verboseTog
                        disp('           Cancelled')
                    end
                end
            else
                disp( 'StepNet:  Specify data file first using `O`')
            end
        case '2' % change stop time for functional net calculation
            if varDataFilenameSet
                if startTimeSet
                    stopTimeLimStr = [ 'start time ' num2str(startTime,'%3.2f')];
                else
                    stopTimeLimStr = num2str(availStTime,'%3.2f');
                end
                if stopTimeSet
                    result = DialogBox(['New stop time  ( > ' stopTimeLimStr ' )'],stopTime,verboseTog);
                else
                    result = DialogBox(['New stop time  ( > ' stopTimeLimStr ' )'],[],verboseTog);
                end
                if ~isempty(result)
                    resultVal = str2num(result);
                    if resultVal > availStTime && resultVal <= availSpTime
                        if ~startTimeSet
                            stopTime = resultVal;
                            stopTimeSet = true;
                            if verboseTog
                                fprintf('           Stop time set to %.2f\n',stopTime);
                            end
                        else
                            if resultVal > startTime
                                stopTime = resultVal;
                                stopTimeSet = true;
                                if verboseTog
                                    fprintf('           Stop time set to %.2f\n',stopTime);
                                end
                            else
                                disp('StepNet:  Invalid stop time -- must be larger than current start time')
                            end
                        end
                    else
                        disp('StepNet:  Invalid stop time -- check available start/stop limits using `=`')
                    end
                    existFNpars = CFGfilenameSet & varDataFilenameSet & nOutSet & startTimeSet & stopTimeSet & dScaleThreshSet;
                else
                    if verboseTog
                        disp('           Cancelled')
                    end
                end
            else
                disp( 'StepNet:  Specify data file first using `O`')
            end
        case '3' % change dominance and time scale thresholds
            if dScaleThreshSet
                result = DialogBox('Dom. scale threshold denominator ( > 1.0 )',dScaleThresh,verboseTog);
            else
                result = DialogBox('Dom. scale threshold denominator ( > 1.0 )',[],verboseTog);
            end
            if ~isempty(result)
                resultVal = str2num(result);
                if resultVal > 1
                    dScaleThresh = resultVal;
                    dScaleThreshSet = true;
                    if verboseTog
                        fprintf('           Dominance scale threshold set to 1/%.4f\n',dScaleThresh);
                    end
                else
                    disp('StepNet:  Invalid scale threshold')
                end
                existFNpars = CFGfilenameSet & varDataFilenameSet & nOutSet & startTimeSet & stopTimeSet & dScaleThreshSet;
            else
                if verboseTog
                    disp('           Cancelled')
                end
            end
            if tScaleThreshSet
                result = DialogBox('New time-scale threshold denominator ( > 1.0 )',tScaleThresh,verboseTog);
            else
                result = DialogBox('New time-scale threshold denominator ( > 1.0 )',[],verboseTog);
            end
            if ~isempty(result)
                resultVal = str2num(result);
                if resultVal > 1
                    tScaleThresh = resultVal;
                    tScaleThreshSet = true;
                    if verboseTog
                        fprintf('           Time scale threshold set to 1/%.4f\n',tScaleThresh);
                    end
                else
                    disp('StepNet:  Invalid scale threshold')
                end
            else
                if verboseTog
                    disp('           Cancelled')
                end
            end
        case '4' % change resampling-step in functional net calculation
            if nOutSet
                result = DialogBox('New nOut re-sampling step',nOut,verboseTog);
            else
                result = DialogBox('New nOut re-sampling step',[],verboseTog);
            end
            if ~isempty(result)
                resultVal = str2num(result);
                if resultVal <= 0 || ( resultVal > 0 && round(resultVal) ~= resultVal )
                    disp('StepNet:  Value out of range. Please enter a positive integer')
                else
                    nOut = resultVal;
                    nOutSet = true;
                    if verboseTog
                        fprintf('           nOut set to %i\n',nOut);
                    end
                end
                existFNpars = CFGfilenameSet & varDataFilenameSet & nOutSet & startTimeSet & stopTimeSet & dScaleThreshSet;
            else
                if verboseTog
                    disp('           Cancelled')
                end
            end
        case '5' % change observable of focus in functional network
            if existNet
                if obsIxSet
                    fprintf('StepNet:  Current observable of focus is #%i, named `%s`\n',obsIx,varnames{obsIx})
                end
                if verboseTog
                    disp('StepNet:  Choose new observable (external variable) of focus from the following:')
                    disp('           Valid variable names:')
                    for vix=1:numExt
                        disp(['             ' varnames{vix}])
                    end
                else
                    disp('StepNet:  Choose new observable (external variable) of focus')
                end
                if obsIxSet
                    defIx = obsIx;
                else
                    defIx = 1;
                end
                result = DialogBox('New observable of focus',varnames{defIx},verboseTog);
                if ~isempty(result)
                    newVarName = result;
                    [p pIx] = ismember( newVarName, varnames );
                    if p
                        if pIx <= numExt
                            obsIx = pIx;
                            if attEstParsSet
                                attEstParsSet = false; % reset, so that defaults (w/ bounds) are correctly recomputed for an AttEst() call
                            end
                            if verboseTog
                                disp(['           Using external variable `' newVarName '`'])
                            end
                        else
                            disp('          Invalid external variable name')
                        end
                    else
                        disp('          Invalid external variable name')
                    end
                end
            else
                disp('StepNet:  No functional network calculated')
            end
        case '6' % change variable in Variable Viewer
            if CFGfilenameSet
                if verboseTog
                    disp(       'StepNet:  Variables that can be viewed:')
                    for vix = 1:numTot
                        fprintf('             %s\n',varnames{vix});
                    end
                end
                result = DialogBox(['New variable to plot' ],varnames{varViewIx},verboseTog);
                if ~isempty(result)
                    newVarName = result;
                    [p pIx] = ismember( newVarName, varnames );
                    if pIx>0
                        varViewIx = pIx;
                        if verboseTog
                            disp(['           Changed variable to view to `' result '`'])
                        end
                        PlotVBox(handVVaxes,seqTimes,vars(:,varViewIx),varViewIx,varnames)
                    else
                        disp('           Invalid variable name')
                    end
                    axes(handFNaxes) % return to FN axes
                else
                    if verboseTog
                        disp('           Cancelled')
                    end
                end
            end
        case '-' % switch between display of functional net for
            % influence strengths (Psis) and for term sizes
            if existNet
                if diagSw == TERMS
                    disp_name = 'influence strengths (Psi values)';
                    diagSw = PSIS;
                else
                    disp_name = 'term sizes';
                    diagSw = TERMS;
                end
                if diagSw == PSIS
                    diag = diag_Psi;
                    Eseq = Eseq_Psi;
                    domGraph = domGraph_Psi;
                    domInfo = domInfo_Psi;
                else
                    diag = diag_term;
                    Eseq = Eseq_term;
                    domGraph = domGraph_term;
                    domInfo = domInfo_term;
                end
                delete(diagSw_handle)
                if verboseTog
                    fprintf('StepNet:  Switched FN diagram to show dominance in terms of %s\n', disp_name)
                end
                diag_state_old = zeros(1, numobjects);  % reset
                refreshDiagState = true; % force update of display
                forceNetRedraw = true;
            end
        case 'U' % run FuncNet if parameters exist
            if existFNpars
                if ~DEqnsCompiled % internal assertion
                    beep
                    disp('StepNet:  Internal error: differential equations not compiled already')
                    return
                end
                if existNet
                    ButtonName = questdlg('Are you sure?','Re-calculate functional network','Yes','No','No');
                    switch ButtonName
                        case 'Yes'
                            proceed = true;
                            disp('StepNet:  Re-calculating functional network')
                        case 'No'
                            proceed = false;
                            if verboseTog
                                disp('StepNet:  Cancelling re-calculation of functional network')
                            end
                    end
                else
                    disp('StepNet:  Calculating functional network')
                    proceed = true;
                end
				%    FuncNet returns this order of params...
				%    seqTimes, diag_states, Eseq, numTot, whosevars, candActsIxMapAbs, ...
				%    candPotsIxMapAbs, allVarNames, var_data, hobj, numSteps, tstep
                %
                % If we have already calc'd these from reading a data file
                % & CFG file, then pass this to FuncNet so it doesn't have
                % to re-read the file -- but don't recycle seqTimes and var_data since these are restricted in time
                % (new start and stop times for this Func Net may be outside the old range) - so use seqTimesAll and varsAll
                % re-nOut'ed using current nOut (which may be different to last time)
                if proceed
                    if exist('actsIxMap','var') && exist('varnames','var') % arbitrarily use these as the indicators
                        if ~isempty(varnames)
                            CFGsetup = cell(1,9);
                            CFGsetup{1} = numInt;
                            CFGsetup{2} = numExt;
							CFGsetup{3} = {varnames{numExt+1:numTot}};
							CFGsetup{4} = {varnames{1:numExt}};
							CFGsetup{5} = varnames;
							CFGsetup{6} = inputsIx;
							CFGsetup{7} = actsIxMap;
							CFGsetup{8} = potsIxMap;
                            if exist('whosevars','var')
                                CFGsetup{9} = whosevars;
                            end
                        else
                            CFGsetup = {};
                        end
                    else
                        CFGsetup = {}; % so that FuncNet can detect that it's not set
                    end
                    figure(DSSRTfigHandle)
                    drawnow
					funcNetData = FuncNet( networkObjectDir, CFGfilename, dScaleThresh, varnames, DEqns_comp, DEpars, ...
                         startTime, stopTime, hobj, DSSRTfigHandle, varBounds, varsAll(1:nOut:length(seqTimesAll),:), ...
                         seqTimesAll(1:nOut:length(seqTimesAll),:), CFGsetup, true );
					if isempty(funcNetData)
                        disp('StepNet:  FuncNet did not return anything ...')
                        disp('           Possible configuration error in setup parameters')
                    elseif ~iscell(funcNetData{1}) && funcNetData{1} == -1
                        % do nothing, user cancelled the operation
                    else
						diag_Psi         = funcNetData{1};
                        diag_term        = funcNetData{2};
						Eseq_Psi         = funcNetData{3};
                        Eseq_term        = funcNetData{4};
						numTot           = funcNetData{5};
						whosevars        = funcNetData{6};
						candActsIxMapAbs = funcNetData{7};
						candPotsIxMapAbs = funcNetData{8};
						varnames         = funcNetData{9};
                        numSteps         = funcNetData{10};
                        timeStep         = funcNetData{11};
                        domGraph_Psi     = funcNetData{12};
                        domGraph_term    = funcNetData{13};
                        domInfo_Psi      = funcNetData{14};
                        domInfo_term     = funcNetData{15};
                        seqTimes         = funcNetData{16}; % just the version sent to it but restricted to start & stop times
                        vars             = funcNetData{17}; % just the version sent to it but restricted to start & stop times
                        if numobjects ~= length(diag_Psi{1})
                            disp('StepNet:  Fatal error running FuncNet:')
                            disp('           Number of graphical objects does not match dimension of diagram')
                            disp('            state variable returned by FuncNet.')
                            disp('           Failed to create a diagram state')
                        else
                            if diagSw == PSIS
                                diag = diag_Psi;
                                Eseq = Eseq_Psi;
                                domGraph = domGraph_Psi;
                                domInfo = domInfo_Psi;
                            else
                                diag = diag_term;
                                Eseq = Eseq_term;
                                domGraph = domGraph_term;
                                domInfo = domInfo_term;
                            end
                            if numVBs ~= length(varnames)
                                disp('StepNet:  Number of vertical bars defined exceeds number of defined variables')
                                disp('            according to size of variables array returned by FuncNet.')
                                disp('           Failed to create a diagram state')
                            else
                                pos = 1; % reset this
%                                 if ~existNet % then need to initialize diagram states for first time
                                % possible bug when 'if ~existNet' retained
                                % around this refresh loop.
                                for obj = 1:numobjects
                                    if hobj{obj}{OBJ_TYPE} == LINK
                                        if hobj{obj}{OBJ_STATE} ~= LINKOFFSTATE
                                            hobj{obj}{OBJ_STATE} = diag{pos}(obj); % refresh initialize diagram states
                                        end % else don't change it
                                    else
                                        hobj{obj}{OBJ_STATE} = diag{pos}(obj); % refresh initialize diagram states
                                    end
                                    hobj{obj}{OBJ_HANDLE} = DrawNet( hobj{obj} );
                                end
                                existNet = true;
%                                 end
                                if markRset
                                    markRset = false;
                                    markRight = 0;
                                    delete(markRhand)
                                    clear markRhand
                                end
                                if markLset
                                    markLset = false;
                                    markLeft = 0;
                                    delete(markLhand)
                                    clear markLhand
                                end
                                if verboseTog
                                    disp('StepNet:  Marks cleared')
                                end
                                beep
                                for vn = 1:numVBs % update vertical bar objects with their variable names
                                    hvb{vn}{OBJ_LABEL} = varnames{vn};
                                end
                                if ~obsIxSet
                                    obsIx = 1; % by default
                                    obsIxSet = true;
                                end
                                % display first variable in VarViewBox unless
                                % one was already defined for a functional
                                % network
                                if ~existNet
                                    varViewIx = 1;
                                end
                                PlotVBox(handVVaxes,seqTimes,vars(:,varViewIx),varViewIx,varnames)
                                % update axes etc. for Time Box and TS box
                                axes(handTBaxes)
                                delete(tBarHand)
                                set(handTBaxes, 'XLim', [ seqTimes(1) seqTimes(numSteps)])
                                tBarHand = Draw_tBar(seqTimes(1),seqTimes(numSteps),'k');
                                axes(handTSaxes)
                                cla
                                if freshTS{1}
                                    if ishandle(freshTS{3})
                                        delete(freshTS{3})
                                    end
                                    freshTS{3} = 0;
                                    freshTS{4} = 0;
                                    freshTS{1} = false;
                                    epochTimesSet = false;
                                end
                                set(handTSaxes, 'XLim', [ seqTimes(1) seqTimes(numSteps)])
                                axes(handFNaxes)

                                refreshDiagState = true; % force update of display
                                forceNetRedraw = true;
                                if verboseTog
                                    beep % possibly make this only happen if t_total > 30 sec
                                end
                                figTitleName = [ CFGfilenameDisp ' : ' varDataFilenameDisp ];
                                disp(['StepNet:  New functional network appearing: `' figTitleName '`'])
                                set(DSSRTfigHandle,'Name',[ figTitle1 figTitleName])
                            end
                        end
                    end
                end % if proceed
            else
                disp('StepNet:  Some parameters for FuncNet.m not yet specified by user')
                disp('           Press `=` to show parameters set')
            end
        case '`' % display variable information
            if existNet
                varFocusId = [];
                if obsIxSet
                    if verboseTog
                        fprintf('StepNet:  External variable in focus is %s ...\n',varnames{obsIx})
                        disp(   '           Restrict variable info to the dependents of this variable? (see dialog box)')
                    end
                    ButtonName = questdlg(['Restrict variable info to focused variable ' varnames{obsIx} '?'],'Variable info...');
                    switch ButtonName
                        case 'Yes'
                            varFocusId = [obsIx];
                            proceed = true;
                        case 'No'
                            proceed = true;
                        otherwise
                            proceed = false;
                    end                    
                end
                if ~isempty(varFocusId) % then add all obsIx-dependents (including internals) to list
                    for dix=1:numTot % allows self-referential variables (i.e. those that directly input to themselves!)
                        if ismember(varFocusId(1),inputsIx{dix})
                            varFocusId = [varFocusId, dix];
                        end
                    end
                    focusedStr = ['focused on variable ' varnames{obsIx}];
                else
                    focusedStr = 'for all variables';
                    varFocusId = 1:numTot;
                end
                if proceed
                    if verboseTog
                        ts_temp_str = '           ';
                    else
                        ts_temp_str = 'StepNet:  ';
                    end
                    disp([ ts_temp_str 'Variable info ' focusedStr ' at time ' num2str(seqTimes(pos),'%.3f') ':'])
                    disp( '             Variable values:')
                    for varIx=varFocusId
                        fprintf('              %s     %4.5f\n', varnames{varIx}, vars(pos,varIx))
                    end
                    disp( '             Variable time-scales (for those governed by differential equations):')
                    focusVarNames = {varnames{varFocusId}};
                    deqnIxMap = zeros(1,numTot); % (only for focused variables) - zero entries will denote no corresponding differential equation present
                    for eqn=1:numEqns
                        % create map entry
                        [ispresent1 pos1] = ismember(DEqns{eqn}{DE_NAMEi},focusVarNames); % pos1 is relative to focusVarNames
                        if ispresent1
                            deqnIxMap(varFocusId(pos1)) = eqn;
                        end
                    end
                    if DEqnsCompiled
                        for varIx=varFocusId
                            DEqnIx = allDEixMap(varIx);
                            if DEqnIx > 0
                                gam1 = DEqns_comp{DEqnIx}{DE_GAMMA1i};
                                sg1 = SumGamma1(gam1,vars(pos,:));
                                tscale = DEqns_comp{DEqnIx}{DE_CFACi} / sg1(1);
                                fprintf('              %s     %4.5f\n', varnames{varIx},tscale)
                            end
                        end
                    else
                        disp('StepNet:  Internal error with compilation of differential equations')
                    end
                    disp(  ['             Reduced model accuracy indicators for primary focused variable ' varnames{obsIx} ':'])
                    %R_abs,g_rat_abs,err_abs,sc_sep

                    % ratio of g's for selected Psi or Term-based reduced
                    % model
                    obsDEix    = allDEixMap(obsIx);
                    currActs   = cell(1,numExt);
                    currActs{obsIx} = Eseq{pos}{ESEQ_ACTIVES_BAR}{obsIx};
                    DEqns_comp = SetActivesSwitch(DEqns_comp, obsIx, allDEixMap, numExt, actsIxMap, currActs);
                    sumGam1_model = SumGamma1(DEqns_comp{ obsDEix }{DE_GAMMA1i}, vars(pos,:), true);
                    sumGam1_full  = SumGamma1(DEqns_comp{ obsDEix }{DE_GAMMA1i}, vars(pos,:), false);
                    g_rat      =  sumGam1_model(1) / sumGam1_full(1);
                    g_rat_abs  = abs( g_rat );

                    % Psi values for |R| ("remainder" between reduced model
                    % and full model's sum of Psi*/p (should = 0))
                    actPsi_star = GetPsiVals(vars(pos,:),DEqns_comp{obsDEix},false); % signed
                    R_sum       = 0;
                    num_terms   = 0;
                    for gam1Ix=1:length(DEqns_comp{ obsDEix }{DE_GAMMA1i})
                        if ~DEqns_comp{ obsDEix }{DE_GAMMA1i}{gam1Ix}{DE_ACTSWi}
                            R_sum = R_sum + actPsi_star(gam1Ix) / DEqns_comp{ obsDEix }{DE_GAMMA1i}{gam1Ix}{DE_GAMVARPOWi};
                            num_terms = num_terms + 1;
                        end
                    end
                    for gam2Ix=1:length(DEqns_comp{ obsDEix }{DE_GAMMA2i})
                        if ~DEqns_comp{ obsDEix }{DE_GAMMA2i}{gam2Ix}{DE_ACTSWi}
                            R_sum = R_sum + actPsi_star(gam2Ix+gam1Ix) / DEqns_comp{ obsDEix }{DE_GAMMA2i}{gam2Ix}{DE_GAMVARPOWi};
                            num_terms = num_terms + 1;
                        end
                    end
                    R_abs = abs( R_sum ) / num_terms;
                    
                    % absolute `error` between V_inf and \tilde{V}_inf
                    err_abs = abs( R_sum / g_rat );

                    % scale separation (~ "spectral gap")
                    Psi_data = zeros(1,numTot);
                    mostDom_Psi = domInfo_Psi{obsIx,pos}{4};
                    ixa_dom_Psi = domInfo_Psi{obsIx,pos}{3};
                    mostDomIx_Psi = candActsIxMapAbs{obsIx}(ixa_dom_Psi(1));
                    Psi_data(mostDomIx_Psi) = mostDom_Psi;
                    ratios_act_Psi = domInfo_Psi{obsIx,pos}{2};
                    lowerIx = min(find(ratios_act_Psi > dScaleThresh));
                    if lowerIx == 1
                        sc_sep = ratios_act_Psi(lowerIx);
                    else
                        sc_sep = ratios_act_Psi(lowerIx)/ratios_act_Psi(lowerIx-1);
                    end
                    
                    % input term sizes
                    varDataLine = vars(pos,:);
                    cur_data_tot = 0;
                    cur_data_act = 0; % for local model active inputs
                    currActs_term = cell(1,numExt);
                    currActs_term{obsIx} = Eseq_term{pos}{ESEQ_ACTIVES_BAR}{obsIx};
                    % set actives switches for term-based reduced model
                    DEqns_term_comp = SetActivesSwitch(DEqns_comp, obsIx, allDEixMap, numExt, actsIxMap, currActs);
                    thisDE = DEqns_term_comp{ obsDEix };
                    for gamIx = 1:length(thisDE{DE_GAMMA1i})
                        g1term = thisDE{DE_GAMMA1i}{gamIx};
                        varix = g1term{DE_GAMVARNAMEi};
                        if varix == 0
                            varix = g1term{DE_TARGETi};
                            varVal = 1;
                        else
                            varVal = varDataLine(varix);
                        end
                        tau_recipVal = GetTauRVal(g1term{DE_ISTAUFFILEi},g1term{DE_TAURECIPi},varVal,g1term{DE_GAMVARPOWi});
                        targetVal = GetTargVal(g1term{DE_ISTGTFFILEi},g1term{DE_TARGETi},varVal,varDataLine);
                        if g1term{DE_INTVARNAMEi} > 0
                            g1termVal = tau_recipVal * varDataLine(g1term{DE_INTVARNAMEi})^g1term{DE_INTVARPOWi} * (targetVal - varDataLine(obsIx));
                        else
                            g1termVal = tau_recipVal * (targetVal - varDataLine(obsIx));
                        end
                        cur_data_tot = cur_data_tot + abs(g1termVal);
                        if g1term{DE_ACTSWi}
                            cur_data_act = cur_data_act + abs(g1termVal);
                        end
                    end
                    for gamIx = 1:length(thisDE{DE_GAMMA2i})
                        g2term = thisDE{DE_GAMMA2i}{gamIx};
                        g2termVal = GetTauRVal(g2term{DE_ISTAUFFILEi},g2term{DE_TAURECIPi},varDataLine(g2term{DE_GAMVARNAMEi}),g2term{DE_GAMVARPOWi});
                        cur_data_tot = cur_data_tot + abs(g2termVal);
                        if g2term{DE_ACTSWi}
                            cur_data_act = cur_data_act + abs(g2termVal);
                        end
                    end
                    % because we're only looking at the ratio, there's no
                    % need to involve the capacitance
                    input_rat = abs( cur_data_act / cur_data_tot );

                    fprintf('               |R|/|Acts''| = %.4f, |target difference| = %.4f, \n',R_abs,err_abs)
                    fprintf('               |g_tot_model / g_tot_full| = %.4f, scale separation at cutoff = %.4f\n', g_rat_abs, sc_sep)
                    fprintf('               |active input terms / tot inputs | = %.4f\n', input_rat)
                else
                    if verboseTog
                        disp('           Cancelled')
                    end
                end
            else
                disp('StepNet:  No functional network held in memory')
            end
        case '!' % plot timescales of all active variables (rel. to focused variable) over a marked interval (log plot)
            if existNet
                proceed = false;
                if markLset && markRset                    
                    ts_str = {'timescale', 'speed'};
                    ts_keystr = {'tau_','speed of'};
                    disp('StepNet:  Plot variables timescales or speeds (= 1/timescales)? (see dialog box)')
                    ButtonName = questdlg('Plot variable timescales or speeds?','Timescales or speeds','Times','Speeds','Cancel','Times');
                    switch ButtonName
                        case 'Times'
                            speedSwitch = 0;
                            proceed = true;
                        case 'Speeds'
                            speedSwitch = 1;
                            proceed = true;
                        otherwise
                            proceed = false;
                    end
                else
                    if verboseTog
                        disp('            Both marks not set')
                    end
                end
                if proceed
                    disp('StepNet:  Plot y-axis in log scale? (see dialog box)')
                    ButtonName = questdlg('Plot y-axis on log scale?','Log scale?','Yes','No','Yes');
                    switch ButtonName
                        case 'Yes'
                            plotLogYScale = true;
                        case 'No'
                            plotLogYScale = false;
                    end

                    if obsIxSet
                        disp(['StepNet:  Showing ' ts_str{speedSwitch+1} 's of variables associated with focused'])
                        fprintf('           variable %s over marked interval\n',varnames{obsIx})
                        disp('           (this may take a few moments -- press <escape> to stop)')
                        focusedSet = [obsIx];
                        % add all obsIx-dependents (including internals) to focusedSet that have associated differential eqns
                        for dix=1:numTot % allows self-referential variables (i.e. those that directly input to themselves!)
                            if ismember(obsIx,inputsIx{dix}) && allDEixMap(dix) > 0
                                focusedSet = [focusedSet, dix];
                            end
                        end
                        proceed = true;
                    else
                        if ~obsIxSet
                            disp('StepNet:  Must set variable of focus first')
                        else
                            disp('StepNet:  Cancelled')
                        end
                        proceed = false;
                    end
                end

                if proceed && DEqnsCompiled
                    stop_ts = false;
                    tPlotInterval = [seqTimes(markLeft), seqTimes(markRight)];
                    tsData = zeros(numTot, markRight-markLeft+1);
                    
                    for iPos = markLeft:markRight
                        if mod(iPos,5) == 0 % refresh every 5 steps
                            drawnow
                        end
                        for varIx=focusedSet
                            sg1 = SumGamma1(DEqns_comp{allDEixMap(varIx)}{DE_GAMMA1i},vars(iPos,:),false);
                            tscale = DEqns_comp{allDEixMap(varIx)}{DE_CFACi} ...
                                / sg1(1);
                            if speedSwitch
                                tsData(varIx,iPos-markLeft+1) = abs(1/tscale);
                            else
                                tsData(varIx,iPos-markLeft+1) = abs(tscale);
                            end
                            if get(DSSRTfigHandle,'CurrentCharacter') == 27
                                disp('StepNet:  Cancelling operation')
                                stop_ts = true;
                                break
                            end
                        end
                        if stop_ts
                            break
                        end
                    end
                    
                    % plot data
                    if ~stop_ts
                        tscaleFigHandle = figure('NumberTitle','off','Name', ...
                            [ upper(ts_str{speedSwitch+1}) 's for ' varnames{obsIx} ' etc. : ' figTitleName ' (eps=1/' num2str(dScaleThresh,'%.3f') ')' ], ...
                            'Position', [150, 200, 540, 440]);
                        tscaleAxes = axes('position',[0.1,0.1,0.8,0.8]);
                        hold on
                        styles = {':','-'};
                        ts_temp = GetColStyleSet(length(focusedSet),[0.7,0.2,0.0],styles,0.5);
                        ts_cols = ts_temp{1};
                        ts_styles = ts_temp{2};
                        for varPos=1:length(focusedSet)
                            varIx = focusedSet(varPos);
                            plot(seqTimes(markLeft:markRight),tsData(varIx,:),styles{ts_styles(varPos)},'Color',ts_cols(varPos,:),'LineWidth',2)
                        end
                        hold off
                        if plotLogYScale
                            set(gca,'YScale','log')
                        end
                        axis tight
                        title([ts_str{speedSwitch+1} ' values for focused variable ' varnames{obsIx} ' and associated variables'], ...
                            'Interpreter', 'none')
                        keyFigHandle = figure('NumberTitle','off','Name', ...
                            [ 'KEY: ' upper(ts_str{speedSwitch+1}) 's for ' varnames{obsIx} ' etc. : ' figTitleName ' (eps=1/' num2str(dScaleThresh,'%.3f') ')' ], ...
                            'Position', [660, 240, 200, 220]);
                        keyAxes = axes('position',[0.05, 0.05, 0.9, 0.9]);
                        hold on
                        axis off
                        for varPos=1:length(focusedSet)
                            varIx = focusedSet(varPos);
                            plot([0 1],[varPos-0.9 varPos-0.9],styles{ts_styles(varPos)},'Color',ts_cols(varPos,:),'LineWidth',2)
                            textHandle = text( 0.23, varPos-0.6, [ts_keystr{speedSwitch+1} ' ' varnames{varIx}], 'Interpreter', 'none');
                            set(textHandle,'Color',ts_cols(varPos,:))
                        end
                        hold off
                        axis([0 1 0 varPos-0.45])
                        disp('StepNet:  Figure displayed (colour key in separate window)')
                    end % if stop_ts
                end % if proceed
            end
        case '~' % do plots vs. time over a marked interval:
            % 1) a variable & its asymptotic 'target' value (original and reduced for epochs)
            % 2) all Psi values and input term sizes for the focused variable (option for log plot)
            if existNet
                if markLset & markRset % & freshTS{3} > 0 & ishandle(freshTS{3})
                    focusedIx = [];
                    if obsIxSet
                        defIx = obsIx;
                    else
                        defIx = 1;
                    end
                    if verboseTog
                        disp('StepNet:  Choose variable to plot information for (see dialog box)')
                        disp('           Valid variable names:')
                        for vix=1:numExt
                            if allDEixMap(vix) > 0
                                disp(['             ' varnames{vix}])
                            end
                        end
                    end
                    result = DialogBox('Plot for variable:',varnames{defIx},verboseTog);
                    proceed = false;
                    if ~isempty(result)
                        newVarName = result;
                        [p pIx] = ismember( newVarName, varnames );
                        if p
                            if pIx <= numExt && allDEixMap(pIx) > 0
                                focusedIx = pIx;
                                proceed = true;
                                disp(['           Using external variable ' newVarName])
                            else
                                disp('          Invalid external variable name')
                            end
                        else
                            disp('          Invalid external variable name')
                        end
                    else
                        disp('          Cancelled')
                    end
                    if proceed
                        disp('StepNet:  Plot y-axis on log scale? (see dialog box)')
                        ButtonName = questdlg('Plot on log scale?','Log scale?','Yes','No','No');
                        switch ButtonName
                            case 'Yes'
                                plotLogYScale = true;
                            case 'No'
                                plotLogYScale = false;
                            otherwise
                                proceed = false;
                                result = []; % to create the correct message at end of this case
                        end
                    end

                    if proceed && ~isempty(focusedIx) && DEqnsCompiled %& ismember(focusedIx,freshTS{5})
                        % freshTS = {true, length(TransSeqs), 0, seqTimes(markLeft), TSfocusId};
                        % ( length(TS) is = the TS id # )
%                         thisTS = TransSeqs{freshTS{2}};
%                         could use epochTimesSet, epochTimes in future
	
                        % use markers to get x-axis time interval for plot
                        tPlotInterval = [seqTimes(markLeft), seqTimes(markRight)];
%                         tPosInterval = [markLeft, markRight];
                        qsfpData = zeros(2, markRight-markLeft+1);
                        focEqnIx = allDEixMap(focusedIx);
                        thisDE = DEqns_comp{focEqnIx}; % do this for compiled eqns only

                        % DON'T DO THIS ANY MORE!
                        % do not include terms w/o associated differential
                        % equations in this 'actives' set focActsAbs (e.g. leak current),
                        % because they are trivially retained in a reduced model anyway.
%                         noDiffEqsList = [];
%                         for gamIx = 1:length(thisDE{DE_GAMMA1i})
%                             if allDEixMap( thisDE{DE_GAMMA1i}{gamIx}{DE_GAMVARNAMEi} ) == 0
%                                 noDiffEqsList = [noDiffEqsList, thisDE{DE_GAMMA1i}{gamIx}{DE_GAMVARNAMEi}];
%                             end
%                         end
%                         for gamIx = 1:length(thisDE{DE_GAMMA2i})
%                             if allDEixMap( thisDE{DE_GAMMA2i}{gamIx}{DE_GAMVARNAMEi} ) == 0
%                                 noDiffEqsList = [noDiffEqsList, thisDE{DE_GAMMA2i}{gamIx}{DE_GAMVARNAMEi}];
%                             end
%                         end

                        for iPos = markLeft:markRight
                            if mod(iPos,5) == 0 % refresh every 5 steps
                                drawnow
                            end
                            % get V_inf data from GetQsfpVal for y-axis plot #1
                            qsfpData(1,iPos-markLeft+1) = GetQsfpVal(DEqns_comp{focEqnIx},vars(iPos,:),false);

                            % step through relevant epochs (setting actives switch accordingly) to get sequence of reduced model V_inf 's
                            % ... set actives for this time epoch
                            focActsIxMap = actsIxMap{focusedIx};
                            focActsAbs = focActsIxMap( Eseq{iPos}{ESEQ_ACTIVES_BAR}{focusedIx} );

                            % ensure terms w/o differential eqns are included for the purposes of the local model
                            % DON'T DO THIS ANY MORE!
%                             for termEntry=noDiffEqsList
%                                 if ~ismember(termEntry,focActsAbs)
%                                     focActsAbs = [focActsAbs, termEntry];
%                                 end
%                             end

                            %%%% NEED TO MAKE THIS MORE EFFICIENT!!
                            % updating this for every step rather than using known transition sequence is inefficient
                            DEqns_comp{focEqnIx} = SetActivesSwitchSingle(DEqns_comp{focEqnIx}, focusedIx, focEqnIx, focActsAbs);

                            % ... connect data together for y-axis plot #2
                            qsfpData(2,iPos-markLeft+1) = GetQsfpVal(DEqns_comp{focEqnIx},vars(iPos,:),true);
                            
                            % include variable data for comparison
                            qsfpData(3,iPos-markLeft+1) = vars(iPos,focusedIx);
                        end

                        % plot 'targets' data
                        asympFigHandle = figure('NumberTitle','off','Name', ...
                            [ 'Orbit, targets for ' varnames{focusedIx} ' : ' figTitleName ' (eps=1/' num2str(dScaleThresh,'%.3f') ')' ]);
                        hold on
                        plot(seqTimes(markLeft:markRight),qsfpData(1,:),'-g')
                        plot(seqTimes(markLeft:markRight),qsfpData(2,:),'-r')
                        plot(seqTimes(markLeft:markRight),qsfpData(3,:),'-b')
                        hold off
                        axis tight
                        if diagSw == PSIS
                            redstr = 'Psis';
                        else
                            redstr = 'terms';
                        end
                        title(['Orbit (b), full (g) and reduced-by-' redstr ' (r) asymptotic values for variable ' varnames{focusedIx}], ...
                            'Interpreter', 'none')
                        disp('StepNet:  Targets figure displayed')
                        
                        % plot Psi data
                        % domInfo{evar,pos} = { E_b, ratios_act, ixa_dom, sort_act_dom, test_act_dom, pot_ix };
                        focusedSet = candActsIxMapAbs{focusedIx};
                        if length(focusedSet) < 2
                            proceed = false; % default value if not enough candidate actives to plot a graph
                        else
                            proceed = true; % default value
                        end
                        if proceed
                            Psi_data = zeros(numTot,markRight-markLeft); % influence strengths
                            term_data = zeros(numTot,markRight-markLeft); % individual additive terms on rhs (currents for neural systems)
                            for iPos = markLeft:markRight
                                if mod(iPos,5) == 0 % refresh every 5 steps
                                    drawnow
                                end
                                if isempty(domInfo_Psi{focusedIx,iPos}{1}) || isempty(domInfo_term{focusedIx,iPos}{1})
                                    disp('StepNet:  Focused variable has no candidate `active` inputs')
                                    proceed = false;
                                    break
                                end
                                % Psi values
                                mostDom_Psi = domInfo_Psi{focusedIx,iPos}{4};
                                ixa_dom_Psi = domInfo_Psi{focusedIx,iPos}{3};
                                mostDomIx_Psi = candActsIxMapAbs{focusedIx}(ixa_dom_Psi(1));
                                Psi_data(mostDomIx_Psi,iPos-markLeft+1) = mostDom_Psi;
                                ratios_act_Psi = domInfo_Psi{focusedIx,iPos}{2};
                                for i=1:length(ixa_dom_Psi)-1
                                    Psi_data(candActsIxMapAbs{focusedIx}(ixa_dom_Psi(i+1)),iPos-markLeft+1) = mostDom_Psi/ratios_act_Psi(i);
                                end
                                % term sizes (e.g. transmembrane "currents")
                                mostDom_term = domInfo_term{focusedIx,iPos}{4};
                                ixa_dom_term = domInfo_term{focusedIx,iPos}{3};
                                mostDomIx_term = candActsIxMapAbs{focusedIx}(ixa_dom_term(1));
                                term_data(mostDomIx_term,iPos-markLeft+1) = mostDom_term;
                                ratios_act_term = domInfo_term{focusedIx,iPos}{2};
                                for i=1:length(ixa_dom_term)-1
                                    term_data(candActsIxMapAbs{focusedIx}(ixa_dom_term(i+1)),iPos-markLeft+1) = mostDom_term/ratios_act_term(i);
                                end
                            end
                        end % if proceed
                        if proceed
                            maxPsi = max( max( Psi_data ));
                            maxTerm = max( max( term_data ));
                            styles = {'-'};
                            psi_cols_temp = GetColStyleSet(length(focusedSet),[0.7,0.2,0.0],styles,0.5);
                            psi_cols = psi_cols_temp{1};
                            psi_styles = psi_cols_temp{2};
                            PsiFigHandle = figure('NumberTitle','off','Name', ...
                                [ 'Influence strengths for ' varnames{focusedIx} ' : ' figTitleName]);
                            if plotLogYScale
                                set(gca,'YScale','log')
                            end
                            hold on
                            for varPos=1:length(focusedSet)
                                plot(seqTimes(markLeft:markRight),Psi_data(focusedSet(varPos),:),'-','Color',psi_cols(varPos,:), 'LineWidth', 2) % styles{psi_styles(varPos)}
                            end
                            hold off
                            axis tight
                            title(['Influence strengths for focused variable ' varnames{focusedIx}], 'Interpreter', 'none')
                            TermFigHandle = figure('NumberTitle','off','Name', ...
                                [ 'Term input sizes for ' varnames{focusedIx} ' : ' figTitleName]);
                            if plotLogYScale
                                set(gca,'YScale','log')
                            end
                            hold on
                            for varPos=1:length(focusedSet)
                                plot(seqTimes(markLeft:markRight),term_data(focusedSet(varPos),:),'-','Color',psi_cols(varPos,:), 'LineWidth', 2)
                            end
                            hold off
                            axis tight
                            title(['Term input sizes for focused variable ' varnames{focusedIx}], 'Interpreter', 'none')

                            % Colour key
                            keyFigHandle = figure('NumberTitle','off','Name', ...
                                [ 'KEY: Influence strengths for ' varnames{obsIx} ' : ' figTitleName ], ...
                                'Position', [660, 240, 200, 220]);
                            keyAxes = axes('position',[0.05, 0.05, 0.9, 0.9]);
                            hold on
                            axis off
                            for varPos=1:length(focusedSet)
                                varIx = focusedSet(varPos);
                                plot([0 1],[varPos-0.9 varPos-0.9],'-','Color',psi_cols(varPos,:),'LineWidth',2) % styles{psi_styles(varPos)}
                                textHandle = text( 0.23, varPos-0.6, varnames{varIx}, 'Interpreter', 'none' );
                                set(textHandle,'Color',psi_cols(varPos,:))
                            end
                            hold off
                            axis([0 1 0 varPos-0.45])
                            disp('StepNet:  Influence strength figure displayed (colour key in separate window)')
                        else
                            disp('StepNet:  Psi data collection not possible: no candidate actives for chosen variable?')
                        end % if proceed
                    else
%                         if ~ismember(focusedIx,freshTS{5}) % then TS does not include the relevant data for this focusedIx
%                             disp('StepNet:  Chosen variable of focus is not included in the focus set of the freshest transition sequence')
                        if ~proceed && ~isempty(result) % then it was a different reason for getting here
                            beep
                            disp('StepNet:  Internal error with compilation of differential equations or variable of focus')
                        elseif verboseTog
                            disp('          Cancelled')
                        end
                    end
                else
                    disp('StepNet:  Marks not set') % or fresh transition sequence not present') % Use of TS not yet implemented
                end
            else
                disp('StepNet:  No functional network held in memory or current observable of focus not set')
            end
        case '@' % find relative scales of influence for all inputs to a given observable at current time position
            if existNet && obsIxSet
                obsName = varnames{obsIx}; % obsIx is an external variable
                if isempty(domInfo{obsIx,pos}{3}) % pick any of the entries to test for emptiness
                    disp('StepNet:  No candidate `active` inputs are specified for the chosen observable')
                else
                    if verboseTog
                        if diagSw == PSIS
                            fprintf('StepNet:  Relative scales of influence for all inputs to observable %s:\n',obsName)
                        else
                            fprintf('StepNet:  Relative sizes of input terms for all inputs to observable %s:\n',obsName)
                        end
                    end
                    if length(domInfo{obsIx,pos}{3}) > 1
                        if diagSw == PSIS
                            typestr = 'influence strength';
                        else
                            typestr = 'input size';
                        end
                        % Show actives
                        if verboseTog
                            fprintf('           Order/scaling of candidate `actives` at pos %3i (time %.4f) is:\n',pos,seqTimes(pos))
                        else
                            fprintf('StepNet:  Order/scaling of candidate `actives` for `%s` at pos %i (time %.4f) is:\n',obsName,pos,seqTimes(pos))
                        end
                        fprintf('           Rank 1 %s value = %.5f\n', typestr, domInfo{obsIx,pos}{4}) 
                        fprintf('           Rank  |  Input Var  |  Fraction of rank #1`s %s\n', typestr)
                        for actIx = 1:length(domInfo{obsIx,pos}{3})
                            actName = varnames{  actsIxMap{obsIx}( domInfo{obsIx,pos}{3}(actIx) )  };
                            if actIx == 1
                                actRatStr = ''; % only used for actIx > 1
                            else
                                actRatStr  = [ '                  1/' num2str(domInfo{obsIx,pos}{2}(actIx-1),'%7.3f') ];
                            end
                            fprintf(['           %3i   :   %s' actRatStr '\n' ],actIx,actName);
                        end
                        % Now show potentials
                        numActs = length(actsIxMap{obsIx});
                        numPots = length(domInfo{obsIx,pos}{6});
                        potNames = cell(1,numPots);
                        fprintf('           %s values of `potentials`:\n', typestr)
                        fprintf(    '           Pots. | Input  ')
                        for var = 1:numActs-1
                            fprintf([ varnames{actsIxMap{obsIx}(var)} '       ']) % this number of spaces for length 3 var names
                        end
                        fprintf([ varnames{actsIxMap{obsIx}(numActs)} '\n'])
                        for varIx = 1:numPots
                            maxDom = [0,0];
                            potNames{varIx} = varnames{ actsIxMap{obsIx}( domInfo{obsIx,pos}{6}(varIx) ) };
                            fprintf('            %s        ', potNames{varIx})
                            for i=1:numActs
                                domVal = domInfo{obsIx,pos}{5}{varIx}(i);
                                if domVal > maxDom(2)
                                    maxDom = [i, domVal];
                                end
                                if strcmp(potNames{varIx}, varnames{actsIxMap{obsIx}(i)})
                                    thisVarDom = domVal; % then this is the Psi value for the current potential active variable
                                end
                                fprintf('%8.4f  ',domVal)
                            end
                            fprintf('-- 1/%.4f of max %s `%s`\n',maxDom(2)/thisVarDom,typestr,varnames{actsIxMap{obsIx}(maxDom(1))})
                        end
                    else
                        fprintf('StepNet:  Input %s is the only active input to the chosen observable at this time\n',varnames{actsIxMap{obsIx}(1)})
                    end                    
                end
            else
                disp('StepNet:  No functional network held in memory or current observable of focus not set')
            end
        case 'W' % write functional net data to file
            if existNet
                if FNdataFnameSet
                    result = DialogBox('WRITE: Func net data file name',EndStripStr(FNdataFilenameDisp,'_FN',true,7),verboseTog);
                else
                    result = DialogBox('WRITE: Func net data file name',[],verboseTog);
                end
                if ~isempty(result)
                    result = EndStripStr(result,'_FN',true,7);
                    FNwrite = true; % initial value
                    if exist([ networkObjectDir '/' result '_FN.mat'],'file') == 2
                        ButtonName = questdlg('File exists. Overwrite?');
                        switch ButtonName
                            case 'Yes'
                                FNwrite = true;
                            otherwise
                                FNwrite = false;
                        end
                    end
                    if FNwrite
                        FNdataFilenameDisp = [result '_FN.mat'];
                        FNdataFname        = [ networkObjectDir '/' FNdataFilenameDisp ];
                        FNdataFnameSet     = true;
                        funcNetDataSave = { diag Eseq numTot whosevars candActsIxMapAbs candPotsIxMapAbs ...
                                varnames numSteps timeStep domGraph domInfo CFGfilename varDataFilename markLeft ...
                                markRight startTime stopTime availStTime availSpTime ...
                                dScaleThresh obsIx varViewIx varBounds seqTimesAll varsAll seqTimes vars ...
                                diag_term Eseq_term domGraph_term domInfo_term};
                        refresh(DSSRTfigHandle)
                        save(FNdataFname, 'funcNetDataSave');
                        clear funcNetDataSave ; % clear this variable, so that can free memory
                                                % AND be able to check if a load created this
                                                % variable successfully
                        disp('StepNet:  FN data file written')
                    else
                        if verboseTog
                            disp('           Cancelled - FN data file not written')
                        end
                    end
                else
                    if verboseTog
                        disp('           Cancelled')
                    end
                end
            else
                disp('StepNet:  No functional network held in memory!')
            end
        case 'E' % read functional net data from file
            readFNokay = false;
            if FNdataFnameSet
                result = DialogBox('READ: Func. net. data file name',EndStripStr(FNdataFilenameDisp,'_FN',true,7),verboseTog);
            else
                result = DialogBox('READ: Func. net. data file name','demo',verboseTog); % TEMPORARY
            end
            if ~isempty(result)
                result = EndStripStr(result,'_FN',true,7);
                if exist([ networkObjectDir '/' result '_FN.mat'],'file') == 2
                    FNdataFilenameDisp = [result '_FN.mat'];
                    FNdataFname    = [ networkObjectDir '/' FNdataFilenameDisp ];
                    FNdataFnameSet = true;
                    refresh(DSSRTfigHandle)
                    load(FNdataFname);
                    if exist('funcNetDataSave','var')~=1 % then file didn't have correct variable name inside
                        fprintf('StepNet:  Unreadable FN data file %s\n',FNdataFilenameDisp)
                        readFNokay = false;
                    else
                        disp('StepNet:  FN data file accessed')
                        readFNokay = true;
                    end
                else
                    fprintf('StepNet:  Func. net. data file name %s_FN.mat does not exist\n',result)
                end
            else
                if verboseTog
                    disp('           Cancelled')
                end
            end
            if readFNokay
                if length(funcNetDataSave) ~= 31
                    disp('StepNet:  Corrupt FN data file. Wrong number of received parameter values')
                else
					diag_Psi         = funcNetDataSave{1};
					Eseq_Psi         = funcNetDataSave{2};
					numTot           = funcNetDataSave{3};
					whosevars        = funcNetDataSave{4};
					candActsIxMapAbs = funcNetDataSave{5};
					candPotsIxMapAbs = funcNetDataSave{6};
					varnames         = funcNetDataSave{7};
                    numSteps         = funcNetDataSave{8};
                    timeStep         = funcNetDataSave{9};
                    domGraph_Psi     = funcNetDataSave{10};
                    domInfo_Psi      = funcNetDataSave{11};
                    % These next entries not present in regular funcNetData returned
                    % from FuncNet.m
                    CFGfilename      = funcNetDataSave{12}; % original name when data created
                    varDataFilename  = funcNetDataSave{13}; % original name when data created
                    CFGfilenameDisp  = RootStripStr(CFGfilename);
                    varDataFilenameDisp = RootStripStr(varDataFilename);
                    CFGfilenameSet      = true;
                    varDataFilenameSet  = true;
                    markLeft         = funcNetDataSave{14}; % recover markers if they existed when saved
                    markRight        = funcNetDataSave{15};
                    startTime        = funcNetDataSave{16};
                    stopTime         = funcNetDataSave{17};
                    startTimeSet     = true;
                    stopTimeSet      = true;
                    availStTime      = funcNetDataSave{18};
                    availSpTime      = funcNetDataSave{19};
                    dScaleThresh      = funcNetDataSave{20};
                    if obsIxSet && obsIx ~= funcNetDataSave{21} % changed obsIx
                        attEstParsSet = false; % reset to get correct defaults for an AttEst() call
                    end
                    obsIx            = funcNetDataSave{21};
                    obsIxSet         = true;
                    varViewIx        = funcNetDataSave{22};
                    varBounds        = funcNetDataSave{23};
                    seqTimesAll      = funcNetDataSave{24};
                    varsAll          = funcNetDataSave{25};
                    seqTimes         = funcNetDataSave{26};
                    vars             = funcNetDataSave{27};
                    diag_term        = funcNetDataSave{28};
					Eseq_term        = funcNetDataSave{29};
                    domGraph_term    = funcNetDataSave{30};
                    domInfo_term     = funcNetDataSave{31};
                    % if you add to these, ensure that you change the # of expected pars in the `if` statement above
                    if diagSw == PSIS
                        diag = diag_Psi;
                        Eseq = Eseq_Psi;
                        domGraph = domGraph_Psi;
                        domInfo = domInfo_Psi;
                    else
                        diag = diag_term;
                        Eseq = Eseq_term;
                        domGraph = domGraph_term;
                        domInfo = domInfo_term;
                    end
                    diag_state_old = zeros(1, numobjects);  % reset
                    clear funcNetDataSave
                    if markLeft > 0
                        markLset = true;
                        markLchanged = true;
                    else
                        markLset = false;
                        markLchanged = true; % to make sure old markers get removed from display
                    end
                    if markRight > 0
                        markRset = true;
                        markRchanged = true;
                    else
                        markRset = false;
                        markRchanged = true; % to make sure old markers get removed from display
                    end
                    existFNpars = CFGfilenameSet & varDataFilenameSet & nOutSet & startTimeSet & stopTimeSet & dScaleThreshSet;
                    if numobjects ~= length(diag{1})
                        disp('StepNet:  Fatal error loading func net data:')
                        disp('          Number of graphical objects does not match dimension of diagram')
                        disp('            state variable found in data file')
                        disp('          Failed to create a diagram state')
                    else % continue to clear workspace
                        attractorEst.exist = false;
                        attractorEst.number = 0;
                        attractorEst.attEstList = {};
                        if numVBs ~= length(varnames)
                            disp('StepNet:  Number of vertical bars defined exceeds number of defined variables')
                            disp('            according to size of variables array found in data file.')
                            disp('          Failed to create a diagram state')
                        else
                            pos = 1; % reset this
                            if ~existNet % then initialize new graphical state
                                for obj = 1:numobjects
                                    if hobj{obj}{OBJ_TYPE} == LINK
                                        if hobj{obj}{OBJ_STATE} ~= LINKOFFSTATE
                                            hobj{obj}{OBJ_STATE} = diag{pos}(obj); % refresh initialize diagram states
                                        end % else don't change it
                                    else
                                        hobj{obj}{OBJ_STATE} = diag{pos}(obj); % refresh initialize diagram states
                                    end
                                    hobj{obj}{OBJ_HANDLE} = DrawNet( hobj{obj} );
								end
                                existNet = true;
                            end
                            for vn = 1:numVBs % update vertical bar objects with their variable names
                                hvb{vn}{8} = varnames{vn}; % entry 8 is the name
                            end
                            % No longer reset varViewIx in V1.5 -- now loaded with FN file
                            %% display first variable in VarViewBox
                            %varViewIx = 1; % reset this
                            PlotVBox(handVVaxes,seqTimes,vars(:,varViewIx),varViewIx,varnames)
                            % update trans. seq. box
                            axes(handTSaxes)
                            cla
                            if freshTS{1}
                                if ishandle(freshTS{3})
                                    delete(freshTS{3})
                                end
                                freshTS{3} = 0;
                                freshTS{4} = 0;
                                freshTS{1} = false;
                                epochTimesSet = false;
                            end
                            set(handTSaxes, 'XLim', [ seqTimes(1) seqTimes(numSteps)])
                            % update axes etc. for Time Box
                            axes(handTBaxes)
                            delete(tBarHand)
                            set(handTBaxes, 'XLim', [ seqTimes(1) seqTimes(numSteps)])
                            tBarHand = Draw_tBar(seqTimes(1),seqTimes(numSteps),'k');
                            axes(handFNaxes) % return to FN axes

                            refreshDiagState = true; % force update of display
                            forceNetRedraw = true;
                            figTitleName = [ CFGfilenameDisp ' : ' varDataFilenameDisp ];
                            disp(['StepNet:  New functional network appearing: `' figTitleName '`'])
                            set(DSSRTfigHandle,'Name',[ figTitle1 figTitleName])
                        end
                    end % if numobjects
                end % if length...
            end % if readFNokay
        case 'I' % redisplay reduced regime info and where possible, plot phase diagrams for reduced dynamical regimes
            if existRegime
                if verboseTog
                    disp('StepNet:  Re-displaying reduced dynamical regime information')
                end
                DisplayRegimes(regimeStruct,verboseTog)
                % plot graphs
                disp('StepNet:  Plotting of phase diagrams not yet available')
            else
                disp('StepNet:  No reduced dynamical regimes available')
            end
        case 'O' % specify data file name
            if varDataFilenameSet
                if existNet
                    disp('StepNet:  WARNING -- current functional network will be lost')
                    beep
                end
                result = DialogBox('Data file name',EndStripStr(varDataFilenameDisp,'.',true),verboseTog);
            else
                result = DialogBox('Data file name',[],verboseTog);
            end
            if ~isempty(result)
                result = EndStripStr(result,'.',true);
                if exist([ networkObjectDir '/' result '.dat'],'file') == 2
                    varDataFilename     = [networkObjectDir '/' result '.dat'];
                    varDataFilenameDisp = [result '.dat'];
                    varDataFilenameSet = true;
                    fprintf('StepNet:  Data file name set to %s.dat\n',result)
                    if existNet
                        clearFNet = true;
                    end
                    if CFGfilenameSet
                        resetView = false;
                    end
                    attractorEst.exist = false;
                    attractorEst.number = 0;
                    attractorEst.attEstList = {};

                    resultVAT = GetVarsAndTimes([DSSRT_ROOTPATH networkObjectDir],varDataFilename,CFGfilename,CFGfilenameSet,verboseTog);
                    if isempty(resultVAT)
                        disp('StepNet:  GetVarsAndTimes returned nothing. No changes were made')
                        proceed = false;
                    else
                        proceed = true;
                    end
                    
                    if proceed
                        if ~isempty(resultVAT{4}) % a way to test that these things got initialized in GetVarsAndTimes()
                            if resultVAT{2} ~= 0 % available stop time == 0 means something messed up
                                availStTime = resultVAT{1};
                                availSpTime = resultVAT{2};
                                varsAll     = resultVAT{5};
                                seqTimesAll = resultVAT{6};
                                vars        = varsAll(1:nOut:length(seqTimesAll),:);
                                seqTimes    = seqTimesAll(1:nOut:length(seqTimesAll));
                            end
                            if resultVAT{2} == 0
                                disp('StepNet:  Warning. Available stop time was zero. New variable data file was not used')
                            end
                        else
                            beep
                            disp('StepNet:  Problem getting initialisations back from GetVarsAndTimes()')
                            disp('           Bad data file?')
                            if varDataFilenameSet
                                varDataFilename = [ networkObjectDir '/' varDataFilenameDisp ]; % way to restore to previous name
                            else
                                varDataFilename = ''; % reset these
                                startTimeSet = false;
                                stopTimeSet  = false;
                            end
                        end
                        if CFGfilenameSet % only provided seqTimes etc. is defined
                            PlotVBox(handVVaxes,seqTimes,vars(:,varViewIx),varViewIx,varnames)
                            axes(handFNaxes)
                        end
                    end % if proceed
                else
                    fprintf('StepNet:  Data file name %s.dat does not exist\n',result)
                end
            else
                if verboseTog
                    disp('           Cancelled')
                end
            end
            existFNpars = CFGfilenameSet & varDataFilenameSet & nOutSet & startTimeSet & stopTimeSet & dScaleThreshSet;
        case 'P' % read parameters file
            result = DialogBox('Parameters file name (for reading)','defaults',verboseTog);
            if ~isempty(result)
                result = EndStripStr(result,'.',true);
                if exist([ networkObjectDir '/' result '.par'],'file') == 2
                    parsFilename = [result '.par'];
                    pars = ReadParamsFile([ networkObjectDir '/' parsFilename ]);
                    if ~isempty(pars)
                        if ~CFGfilenameSet % don't change CFGfilename from params file if already set
                            use_new_cfg = true; % but use new hobj and hvb objects created below (for the first time)
                            CFGfilename        = [ networkObjectDir '/' pars{1} ];
                            CFGfilenameDisp    = pars{1};
                            CFGfilenameSet     = true;
                        else
                            use_new_cfg = false;
                            disp('StepNet:  CFG filename already set. This is not updated by reading')
                            disp('            a parameters .par file')
                        end
                        newVDfilename = [ networkObjectDir '/' pars{2} ];
                        if varDataFilenameSet
                            if strcmp( varDataFilename, newVDfilename )
                                use_new_vars = false;
                            else
                                varDataFilename = newVDfilename;
                                use_new_vars = true;
                            end
                        else
                            varDataFilename = newVDfilename;
                            use_new_vars = true;
                        end
                        nOut           = pars{6};
                        nOutSet        = true;
                        attEstParams   = pars{11};
                        attEstParsSet  = true;
                        % Not sure how the following helps! So I took it back out...
                        %  (It creates an error when loading pars file from an empty memory)
                        % tempCFGSet = false; % force this so that vars etc. don't get overwritten if CFGfilenameSet == true
                            % doesn't matter anyway if CFGfilenameSet really is false already.
%                         resultVAT = GetVarsAndTimes([DSSRT_ROOTPATH networkObjectDir],varDataFilename,CFGfilename,tempCFGSet,verboseTog);
                        resultVAT = GetVarsAndTimes([DSSRT_ROOTPATH networkObjectDir],varDataFilename,CFGfilename,CFGfilenameSet,verboseTog);
                        if isempty(resultVAT)
                            disp('StepNet:  GetVarsAndTimes returned nothing. No changes were made')
                            proceed = false;
                        else
                            proceed = true;
                        end

                        if proceed
                            if ~isempty(resultVAT{4}) % a way to test that these things got initialized in GetVarsAndTimes()
                                if use_new_vars && resultVAT{2} ~= 0 % available stop time == 0 means something messed up
                                    availStTime = resultVAT{1};
                                    availSpTime = resultVAT{2};
                                    varsAll     = resultVAT{5};
                                    seqTimesAll = resultVAT{6};
                                    vars        = varsAll(1:nOut:length(seqTimesAll),:);
                                    seqTimes    = seqTimesAll(1:nOut:length(seqTimesAll));
                                    % don't reset this until the end (when VVbox may be redisplayed)
                                end
                                if resultVAT{2} == 0
                                    disp('StepNet:  Warning. Available stop time was zero. New variable data file was not used')
                                end
                                if use_new_cfg
                                    numTot      = resultVAT{3};
                                    varnames    = resultVAT{4};
                                    numInt      = resultVAT{7};
                                    numExt      = resultVAT{8};
                                    actsIxMap   = resultVAT{9};
                                    potsIxMap   = resultVAT{10};
                                    inputsIx    = resultVAT{11};
                                    hobj        = resultVAT{12};
                                    hvb         = resultVAT{13};
                                    DEqns       = resultVAT{14};
                                    DEpars      = resultVAT{15};
                                    varBounds   = resultVAT{16};
                                    use_new_cfg = false; % reset
                                end
                            else
                                beep
                                disp('StepNet:  Problem getting initialisations back from GetVarsAndTimes()')
                                disp('           Bad data file?')
                                if varDataFilenameSet
                                    varDataFilename = [ networkObjectDir '/' varDataFilenameDisp ]; % way to restore to previous name
                                else
                                    varDataFilename = ''; % reset these
                                    startTimeSet = false;
                                    stopTimeSet  = false;
                                end
                                use_new_vars = false; % reset
                            end
                            varDataFilenameDisp = pars{2};
                            varDataFilenameSet  = true;
                            startTime           = pars{3};
                            stopTime            = pars{4};
                            if ~ ( startTime >= availStTime && startTime < min(stopTime,availSpTime) )
                                startTime = availStTime;
                            end
                            startTimeSet        = true;
                            if ~ ( stopTime <= availSpTime && stopTime > max(startTime, availStTime) )
                                stopTime = availSpTime;
                            end
                            stopTimeSet         = true;
                            dScaleThresh        = pars{5};
                            dScaleThreshSet     = true;
                            ignoreActsPotsOrder = pars{7};
                            guessPeriod         = pars{8};
                            marginPC            = pars{9};
                            tScaleThresh        = pars{10};
                            tScaleThreshSet     = true;
                            % If using new vars, re-draw current variable's data (don't need to
                            %   reset to 1 since haven't changed system)
                            if CFGfilenameSet && use_new_vars % only provided seqTimes etc. is defined
                                use_new_vars = false; % reset
                                PlotVBox(handVVaxes,seqTimes,vars(:,varViewIx),varViewIx,varnames)
                                axes(handFNaxes)
                            end
                        else
                            use_new_vars = false; % reset
                        end
                    end
                else
                    fprintf('StepNet:  Parameters file name %s.par does not exist\n',result)
                end
            else
                if verboseTog
                    disp('           Cancelled')
                end
            end
            existFNpars = CFGfilenameSet & varDataFilenameSet & nOutSet & startTimeSet & stopTimeSet & dScaleThreshSet;
        case '=' % list parameters specified for FuncNet, periodic cycle shooting, and continuation
            notSetStr = 'NOT SET';
            if CFGfilenameSet,   CFGstr = ['`' CFGfilenameDisp '`'];        else CFGstr = notSetStr; end
            if varDataFilenameSet,  datStr = ['`' varDataFilenameDisp '`']; else datStr = notSetStr; end
            if nOutSet,          nOutStr = num2str(nOut);               else nOutStr = notSetStr; end
            if availStTime>=0,   aStTstr = num2str(availStTime);        else aStTstr = notSetStr; end
            if availSpTime>0,    aSpTstr = num2str(availSpTime);        else aSpTstr = notSetStr; end
            if startTimeSet,     StTmStr = num2str(startTime,'%.2f');   else StTmStr = notSetStr; end
            if stopTimeSet,      SpTmStr = num2str(stopTime,'%.2f');    else SpTmStr = notSetStr; end
            if dScaleThreshSet,   dsThStr = num2str(dScaleThresh,'%.4f'); else dsThStr = notSetStr; end
            if tScaleThreshSet,   tsThStr = num2str(tScaleThresh,'%.4f'); else tsThStr = notSetStr; end
            if attEstParsSet
                attDelStr = num2str(attEstParams.delta,'%.4f');
                attdVStr = num2str(attEstParams.dVmax,'%.4f');
                attVsStr = num2str(attEstParams.Vinres,'%.4f');
                attVpStr = num2str(attEstParams.Vpert,'%.4f');
                attsRStr = num2str(attEstParams.searchRes,'%2i');
                attlrMStr = num2str(attEstParams.lowResMultiple,'%2i');
                attdTStr = num2str(attEstParams.derivThresh,'%.4f');
            end
            if existNet
                if obsIxSet,  ObsIStr = ['#' num2str(obsIx) ' = `' varnames{obsIx} '`']; else ObsIStr = notSetStr; end
                VarVStr = ['#' num2str(varViewIx) ' = `' varnames{varViewIx} '`'];
            end
            disp(       'StepNet:  Main parameters:')
            spacesStr = '          ';
            fprintf([spacesStr '  CFG file name               = %s\n'], CFGstr);
            fprintf([spacesStr '  data file name              = %s\n'], datStr);
            fprintf([spacesStr '  available start time        = %s\n'], aStTstr);
            fprintf([spacesStr '  available stop time         = %s\n'], aSpTstr);
            fprintf([spacesStr '  start time                  = %s\n'], StTmStr);
            fprintf([spacesStr '  stop time                   = %s\n'], SpTmStr);
            fprintf([spacesStr '  dominance scale threshold   = 1/%s\n'], dsThStr);
            fprintf([spacesStr '  time scale threshold        = 1/%s\n'], tsThStr);
            fprintf([spacesStr '  nOut                        = %s\n'], nOutStr);
            disp([spacesStr 'Parameters for periodic state shooting:'])
            if ignoreActsPotsOrder,  iAO_str = 'ON';  else  iAO_str = 'OFF';  end
            fprintf([spacesStr '  ignore actives/pots order   = %s\n'], iAO_str);
            fprintf([spacesStr '  guess period                = %4.4f\n'], guessPeriod);
            fprintf([spacesStr '  margin p.c. of period       = %4.2f\n'], marginPC);
            if existNet
                disp([spacesStr 'Parameters for running AttEst:'])
                if attEstParsSet
                    fprintf([spacesStr '  delta                  = ' attDelStr ' (currently unused)\n'])
                    fprintf([spacesStr '  dVmax                  = ' attdVStr '\n'])
                    fprintf([spacesStr '  Vinres                 = ' attVsStr '\n'])
                    fprintf([spacesStr '  Vpert                  = ' attVpStr '\n'])
                    fprintf([spacesStr '  search resolution      = ' attsRStr '\n'])
                    fprintf([spacesStr '  low res multiple       = ' attlrMStr '\n'])
                    fprintf([spacesStr '  derivative threshold   = ' attdTStr '\n'])
                else
                    fprintf([spacesStr '  NONE SET\n'])
                end
                disp([spacesStr 'Functional Network parameters (saved in _FN files)'])
                disp([spacesStr ' Current observable of focus:   ' ObsIStr])
                disp([spacesStr ' Current variable data to view: ' VarVStr])
            end
        case '+' % write out parameters to file, provided all have been set!
            if CFGfilenameSet && varDataFilenameSet && nOutSet && startTimeSet && stopTimeSet ...
                    && dScaleThreshSet && tScaleThreshSet && attEstParsSet
                result = DialogBox('Parameters file name (for writing)',[],verboseTog);
                if ~isempty(result)
                    result = EndStripStr(result,'.',true);
                    parFilename = [ networkObjectDir '/' result '.par'];
                    if exist(parFilename,'file')==2
                        ButtonName = questdlg('File exists. Overwrite?');
                        switch ButtonName
                            case 'Yes'
                                parWrite = true;
                            otherwise
                                parWrite = false;
                        end
                    else
                        parWrite = true;
                    end
                    if parWrite
                        wrotePars = WriteParamsFile( parFilename, ...
                            { CFGfilenameDisp, varDataFilenameDisp, startTime, stopTime, dScaleThresh, nOut, ...
                              ignoreActsPotsOrder, guessPeriod, marginPC, tScaleThresh, attEstParams });
                    else
                        wrotePars = 0;
                    end
                else
                    wrotePars = 0;
                    if verboseTog
                        disp('           Cancelled')
                    end
                end
                if wrotePars
                    fprintf('StepNet:  Wrote parameter file %s.par\n',result)
                else
                    disp('           Did not write parameter file')
                end
            else
                disp('StepNet:  Cannot write to parameters file. Some parameters not')
                disp('            yet specified. Press `=` to see which')
            end
        case ';' % jump to left marker
            if existNet
                if markLset
                    if pos ~= markLeft
                        pos = markLeft;
                        refreshDiagState = true;
                        if verboseTog
                            disp('StepNet:  Jumped to left marker')
                        end
                    end
                else
                    if verboseTog
                        disp('StepNet:  Cannot jump -- no left mark set')
                    end
                end
            else
                disp('StepNet:  No functional network calculated!')
            end
        case ':' % jump to right marker
            if existNet
                if markRset
                    if pos ~= markRight
                        pos = markRight;
                        refreshDiagState = true;
                        if verboseTog
                            disp('StepNet:  Jumped to right marker')
                        end
                    end
                else
                    if verboseTog
                        disp('StepNet:  Cannot jump -- no right mark set')
                    end
                end
            else
                disp('StepNet:  No functional network calculated!')
            end
        case '<' % jump to beginning
            if existNet
                if pos ~= 1
                    pos = 1;
                    refreshDiagState = true;
                    if verboseTog
                        disp('StepNet:  Jumped to beginning position')
                    end
                end
            else
                disp('StepNet:  No functional network calculated!')
            end
        case '>' % jump to end
            if existNet
                if pos ~= numSteps
                    pos = numSteps;
                    refreshDiagState = true;
                    if verboseTog
                        disp('StepNet:  Jumped to end position')
                    end
                end
            else
                disp('StepNet:  No functional network calculated!')
            end
        case '%' % draw graph of dominances for current observable
            if existNet
                if obsIxSet
                    if sum(sum(domGraph(obsIx,:,:)))~=0 % then non-zero entries => proceed
                        graphLen = length(domGraph(obsIx,:,:)); % number of time-steps present
                        DG_TIME = 1; % indices
                        DG_DOM  = 2;
                        DG_VAR  = 3;
                        dGYlim = 7;
                        domVar = domGraph(obsIx,1,DG_VAR); % look at dominant variable at time step 1
                        domVarChanged = true;
                        graphDataX = zeros(1,graphLen);
                        graphDataY = zeros(1,graphLen);
                        textXY = {};
                        domFig = figure('NumberTitle','off','Name',[ figTitleName ': "Dominance factor" of most dominant inputs to observable ' varnames{obsIx}]);
                        hold on
                        for x = 1:graphLen
                            graphDataX(x) = domGraph(obsIx,x,DG_TIME);
                            graphDataY(x) = domGraph(obsIx,x,DG_DOM);
                            if domGraph(obsIx,x,DG_VAR) ~= domVar
                                domVarChanged = true; % so display new dominant variable's name on plot
                                domVar = domGraph(obsIx,x,DG_VAR);
                            end
                            if domVarChanged
                                domVarChanged = false;
                                textXY = { textXY{:}, [graphDataX(x), x, domGraph(obsIx,x,DG_VAR)] }; % second entry used later when plotting text!
                            end
                        end
                        plot(graphDataX,graphDataY,'b-')
                        for textIx = 1:length(textXY)
                            text( textXY{textIx}(1), min(graphDataY(textXY{textIx}(2)+4)+2*rand,0.5*dGYlim), ...
                                varnames{textXY{textIx}(3)}, 'Color', 'r', 'Interpreter', 'none')
                        end
                        axis([domGraph(obsIx,1,DG_TIME) domGraph(obsIx,graphLen,DG_TIME) 1 dGYlim])
                        title(['Dominance graph for observable ' varnames{obsIx}], 'Interpreter', 'none')
                        if verboseTog
                            disp('StepNet:  y-axis limits for the dominance graph can be set manually using the')
                            disp('            select tool in the graph`s menu bar')
                        end
                    else
                        disp( 'StepNet:  A dominance graph cannot be computed for observables with one or fewer candidate actives')
                    end
                else
                    disp('StepNet:  You must specify an observable first')
                end
            else
                disp('StepNet:  No functional network calculated!')
            end
        case 'Q' % DUMMY
            % ignore here (need to catch so that not trapped by `otherwise` statement)
        case 'a' % DUMMY used for initial value of keypress on DSSRT loading
            % ignore here
        otherwise
            if verboseTog
                if ischar(keypress)
                    fprintf('StepNet:  Command `%s` not recognised. Press `H` for help or `Q` to quit!\n', keypress)
                else
                    fprintf('StepNet:  Key code #%d not recognised. Press `H` for help or `Q` to quit!\n', keypress)
                end
            end
    end % switch keypress

    %%%% BEGIN { CLEAR NET }
    if clearFNet && existNet
        disp('StepNet:  Clearing functional network memory')
        existNet = false;
        clearFNet = false;
        axes(handTSaxes)
        cla
        if freshTS{1}
            freshTS{1} = false; % retain previous TS info for possible future use
            if ishandle(freshTS{3})
                delete(freshTS{3})
            end
            freshTS{3} = 0; % reset TSbox points plot handle
            freshTS{4} = 0;
        end
        epochTimesSet = false;
        % clear FN vars
        % except: seqTimes, numTot, varnames & vars (these already refreshed or same as before, since we got here from re-loading an CFG or data file)
        clear diag Eseq whosevars candActsIxMapAbs candPotsIxMapAbs numSteps timeStep
        % reset diagram
        axes(handFNaxes)
        cla

        % Get new object handles and re-initialise diagram etc.
        old_handles = { handVVaxes, handTSaxes, handTBaxes, tBarHand, tTickHand, handFNaxes };
		initFNresult = InitFuncNet(DSSRTfigHandle,true, old_handles);
		if isempty(initFNresult)
            disp('StepNet:  Cannot proceed')
            disp('           Perhaps check that a directory for network object specs is present and correct')
            return
		end
		
		handVVaxes = initFNresult{1}; % Variable viewer subplot
        handTSaxes = initFNresult{2}; % Transition Sequence subplot
		handTBaxes = initFNresult{3}; % Time box subplot
		handFNaxes = initFNresult{4}; % Functional net subplot
		tBarHand   = initFNresult{5}{1};
		tTickHand  = initFNresult{5}{2};

        % marks
        if markLset
            delete(markLhand)
            clear markLhand
            markLset = false;
            markLeft = 0;
        end
        if markRset
            delete(markRhand)
            clear markRhand
            markRset = false;
            markRight = 0;
        end

        if exist('total_actives','var')
            clear total_actives
        end
        
        if obsIxSet
            obsIx = 1;
        end

        if CFGfilenameSet
            if resetView
                varViewIx = 1; % reset this
            end
            if varDataFilenameSet % refresh to full availability of time data by resetting vars, seqTimes
                if ~nOutSet
                    beep
                    disp('StepNet:  Internal error! nOut not set before GetVarsAndTimes called during ClearNet')
                    disp('           (varDataFilenameSet == true)')
                    disp('           Press any key to continue with nOut = 2')
                    nOut = 2;
                    nOutSet = true;
                end
                resultVAT = GetVarsAndTimes([DSSRT_ROOTPATH networkObjectDir],varDataFilename,CFGfilename,CFGfilenameSet,verboseTog);
                if isempty(resultVAT)
                    disp('StepNet:  Fatal error. GetVarsAndTimes returned nothing')
                    return
                end
                availStTime = resultVAT{1};
                availSpTime = resultVAT{2};
                if ~isempty(resultVAT{4})
                    numTot      = resultVAT{3};
                    varnames    = resultVAT{4};
                    varsAll     = resultVAT{5};
                    seqTimesAll = resultVAT{6};
                    vars        = varsAll(1:nOut:length(seqTimesAll),:);
                    seqTimes    = seqTimesAll(1:nOut:length(seqTimesAll));
                    numInt      = resultVAT{7};
                    numExt      = resultVAT{8};
                    actsIxMap   = resultVAT{9};
                    potsIxMap   = resultVAT{10};
                    inputsIx    = resultVAT{11};
                    hobj        = resultVAT{12};
                    hvb         = resultVAT{13};
                    DEqns       = resultVAT{14};
                    DEpars      = resultVAT{15};
                    varBounds   = resultVAT{16};
                else
                    beep
                    disp('StepNet:  Problem getting initialisations back from GetVarsAndTimes()')
                    disp('           Bad data file?')
                end
                if startTime < availStTime
                    startTime = availStTime;
                    startTimeSet = true;
                end
                if stopTime > availSpTime
                    stopTime = availSpTime;
                    stopTimeSet = true;
                end
                PlotVBox(handVVaxes,seqTimes,vars(:,varViewIx),varViewIx,varnames)
                % update axes etc. for Time Box and TS box
                axes(handTBaxes)
                delete(tBarHand)
                set(handTBaxes, 'XLim', [ seqTimes(1) seqTimes(length(seqTimes))])
                tBarHand = Draw_tBar(seqTimes(1),seqTimes(length(seqTimes)),'k');
                axes(handTSaxes)
                set(handTSaxes, 'XLim', [ seqTimes(1) seqTimes(length(seqTimes))])
                axes(handFNaxes) % return to FN axes
            else % varDataFilenameSet not set
                if ~nOutSet
                    beep
                    disp('StepNet:  Internal error! nOut not set before GetVarsAndTimes called during ClearNet')
                    disp('           (varDataFilenameSet == false)')
                    disp('           Press any key to continue with nOut = 2')
                    nOut = 2;
                    nOutSet = true;
                end
                resultVAT = GetVarsAndTimes([DSSRT_ROOTPATH networkObjectDir],'',CFGfilename,CFGfilenameSet,verboseTog);
                if isempty(resultVAT)
                    disp('StepNet:  Fatal error. GetVarsAndTimes returned nothing')
                    return
                end
                availStTime = resultVAT{1};
                availSpTime = resultVAT{2};
                if ~isempty(resultVAT{4})
                    numTot      = resultVAT{3};
                    varnames    = resultVAT{4};
                    % These two come from the data file, which we don't
                    % have a specification for yet, so ignore...
                    % varsAll   = resultVAT{5};
                    % seqTimesAll = resultVAT{6};
                    numInt      = resultVAT{7};
                    numExt      = resultVAT{8};
                    actsIxMap   = resultVAT{9};
                    potsIxMap   = resultVAT{10};
                    inputsIx    = resultVAT{11};
                    hobj        = resultVAT{12};
                    hvb         = resultVAT{13};
                    DEqns       = resultVAT{14};
                    DEpars      = resultVAT{15};
                    varBounds   = resultVAT{16};
                    if startTime < availStTime
                        startTime = availStTime;
                        startTimeSet = true;
                    end
                    if stopTime > availSpTime
                        stopTime = availSpTime;
                        stopTimeSet = true;
                    end
                else
                    beep
                    disp('StepNet:  Problem getting initialisations back from GetVarsAndTimes()')
                    disp('           Bad data file?')
                end
            end
		end

		numVBs = length(hvb); % didn't need re-assigning really
        %%%% set up initial stuff on screen...
		% display text for relevant objects, using second argument (initialize) == true
		axes(handFNaxes) % set to current axes
		for obj = 1:numobjects
            DrawNet(hobj{obj}, true); % discard returned (null) handle for initialization
            hobj{obj}{OBJ_HANDLE} = DrawNet( hobj{obj} ); % draw inactive state
		end
        
        for vb = 1:numVBs
            hvbc = hvb{vb};
            xb = hvbc{3};
            yb = hvbc{4};
            yh = hvbc{7};
            hvb{vb}{1} = plot([xb xb],[yb yb+yh],'b-');
            hvb{vb}{2} = plot([xb-0.01 xb+0.01], [yb yb],'k-','LineWidth',2);
		end

        refreshDiagState = true; % this line might be redundant!
    end
    %%%% END { CLEAR NET }

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%     UPDATE DIAGRAM IF SOMETHING HAPPENED!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if existNet && freshTS{1} % then update transition sequence display in main window
        freshTS{1} = false; % reset this flag, but retain other entries for possible future use
        figure(DSSRTfigHandle) % ensure our figure is current
        axes(handTSaxes)
        plotTS = TransSeqs{freshTS{2}}{TS_TSEQ};
        lenTS = length(plotTS);
        epochTimes = [];
        for entry=1:lenTS
            % collect all epoch time coords:
            % second entry in each is the relative time of epoch
            epochTimes = [epochTimes, plotTS{entry}{2}];
        end
        epochTimesSet = true;
        % absolute time given by offset of markLeft position (4th entry in freshTS)
        freshTS{3} = Draw_points(epochTimes+freshTS{4});
        axes(handFNaxes) % return to FN axes
    end
    
    if existNet && refreshDiagState
        if diagSw_handle == -1 || forceNetRedraw
            axes(handFNaxes)
            if diagSw == TERMS
                diagSw_str = 'Term';
                diagSw_fontsize = 20;
                diagSw_xval = -0.05;
            else
                diagSw_str = '\Psi';
                diagSw_fontsize = 25;
                diagSw_xval = -0.05;
            end
            diagSw_handle = text(diagSw_xval, -0.03, diagSw_str, 'FontSize', diagSw_fontsize);
        end
        %%%% Bounds checking for state position index
        if pos < 1
            pos = 1;
        elseif pos > numSteps
            pos = numSteps;
        end
        refreshDiagState = false; % reset this

        diag_state_new = diag{pos};
        figure(DSSRTfigHandle); % ensure our figure is current
        if pos == markLeft
            markStr = ' [L]';
        elseif pos == markRight
            markStr = ' [R]';
        else
            markStr = '     ';
        end
        axes(handFNaxes)
		titFNhand = title(['Position ' num2str(pos,'%-4d') ' @ time ' num2str(seqTimes(pos),'%-5.2f') ' (rel. + ' num2str(seqTimes(pos)-seqTimes(1),'%-5.2f') ') ' markStr], ...
            'FontSize', 12); % update time appearing in title

        %%%% BEGIN { UPDATE VERTICAL BARS }
		for npv = 1:numTot
            ext_var = whosevars(npv); % list of external variables this variable is used by
            hvbc = hvb{npv}; % current handle of vertical bar

            if ext_var ~= 0 % then assume this object outputs to no object, so we cannot colour the vertical bar!
                allacts = Eseq{pos}{ESEQ_ACTIVES_BAR}{ext_var}; % e.g. [1 4]
                allpots = Eseq{pos}{ESEQ_POTENTS_BAR}{ext_var};
                alla_len = length(allacts);
                allp_len = length(allpots);
                if alla_len > 0 % map the indices for actives (local to each ext_var in Eseq) onto the absolute indices, 1:numTot
                    acts_ix = zeros(1,alla_len);
                    for i=1:alla_len
                        acts_ix(i) = candActsIxMapAbs{ext_var}(allacts(i)); % absolute index
                    end
                else
                    acts_ix = [];
                end
                if allp_len > 0 % map the indices for potentials (local to each ext_var in Eseq) onto the absolute indices, 1:numTot
                    pots_ix = zeros(1,allp_len);
                    for i=1:allp_len
                        pots_ix(i) = candActsIxMapAbs{ext_var}(allpots(i)); % was PotsIx (allpots is defined in terms
                        % of the ActsIx indices, although this does not mean we'll get an active erroneously in the pots list!)
                    end
                else
                    pots_ix = [];
                end
            
                if ismember(npv,acts_ix) % then it's active
                    vBarState = 2;
                elseif ismember(npv,pots_ix) % then it's potentiated
                    vBarState = 1;
                else % inactive, unpotentiated
                    vBarState = 0;
                end
            else
                vBarState = 0; % i.e. unused for output-less nodes
            end
            
            hvb{npv} = RedrawVBar(npv, hvbc, vBarState, vars(pos,npv), toggle_ls, log_scale, lmag);
        end % for npv
        toggle_ls = false; % ensure reset log_scale's 'updated' toggle (whether it was false or not)
        %%%% END { UPDATE VERTICAL BARS }
	
        %%%% BEGIN { UPDATE NETWORK OBJECTS }
        % find out which states have changed since last time step, and update diagram
        if forceNetRedraw
            test=ones(1,numobjects);
            forceNetRedraw = false;
        else
            test=zeros(1,numobjects);
            for testnum = 1:numobjects
                test(testnum) = ( diag_state_new(testnum)~=diag_state_old(testnum) );
            end
        end
		if sum(test) > 0 % then something was different, so update diagram according to what changed
            figure(DSSRTfigHandle)
            for obj = 1:numobjects
				if test(obj) % don't need to check for LINKOFFSTATE for link objects here, since those states will never
                    % have changed, so test(obj) will always be false anyway for those objects
                    hobj{obj}{OBJ_STATE} = diag_state_new(obj);
					hobj{obj}{OBJ_HANDLE} = DrawNet( hobj{obj} ); % redraw changed object and update object graphics handle
				end
			end
            diag_state_old = diag_state_new; % only copying this if something was new this step
        end % else don't do anything more
        %%%% END { UPDATE NETWORK OBJECTS }
        
        %%%% BEGIN { UPDATE TIME BAR BOX }
        axes(handTBaxes)
        set(tTickHand, 'XData', [ seqTimes(pos) seqTimes(pos) ], 'Color', 'b');
        if markLchanged
            markLchanged = false;
            if ~markLset
                if exist('markLhand','var')
                    delete(markLhand)
                    clear markLhand
                end
            else
                if exist('markLhand','var')
                    set(markLhand, 'XData', [seqTimes(markLeft) seqTimes(markLeft)]);
                else
                    markLhand = Draw_mark(seqTimes(markLeft));
                end
            end
        end
        if markRchanged
            markRchanged = false;
            if ~markRset
                if exist('markRhand','var')
                    delete(markRhand)
                    clear markRhand
                end
            else
                if exist('markRhand','var')
                    set(markRhand, 'XData', [seqTimes(markRight) seqTimes(markRight)]);
                else
                    markRhand = Draw_mark(seqTimes(markRight));
                end
            end
        end

        axes(handFNaxes)
        %%%% END  { UPDATE TIME BAR BOX }
        
    end % if refreshDiagState
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% Wait for user command over figure
    keyokay = false;
    while ~keyokay
        w = 0;
        while w~=1
            w = waitforbuttonpress;
        end
        keypress = get(DSSRTfigHandle,'CurrentCharacter');
        if ischar(keypress) & ~ismember(keypress, [27, 8, 13, 14, 15])
            keypress = upper(keypress);
            keyokay = true;
        end
	end
    if keypress == 'Q'
        disp('StepNet:  Really quit? (see dialog box)')
        ButtonName = questdlg('Really quit?','Uh-oh','Yes','No','No');
        switch ButtonName
            case 'Yes'
                doQuit = true;
                disp('          Quitting ...')
                disp(' ')
            otherwise
                disp('           Quit cancelled')
        end
    end

end % main while loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% EXIT SETTINGS
% Either 'Q' has been pressed, or there's been an untrapped StepNet error.
if errStatus
    fprintf('\nStepNet:  Nasty untrapped StepNet error must have occurred.')
    disp('          Quitting in disgrace ...')
    beep
    beep
end
% Also remove figure window
if DSSRTfigHandle ~= 0 & ishandle(DSSRTfigHandle)
   delete(DSSRTfigHandle)
end

% remove networkObjectDir from path (which is global)
rmpath([DSSRT_ROOTPATH networkObjectDir])

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                Secondary functions only used by StepNet
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% INFO: fields for regimeStruct
%regimeStruct.regimes
%regimeStruct.numEps
%regimeStruct.numRegs
%regimeStruct.transSeq
%regimeStruct.isCycle
%regimeStruct.varnames
%regimeStruct.period
%regimeStruct.startTime    % not necessarily the same as the value passed to RegimeDet
%regimeStruct.tScaleThresh
%regimeStruct.dScaleThresh
%regimeStruct.algOptions
%regimeStruct.focVarIx
function DisplayRegimes(regimeStruct,verboseTog)
global TS_NAME
if verboseTog
    disp(' ')
end
regimes = regimeStruct.regimes;
disp(   'StepNet:  Key:  (F) = variable has fast response,  (S) = variable has slow response')
disp(   '                (C) = variable is near-constant,   (X) = variable is a cross-multiplier in an input')
disp(   '                (I) = variable is inactive (not in set of actives)')
disp(   '              Dim_A = all actives with assoc. DE`s')
disp(   '              Dim_B = Dim_A - fast and slow vars only')
disp(   '              Dim_C = Dim_A - fast, slow & near-const vars')
fprintf('StepNet:  Determined %i reduced regimes within %i epochs\n',regimeStruct.numRegs,regimeStruct.numEps)
fprintf('           from transition sequence `%s`\n',regimeStruct.transSeq{TS_NAME})
fprintf('           using time scale threshold 1/%.4f\n',regimeStruct.tScaleThresh)
algStrs = {'OFF','ON'};
sfStrs = {'','(not strict)','(strict)'};
disp(   '          Algorithm options used:')
fprintf('           speed fussiness:       %s %s\n',algStrs{regimeStruct.algOptions.speedFussy+1}, sfStrs{regimeStruct.algOptions.speedFussyLevel+1})
fprintf('           fast leaver fussiness: %s\n',algStrs{regimeStruct.algOptions.fastLeaveFussy+1})
fprintf('           long epoch fussiness:  %s\n',algStrs{regimeStruct.algOptions.longFussy+1})
if regimeStruct.isCycle && regimeStruct.period > 0
    fprintf('           Cycle period = %.4f\n',regimeStruct.period)
end
maxDim = [1 1 1];
for reg=1:regimeStruct.numRegs
    fprintf('           %i : t = [%.4f, %.4f),  rel. + [%.4f, %.4f) to first regime\n',reg, ...
        regimes(reg).timeInt(1),regimes(reg).timeInt(2),regimes(reg).timeInt(1)-regimeStruct.startTime, ...
        regimes(reg).timeInt(2)-regimeStruct.startTime)
    fprintf('             containing epochs:        ')
    for epochIx=regimes(reg).epochs
        fprintf('%i  ',epochIx)
    end
    dims = [ regimes(reg).dimensionA regimes(reg).dimensionB regimes(reg).dimensionC ];
    for d=1:3
        if dims(d) > maxDim(d)
            maxDim(d) = dims(d);
        end
    end
    fprintf('\n')
    fprintf('             tot. dynamic dimension:   %i (A), %i (B), %i (C)\n',dims(1), dims(2), dims(3))
    fprintf('             dynamic vars:             ')
    for var=union(regimes(reg).dynVars, regimeStruct.focVarIx)
        if ismember(var,regimes(reg).fastVars)
            tstring = '(F)';
        elseif ismember(var,regimes(reg).slowVars)
            tstring = '(S)';
        else
            tstring = '';
        end
        if ismember(var,regimes(reg).constVars)
            tstring = [tstring '(C)'];
        end
        if ~ismember(var,regimes(reg).qsPars)
            fprintf('%s%s  ',regimeStruct.varnames{var},tstring)
        end
    end
    for var=regimes(reg).dynCross
        if ismember(var,regimes(reg).fastVars)
            tstring = '(F)';
        elseif ismember(var,regimes(reg).slowVars)
            tstring = '(S)';
        else
            tstring = '';
        end
        if ismember(var,regimes(reg).constVars)
            tstring = [tstring '(C)'];
        end
        if ~ismember(var,regimes(reg).qsPars)
            fprintf('%s(X)%s  ',regimeStruct.varnames{var},tstring)
        end
    end
    fprintf('\n')
    fprintf('             non-dynamic vars:         ')
    for var=regimes(reg).nonDynVars
        fprintf('%s  ',regimeStruct.varnames{var})
    end
    fprintf('\n')
    fprintf('             quasi-static bif`n pars:  ')
    for var=regimes(reg).qsPars
        if ismember(var,regimes(reg).fastVars)
            tstring = '(F)';
        elseif ismember(var,regimes(reg).slowVars)
            tstring = '(S)';
        else
            tstring = '(I)'; % inactive (not in set of actives) is the only other possibility
        end
        fprintf('%s%s  ',regimeStruct.varnames{var},tstring)
    end
    fprintf('\n\n')
end
maxDimRegs = {{'';' '}, {'';' '}, {'';' '}};
for reg=1:regimeStruct.numRegs % go around again to display info about maximum dimension of regime models
    dims = [ regimes(reg).dimensionA regimes(reg).dimensionB regimes(reg).dimensionC ];
    for d=1:3
        if dims(d) == maxDim(d)
            maxDimRegs{d}{1} = [maxDimRegs{d}{1} num2str(reg,'%3i') ' ' ];
        end
    end
end
for d=1:3
    if length(maxDimRegs{d}{1})>2
        maxDimRegs{d}{2} = 's';
    end
end
fprintf('           Maximum regime dimensions: (A) %i -- for regime%s %s\n',maxDim(1),maxDimRegs{1}{2},maxDimRegs{1}{1})
fprintf('                                      (B) %i -- for regime%s %s\n',maxDim(2),maxDimRegs{2}{2},maxDimRegs{2}{1})
fprintf('                                      (C) %i -- for regime%s %s\n\n',maxDim(3),maxDimRegs{3}{2},maxDimRegs{3}{1})
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newPos = FindPosFromTime(t,seqTimes,numSteps)
newPos = 1;
if t < seqTimes(1)
    disp('StepNet:  Time specified is too small for this functional network. Using lowest time in set')
    return % not strictly necessary as cannot pass through next while statement anyway
end
while seqTimes(newPos) < t % find position associated with this time
    newPos = newPos + 1;
    if newPos == numSteps
        disp('StepNet:  Time specified is too large for this functional net. Using largest time in set')
        break
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% copy of this function also exists in AttEst_*.m and PlotPhaseDiagrams.m
function DEqns = SetActivesSwitch(DEqns, focusSet, deqnIxMap, numExt, actsIxMap, allActs)
% do this only on COMPILED Gamma sets, for storing in either compiled or uncompiled equations
global DE_NAMEi DE_GAMMA1i DE_GAMMA2i DE_ACTSWi DE_ISTAUFFILEi DE_TAURECIPi DE_GAMVARNAMEi ...
    DE_GAMVARPOWi DE_ISTGTFFILEi DE_TARGETi DE_INTVARNAMEi DE_INTVARPOWi

for thisAbsIx = focusSet
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function new_hvbc = RedrawVBar(npv, hvbc, state, val, toggle_ls, log_scale, lmag)
switch state
    case 0 % inactive
        col = 'b';
    case 1 % potentiated
        col = 'g';
%         col = 'b'; % swap commenting here to ignore potentials in func net figure
    case 2 % active
        col = 'r';
end

new_hvbc = hvbc; % initial value

if hvbc{9} > 0 % then this object can be log-scaled
    if log_scale % only use log scale for those vars with that option enabled
        % Currently assumes this is only applied to gating vars running
        % from 0..1  ... since it doesn't adjust for other max and min of possible values
        % (see maple file FuncNet_logscale.mws for details of how to
        % re-scale for other variables)
        if hvbc{9} == 2 % then 'reverse' log to magnify top range of the variable
            yc = hvbc{4}+hvbc{7}*(1-log(lmag*(1-val)+1)/log(lmag+1)-hvbc{5})/hvbc{6};
            if toggle_ls
                new_hvbc{10}=plot(hvbc{3}, hvbc{4}+hvbc{7}+0.01, 'b.'); % show which vars are in log scale
            end
        else % assume == 1
            yc = hvbc{4}+hvbc{7}*(log(lmag*val+1)/log(lmag+1)-hvbc{5})/hvbc{6};
            if toggle_ls
                new_hvbc{10}=plot(hvbc{3}, hvbc{4}-0.01, 'b.'); % show which vars are in log scale
            end
        end
    else
        yc = hvbc{4}+hvbc{7}*(val-hvbc{5})/hvbc{6};
        if hvbc{10} ~= 0 && toggle_ls
            delete(hvbc{10})
            new_hvbc{10} = 0; % reset handle in bar's data structure
        end
    end
else % always linear scale
    yc = hvbc{4}+hvbc{7}*(val-hvbc{5})/hvbc{6};
end
if yc < hvbc{4},           yc=hvbc{4}; end
if yc > hvbc{4}+hvbc{7},   yc=hvbc{4}+hvbc{7}; end
set(hvbc{2}, 'YData', [yc yc], 'Color', col)

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function result = InitFuncNet(fig,restarting,handles)
result = {}; % default value
switch nargin
    case 1
        disp('InitFuncNet:  Initializing for first time')
        restarting = false; % default value
        handles = {};
    case 3
        if restarting && isempty(handles)
            disp('InitFuncNet:  Cannot restart with empty `handles` parameter')
            return
        end
    otherwise
        disp('InitFuncNet:  Wrong number of arguments passed')
        return
end

if fig ~= 0
    figure(fig);
else
    disp('InitFuncNet:  Figure handle must be provided to InitFuncNet in order to initialize')
    return
end

if ~restarting
	% variable view box
	handVarVBox = subplot('Position', [0.1 0.885 0.85 0.07]); % NOTE - this position is also stored in refresh command '0'
	title('variable unspecified')
	axis(handVarVBox,'on')
    axis tight
	set(handVarVBox,'XTickLabelMode','manual')
	set(handVarVBox,'YTickLabelMode','manual')
	set(handVarVBox,'XTickLabel','')
	set(handVarVBox,'YTickLabel','')

    % transition sequence box
    handTSeqBox = subplot('Position', [0.1 0.84 0.85 0.01]); % NOTE - this position is also stored in refresh command '0'
    hold on
    axis(handTSeqBox,'off')

	% time box
	handTimeBox = subplot('Position', [0.1 0.8 0.85 0.04]); % NOTE - this position is also stored in refresh command '0'
	hold on
	tBar = Draw_tBar(0,1,'k');
	tTick = Draw_tTick( 0, 'k' );
    set(handTimeBox,'DrawMode','fast')
	axis(handTimeBox,'off')
	
	% func net box
	handFNetBox = subplot('Position', [0.05 0.04 0.9 0.71] ); % NOTE - this position is also stored in refresh command '0'
	set(handFNetBox,'PlotBoxAspectRatioMode','manual')
	set(handFNetBox,'XLimMode','manual')
	set(handFNetBox,'YLimMode','manual')
	set(handFNetBox,'DrawMode','fast')
	axis(handFNetBox,'off')
	title('No functional network. Press `H` for help, `Q` to quit')
	hold on
else % re-start
    handVarVBox = handles{1};
    axes(handVarVBox)
    cla
    title('variable unspecified')
    set(handVarVBox,'XTickLabelMode','manual')
	set(handVarVBox,'YTickLabelMode','manual')
	set(handVarVBox,'XTickLabel','')
	set(handVarVBox,'YTickLabel','')

    handTSeqBox = handles{2};
    axes(handTSeqBox)
    cla
    
    handTimeBox = handles{3};
    axes(handTimeBox)
    tBar = handles{4};
    tTick = handles{5};
    delete(tBar)
    delete(tTick)
    tBar = Draw_tBar(0,1,'k');
	tTick = Draw_tTick( 0, 'k' );

    handFNetBox = handles{6};
    axes(handFNetBox)
    cla % reset % may have to use `reset` command to delete any hidden objects (like vbars)
              % ... to cover future implementation of vbar hiding
	set(handFNetBox,'PlotBoxAspectRatioMode','manual')
	set(handFNetBox,'XLimMode','manual')
	set(handFNetBox,'YLimMode','manual')
	set(handFNetBox,'DrawMode','fast')
	axis(handFNetBox,'off')
	title('No functional network. Press `H` for help, `Q` to quit')
	hold on
end % if ~ restarting

result = {handVarVBox, handTSeqBox, handTimeBox, handFNetBox, {tBar, tTick}};
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = DialogBox(promptTitle,defVal,verboseTog)
% This function assumes that negative numbers will not be legal
%  values for any parameters to be got from user input!
if nargin ~=3
    disp('DialogBox:  Wrong number of arguments passed - expected 3')
    return
end
result = '';
if verboseTog
    comLineStr1 = 'StepNet:  Enter ';
    comLineStr2 = ' in dialog box';
    disp([ comLineStr1 lower(promptTitle) comLineStr2 ]);
end

prompt  = [ promptTitle '?' ];
ptitle  = 'StepNet parameter entry:';
lines   = 1;

if isempty(defVal)
    def = {''};
else
    if isnumeric(defVal)
        def = {num2str(defVal)};
	elseif ischar(defVal)
        def = {defVal};
	else
        disp('DialogBox:  Internal error. `defVal` must be a numeric or string')
        return
	end
end

answer  = inputdlg(prompt,ptitle,lines,def,'off');
if ~isempty(answer)
    if ~isempty(answer{1})
        result = answer{1};
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = ReadParamsFile(filename)
global NUMPARS
result = {}; % initial and default value (signifies error reading file)
errStatus = false;
errLine   = 0;
fileLines = cell(1,NUMPARS);
fid = fopen(filename,'r');
if fid ~= -1
    for i = 1:NUMPARS
        fileLines{i} = fgetl(fid);
        if fileLines{i} == -1
            errStatus = true;
            break
        end
    end
    if ~errStatus
        CFGfilename    = sscanf(fileLines{1},'%s');
        varDataFilename = sscanf(fileLines{2},'%s');
        startTime      = sscanf(fileLines{3},'%f');
        stopTime       = sscanf(fileLines{4},'%f');
        dScaleThresh   = sscanf(fileLines{5},'%f');
        nOut           = sscanf(fileLines{6},'%i');
        ignoreActsSt   = sscanf(fileLines{7},'%i');
        guessPeriod    = sscanf(fileLines{8},'%f');
        marginPC       = sscanf(fileLines{9},'%f');
        tScaleThresh   = sscanf(fileLines{10},'%f');
        attEstP.delta  = sscanf(fileLines{11},'%f'); % no bounds checking for this (unused parameter)
        attEstP.dVmax  = sscanf(fileLines{12},'%f');
        attEstP.Vinres = sscanf(fileLines{13},'%f');
        attEstP.Vpert  = sscanf(fileLines{14},'%f');
        attEstP.searchRes       = sscanf(fileLines{15},'%i');
        attEstP.lowResMultiple  = sscanf(fileLines{16},'%i');
        attEstP.derivThresh     = sscanf(fileLines{17},'%f');

        % perform basic bounds checking
        if startTime < 0,                          errStatus = true; errLine = 3; end
        if stopTime <= 0,                          errStatus = true; errLine = 4; end
        if dScaleThresh < 1,                       errStatus = true; errLine = 5; end
        if nOut <= 0 || (nOut ~= round(nOut)),     errStatus = true; errLine = 6; end % not a positive integer
        if ignoreActsSt ~= 0 && ignoreActsSt ~= 1, errStatus = true; errLine = 7; end % not a boolean
        if guessPeriod <= 0,                       errStatus = true; errLine = 8; end % not a valid time
        if marginPC <=0 || marginPC >= 1,          errStatus = true; errLine = 9; end % not a valid non-zero percentage
        if tScaleThresh < 1,                       errStatus = true; errLine = 10; end
        if attEstP.dVmax <= 0,                                     errStatus = true; errLine = 12; end % checking against variable's bounds done in command 'R' code
        if attEstP.Vpert <= 0 || attEstP.Vpert >= attEstP.dVmax,   errStatus = true; errLine = 14; end
        if attEstP.Vinres <= 0 || attEstP.Vinres >= attEstP.dVmax, errStatus = true; errLine = 13; end
        if attEstP.searchRes <= 0 || attEstP.searchRes > 500 || (round(attEstP.searchRes) ~= attEstP.searchRes)
            errStatus = true; errLine = 15;
        end
        if attEstP.lowResMultiple <= 0 || attEstP.lowResMultiple > 100 || (round(attEstP.lowResMultiple) ~= attEstP.lowResMultiple)
            errStatus = true; errLine = 16;
        end
        if attEstP.derivThresh <= 0, errStatus = true; errLine = 17; end
    end % else wrong number of lines in param file (corrupted), for instance
else
    errStatus = true;
    fprintf('ReadParamsFile:  Error opening parameters file %s\n',filename)
end

if ~errStatus
    result{1} = CFGfilename;
    result{2} = varDataFilename;
    result{3} = startTime;
    result{4} = stopTime;
    result{5} = dScaleThresh;
    result{6} = nOut;
    result{7} = ignoreActsSt;
    if ignoreActsSt,  iAO_str = 'ON';  else  iAO_str = 'OFF';  end
    result{8} = guessPeriod;
    result{9} = marginPC;
    result{10} = tScaleThresh;
    result{11} = attEstP;

    attDelStr = num2str(attEstP.delta,'%.4f');
    attdVStr = num2str(attEstP.dVmax,'%.4f');
    attVsStr = num2str(attEstP.Vinres,'%.4f');
    attVpStr = num2str(attEstP.Vpert,'%.4f');
    attsRStr = num2str(attEstP.searchRes,'%3i');
    attlrMStr = num2str(attEstP.lowResMultiple,'%3i');
    attdTStr = num2str(attEstP.derivThresh,'%.4f');

    fprintf('           CFG filename                 = `%s`\n',CFGfilename)
    fprintf('           data filename                = `%s`\n',varDataFilename)
    fprintf('           start time                   = %4.2f\n',startTime)
    fprintf('           stop time                    = %4.2f\n',stopTime)
    fprintf('           dominance scale threshold    = 1/%.4f\n',dScaleThresh)
    fprintf('           time scale threshold         = 1/%.4f\n',tScaleThresh)
    fprintf('           nOut                         = %i\n',nOut)
    fprintf('           ignore actives/pots order    = %s\n',iAO_str)
    fprintf('           guess period                 = %.2f\n',guessPeriod)
    fprintf('           margin p.c. of period        = %.2f\n',marginPC)
    fprintf('           Att Est delta                = %s (currently unused)\n',attDelStr)
    fprintf('           Att Est dVmax                = %s\n',attdVStr)
    fprintf('           Att Est Vinres               = %s\n',attVsStr)
    fprintf('           Att Est Vpert                = %s\n',attVpStr)
    fprintf('           Att Est search resolution    = %s\n',attsRStr)
    fprintf('           Att Est low res multiple     = %s\n',attlrMStr)
    fprintf('           Att Est derivative threshold = %s\n',attdTStr)
else
    fprintf('ReadParamsFile:  Error (e.g. violation of parameter bounds) in parameters file %s\n',filename)
    fprintf('                 on line %i\n',errLine)
end
fclose(fid);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = WriteParamsFile(filename,pars)
global NUMPARS
result = false;
if length(pars) + 6 ~= NUMPARS % + 6 for AttEst parameters
    fprintf('WriteParamsFile:  Internal error! Incorrect number of parameter values passed\n')
    return
end
attEstP = pars{11};
aep1 = attEstP.delta;
aep2 = attEstP.dVmax;
aep3 = attEstP.Vinres;
aep4 = attEstP.Vpert;
aep5 = attEstP.searchRes;
aep6 = attEstP.lowResMultiple;
aep7 = attEstP.derivThresh;
fid = fopen(filename,'w');
if fid ~= -1
    fprintf(fid, strcat( [ pars{1} '\n' pars{2} '\n' num2str(pars{3},'%5.4f') '\n'],...
        [num2str(pars{4},'%5.4f') '\n' num2str(pars{5},'%2.4f') '\n' num2str(pars{6},'%2i') '\n' ], ...
        [num2str(pars{7},'%2i') '\n' num2str(pars{8},'%4.4f') '\n' num2str(pars{9},'%4.2f') '\n' ], ...
        [num2str(pars{10},'%2.4f') '\n'], [num2str(aep1,'%2.4f') '\n'], [num2str(aep2,'%2.4f') '\n'], ...
        [num2str(aep3,'%2.4f') '\n'], [num2str(aep4,'%2.4f') '\n'], ...
        [num2str(aep5,'%2i') '\n'], [num2str(aep6,'%2i') '\n'], [num2str(aep7,'%2.4f') '\n']) );
    result = true;
else
    fprintf('WriteParamsFile:  Error opening parameters file %s\n',filename)
end
fclose(fid);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = GetVarsAndTimes(networkObjectDir,varDataFilename,CFGfilename,CFGset,verboseTog)
% networkObjectDir includes full root path
result = {}; % initial value (and error status)
setupDone = false;
varsTimesDone = false;
VTfilename = [ EndStripStr(varDataFilename,'.',true) '_VT.mat' ];
if exist( VTfilename, 'file') == 2 % then file of Vars, Times, and other setup info exists
    if verboseTog
        disp('GetVarsAndTimes:  Attempting to load setup info from VT file')
    end
    load( VTfilename );
    if exist('setupCFG','var') && exist('seqTimes','var')
        setupDone = true; % load was successful
        if verboseTog
            disp('                  Load successful.')
        end
    else
        disp('GetVarsAndTimes:  Error loading VT file. Re-determining setup info')
    end
end
    
if CFGset
    if verboseTog
        disp('GetVarsAndTimes:  Reading CFG file for setup info')
    end
    if ~setupDone % then need to read from file
        setupCFG = ReadConfigFile(CFGfilename,networkObjectDir);
    end
    if ~isempty(setupCFG)
        numInt    = setupCFG{1};
        numExt    = setupCFG{2};
        numTot    = numInt+numExt;
        varnames  = setupCFG{5};
        inputsIx  = setupCFG{6};
        ActsIxMap = setupCFG{7};
        PotsIxMap = setupCFG{8};
        diagObj   = setupCFG{9};
        vBobj     = setupCFG{10};
        DEqns     = setupCFG{11};
        DEpars    = setupCFG{12};
        varBounds = setupCFG{13};
        if verboseTog
            disp('                  File read successfully')
        end
    else
        disp('GetVarsAndTimes:  Error reading CFG file')
        return
    end
end

if ~setupDone
    if verboseTog
        disp('GetVarsAndTimes:  Reading variable data file')
    end
	if ~isempty(varDataFilename)
        if exist( varDataFilename, 'file') == 2
            allVarDat = load(varDataFilename);
        else
            fprintf('GetVarsAndTimes:  Variable data file %s not found\n',varDataFilename)
            return
        end
        varColLen = length(allVarDat(1,:));
        varRowLen = length(allVarDat(:,1));
        if varColLen < numTot + 1
            disp('GetVarsAndTimes:  Data file had fewer columns than number of declared variables!')
            return
        end
        if varRowLen == 1
            disp('GetVarsAndTimes:  Data file had only one time entry!')
            return
        end
        seqTimes = allVarDat(:,1);
        var_data = allVarDat(:,2:1+numTot);
        availStTime = seqTimes(1);
        availSpTime = seqTimes(length(seqTimes));
        if availSpTime > availStTime
            if verboseTog
                disp('                  Variable data file read successfully')
            end
		else
            disp('                  Data file error. No variable and time info retrieved')
            return
		end
	end % if ~isempty(varDataFilename)
end

if ~setupDone % then didn't read everything from a _VT file, so make a new one for future use
    save( VTfilename, 'setupCFG', 'availStTime', 'availSpTime', 'var_data', 'seqTimes' );
    disp('GetVarsAndTimes:  VT file saved for future use')
end

result = {availStTime, availSpTime, numTot, varnames, var_data, seqTimes, numInt, numExt, ...
          ActsIxMap, PotsIxMap, inputsIx, diagObj, vBobj, DEqns, DEpars, varBounds };
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = FindIdentState(Eseq,startPos,EseqIx,direction,no_order,guessPeriod,marginPC,evs,verboseTog)
global ESEQ_TIMES % index in Eseq
result = -1; % initial and default value
switch nargin
    case 5 % set defaults
        marginPC = 0.9;
        no_order = 1;
        evs = [];
        verboseTog = false;
    case 6 % set defaults
        no_order = 1;
        evs = [];
        verboseTog = false;
    case 7 % set defaults
        evs = [];
        verboseTog = false;
    case 8 % set defaults
        verboseTog = false;
    case 9 % ignore
    otherwise
        disp('FindIdentState:  Wrong number of arguments passed')
        return
end

if EseqIx < 1 || EseqIx > 4
    disp('FindIdentState:  Index into Eseq out of range [1,4]')
    return
end

maxPos = length(Eseq);
timeStep = Eseq{2}{ESEQ_TIMES} - Eseq{1}{ESEQ_TIMES};
startPdLen = round(guessPeriod*marginPC/timeStep); % start looking within marginPC% of guessPeriod
% initial position to start search
pos = startPos + direction*startPdLen;
if verboseTog
    if direction == 1
        dirnStr = 'forwards';
    elseif direction == -1
        dirnStr = 'backwards';
    else
        disp('FindIdentState:  Bad direction parameter')
        return
    end
    fprintf(['FindIdentState:  Starting search %i steps ' dirnStr ' from %i\n'], startPdLen, startPos)
end

if pos < 1 || pos > maxPos
    disp('FindIdentState:  The value of `guess period x margin%` is out of time range for this Eseq')
    disp('                 Cannot start search!')
    return
end

foundID = false;
init_elt = Eseq{startPos}{EseqIx};
Nv = length(Eseq{1}{EseqIx}); % this is constant for all positions in the sequence
                              % and is the number of observables (external variables) recorded
test = zeros(1,Nv); % only need to initialize once, always gets overwritten
if isempty(evs)
    evs = 1:Nv
end
while ~foundID
    pos = pos + direction;
    if pos > maxPos || pos < 1,  break, end % nothing found 
    for ev=evs
        cand_elt = Eseq{pos}{EseqIx};
        len_init = length(init_elt{ev});
        len_cand = length(cand_elt{ev});
        if len_cand == len_init
            if no_order
                test(ev) = (len_init ~= length(intersect(init_elt{ev},cand_elt{ev})));
            else
                test(ev) = sum(init_elt{ev} ~= cand_elt{ev}); % ~= returns 0 if sets equal
            end
        else
            test(ev) = 1; % not equal length so sets not equal
        end
    end
    if sum(test) == 0 % compare succeeded -- this is an identical state
        foundID = true;
    end
end

if foundID
    result = pos; % pos has been increased/decreased by one too many in search
                  % loop to correspond to the actual identical position
else
    result = -1;
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function DEqn = SetActivesSwitchSingle(DEqn, thisAbsIx, thisDEix, thisAbsActs)
% This version (cf. SetActivesSwitch()) takes only single focused vars and presumes that DEqnIx is correct
%  and that actives set has already been determined
% Us this only on COMPILED Gamma sets, for storing in either compiled or uncompiled equations
global DE_NAMEi DE_GAMMA1i DE_GAMMA2i DE_ACTSWi DE_ISTAUFFILEi DE_TAURECIPi DE_GAMVARNAMEi ...
    DE_GAMVARPOWi DE_ISTGTFFILEi DE_TARGETi DE_INTVARNAMEi DE_INTVARPOWi

thisGamma1 = DEqn{DE_GAMMA1i};
thisGamma2 = DEqn{DE_GAMMA2i};
lenG1 = length(thisGamma1);
lenG2 = length(thisGamma2);
for g1t=1:lenG1
    DEqn{DE_GAMMA1i}{g1t}{DE_ACTSWi} = ismember( thisGamma1{g1t}{DE_GAMVARNAMEi}, thisAbsActs );
end
for g2t=1:lenG2
    DEqn{DE_GAMMA2i}{g2t}{DE_ACTSWi} = ismember( thisGamma2{g2t}{DE_GAMVARNAMEi}, thisAbsActs );
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


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
% A copy of this function is also present inside AttEst.m (it's different
% in FuncNet.m)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function qsfp = GetQsfpVal(thisDE,varDataLine,onlyActives)
% only accepts "compiled" thisDE
global DE_NAMEi DE_GAMMA1i DE_GAMMA2i DE_ACTSWi DE_ISTAUFFILEi DE_TAURECIPi DE_GAMVARNAMEi ...
    DE_GAMVARPOWi DE_ISTGTFFILEi DE_TARGETi DE_INTVARNAMEi DE_INTVARPOWi LARGEBOUND
% V0eqn in terms of conductances and currents (for example of voltage equations):
%  capacitance is irrelevant here so it's ignored
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
%     disp('StepNet:  Warning. Division by zero in GetQsfpVal()')
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
    if ~onlyActives || ( onlyActives && g1term{DE_ACTSWi} )
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
    if ~onlyActives || ( onlyActives && g2term{DE_ACTSWi} )
        tau_recipVal = GetTauRVal(g2term{DE_ISTAUFFILEi},g2term{DE_TAURECIPi},...
            varDataLine(g2term{DE_GAMVARNAMEi}),g2term{DE_GAMVARPOWi});
        sum = sum + tau_recipVal;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function thisDE_new = CompileGammas(thisDE_old, DEpars, varnames)
% gam1term order: actSw, filefuncflag, taurecip, var, power, filefuncflag, target [, intvar, power]
% gam2term order: actSw, filefuncflag, taurecip, var, power
global DE_NAMEi DE_GAMMA1i DE_GAMMA2i DE_ACTSWi DE_ISTAUFFILEi DE_TAURECIPi DE_GAMVARNAMEi ...
    DE_GAMVARPOWi DE_ISTGTFFILEi DE_TARGETi DE_INTVARNAMEi DE_INTVARPOWi

thisDE_new = {};

Gamma1 = thisDE_old{DE_GAMMA1i}; % to become the index- and par- ready versions
Gamma2 = thisDE_old{DE_GAMMA2i};
lenG1 = length(Gamma1);
lenG2 = length(Gamma2);

if lenG1 + lenG2 == 0
    beep
    disp('CompileGammas:  FATAL ERROR - no terms in focused equation!')
    fprintf('  (Equation for %s)\n', thisDE_old{DE_NAMEi})
    return
end

for g1t=1:lenG1
    g1term = Gamma1{g1t};
    if ~g1term{DE_ISTAUFFILEi}
        result = LookupDEpar( g1term{DE_TAURECIPi}, DEpars );
        if result{2}
            Gamma1{g1t}{DE_TAURECIPi} = result{1};
        else
            fprintf('CompileGammas:  Bad DEpar lookup for %s\n',g1term{DE_TAURECIPi})
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

for g2t=1:lenG2
    g2term = Gamma2{g2t};
    if ~g2term{DE_ISTAUFFILEi}
        result = LookupDEpar( g2term{DE_TAURECIPi}, DEpars);
        if result{2}
            Gamma2{g2t}{DE_TAURECIPi} = result{1};
        else
            fprintf('CompileGammas:  Bad DEpar lookup for %s\n',g2term{DE_TAURECIPi})
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% strip leading path from a filename
function stripStr = RootStripStr(inputStr)
ix1 = max(findstr('\',inputStr)); % cover us for dos filenames
ix2 = max(findstr('/',inputStr)); % cover us for unix filenames
if isempty(ix1)
    ix = ix2;
elseif isempty(ix2)
    ix = ix1;
else
    ix = max(ix1,ix2);
end
if ix > 0
    if ix<length(inputStr)
        stripStr = inputStr(ix+1:length(inputStr)); % start from position just after *final* '\' (or '/') char in case of extra sub-directories in path
    else
        stripStr = inputStr;
    end
else
    stripStr = inputStr;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stripStr = EndStripStr(inputStr,findChar,lastDelim,lastScope)
% strip trailing extension in a filename, using findChar to mark the beginning of the extension
if nargin==2
    lastDelim = false;
    lastScope = 5;
end
if nargin==3
    lastScope = 5;
end
if lastScope > length(inputStr)-1
    lastScope = length(inputStr)-1;
end
ix = findstr(findChar,inputStr);
if ~isempty(ix)
    if lastDelim
        if max(ix) > length(inputStr)-lastScope % filter out all except the 'lastScope' number of possibilities
                                                % (e.g. for filename extensions of known size)
            stripStr = inputStr(1:max(ix)-1); % keep all other occurrences of findChar except the final one
        else
            stripStr = inputStr;
        end
    else
        if min(ix)>1
            stripStr = inputStr(1:min(ix)-1); % keep only up to the first occurrence of findChar
        else
            stripStr = inputStr;
        end
    end
else
    stripStr = inputStr;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stripStr = EndRemainStr(inputStr,findChar,lastDelim)
% this function returns the remainder (e.g. the extension) of a string, without the delimiter findChar 
if nargin==2
    lastDelim = false;
end
ix = findstr(findChar,inputStr);
if ~isempty(ix)
    if lastDelim
        if max(ix)>1
            stripStr = inputStr(max(ix)+1:length(inputStr));
        else
            stripStr = [];
        end
    else
        if min(ix)>=1
            stripStr = inputStr(min(ix)+1:length(inputStr));
        else
            stripStr = [];
        end
    end
else
    stripStr = [];
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function filename_write = GetNextFilename(filename,extStr)
% (c) Robert Clewley, BU, 2004
% Can only produce filenames for up to 999 different filenames with this extension in the local directory.
% extStr is the extension used for the filename. While this should NOT contain the '.' character,
%   this function will remove it if present in the first position of the string.
findDot = strfind(extStr,'.');
if ~isempty(findDot)
    if length(findDot) > 1
        disp('GetNextFilename error: Too many `.` sub-extensions in extStr. These are not supported in this function')
        filename_write = '';
        return
    else
        if findDot ~=1
            disp('GetNextFilename error: `.` was found inside extStr argument.')
            filename_write = '';
            return
        else % the dot is the first character
            extStr = extStr(2:length(extStr));
        end
    end
end

filename_write = filename; % initial value
filename_extstrip = EndStripStr(filename,'.',true);
trailingStr = EndRemainStr(filename_extstrip,'_');
foundFreeNum = false; % initial value
numCounter = 1;
while ~foundFreeNum && numCounter <= 999
	if ~isempty(trailingStr) && isNum(trailingStr)
        newNumStr = num2str(str2num(trailingStr)+numCounter,'%3.3i');
	else
        newNumStr = num2str(numCounter,'%3.3i');
	end
    filename_write = [EndStripStr(filename_extstrip,'_',true) '_' newNumStr '.' extStr];
    if ~exist(filename_write,'file')
        foundFreeNum = true;
    else
        numCounter = numCounter + 1;
    end
end

if numCounter > 999
    beep
    disp(['GetNextFilename error. Too many `.' extStr '` files (more than 999)! Cannot continue'])
    filename_write = '';
    return
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function index = VarnameIxMap(name, varnames)
% assumes varnames is a cell array, containing no repetitions!
[present index] = ismember(name,varnames);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% useful as a dummy argument for KeyPressFcn callback, to prevent
% echo of key presses to command window
function DoNothing()
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fignum = AttFig(attData)
% Routine to produce Attractor Estimate figure
% plotData is a cell array of these entries:
%  {epochLabel, epochTlength, epochTdata, epochVdata, epochActNames, epochAttInterval}
if nargin > 2 || nargin == 0
    beep
    disp('AttFig:  Internal error! Wrong number of parameters passed. Expected 1 or 2')
    fignum = -1; % error status
    return
end
if nargin == 1
    title_str = '';
end

fignum = figure('NumberTitle', 'off', 'Name', attData.nameStr, 'Position', [280, 200, 420, 260]);
figure(fignum)
title(attData.figTitle)
axis manual

minY = 1000;
maxY = -1000;

gcax = gca;
set(gcax,'Position',[0.1 0.1 0.85 0.8])
hold on

drawnow
figure(fignum); % to try to stop Matlab getting stuck on a different figure before drawing this

% These time values have been hacked to counter mismatch in figure -- but don't know source!!
% I'm assuming it really is dt ... and why doesn't it apply to the original V data in blue?
dt = attData.dt; % more readable to use a smaller var name in the plotting below
for i=1 : attData.numEpochs    
    pData = attData.plotData{i};
    interval = pData{6};
    if min(interval) < minY
        minY = min(interval);
    end
    if max(interval) > maxY
        maxY = max(interval);
    end
    
    ni = 1+i; % this is allowed in the sense that the extra epoch is for padding,
              % and only contains the initial vertical interval
    % ni is defined explicitly here, in case in future, for cycles, we use a modulo
    % scheme to refer back to first epoch

    if pData{7} % then original epoch boundary
        plot([pData{3}(1)-dt,pData{3}(1)-dt],interval,'-r') % cross-section of attractor estimate
    end % don't plot the non-original vertical cross-sections
    
    % linear interpolation lines between cross-sections
    if i ~= attData.numEpochs || ( i == attData.numEpochs && length(attData.plotData) == ni ) % then we can safely access an entry for the next (dummy) epoch
        plot([pData{3}(1)-dt,pData{3}(1)+pData{2}-dt],[interval(1), attData.plotData{ni}{6}(1)],'-g')
        plot([pData{3}(1)-dt,pData{3}(1)+pData{2}-dt],[interval(2), attData.plotData{ni}{6}(2)],'-g')
    end
    plot(pData{3},pData{4}','-b') % original orbit
%     drawnow % experiment to see if this stops writing to main window
end
hold off

extra = 0.03*abs(minY-maxY); % add an extra 3% to limits
axis([ 0 attData.Tperiod minY-extra maxY+extra ])
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                 TRANSITION SEQUENCE ROUTINES                  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function result = GenTransSeq(Eseq, start_pos, exp_rep_pos, checkcycle, no_order, ...
    focusId, caIxMap, ignoreVars, verboseTog)
% GenTransSeq takes curly E set data (sets of actives per observable, over time)
%  that defines a functional network, throws away the time data, and only records
%  changing states (the sets of actives) in the order they appear in the curly E set.
%
% (c) Robert Clewley, Center for Biodynamics, 2002 - 2004

% Update May 2004: added absolute position to cell array
% Update April 2004: fix of possible bug: added filtering of initial set EtransSeq
%   and TS set throughout orbit by `ignoreVars` variables (e.g. passive vars in model)
% Update November 2005: added quick fix for focusId optional restriction to focused vars
%   (essentially the same as ignoreVars -- perhaps should unify them)

%%%%% Other notes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For checkcyle == true, the expected end position of cycle is checked for being a true
%   repeating state vs. that at starting position. No adaptive end matching implemented
%   for checkcycle option, which might be defunct now, since StepNet can now shoot for
%   distant identical states in Eseq. This routine assumes user has a good idea of where
%   repeating state falls in Eseq (since state could conceivably appear multiple times
%   in the sequence but not near the ends of the candidate limit cycle, as starting at
%   start_pos).
% no_order option makes e.g. the ordered sets of actives [1 4] and [4 1], from E_b,
%   equivalent.
% focusId option makes the function return a TS only according to changes in the diagram
%   states for whichever variables are passed (intended to be an observable variable and
%   its immediate dependents).
% ignoreVars set (if non-empty) tells this function to ignore changes in these variables
%   but to include them always in the TS (uses the candidate actives Ix map).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result = {}; % default value
if verboseTog
    disp(' ')
end
switch nargin
% only if pass in varnames as first argument, for debugging purposes
%     case 10
%         % TEMP
    case 9
        % do nothing
    case 8
        verboseTog = false;
    case 7
        verboseTog = false;
        ignoreVars = [];
    case 6
        verboseTog = false;
        ignoreVars = [];
        caIxMap = {};
    case 5
        verboseTog = false;
        focusId = [];
        ignoreVars = [];
        caIxMap = {};
    case 4
        verboseTog = false;
        no_order = true;
        focusId = [];
        ignoreVars = [];
        caIxMap = {};
    case 3
        verboseTog = false;
        focusId = [];
        checkcycle = 0;
        caIxMap = {};
        no_order = true;
        ignoreVars = [];
    otherwise
        disp('GenTransSeq:  Function requires between 3 and 9 arguments')
        return
end

if no_order
    ignoreStr = '';
else
    ignoreStr = 'not ';
end
if verboseTog
    disp( 'GenTransSeq:  Generating transition sequence ...')
    disp(['               ' ignoreStr 'ignoring order of dominant variables at each step'])
end

maxpos = length(Eseq);
if start_pos >= maxpos || exp_rep_pos > maxpos
    disp('GenTransSeq:  Start and end positions larger than size of sequence')
    return
end
if start_pos >= exp_rep_pos
    disp('GenTransSeq:  Start position must be less than end position')
    return
end
if checkcycle
    disp('GenTransSeq:  Checking for cycle repeat at specified end position')
end

Nv = length(Eseq{1}{2}); % this is constant for all positions in the sequence
   % and is the number of observables (external variables) recorded
first_elt = Eseq{start_pos}{2}; % need this later when checking tie up with end of cycle
first_elt_t = Eseq{start_pos}{1}; % initial ABSOLUTE real time of sequence start
emptyIV = isempty(ignoreVars);
% OLD: filter the initial Eseq according to ignoreVars (e.g. to eliminate passive variables that were omitted from focusId)
% NEW: keep the Eseq's intact, do all the filtering only on temporary sets for the transition decisions in the loop
% if ~emptyIV
%     for ev=1:Nv
%         focusedInitial_elt = first_elt{ev};
%         temp = [];
%         for ix=focusedInitial_elt
%             if ~ismember(caIxMap{ev}(ix),ignoreVars)
%                 temp = [temp, ix];
%             end
%         end
%         first_elt{ev} = temp;
%     end
% end
EtransSeq{1} = {first_elt, 0, start_pos}; % second elt is RELATIVE real time of sequence start
Etlen = 1; % initial value of end position of the transition sequence
pos = start_pos;
done = false;
succeed = true; % updated if cycle end matching fails
emptyF = isempty(focusId);

while ~done
    pos=pos+1; % position in Eseq

    % get candidate network state (cell array of Nv cell arrays, each holding lists of active variables)
    Ets_cand = Eseq{pos}{2}; % 2 selects actives from all the E_b sets
    Ets_cand_t = Eseq{pos}{1} - first_elt_t; % RELATIVE (REAL) TIME in sequence [arg. 1 selects the current real time]
    
    % test to see if candidate is different from last recorded distinct network state
    test = zeros(1,Nv);
    for ev=1:Nv
        if emptyF || ismember(ev, focusId)
            Ecand = Ets_cand{ev};
            Elast = EtransSeq{Etlen}{1}{ev};
            if ~emptyIV
                temp_c = [];
                for ix=Ecand
                    if ~ismember(caIxMap{ev}(ix),ignoreVars)
                        temp_c = [temp_c, ix];
                    end
                end
                temp_l = [];
                for ix=Elast
                    if ~ismember(caIxMap{ev}(ix),ignoreVars)
                        temp_l = [temp_l, ix];
                    end
                end
                if no_order
                    Ecand = sort(temp_c);
                    Elast = sort(temp_l);
                else
                    Ecand = temp_c;
                    Elast = temp_l;
                end
%                 if ~isempty(Ecand)
%                     Ets_cand{ev} = Ecand;
%                 end
            end
%             fprintf('\nGenTS: %s for t = %.5f', varnames{ev}, Ets_cand_t)
%             for ix=Ecand
%                 fprintf('   new cand = %s\n',varnames{caIxMap{ev}(ix)})
%             end
%             for ix=Elast
%                 fprintf('   old cand = %s\n',varnames{caIxMap{ev}(ix)})
%             end
            if length(Ecand) == length(Elast)
                if no_order && emptyIV
                    test(ev) = sum(sort(Ecand) ~= sort(Elast));
                else % or already sorted in `if ~emptyIV` above, or no_order == 0
                    test(ev) = sum(Ecand ~= Elast); % ~= returns 0 if sets equal: sum for vector style return type when '==' with vectors
                end
            else
                test(ev) = 1; % failure if sets not same length ... they can't be equal
            end
        end % else don't test this variable, so force it false (not a different state)
    end

    if sum(test) > 0 % compare failed -- therefore next distinct member of sequence found!
        EtransSeq = {EtransSeq{:}, {Ets_cand, Ets_cand_t, pos}};
        Etlen = Etlen + 1;
    end

    % end condition for checkcycle == true (but when not adaptively matching cycle ends)
    if pos == exp_rep_pos % should be 1 + cycle length if checkcycle
        if checkcycle
            test = zeros(1,Nv);
            for ev=1:Nv
                len_cand = length(first_elt{ev});
                len_curr = length(EtransSeq{Etlen}{1}{ev});
                if len_cand == len_curr
                    if no_order
                        test(ev) = (len_cand ~= length(intersect(first_elt{ev},EtransSeq{Etlen}{1}{ev})));
                    else
                        test(ev) = sum(first_elt{ev} ~= EtransSeq{Etlen}{1}{ev}); % ~= returns 0 if sets equal
                    end
                else
                    test(ev) = true;
                end
            end
            if sum(test) > 0 % compare failed -- this isn't the end of the cycle
                fprintf('GenTransSeq:  exp_rep_pos %i is not end of a cycle\n',exp_rep_pos)
                succeed = false;
            end
        end % else just end at specified position!
        done = true;
    end
end

if succeed
    if verboseTog
        if checkcycle
            fprintf('GenTransSeq:  Sequence length (# events) over periodic cycle candidate is %i\n\n',Etlen)
        else
            fprintf('GenTransSeq:  Sequence length (# events) is %i\n\n',Etlen)
        end
    end
    % pos-1 == (exp_rep_end - start_pos) if no adaptive end matching & ends match correctly
    result = EtransSeq; % could argue that shouldn't include final state as it is for an epoch beyond the right mark
    % but it does contain the end time of the TS (i.e. the end time of the last epoch)
else
    if checkcycle
        disp('GenTransSeq:  End point matching failed for candidate cycle. Returned partial sequence ...')
    else
        disp('GenTransSeq:  Unknown problem. Returning partial sequence ...')
    end
    result = EtransSeq;
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function result = CompareTransSeq(Etran1, Etran2, no_order, checkall, padSeq, ...
                                  focusId, verbose_level)
result = {}; % initial and default value
switch nargin
    case 2
        no_order=true;
        checkall=true; % this goes through and notes all differences in the sets
        padSeq=[false, 2, true]; % intelligently pad out shorter sequence to best match the other,
                           % up to # times given by second value, with
                           % mismatched neighbouring-pair swapping
                           % activated by third value
        verbose_level=0; % if this is on then shows the set contents
        focusId = [];
    case 3
        checkall=true;
        padSeq=[false, 2, true];
        verbose_level=0;
        focusId = [];
    case 4
        padSeq=[false, 2, true];
        verbose_level=0;
        focusId = [];
    case 5
        verbose_level=0;
        focusId = [];
    case 6
        verbose_level=0;
    case 7 % ignore
    otherwise
        disp('CompTransSeq:  This function requires between 2 and 7 arguments')
        return
end

lenE1 = length(Etran1);
lenE2 = length(Etran2);
if lenE1 <=1 || lenE2 <= 1
    disp('CompTransSeq:  Only pass sequences with length > 1')
    return
end
if lenE1 ~= lenE2 && ~padSeq(1)
    disp('CompTransSeq:  Only pass equal length sequences unless you set padSeq(1) == TRUE')
    return
end
if no_order && verbose_level >= 1
    disp('CompTransSeq:  Order of dominant variables being ignored')
end

if verbose_level > 0
    disp('CompTransSeq:  Comparing transition sequences ...')
end

Nv = length(Etran1{1}); % this is constant for all positions in the sequence
   % and is the number of observables (external variables) recorded

if ~padSeq(1)
    compResult = DoCompTS(Etran1,Etran2,lenE1,lenE2,Nv,checkall,no_order,focusId,verbose_level);
else
    if lenE1 == lenE2 % needn't have set padSeq(1) to TRUE
        compResult = DoCompTS(Etran1,Etran2,lenE1,lenE2,Nv,checkall,no_order,focusId,verbose_level); % just compare them, since same length
        padSeq(1) = false; % option wasn't needed, and this enables correct 2nd return argument in `result`
	else % cool case
        checkall = true; % force it on
        if lenE2 < lenE1
            sEt  = Etran2; % short
            lEt  = Etran1; % long
            lenS = lenE2;
            lenL = lenE1;
        else
            sEt  = Etran1;
            lEt  = Etran2;
            lenS = lenE1;
            lenL = lenE2;
        end
        padAmount = abs(lenE2-lenE1);
        % padTimes  = padSeq(2);
        % temporary diagnostic choice
        if verbose_level > 0
            disp('CompTransSeq:  Algorithm in development. Will pad up to 2 times only')
        end
        padTimes  = 2;
        padAms    = zeros(1,padTimes);
        padPos    = zeros(1,padTimes);
        bestCandSeq  = { sEt{:} lEt{lenS+1:lenL} }; % initial best candidate padded sequence
        bestNumErrs  = 1000000; % initial value
        bestErrPos   = [];
        bestPads     = {};

        Cresult      = DoCompTS(bestCandSeq, lEt, lenL, lenL, Nv, true, no_order, focusId, 0);
        %%%%% DIAGNOSTIC
        if verbose_level >= 1
            fprintf('CompTransSeq:  DIAGNOSTIC -- pre-padding # errors = %i, at positions\n', Cresult(1));
            fprintf('               ')
            for errnum = 1:Cresult(1)
                fprintf('%4i ',Cresult(errnum+1));
            end
            fprintf('\n')
            fprintf('CompTransSeq:  DIAGNOSTIC -- Padding once...\n')
        end
        %%%%%%%%%%%%%%%%
        % don't have to worry about padPos(1) being out of range of sEt,
        % since the pre-padded bestCandSeq is all lEt from lenS onwards,
        % and so won't have any error positions past the index lenS
        padPos(1) = FindLargestErrRun(Cresult); % was Cresult(2); % ... merely the first error position

        % first time pad out with one padding occurrence
        if padPos(1) + padAmount -1 > lenL
            padAms(1) = lenL + 1 - padPos(1);
        else
            padAms(1) = padAmount;
        end
        if Cresult(1) < bestNumErrs
            bestNumErrs = Cresult(1);
            bestErrPos = Cresult(2:length(Cresult));
            bestPads   = {padPos,padAms};
            if padPos(1) > 1
                newTS  = { sEt{1:padPos(1)-1} lEt{padPos(1):padPos(1)+padAms(1)-1} sEt{padPos(1):lenS} };
            else
                newTS  = { lEt{padPos(1):padPos(1)+padAms(1)-1} sEt{padPos(1):lenS} };
            end
        end
        Cresult = DoCompTS(newTS, lEt, lenL, lenL, Nv, true, no_order, focusId, 0);
        if Cresult(1) < bestNumErrs
            bestNumErrs = Cresult(1);
            bestErrPos = Cresult(2:length(Cresult));
            bestPads   = {padPos,padAms};
            bestCandSeq = newTS;
            %%%%% DIAGNOSTIC
            if verbose_level >= 1
                fprintf('                     Pad #1 inserted lEt{%i:%i} at pos %i\n',padPos(1),padPos(1)+padAms(1)-1,padPos(1));
            end
            %%%%%%%%%%%%%%%%
        end
        if padTimes > 1 && padAmount >= padTimes
        
            %%%%% DIAGNOSTIC
            if verbose_level >= 1
                fprintf('CompTransSeq:  DIAGNOSTIC -- Post-padding #1... # errors = %i, at positions\n', bestNumErrs);
                fprintf('               ')
                for errnum = 1:bestNumErrs
                    fprintf('%4i ',bestErrPos(errnum));
                end
                fprintf('\n')
                disp('CompTransSeq:  DIAGNOSTIC -- Looking at second padding ...')
            end
            %%%%%%%%%%%%%%%%
            padPos(2) = FindLargestErrRun(Cresult); % was Cresult(2); %... merely the first of remaining error positions
			for shiftAm = 1 : padAmount-1;
                if padPos(2) + shiftAm -1 > lenL
                    actualShift = lenL + 1 - padPos(2);
                    padAms(2) = actualShift;
                else
                    actualShift = shiftAm;
                    padAms(2) = actualShift;
                end
                padAms(1) = padAmount - actualShift;
                if padPos(1) > 1
                    newTS = { sEt{1:padPos(1)-1} };
                else
                    newTS = {};
                end
                for pos = 1:padTimes
                    newTS = { newTS{:} lEt{padPos(pos):padPos(pos)+padAms(pos)-1} };
                    if pos < padTimes
                        newTS = { newTS{:} sEt{padPos(pos):padPos(pos+1)-1} };
                    else
                        newTS = { newTS{:} sEt{padPos(pos):lenS} };
                    end
                end
                Cresult = DoCompTS(newTS, lEt, lenL, lenL, Nv, true, no_order, focusId, 0);
                if Cresult(1) < bestNumErrs
                    bestNumErrs = Cresult(1);
                    bestErrPos = Cresult(2:bestNumErrs+1);
                    bestPads = {padPos,padAms};
                    bestCandSeq = newTS;
                    %%%%% DIAGNOSTIC
                    if verbose_level >= 1
                        fprintf('CompTransSeq:  DIAGNOSTIC -- Pad #2: inserted lEt{%i:%i} @ %i\n',padPos(1),padPos(1)+padAms(1)-1,padPos(1))
                        fprintf('                                         and lEt{%i:%i} @ %i\n',padPos(2),padPos(2)+padAms(2)-1,padPos(2))
                    end
                    %%%%%%%%%%%%%%%%
                end
            end % for 

        end % if padTimes > 1 & ...
            
        %%%%% DIAGNOSTIC
        if verbose_level >= 1
            fprintf('CompTransSeq:  DIAGNOSTIC -- post-padding #2... # errors = %i, at positions\n', bestNumErrs);
            fprintf('               ')
            for errnum = 1:bestNumErrs
                fprintf('%4i ',bestErrPos(errnum));
            end
            fprintf('\n')
        end
        %%%%%%%%%%%%%%%%
            
        if padSeq(3) % neighbouring-pair switching is ON, so go further ...
            tryAgain = true;
            while tryAgain
                if bestNumErrs > 1
                    for errnum = 1:bestNumErrs-1
                        switchPosCand = bestErrPos(errnum); % candidate position for switching values around
                        if bestErrPos(errnum+1) - switchPosCand == 1 % then the errors are at neighbouring positions
                            newTS = { bestCandSeq{1:switchPosCand-1} bestCandSeq{switchPosCand+1} ...
                                      bestCandSeq{switchPosCand} bestCandSeq{switchPosCand+2:lenL} }; % switched pair!
                            Cresult = DoCompTS(newTS, lEt, lenL, lenL, Nv, true, no_order, focusId, 0);
                            if Cresult(1) < bestNumErrs
                                %%%%% DIAGNOSTIC
                                if verbose_level >= 1
                                    fprintf('CompTransSeq:  DIAGNOSTIC -- Successfully switched positions %i and %i\n',switchPosCand,switchPosCand+1)
                                    fprintf('CompTransSeq:  DIAGNOSTIC -- After switching, # errors = %i, at positions\n', Cresult(1));
                                    fprintf('               ')
                                    for errnum = 1:Cresult(1)
                                        fprintf('%4i ',Cresult(errnum+1));
                                    end
                                    fprintf('\n')
                                end
                                %%%%%%%%%%%%%%%%
                                bestNumErrs = Cresult(1);
                                bestErrPos = Cresult(2:length(Cresult));
                                % bestPads = {padPos,padAms};
                                bestCandSeq = newTS;
                                tryAgain = true;
                                break
                            end % if Cresult(1)...
                        end % if bestErrPos...
                    end % for errnum
                    tryAgain = false; % no switching helped
                else
                    %%%%% DIAGNOSTIC
                    if verbose_level >= 1
                        disp('CompTransSeq:  DIAGNOSTIC -- Neighbouring-pair switching cannot be invoked... already only 1 or fewer errors')
                    end
                    %%%%%%%%%%%%%%%%
                    tryAgain = false;
                end % if
            end % while
        end % if padSeq(3)

        %for padTime = 2:padTimes
        %    for shiftAm = 1 : padAmount-1;
        %        padAms(padTime) = shiftAm;
        %        for i = padTime:2
        %            padAms(i-1) = padAms(i-1) - (shiftAm-SOMTHING);
        %        end
        %         padBit = 
        %    end
        %end
            
    end
end

if padSeq(1)
    result = {[bestNumErrs, bestErrPos], bestCandSeq, bestPads}; % bestPads is {padPos,padAms}
else
    result = {compResult, {}, {}}; % compResult is [errors errpos]
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Secondary functions to CompareTransSeq()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function result = DoCompTS(Etran1,Etran2,lenE1,lenE2,Nv,checkall,no_order,foc,verbose_level)
errors = 0;
errpos = [];
result = [];
actualLenE1 = length(Etran1);
actualLenE2 = length(Etran2);
if actualLenE1 ~= lenE1 || actualLenE2 ~= lenE2
	disp('DoCompTS:  Internal error -- mismatch between expected and actual sequence lengths passed')
end
if lenE1 ~= lenE2
    disp('DoCompTS:  Internal error -- transition sequences passed here must have same length')
    return
end
emptyfoc = isempty(foc);
for p = 1:lenE1
	E1p = Etran1{p}; % Nv cell arrays in these two objects
    E2p = Etran2{p};
    test = zeros(1,Nv);
    for ev=1:Nv
        if emptyfoc || ismember(ev, foc)
            len1 = length(E1p{1}{ev});
            len2 = length(E2p{1}{ev});
            if len1 == len2
                if no_order
                    test(ev) = (len1 ~= length(intersect(E1p{1}{ev},E2p{1}{ev})));
                else
                    test(ev) = sum(E1p{1}{ev} ~= E2p{1}{ev}); % ~= returns 0 if sets equal: sum since ~= returns vector when comparing vectors
                end
            else
                test(ev) = 1; % failure if sets not equal length
            end
        end  % else remains 0
    end
    if sum(test) > 0
        if verbose_level == 2
            fprintf('CompTransSeq:  DIAGNOSTIC -- Compare failed at position %i\n',p)
            disp(   '               sequence values at this position:')
        end
        errors = errors + 1;
        errpos = [errpos, p];
        if checkall
            if verbose_level == 2
                celldisp(E1p{1})
                celldisp(E2p{1})
                disp('******************************************')
            end
        else
            break % as soon as 1 error found if ~checkall
        end
    end
end
result = [errors, errpos];
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function epos = FindLargestErrRun(Cresult,nth)
if nargin == 1
    nth = 1; % i.e. 1st largest (can get 2nd, 3rd etc.)
end
numerrs = Cresult(1);
errs = Cresult(2:numerrs+1);
longestRunLen = 0;
longestRunPos = 1;
lastErrPos = -1; % initial value trick to avoid additional check in loop when ep==1 and errs(1) might be 1
allRunPos = [];
currRunLen = 0;
currRunPos = 0;
for ep = 1:numerrs
    if errs(ep) ~= lastErrPos + 1
        % check if last run (that has just ended unless ep==1) was longest so far
        if currRunLen > longestRunLen
            allRunPos = [currRunPos allRunPos];
            longestRunLen = currRunLen;
            longestRunPos = currRunPos;
        end
        % start new run
        currRunPos = errs(ep);
        currRunLen = 1;
    else
        currRunLen = currRunLen + 1;
        % now catch case: if longest run goes all the way to the end,
		% then update longest
        if ep == numerrs
            if currRunLen > longestRunLen
                allRunPos = [currRunPos allRunPos];
                longestRunLen = currRunLen;
                longestRunPos = currRunPos;
            end
        end
    end
    lastErrPos = errs(ep);
end
numruns = length(allRunPos);
if nth > numruns
    epos = allRunPos(numruns);
else
    epos = allRunPos(nth);
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function ViewTransSeq( TSfigHandle, TSlist, varnames, ignoreOrder, numInt, numExt, ...
    actsIxMap, potsIxMap )
% ViewTransSeq displays the contents of transition sequences graphically.
%
% Key to horizontal `bar` display of viewed transition sequences:
% - Contiguous blue lines indicate original portions of the sequence.
% - Broken lines in black indicate portions that were padded by the
%   sequence matching algorithm in order to minimise the number of
%   errors between sequence members over the length of the sequences, and
%   to make the sequences the same length.
% - Remaining errors in sequences generated by the matching algorithm
%   (see StepNet) for a sequence are shown by red markers at their
%   sequence positions.
%
% (c) Robert Clewley, Center for Biodynamics, Boston University, Feb 2003

% Last update: November 2005 -- arrange var names according to TS focus set

fprintf('\n')
if nargin ~= 9
    beep
    disp('ViewTransSeq:  Internal error! Wrong number of arguments passed. Expected 9')
    disp('                 Exiting ViewTransSeq without displaying anything')
    return
end

%%%% CONSTANTS
global ESEQ_TIMES ESEQ_ACTIVES_BAR ESEQ_ACTIVES ESEQ_POTENTS_BAR ESEQ_POTENTS
global TS_TSEQ TS_NAME TS_VALX TS_VALY TS_ERRS TS_PADS TS_FOCUS

TS_NAME_LIMIT = 14; % # of characters to fit in name field of Trans Seq window
EV_NAME_LIMIT = 5; % # of characters to fit in external var name field in Eseq window
EXTVAR_LIMIT = 5;  % Max # of external variables that can be viewed at once in the Eseq window
axesXadjust = 0.02; % extra space at the ends of X axes for Time Box and Trans Seqs box

%%%% Some initializations
numTS    = length(TSlist);
TSlength = length(TSlist{1}{TS_TSEQ}); % same for all Trans Seqs
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% BEGIN { Initialize figure window }
handTSaxes = subplot('Position', [0.15 0.7 0.8 0.27]);
% axis(handTSaxes,'off')
hold on
set(handTSaxes,'Color', [0.8 0.8 0.8])
set(handTSaxes,'XLimMode','manual')
set(handTSaxes,'YLimMode','manual')
set(handTSaxes,'YLim', [ 0.5 numTS+0.5 ])
set(handTSaxes,'XLim', [ 1-axesXadjust TSlength+axesXadjust ])
set(handTSaxes,'XTickLabelMode','manual')
set(handTSaxes,'YTickLabelMode','manual')
set(handTSaxes,'YTick',[1:numTS])
set(handTSaxes,'YDir','reverse')
% show only some fraction of ticks (minimum every 2) for every 20 trans seq entries -- in multiples of 2
tickFrac = 2*ceil(TSlength/20);
numXticks = ceil(TSlength/tickFrac)+1;
TS_XtickNames = cell(1,numXticks);
TS_XtickVals(1) = 1; % ensure '1' is first tick (MATLAB's auto mode seems to always drop this off)
TS_XtickNames{1} = '1';
for x = 1:numXticks
    TS_XtickVals(x+1) = x*tickFrac;
    TS_XtickNames{x+1} = num2str(x*tickFrac);
end
% ensure final value is the last number in the sequence, so overwrite last
% entry computed in previous loop
TS_XtickVals(numXticks) = TSlength;
TS_XtickNames{numXticks} = num2str(TSlength);
set(handTSaxes,'XTick',TS_XtickVals);
set(handTSaxes,'XTickLabel',TS_XtickNames)

% time box
handTBaxes = subplot('Position', [0.15 0.62 0.8 0.03]);
hold on
tBarHand = Draw_tBarTSview(1-axesXadjust,TSlength+axesXadjust,'k');
tTickHand = Draw_tTickTSview( 1, 'k' );
axis(handTBaxes,'off')
axis tight

% Eseq active/potentials sets box
handESaxes = subplot('Position', [0.05 0.05 0.9 0.5] );
% set(handESaxes,'PlotBoxAspectRatioMode','manual') %  NO! don't do this!
set(handESaxes,'XLimMode','manual')
set(handESaxes,'YLimMode','manual')
set(handESaxes,'DrawMode','fast')
axis(handESaxes,'off')
set(handESaxes,'YDir','reverse')
hold on
plot( [0 0], [0 1], 'k-')
plot( [1 1], [0 1], 'k-')
plot( [0 1], [1 1], 'k-')
plot( [0 1], [0 0], 'k-')
axis tight

% Double buffering setup for smoother graphics
set(TSfigHandle,'RendererMode','manual')
set(TSfigHandle,'Renderer','painters')
set(TSfigHandle,'DoubleBuffer','on')
%%%% END { Initialize figure window }

%%%% BEGIN { Initialize Trans Seqs sub-plot }
TSnameStrs = cell(1,numTS);
handTSnames = zeros(1,numTS);
for ts = 1:numTS
    TSname = TSlist{ts}{TS_NAME};
    if length(TSname) > TS_NAME_LIMIT
        TSnameStrs{ts} = TSname(1:TS_NAME_LIMIT);
        disp('ViewTransSeq:  Trans sequence names truncated in display')
    else
        TSnameStrs{ts} = TSname;
    end
    TSpads = TSlist{ts}{TS_PADS};
    if ~isempty(TSpads)
        tempNumPads = length( TSpads{1} );
        numPads = tempNumPads; % initial value
        for padPos = 1:tempNumPads
            if TSlist{ts}{TS_PADS}{1}(padPos) == 0
                numPads = numPads - 1; % substract one for every unused padding
            end
        end
    else
        numPads = 0;
    end
    if ~isempty(TSpads)
        if TSpads{1}(1)==1
            numParts = 2*numPads;
        else
            numParts = 2*numPads + 1; % i.e. normal part, padded part, normal part, etc. (NB may start with a padded part!)
        end
    else
        numParts = 1;
    end
    TSparts = cell(1,numParts);
    TStypeList = zeros(1,TSlength);
	inPad = 0;
    if numPads > 0
		for TSpos = 1:TSlength
            for padPos = 0:numPads
                if padPos == 0
                    if TSpos < TSpads{1}(padPos+1)
                        TStypeList(TSpos) = not(inPad);
                    end
                elseif padPos > 0 && padPos < numPads
                    if TSpos < TSpads{1}(padPos+1) && TSpos > TSpads{1}(padPos) + TSpads{2}(padPos) - 1
                        TStypeList(TSpos) = not(inPad);
                    end
                else % padPos == numPads
                    if TSpos > TSpads{1}(padPos) + TSpads{2}(padPos) - 1
                        TStypeList(TSpos) = not(inPad);
                    end
                end
            end % for padPos
		end % for TSpos
    else
        TStypeList = ones(1,TSlength);
    end % if numPads > 0
    TSpos = 1;
	for part = 1:numParts
        oldTSpos = TSpos;
        if TSpos == TSlength
            break
        end
        breakWhile = false;
        while TStypeList(TSpos+1) == TStypeList(oldTSpos) % then ~inPad
            TSpos = TSpos + 1;
            if TSpos == TSlength
                breakWhile = true;
                break % out of while loop
            end
        end
        TSposBeg = oldTSpos;
        TSposEnd = TSpos;
        TSparts{part} = [ TStypeList(oldTSpos), TSposBeg, TSposEnd ];
        TSpos = TSpos + 1;
        if breakWhile
            break % out of for loop
        end
	end % for part
    % errors
    numErrs = 0; % default
    if ~isempty(TSlist{ts}{TS_ERRS})
        numErrs = length(TSlist{ts}{TS_ERRS});
        TSerrs = zeros(1,numErrs);
        for errnum = 1:numErrs
            TSerrs(errnum) = TSlist{ts}{TS_ERRS}(errnum);
        end
    end
    axes(handTSaxes)
    handTSbars  = Draw_TSBar( TSparts, numParts,ts ); % this returns vectors (of length 2 * numParts)
    if numErrs > 0
        handTSerrs  = Draw_TSErr( TSerrs, ts, numTS );
    end
end % for ts
set(handTSaxes,'YTickLabel',TSnameStrs)
%%%% END { Initialize Trans Seqs sub-plot }

%%%% BEGIN { Initialize Eseq sub-plot }
% organize varnames according to TS_FOCUS set, followed by remainder
% -- use the largest TS_FOCUS set present in the TS's given
% if all TS_FOCUS sets are empty this will fall through to making
% varnames_ordered = varnames
lenfoc = 0;
TSfocus = [];
for ts = 1:numTS
    if length(TSlist{ts}{TS_FOCUS}) > lenfoc
        TSfocus = TSlist{ts}{TS_FOCUS};
        lenfoc = length(TSfocus);
    end
end
if isempty(TSfocus)
    varnames_ordered = varnames;
else
    varnames_ordered = varnames{TSfocus};
    % fill in rest of varnames_ordered from remaining varnames
    for vi = 1:length(varnames)
        v = varnames{vi};
        if ~ismember(v, varnames_ordered)
            varnames_ordered = {varnames_ordered{:}, v};
        end
    end
end
if numExt > EXTVAR_LIMIT
%     fprintf('ViewTransSeq:  Can only view up to %i external variables at a time\n',EXTVAR_LIMIT)
%     disp(   '               Enter a space separated list of up to this number of')
%     disp(   '                 external variable indices in dialog box (see list below)')
%     disp(   '               Variable names of indices appearing in the transition sequences:')
%     for i=1:length(varnames)
%         fprintf('                   %3i  %s\n',i,varnames_ordered{i})
%     end
%     %%%% TEMPORARY
%     disp(   '               Dialog box not yet implemented ... truncating to max # external vars in display')
%     beep
    numExt = EXTVAR_LIMIT;
end
shownTruncateEV = false; % flag so that only displays that truncated names in Eseq window once
divEvX = InitESbox(handESaxes,numInt,numExt,varnames_ordered,EV_NAME_LIMIT,numTS);
%%%% END { Initialize Eseq sub-plot }

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%         Main control loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ignoreOrder
    iO_str = 'ON';
else
    iO_str = 'OFF';
end
disp(['ViewTransSeq:  `ignore order of variables` switch is ' iO_str])
disp('ViewTransSeq:  Error positions in compared sequences (for the selected focused variable) are shown by red markers in TS viewbox')
keypress = 'H'; % initialized to 'help' to display commands info at beginning
pos = 1; % initial value
doQuit = false;
while ~doQuit
    switch keypress
        case 'H' % help
            fprintf('\n')
            disp(' **********************************************************************************')
			disp('   ViewTransSeq single key commands  (active only over ViewTransSeq figure window)')
            disp(' *************************************************************************************')
			disp(' Epoch change:    N(ext)  B(ack)   , (+10)  . (-10)  SPACEBAR (+1)')
            disp('                  < (go to start)  > (go to end)')
			disp(' General:         Q (Return to functional net)   S(ave figure)   H(elp)')
            disp('                  F(etch names of indices)       A (toggle ignore var order)' )
            disp('                  X (show sequence info)         * (colour key & about ViewTransSeq)')
            disp(' *************************************************************************************')
            fprintf('\n')
        case '*'
            fprintf('\n')
            help ViewTransSeq
        case {'N', ' '} % forward 1
            pos = pos + 1;
        case 'B' % back 1
            pos = pos - 1;
        case '.' % back 10
            pos = pos + 10;
        case ',' % forward 10
            pos = pos - 10;
        case '<' % to beginning
            pos = 1;
        case '>' % to end
            pos = TSlength;
        case 'S' % save figure 'screenshot'
            saveas(TSfigHandle,strcat('Figures/TransSeq_Step',num2str(pos),'.fig'))
            disp('ViewTransSeq:  Figure saved in DSSRT root directory')
        case 'Q' % return control to StepNet
            % ignore (dealt with elsewhere)
        case 'X'
            fprintf('ViewTransSeq:  Transition sequences in memory 1 -> %i:\n',numTS)
            for ts = 1:numTS
                currTS = TSlist{ts};
                fprintf(    '           TransSeq #%i: %s\n', ts, currTS{TS_NAME});
                fprintf(    '                 Variable X value: %3.4f\n', currTS{TS_VALX});
                fprintf(    '                 Variable Y value: %3i\n', currTS{TS_VALY});
                fprintf(    '                 Length:           %3i\n',length(currTS{TS_TSEQ}));
                if ~isempty( currTS{TS_ERRS} )
                    fprintf('                 Error positions vs. original compared sequence:\n')
                    fprintf('                  ')
                    for errnum = 1:currTS{TS_VALY}
                        fprintf('%4i ',currTS{TS_ERRS}(errnum)); % error positions in the remainder of array result{1}
                    end
                    fprintf('\n')
                end
                if ~isempty( currTS{TS_PADS} )
                    fprintf('                 Padded at positions:\n')
                    fprintf('                  ')
                    for padPos = 1:length(currTS{TS_PADS}{1})
                        if currTS{TS_PADS}{1}(padPos) > 0 % unused paddings entered as 0
                            fprintf('%3i x%3i,',currTS{TS_PADS}{1}(padPos),currTS{TS_PADS}{2}(padPos));
                        end
                    end
                    fprintf('\n')
                end
            end
        case 'F' % fetch names of set indices
            disp(       'ViewTransSeq:  Variable names of indices appearing in the transition sequences:')
            for i=1:length(varnames)
                fprintf('                 %3i  =  %s\n',i,varnames_ordered{i})
            end
            fprintf('\n')
        case 'A'
            ignoreOrder = not(ignoreOrder);
            if ignoreOrder
                iO_str = 'ON';
			else
                iO_str = 'OFF';
			end
			disp(['ViewTransSeq:  `ignore order of variables` switch is ' iO_str])
        otherwise
            disp('ViewTransSeq:  Command not recognised')
    end

    if pos < 1
        pos = 1;
    elseif pos > TSlength
        pos = TSlength;
    end

    figure(TSfigHandle); % ensure our figure is current
    drawnow

    %%%% BEGIN { UPDATE Eseqs BOX }
    axes(handESaxes)
    cla
    % divEvX doesn't need updating, so ignore returned val
    InitESbox(handESaxes,numInt,numExt,varnames_ordered,EV_NAME_LIMIT,numTS);
    title(['Dominant variables at sequence position ' num2str(pos,'%-4i')]) % update step position in title
    currTSeqs = cell(1,numTS);
    for ts = 1: numTS
        tempTS = TSlist{ts}{TS_TSEQ}{pos};
        if ignoreOrder
            for ev = 1:numExt
                tempTS{1}{ev} = sort(tempTS{1}{ev});
            end
        end
        currTSeqs{ts} = tempTS;
    end
    for ev = 1:numExt
        for ts = 1:numTS
            for entry = 1:length( currTSeqs{ts}{1}{ev} )
                nameVar = varnames{ actsIxMap{ev}(currTSeqs{ts}{1}{ev}(entry)) }; % was varnames( ) but hadn't flagged error?
                if length(nameVar) > EV_NAME_LIMIT
                    nameVar = nameVar(1:EV_NAME_LIMIT);
                    if ~shownTruncateEV
                        disp('ViewTransSeq:  ext. var. names truncated in display')
                        shownTruncateEV = true;
                    end
                end
                text( GetESxpos(divEvX,numTS,ev,ts), 0.11+entry*0.11, nameVar )
            end
        end
    end
    text( 0.6, 0.9, ['(relative) real time = ' num2str(tempTS{2}) ], 'Color', 'k', 'FontWeight','demi' )
    %%%% END { UPDATE Eseqs BOX }

    %%%% BEGIN { UPDATE TIME BAR BOX }
    axes(handTBaxes)
    set(tTickHand, 'XData', [ pos pos ], 'Color', 'b');
    axes(handESaxes)
    %%%% END  { UPDATE TIME BAR BOX }

    %%%% wait for user command over viewTS figure
    keyokay = false;
    while ~keyokay
        w = 0;
        while w~=1
            w = waitforbuttonpress;
        end
        if ~ishandle(TSfigHandle)
            doQuit = true;
            break
        else
            keypress = get(TSfigHandle,'CurrentCharacter');
            if ischar(keypress) & ~ismember(keypress, [27, 8, 13, 14, 15])
                keypress = upper(keypress);
                keyokay = true;
            end
        end
	end
    if keypress == 'Q'
        disp('ViewTransSeq:  Really return to StepNet? (see dialog box)')
        ButtonName = questdlg('Really return to StepNet?','Return?','Yes','No','No');
        switch ButtonName
            case 'Yes'
                doQuit = true;
                disp('               Returning ...')
            otherwise
                disp('               Return cancelled')
        end
    end
    
end % while

fprintf('\n')
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                  Secondary functions to ViewTransSeq()                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handleTbar = Draw_tBarTSview(x1,x2,colStr)
handleTbar = plot( [x1 x2], [0.01 0.01], [ colStr '-' ], 'LineWidth', 2 );
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handleTtick = Draw_tTickTSview(xpos,colStr)
% xpos in units of real time (or [0,1] if a variable is undefined)
handleTtick = plot( [xpos xpos], [0.008 0.012], [ colStr '-' ], 'LineWidth', 2);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = Draw_TSErr( TSerrs, ypos, numTS )
handles = zeros(length(TSerrs),1); % initial value
for errnum = 1:length(TSerrs)
    handles(errnum) = plot( [ TSerrs(errnum) TSerrs(errnum) ], [ypos-0.03*numTS ypos+0.03*numTS], 'r-');
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = Draw_TSBar( TSparts, numParts, ypos )
Xparts = zeros(numParts,2);
Yparts = zeros(numParts,2);
handles = zeros(numParts*2);
for p = 1:numParts
    partType = TSparts{p}(1);
    if partType == 0
        styleStr = 'k:';
    else
        styleStr = 'b-';
    end
    handles(p) = plot( [TSparts{p}(2) TSparts{p}(3)], [ypos ypos], styleStr, 'LineWidth', 2 );
    handles(p+numParts) = plot( [TSparts{p}(2) TSparts{p}(3)], [ypos ypos], ['.' styleStr(1)], 'LineWidth', 2 );
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function divEvX = InitESbox(handESaxes,numInt,numExt,varnames,EV_NAME_LIMIT,numTS)
global TS_NAME
axes(handESaxes)
text(0.001,0.07,'ext. var.','FontWeight','demi')
text(0.001,0.13,'t.s. #','FontWeight','demi')
plot([0.001 0.998],[0.16 0.16], 'k-', 'LineWidth', 2)
divEvX = (0.93-0.07)/numExt;
numExtStrs = cell(1,numExt);
for ev = 1:numExt
    EVname = varnames{ev};
    if length(EVname) > EV_NAME_LIMIT
        numExtStrs{ev} = EVname(1:EV_NAME_LIMIT);
    else
        numExtStrs{ev} = EVname;
    end
    xpos = 0.11 + 1.16*divEvX*(ev-1);
    for ts = 1:numTS
        text( GetESxpos(divEvX,numTS,ev,ts), 0.13, num2str(ts,'%3i'), 'Color', 'b' )
    end
    text(xpos, 0.07, numExtStrs{ev}, 'Color', 'k', 'FontWeight','demi')
    plot([xpos+0.6*1.16*divEvX xpos+0.6*1.16*divEvX], [0.04 0.8], 'k-')
end
return

%%%%%%%%%%%%%%%%%%%%%%%

function xpos = GetESxpos(divEvX,numTS,ev,ts)
xpos = 0.09 + 0.07*divEvX*numTS + 1.18*divEvX*(ev-1) + 0.31*divEvX*(ts - 1 -(numTS-1)/2);
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%               StepNet  DRAWING ROUTINES                   %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function handleTbar = Draw_tBar(x1,x2,colStr)
handleTbar = plot( [x1 x2], [0.01 0.01], [ colStr '-' ] );
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handleTtick = Draw_tTick(xpos,colStr)
% xpos in units of real time (or [0,1] if a variable is undefined)
handleTtick = plot( [xpos xpos], [0.008 0.012], [ colStr '-' ]);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handleMark = Draw_mark(xpos)
% for markers
handleMark = plot( [xpos xpos], [0.009 0.011], 'r', 'LineWidth', 2);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handlePoints = Draw_points(xposList)
% for TS list
handlePoints = plot( xposList, 0.001, 'r.');
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PlotVBox(handVVaxes,seqTimes,var,VIx,varnames)
axes(handVVaxes)
plot(seqTimes,var)
xLims = get(handVVaxes, 'XLim');
xLims(1) = floor(xLims(1));
xLims(2) = ceil(xLims(2));
set(handVVaxes, 'XLim', xLims)
axis tight
yLims = get(handVVaxes, 'YLim');
if yLims(2)-yLims(1) > 0.2
    yLims(1) = floor(yLims(1));
    yLims(2) = ceil(yLims(2));
else % got really small changes, so zoom in just a little
    yLims(2) = 0.2 + yLims(1);
end
set(handVVaxes, 'YLim', yLims)
set(handVVaxes, 'YTickMode','manual')
set(handVVaxes, 'YTick', [yLims(1); yLims(1)+(yLims(2)-yLims(1))/2; yLims(2)])
if ~isempty(varnames)
	if ~isempty(varnames{VIx})
        title(['Variable ' varnames{VIx} ], 'Interpreter', 'none')
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = DrawNet( obj, init )
global NODE LINK OBJ_TYPE
if nargin == 1
    init = false;
end
if obj{OBJ_TYPE} == NODE
    h = DrawNode( obj, init );
else
    h = DrawLink( obj, init );
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = DrawLink(obj,init)
global OBJ_STATE OBJ_HANDLE OBJ_TYPE OBJ_DEPVAR OBJ_ACTSMAP OBJ_POTSMAP OBJ_STATELT OBJ_LABEL
global OBJ_OCOORD OBJ_LCOORD OBJ_OBJVAR LINKOFFSTATE
% OBJ_OBJVAR for link objects only
if nargin < 2
    init = 0;
end
state     = obj{OBJ_STATE};
objCoords = obj{OBJ_OCOORD};
h         = obj{OBJ_HANDLE};
if init
    return % no initialization for link objects
end
if h == 0 % no figure object handle set
    if state ~= LINKOFFSTATE
        h=draw_arr(objCoords(1), objCoords(2), objCoords(3), objCoords(4), 'g-');
        set(h,'Visible','off')
    end % else leave alone
end
switch state
    case 0
        set(h,'Visible','off')
    case 1
        set(h, 'Color', 'g');
        set(h,'Visible','on')
%         set(h,'Visible','off') % swap commenting here to suppress potentials in func net figure
    case 2
        set(h, 'Color', 'r');
        set(h,'Visible','on')
    otherwise
        % e.g. when state == LINKOFFSTATE
        h=0;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = DrawNode(obj,init)
global OBJ_STATE OBJ_HANDLE OBJ_TYPE OBJ_DEPVAR OBJ_ACTSMAP OBJ_POTSMAP OBJ_STATELT OBJ_LABEL
global OBJ_OCOORD OBJ_LCOORD OBJ_OBJVAR LINKOFFSTATE
% OBJ_OBJVAR for link objects only
if nargin == 1
    init = 0;
end
state       = obj{OBJ_STATE};
h           = obj{OBJ_HANDLE};
objCoords   = obj{OBJ_OCOORD};
objStLT     = obj{OBJ_STATELT};
labelStr    = obj{OBJ_LABEL};
labelCoords = obj{OBJ_LCOORD};
numStates   = length(objStLT)+1; % total number of possible object states
if init
    draw_text(labelCoords(1), labelCoords(2), labelStr, labelCoords(3), 0, 'b-', 2);
    return
end
if h==0 % then drawing initial state (init not set)
    h=draw_circ(objCoords(1), objCoords(2), objCoords(3), -180, 180, 'b:', 2);
else
    for s = 0:numStates-1
    	if state == s
            if s == 0
                set(h,'Color', 'b', 'LineStyle', ':');
            else
                set(h,'Color', objStLT{s}(1), 'LineStyle', objStLT{s}(2));
            end
        end
	end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following drawing routines by F. Auger, March 1999 with contributions as noted.
% Send any error you find or any comment via e-mail to : f.auger@ieee.org
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function harr=draw_arr(x1,y1,x2,y2,style)
xnorm=[0 1 0.9 1  0.9];
ynorm=[0 0 0.08 0 -0.08];

teta=atan2(y2-y1,x2-x1);
r=sqrt((x2-x1)^2+(y2-y1)^2);        

xcoord=xnorm*cos(teta)-ynorm*sin(teta);
ycoord=xnorm*sin(teta)+ynorm*cos(teta);
harr=plot(x1+xcoord*r, y1+ycoord*r,style,'LineWidth',2);
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h=draw_circ(xc, yc, r, theta1, theta2, style, thick)
% Updated to include line thickness by R. Clewley, Dec 2002
% set N=12 for fewer polygon sides - to improve plotting speed
if nargin == 6
    thick=1;
elseif nargin == 3
    thick=1;
    theta1=-180;
    theta2=180;
    style='b-';
elseif nargin < 3 || nargin == 4 || nargin == 5
        disp('Not enough paramaters to draw_circ')
        return
end

N=14; theta=pi*linspace(0,theta2-theta1,N)/180.0;
xnorm=r*cos(theta);
ynorm=r*sin(theta);
xcoord=xnorm*cos(pi*theta1/180)-ynorm*sin(pi*theta1/180);
ycoord=xnorm*sin(pi*theta1/180)+ynorm*cos(pi*theta1/180);
h=plot(xc+xcoord, yc+ycoord, style ,'LineWidth',thick);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_char(x,y,TheChar,scale,angle,style,thick)
% drawchar(x,y,TheChar,scale,angle,style,thick) draws a character 
% at the (x,y) point, with scale as desired scaling factor 
% and angle (in degrees) as orientation. 
% style is the third plot parameter (see the help of plot).
%
% Includes also several contributions of Fred M. Staudaher (freds@packet.net)
% Updated to include line thickness by R. Clewley, Dec 2002

if nargin==0, help draw_char; return; end;
if (nargin<3), error('At least 3 parameters required.');
elseif (nargin==3), scale=1; angle=0; style='-'; thick=1;
elseif (nargin==4), angle=0; style='-'; thick=1;
elseif (nargin==5), style='-'; thick=1;
elseif (nargin==6), thick=1;
end;

TheChar=abs(TheChar);
teta=pi*angle/180.0;
if (TheChar==33),                                     % !
 % contributed by Fred M. Staudaher (freds@packet.net)
 xpoints=[0.40 0.50 0.50 0.40 0.40];
 ypoints=[0.05 0.05 0.15 0.15 0.05];
 xcoord=xpoints*cos(teta)-ypoints*sin(teta);
 ycoord=xpoints*sin(teta)+ypoints*cos(teta);
 plot(x+xcoord*scale, y+ycoord*scale,style,'LineWidth',thick);
 xpoints=[0.45 0.45]; 
 ypoints=[0.4 1];
elseif (TheChar==34),                                 % "
 % contributed by Fred M. Staudaher (freds@packet.net)
 xpoints=[0.4 0.4];
 ypoints=[0.75 0.95];
 xcoord=xpoints*cos(teta)-ypoints*sin(teta);
 ycoord=xpoints*sin(teta)+ypoints*cos(teta);
 plot(x+xcoord*scale, y+ycoord*scale,style,'LineWidth',thick);
 xpoints=[0.6 0.6];
 ypoints=[0.75 0.95];
elseif (TheChar==35),                                 % #
 xpoints=[0.2 0.8];
 ypoints=[0.6 0.6];
 xcoord=xpoints*cos(teta)-ypoints*sin(teta);
 ycoord=xpoints*sin(teta)+ypoints*cos(teta);
 plot(x+xcoord*scale, y+ycoord*scale,style,'LineWidth',thick);
 xpoints=[0.2 0.8]; 
 ypoints=[0.4 0.4];
 xcoord=xpoints*cos(teta)-ypoints*sin(teta);
 ycoord=xpoints*sin(teta)+ypoints*cos(teta);
 plot(x+xcoord*scale, y+ycoord*scale,style,'LineWidth',thick);
 xpoints=[0.3 0.5]; 
 ypoints=[0.1 0.9];
 xcoord=xpoints*cos(teta)-ypoints*sin(teta);
 ycoord=xpoints*sin(teta)+ypoints*cos(teta);
 plot(x+xcoord*scale, y+ycoord*scale,style,'LineWidth',thick);

 xpoints=[0.5 0.7]; 
 ypoints=[0.1 0.9];
elseif (TheChar==36),                                 % $
 % contributed by Fred M. Staudaher (freds@packet.net)
 xpoints=[0.45 0.45]; 
 ypoints=[0.1 1.0];
 xcoord=xpoints*cos(teta)-ypoints*sin(teta);
 ycoord=xpoints*sin(teta)+ypoints*cos(teta);
 plot(x+xcoord*scale, y+ycoord*scale,style,'LineWidth',thick);
 xpoints=[0.55 0.55]; 
 ypoints=[0.1 1.0];
 xcoord=xpoints*cos(teta)-ypoints*sin(teta);
 ycoord=xpoints*sin(teta)+ypoints*cos(teta);
 plot(x+xcoord*scale, y+ycoord*scale,style,'LineWidth',thick);
 xpoints=    [0.9 0.7 0.3 0.1 0.1 0.3 0.7 0.9 0.9 0.7 0.3 0.1]; 
 ypoints=0.7*[0.8 1.0 1.0 0.8 0.7 0.5 0.5 0.4 0.2 0.1 0.1 0.2]+0.15;
elseif (TheChar==37),                                 % %
 % contributed by Fred M. Staudaher (freds@packet.net)
 draw_circ(x+.2*scale,y+.8*scale,0.1*scale);
 draw_circ(x+.8*scale,y+.2*scale,0.1*scale);
 xpoints=[0.2 0.8]; 
 ypoints=[0.1 0.9];
elseif (TheChar==38),                                 % &
 xpoints=[0.7 0.5 0.3 0.2 0.2 0.6 0.4 0.2 0.7];
 ypoints=[0.5 0.1 0.1 0.2 0.4 0.8 1.0 0.8 0.1];
elseif (TheChar==39),                                 % '
 xpoints=[0.5 0.5];
 ypoints=[0.8 0.9];
elseif (TheChar==40),                                 % (
 xpoints=[0.6 0.45 0.35 0.45 0.6];
 ypoints=[0.1 0.20 0.55 0.90 1.0];
elseif (TheChar==41),                                 % )
 xpoints=[0.05 0.2 0.30 0.2 0.05];
 ypoints=[0.1  0.2 0.55 0.9 1.00];
elseif (TheChar==42),                                 % *
 xpoints=[0.2 0.8];
 ypoints=[0.5 0.5];
 xcoord=xpoints*cos(teta)-ypoints*sin(teta);
 ycoord=xpoints*sin(teta)+ypoints*cos(teta);
 plot(x+xcoord*scale, y+ycoord*scale,style,'LineWidth',thick);
 xpoints=[0.3 0.7];
 ypoints=[0.3 0.7];
 xcoord=xpoints*cos(teta)-ypoints*sin(teta);
 ycoord=xpoints*sin(teta)+ypoints*cos(teta);
 plot(x+xcoord*scale, y+ycoord*scale,style,'LineWidth',thick);
 xpoints=[0.3 0.7];
 ypoints=[0.7 0.3];
 xcoord=xpoints*cos(teta)-ypoints*sin(teta);
 ycoord=xpoints*sin(teta)+ypoints*cos(teta);
 plot(x+xcoord*scale, y+ycoord*scale,style,'LineWidth',thick);
 xpoints=[0.5 0.5];
 ypoints=[0.2 0.8];
elseif (TheChar==43),                                 % +
 xpoints=[0.2 0.8];
 ypoints=[0.5 0.5];
 xcoord=xpoints*cos(teta)-ypoints*sin(teta);
 ycoord=xpoints*sin(teta)+ypoints*cos(teta);
 plot(x+xcoord*scale, y+ycoord*scale,style,'LineWidth',thick);
 xpoints=[0.5 0.5];
 ypoints=[0.2 0.8];
elseif (TheChar==44),                                 % ,
 xpoints=[0.50 0.40 0.40 0.50 0.50 0.50 0.40];
 ypoints=[0.05 0.05 0.15 0.15 0.05 0   -0.1];
elseif (TheChar==45),                                 % -
 xpoints=[0.2 0.8];
 ypoints=[0.5 0.5];
elseif (TheChar==46),                                 % .
 xpoints=[0.40 0.50 0.50 0.40 0.40];
 ypoints=[0.05 0.05 0.15 0.15 0.05];
elseif (TheChar==47),                                 % /
 xpoints=[0.2 0.8];
 ypoints=[0.2 0.8];
elseif (TheChar==48),                                 % 0
 xpoints=[0.1 0.9 33/170 0.1 0.3 0.7 0.9 0.9 0.7 0.3 0.1 0.1];
 ypoints=[0.1 1.0 35/170 0.3 0.1 0.1 0.3 0.8 1.0 1.0 0.8 0.3];
elseif (TheChar==49),                                 % 1
 xpoints=[0.3 0.5 0.5 0.3 0.7];
 ypoints=[0.8 1.0 0.1 0.1 0.1];
elseif (TheChar==50),                                 % 2
 xpoints=[0.1 0.2 0.9 0.9  0.2  0.1  0.1 0.8 0.9]; 
 ypoints=[0.9 1.0 1.0 0.55 0.55 0.45 0.1 0.1 0.2];
elseif (TheChar==51),                                 % 3
 xpoints=[0.1 0.2 0.8 0.9 0.9  0.8  0.3  0.8  0.9  0.9 0.8 0.2 0.1]; 
 ypoints=[0.9 1.0 1.0 0.9 0.65 0.55 0.55 0.55 0.45 0.2 0.1 0.1 0.2];
elseif (TheChar==52),                                 % 4
 xpoints=[0.8  0.1  0.65 0.65]; 
 ypoints=[0.45 0.45 1.0  0.1];
elseif (TheChar==53),                                 % 5
 xpoints=[0.9 0.1 0.1  0.8  0.9  0.9 0.8 0.2 0.1]; 
 ypoints=[1.0 1.0 0.55 0.55 0.45 0.2 0.1 0.1 0.2];
elseif (TheChar==54),                                 % 6
 xpoints=[0.9 0.8 0.2 0.1 0.1  0.8  0.9  0.9 0.8 0.2 0.1 0.1]; 
 ypoints=[0.9 1.0 1.0 0.9 0.55 0.55 0.45 0.2 0.1 0.1 0.2 0.55];
elseif (TheChar==55),                                 % 7
 xpoints=[0.1 0.9 0.9 0.5  0.5]; 
 ypoints=[1.0 1.0 0.8 0.55 0.1];
elseif (TheChar==56),                                 % 8
 xpoints=[0.2  0.1  0.1 0.2 0.8 0.9 0.9  0.8  0.2  0.8  0.9  0.9 0.8 0.2 0.1 0.1 0.2]; 
 ypoints=[0.55 0.65 0.9 1.0 1.0 0.9 0.65 0.55 0.55 0.55 0.45 0.2 0.1 0.1 0.2 0.45 0.55];
elseif (TheChar==57),                                 % 9
 xpoints=[0.2  0.1  0.1 0.2 0.8 0.9 0.9  0.2  0.9  0.9 0.8 0.2 0.1]; 
 ypoints=[0.55 0.65 0.9 1.0 1.0 0.9 0.55 0.55 0.55 0.2 0.1 0.1 0.2];
elseif (TheChar==58),                                 % :
 xpoints=[0.40 0.50 0.50 0.40 0.40];
 ypoints=[0.65 0.65 0.75 0.75 0.65];
 xcoord=xpoints*cos(teta)-ypoints*sin(teta);
 ycoord=xpoints*sin(teta)+ypoints*cos(teta);
 plot(x+xcoord*scale, y+ycoord*scale,style,'LineWidth',thick);
 xpoints=[0.40 0.50 0.50 0.40 0.40];
 ypoints=[0.05 0.05 0.15 0.15 0.05];
elseif (TheChar==59),                                 % ;
 xpoints=[0.40 0.50 0.50 0.40 0.40];
 ypoints=[0.65 0.65 0.75 0.75 0.65];
 xcoord=xpoints*cos(teta)-ypoints*sin(teta);
 ycoord=xpoints*sin(teta)+ypoints*cos(teta);
 plot(x+xcoord*scale, y+ycoord*scale,style,'LineWidth',thick);
 xpoints=[0.50 0.40 0.40 0.50 0.50 0.50 0.40];
 ypoints=[0.05 0.05 0.15 0.15 0.05 0   -0.1];
elseif (TheChar==60),                                 % <
 xpoints=[0.9 0.1 0.9];
 ypoints=[0.1 0.5 0.9];
elseif (TheChar==61),                                 % =
 xpoints=[0.2 0.8];
 ypoints=[0.6 0.6];
 xcoord=xpoints*cos(teta)-ypoints*sin(teta);
 ycoord=xpoints*sin(teta)+ypoints*cos(teta);
 plot(x+xcoord*scale, y+ycoord*scale,style,'LineWidth',thick);
 xpoints=[0.2 0.8]; 
 ypoints=[0.4 0.4];
elseif (TheChar==62),                                 % >
 xpoints=[0.1 0.9 0.1];
 ypoints=[0.1 0.5 0.9];
elseif (TheChar==63),                                 % ?
 xpoints=[0.45 0.55 0.55 0.45 0.45];
 ypoints=[0.05 0.05 0.15 0.15 0.05];
 xcoord=xpoints*cos(teta)-ypoints*sin(teta);
 ycoord=xpoints*sin(teta)+ypoints*cos(teta);
 plot(x+xcoord*scale, y+ycoord*scale,style,'LineWidth',thick);
 xpoints=[0.3 0.2 0.3 0.7 0.8 0.5 0.50];
 ypoints=[0.6 0.8 1.0 1.0 0.8 0.4 0.25];
elseif (TheChar==64),                                 % @
 % contributed by Fred M. Staudaher (freds@packet.net)
 xpoints=[0.9 0.7 0.5 0.35 0.35 0.5 0.7 0.9 0.9 0.7 0.3 0.1 0.1 0.35 0.9]; 
 ypoints=[0.6 0.7 0.7 0.60 0.40 0.3 0.3 0.4 0.8 0.9 0.9 0.7 0.3 0.10 0.1];
elseif (TheChar==65),                                 % A
 xpoints=[0.1 0.1 0.3 0.7 0.9 0.9 0.9 0.1];
 ypoints=[0.1 0.8 1.0 1.0 0.8 0.1 0.55 0.55];
elseif (TheChar==66),                                 % B
 xpoints=[0.1 0.8 0.9 0.9 0.8 0.2 0.8 0.9 0.9 0.8 0.1 0.2 0.2];
 ypoints=[0.1 0.1 0.1 0.4 0.55 0.55 0.55 0.6 0.9 1.0 1   1   0.1];
elseif (TheChar==67),                                 % C
 xpoints=[0.9 0.8 0.2 0.1 0.1 0.2 0.8 0.9]; 
 ypoints=[0.2 0.1 0.1 0.2 0.9 1.0 1.0 0.9];
elseif (TheChar==68),                                 % D
 xpoints=[0.1 0.8 0.9 0.9 0.8 0.1 0.2 0.2]; 
 ypoints=[0.1 0.1 0.2 0.9 1.0 1.0 1.0 0.1];
elseif (TheChar==69),                                 % E
 xpoints=[0.1 0.8 0.1 0.1 0.8 0.1 0.1 0.8]; 
 ypoints=[0.1 0.1 0.1 0.55 0.55 0.55 1   1];
elseif (TheChar==70),                                 % F
 xpoints=[0.1 0.1 0.8 0.1 0.1 0.8]; 
 ypoints=[0.1 0.55 0.55 0.55 1   1];
elseif (TheChar==71),                                 % G
 xpoints=[0.5 0.9 0.9 0.8 0.2 0.1 0.1 0.1 0.8 0.9]; 
 ypoints=[0.55 0.55 0.1 0.1 0.1 0.2 0.9 1.0 1.0 0.9];
elseif (TheChar==72),                                 % H
 xpoints=[0.15 0.15 0.15 0.85 0.85 0.85]; 
 ypoints=[0.10 1.00 0.55 0.55 1.00 0.10];
elseif (TheChar==73),                                 % I
 xpoints=[0.4 0.6 0.5 0.5 0.6 0.4]; 
 ypoints=[0.1 0.1 0.1 1.0 1.0 1.0];
elseif (TheChar==74),                                 % J
 xpoints=[0.1 0.1 0.3 0.7 0.8 0.8]; 
 ypoints=[0.3 0.2 0.1 0.1 0.2 1.0];
elseif (TheChar==75),                                 % K
 xpoints=[0.9 0.1 0.9 0.1 0.1 0.1]; 
 ypoints=[1.0 0.55 0.1 0.55 1.0 0.1];
elseif (TheChar==76),                                 % L
 xpoints=[0.1 0.1 0.8 0.9]; 
 ypoints=[1.0 0.1 0.1 0.2];
elseif (TheChar==77),                                 % M
 xpoints=[0.1 0.1 0.5 0.9 0.9]; 
 ypoints=[0.1 1.0 0.5 1.0 0.1];
elseif (TheChar==78),                                 % N
 xpoints=[0.1 0.1 0.1 0.9 0.9 0.9]; 
 ypoints=[0.1 1.0 0.8 0.3 0.1 1.0];
elseif (TheChar==79),                                 % O
 xpoints=[0.1 0.3 0.7 0.9 0.9 0.7 0.3 0.1 0.1]; 
 ypoints=[0.3 0.1 0.1 0.3 0.8 1.0 1.0 0.8 0.3];
elseif (TheChar==80),                                 % P
 xpoints=[0.1 0.1 0.8 0.9 0.9 0.8 0.1]; 
 ypoints=[0.1 1.0 1.0 0.9 0.6 0.5 0.5];
elseif (TheChar==81),                                 % Q
 xpoints=[0.1 0.3 0.7 0.8 0.6 0.7 0.9 0.9 0.7 0.3 0.1 0.1]; 
 ypoints=[0.3 0.1 0.1 0.0 0.2 0.1 0.3 0.8 1.0 1.0 0.8 0.3];
elseif (TheChar==82),                                 % R
 xpoints=[0.9 0.5 0.2 0.8 0.9 0.9 0.8 0.1 0.2 0.2];
 ypoints=[0.1 0.5 0.5 0.5 0.6 0.9 1.0 1   1   0.1];
elseif (TheChar==83),                                 % S
 xpoints=[0.9 0.7 0.3 0.1 0.1 0.3 0.7 0.9 0.9 0.7 0.3 0.1]; 
 ypoints=[0.8 1.0 1.0 0.8 0.7 0.5 0.5 0.4 0.2 0.1 0.1 0.2];
elseif (TheChar==84),                                 % T
 xpoints=[0.5 0.5 0.1 0.9]; 
 ypoints=[0.1 1.0 1.0 1.0];
elseif (TheChar==85),                                 % U
 xpoints=[0.1 0.1 0.3 0.7 0.9 0.9]; 
 ypoints=[1.0 0.3 0.1 0.1 0.3 1.0];
elseif (TheChar==86),                                 % V
 xpoints=[0.1 0.5 0.9]; 
 ypoints=[1.0 0.1 1.0];
elseif (TheChar==87),                                 % W
 xpoints=[0.1 0.3 0.5 0.7 0.9];
 ypoints=[1.0 0.1 0.45 0.1 1.0];
elseif (TheChar==88),                                 % X
 xpoints=[0.1 0.9 0.5  0.1 0.9]; 
 ypoints=[0.1 1.0 0.55 1.0 0.1];
elseif (TheChar==89),                                 % Y
 xpoints=[0.1 0.5 0.9 0.5 0.5]; 
 ypoints=[1.0 0.5 1.0 0.5 0.1];
elseif (TheChar==90),                                 % Z
 xpoints=[0.1 0.3 0.9 0.1 0.7 0.9]; 
 ypoints=[0.8 1.0 1.0 0.1 0.1 0.3];
elseif (TheChar==91),                                 % [
 xpoints=[0.6 0.35 0.35 0.6];
 ypoints=[0.1 0.10 1.00 1.0];
elseif (TheChar==92),                                 % \
 xpoints=[0.2 0.8];
 ypoints=1-[0.2 0.8];
elseif (TheChar==93),                                 % ]
 xpoints=[0.4 0.65 0.65 0.4];
 ypoints=[0.1 0.10 1.00 1.0];
elseif (TheChar==94),                                 % ^
 xpoints=[0.2 0.55 0.9];
 ypoints=[0.7 1.00 0.7];
elseif (TheChar==95),                                 % _
 xpoints=[0.1 1.0];
 ypoints=[0.1 0.1];
elseif (TheChar==97),                                 % a
 xpoints=[1.0 0.9 0.9 0.7 0.3 0.1 0.1 0.3 0.7 0.9]; 
 ypoints=[0.1 0.3 0.5 0.7 0.7 0.5 0.3 0.1 0.1 0.3];
elseif (TheChar==98),                                 % b
 xpoints=[0.9 0.9 0.7 0.3 0.1 0.1 0.1 0.3 0.7 0.9]; 
 ypoints=[0.3 0.5 0.7 0.7 0.5 1.0 0.3 0.1 0.1 0.3];
elseif (TheChar==99),                                 % c
 xpoints=[0.9 0.8 0.3 0.1 0.1 0.3 0.8 0.9]; 
 ypoints=[0.6 0.7 0.7 0.5 0.3 0.1 0.1 0.2];
elseif (TheChar==100),                                % d
 xpoints=[0.9 0.9 0.9 0.9 0.7 0.3 0.1 0.1 0.3 0.7 0.9]; 
 ypoints=[0.1 1.0 0.3 0.5 0.7 0.7 0.5 0.3 0.1 0.1 0.3];
elseif (TheChar==101),                                % e
 xpoints=[0.1 0.8 0.9 0.9 0.8 0.3 0.1 0.1 0.3 0.8 0.9]; 
 ypoints=[0.4 0.4 0.5 0.6 0.7 0.7 0.5 0.3 0.1 0.1 0.2];
elseif (TheChar==102),                                % f
 xpoints=[0.2 0.20 0.65 0.20 0.2 0.4 0.8]; 
 ypoints=[0.1 0.55 0.55 0.55 0.9 1   1.0];
elseif (TheChar==103),                                % g
 xpoints=[0.1 0.3 0.7 0.9 0.9 0.9 0.7 0.3 0.1 0.1 0.3 0.7 0.9]; 
 ypoints=[0 -0.1 -0.1 0 0.4 0.6 0.8 0.8 0.6 0.4 0.2 0.2 0.4]-0.1;
elseif (TheChar==104),                                % h
 xpoints=[0.2 0.20 0.20 0.65 0.8 0.8]; 
 ypoints=[0.1 0.9  0.55 0.55 0.45 0.1];
elseif (TheChar==105),                                % i
 xpoints=[0.45 0.55 0.55 0.45 0.45];
 ypoints=[0.87 0.87 0.97 0.97 0.87];
 xcoord=xpoints*cos(teta)-ypoints*sin(teta);
 ycoord=xpoints*sin(teta)+ypoints*cos(teta);
 plot(x+xcoord*scale, y+ycoord*scale,style,'LineWidth',thick);
 xpoints=[0.3 0.7 0.5 0.5 0.7 0.3]; 
 ypoints=[0.1 0.1 0.1 0.7 0.7 0.7];
elseif (TheChar==106),                                % j
 xpoints=[0.45 0.55 0.55 0.45 0.45];
 ypoints=[0.87 0.87 0.97 0.97 0.87];
 xcoord=xpoints*cos(teta)-ypoints*sin(teta);
 ycoord=xpoints*sin(teta)+ypoints*cos(teta);
 plot(x+xcoord*scale, y+ycoord*scale,style,'LineWidth',thick);
 xpoints=[0.2  0.3  0.5 0.5 0.6 0.4];
 ypoints=[0.1 -0.1 -0.1 0.7 0.7 0.7];
elseif (TheChar==107),                                % k
 xpoints=[0.2 0.20 0.2 0.8 0.2 0.8];
 ypoints=[0.1 0.9  0.4 0.7 0.4 0.1];
elseif (TheChar==108),                                % l
 xpoints=[0.3 0.8 0.8 0.7 0.5 0.4 0.4 0.6 0.7 0.8];
 ypoints=[0.3 0.5 0.9 1.0 1.0 0.9 0.2 0.1 0.1 0.2];
elseif (TheChar==109),                                % m
 xpoints=[0.2 0.2 0.5 0.5 0.5 0.7 0.8 0.8]; 
 ypoints=[0.1 0.7 0.7 0.1 0.7 0.7 0.6 0.1];
elseif (TheChar==110),                                % n
 xpoints=[0.25 0.35];
 ypoints=[0.10 0.10];
 xcoord=xpoints*cos(teta)-ypoints*sin(teta);
 ycoord=xpoints*sin(teta)+ypoints*cos(teta);
 plot(x+xcoord*scale, y+ycoord*scale,style,'LineWidth',thick);
 xpoints=[0.75 0.85];
 ypoints=[0.10 0.10];
 xcoord=xpoints*cos(teta)-ypoints*sin(teta);
 ycoord=xpoints*sin(teta)+ypoints*cos(teta);
 plot(x+xcoord*scale, y+ycoord*scale,style,'LineWidth',thick);
 xpoints=[0.3 0.3 0.2 0.7 0.8 0.8]; 
 ypoints=[0.1 0.7 0.7 0.7 0.6 0.1];
elseif (TheChar==111),                                % o
 xpoints=[0.9 0.9 0.7 0.3 0.1 0.1 0.3 0.7 0.9]; 
 ypoints=[0.3 0.5 0.7 0.7 0.5 0.3 0.1 0.1 0.3];
elseif (TheChar==112),                                % p
 xpoints=[0.9 0.9 0.7 0.3 0.1  0.1 0.1 0.3 0.7 0.9]; 
 ypoints=[0.3 0.5 0.7 0.7 0.5 -0.2 0.3 0.1 0.1 0.3];
elseif (TheChar==113),                                % q
 xpoints=[0.9 0.9 0.7 0.3 0.1 0.1 0.3 0.7 0.9]; 
 ypoints=[-0.1 0.5 0.7 0.7 0.5 0.3 0.1 0.1 0.3];
elseif (TheChar==114),                                % r
 xpoints=[0.3 0.3 0.2 0.3 0.3 0.5 0.7 0.8 0.8]; 
 ypoints=[0.1 0.7 0.7 0.7 0.5 0.7 0.7 0.6 0.5];
elseif (TheChar==115),                                % s
 xpoints=[0.9 0.8 0.3 0.1 0.1 0.3 0.8 0.9 0.9 0.8 0.3 0.1];
 ypoints=[0.6 0.7 0.7 0.6 0.5 0.4 0.4 0.3 0.2 0.1 0.1 0.2];
elseif (TheChar==116),                                % t
 xpoints=[0.3 0.7 0.4 0.4 0.4 0.6 0.7 0.8];
 ypoints=[0.7 0.7 0.7 1.0 0.2 0.1 0.1 0.2];
elseif (TheChar==117),                                % u
 xpoints=[0.9 0.9 0.9 0.7 0.3 0.1 0.1]; 
 ypoints=[0.7 0.1 0.3 0.1 0.1 0.3 0.7];
elseif (TheChar==118),                                % v
 xpoints=[0.9 0.9 0.5 0.1 0.1]; 
 ypoints=[0.7 0.6 0.1 0.6 0.7];
elseif (TheChar==119),                                % w
 xpoints=[0.9 0.9 0.7 0.5 0.5 0.5 0.3 0.1 0.1]; 
 ypoints=[0.7 0.4 0.1 0.4 0.7 0.4 0.1 0.4 0.7];
elseif (TheChar==120),                                % x
 xpoints=[0.1 0.9 0.5 0.1 0.9]; 
 ypoints=[0.7 0.1 0.4 0.1 0.7];
elseif (TheChar==121),                                % y
 xpoints=[ 0.1  0.3  0.7  0.8 0.8 0.8 0.7 0.3 0.1 0.1]; 
 ypoints=[-0.1 -0.2 -0.2 -0.1 0.7 0.2 0.1 0.1 0.2 0.7];
elseif (TheChar==122),                                % z
 xpoints=[0.1 0.2 0.9 0.1 0.8 0.9]; 
 ypoints=[0.6 0.7 0.7 0.1 0.1 0.2];
elseif (TheChar==123),                                % {
 xpoints=[0.6 0.45 0.40 0.35 0.30 0.35 0.40 0.45 0.6];
 ypoints=[0.1 0.30 0.55 0.65 0.55 0.45 0.55 0.80 1.0];
elseif (TheChar==124),                                % |
 % contributed by Fred M. Staudaher (freds@packet.net)
 xpoints=[0.5 0.5]; 
 ypoints=[0.1 0.45];
 xcoord=xpoints*cos(teta)-ypoints*sin(teta);
 ycoord=xpoints*sin(teta)+ypoints*cos(teta);
 plot(x+xcoord*scale, y+ycoord*scale,style,'LineWidth',thick);
 xpoints=[0.5 0.5]; 
 ypoints=[0.55 0.9];
elseif (TheChar==125),                                % }
 xpoints=1-[0.6 0.45 0.40 0.35 0.30 0.35 0.40 0.45 0.6];
 ypoints=[0.1 0.30 0.55 0.65 0.55 0.45 0.55 0.80 1.0];
elseif (TheChar==126),                                % ~
 % contributed by Fred M. Staudaher (freds@packet.net)
 xpoints=[0.1 0.20 0.3 0.40 0.5 0.60 0.7 0.80 0.9]; 
 ypoints=[0.9 0.98 1.0 0.98 0.9 0.82 0.8 0.82 0.9]-0.1;
elseif (TheChar==127),                                % space
 return 
else
 xpoints=0; ypoints=0;
 fprintf('unknown char %i\n', TheChar);
end;

xcoord=xpoints*cos(teta)-ypoints*sin(teta);
ycoord=xpoints*sin(teta)+ypoints*cos(teta);
plot(x+xcoord*scale, y+ycoord*scale,style,'LineWidth',thick);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xout,yout]=draw_text(x,y,TheText,scale,angle,style,thick)
% draw_text(x,y,TheText,scale,angle,style,thick) draws a text string 
% at the (x,y) point, with scale as desired scaling factor 
% and angle (in degrees) as orientation. 
% If TheText is a vector, the second row is the subscript, 
% and the third is the exponent.
% style is the third plot parameter (see the help of plot).

% Updated to include line thickness by R. Clewley, Dec 2002

if nargin==0
 help draw_text;
 disp(' **  draw_text() routine by auger@ge44.univ-nantes.fr')
 disp('       and modified by rclewley@bu.edu')
 return
elseif (nargin<3), error('At least 3 parameters required.');
elseif (nargin==3), scale=1; angle=0; style='b-'; thick=1;
elseif (nargin==4), angle=0; style='b-'; thick=1;
elseif (nargin==5), style='b-'; thick=1;
end;

[TextRow,TextCol]=size(TheText);
if (TextRow>1),
 error('TheText must have one row');
end;
TheText=abs(TheText);

teta=pi*angle/180.0; savedscale=scale;
xnorm=1; ynorm=0;
xcoord=xnorm*cos(teta)-ynorm*sin(teta);
ycoord=xnorm*sin(teta)+ynorm*cos(teta);

i=1;
while (i<=TextCol),
 if (TheText(i)==92)&(i<TextCol),                     % LaTeX symbols begining with \
  if (i+2<=TextCol)
   if (TheText(i+1)==34),                             % \"
    if (TheText(i+2)>=97)&(TheText(i+2)<=122),        % lowercase letters
     xnorm=[0.35 0.45 0.45 0.35 0.35];
     ynorm=[0.85 0.85 0.95 0.95 0.85];
    elseif (TheText(i+2)>=65)&(TheText(i+2)<=90),     % uppercase letters
     xnorm=[0.30 0.40 0.40 0.30 0.30];
     ynorm=[1.10 1.10 1.20 1.20 1.10];
    end
    xcoord2=xnorm*cos(teta)-ynorm*sin(teta);
    ycoord2=xnorm*sin(teta)+ynorm*cos(teta);
    plot(x+xcoord2*scale, y+ycoord2*scale,'LineWidth',thick);

    if (TheText(i+2)>=97)&(TheText(i+2)<=122),        % lowercase letters
     xnorm=[0.55 0.65 0.65 0.55 0.55];
     ynorm=[0.85 0.85 0.95 0.95 0.85];
    elseif (TheText(i+2)>=65)&(TheText(i+2)<=90),     % uppercase letters
     xnorm=[0.50 0.60 0.60 0.50 0.50];
     ynorm=[1.10 1.10 1.20 1.20 1.10];
    end
    xcoord2=xnorm*cos(teta)-ynorm*sin(teta);
    ycoord2=xnorm*sin(teta)+ynorm*cos(teta);
    plot(x+xcoord2*scale, y+ycoord2*scale,'LineWidth',thick);
    i=i+2;
    if (TheText(i)==105),                             % special case for the i
     xnorm=[0.3 0.7 0.5 0.5 0.7 0.3]; 
     ynorm=[0.1 0.1 0.1 0.7 0.7 0.7];
     xcoord2=xnorm*cos(teta)-ynorm*sin(teta);
     ycoord2=xnorm*sin(teta)+ynorm*cos(teta);
     plot(x+xcoord2*scale, y+ycoord2*scale,style,'LineWidth',thick);
     i=i+1; x=x+scale*xcoord; y=y+scale*ycoord;
    end
   elseif (TheText(i+1)==39)                          % \'
    if (TheText(i+2)>=97)&(TheText(i+2)<=122),        % lowercase letters
     xnorm=[0.45 0.65];
     ynorm=[0.80 0.90];
    elseif (TheText(i+2)>=65)&(TheText(i+2)<=90),     % uppercase letters
     xnorm=[0.40 0.60];
     ynorm=[1.10 1.20];
    end
    xcoord2=xnorm*cos(teta)-ynorm*sin(teta);
    ycoord2=xnorm*sin(teta)+ynorm*cos(teta);
    plot(x+xcoord2*scale, y+ycoord2*scale,'LineWidth',thick);
    i=i+2;
    if (TheText(i)==105),                             % special case for the i
     xnorm=[0.3 0.7 0.5 0.5 0.7 0.3]; 
     ynorm=[0.1 0.1 0.1 0.7 0.7 0.7];
     xcoord2=xnorm*cos(teta)-ynorm*sin(teta);
     ycoord2=xnorm*sin(teta)+ynorm*cos(teta);
     plot(x+xcoord2*scale, y+ycoord2*scale,style,'LineWidth',thick);
     i=i+1; x=x+scale*xcoord; y=y+scale*ycoord;
    end
   elseif (TheText(i+1)==94)                          % \^
    if (TheText(i+2)>=97)&(TheText(i+2)<=122),        % lowercase letters
     xnorm=[0.40 0.50 0.60];
     ynorm=[0.80 0.90 0.80];
    elseif (TheText(i+2)>=65)&(TheText(i+2)<=90),     % uppercase letters
     xnorm=[0.40 0.50 0.60];
     ynorm=[1.10 1.20 1.10];
    end
    xcoord2=xnorm*cos(teta)-ynorm*sin(teta);
    ycoord2=xnorm*sin(teta)+ynorm*cos(teta);
    plot(x+xcoord2*scale, y+ycoord2*scale,'LineWidth',thick);
    i=i+2;
    if (TheText(i)==105),                             % special case for the i
     xnorm=[0.3 0.7 0.5 0.5 0.7 0.3]; 
     ynorm=[0.1 0.1 0.1 0.7 0.7 0.7];
     xcoord2=xnorm*cos(teta)-ynorm*sin(teta);
     ycoord2=xnorm*sin(teta)+ynorm*cos(teta);
     plot(x+xcoord2*scale, y+ycoord2*scale,style,'LineWidth',thick);
     i=i+1; x=x+scale*xcoord; y=y+scale*ycoord; 
    end;
   elseif (TheText(i+1)==96)                          % \`
    if (TheText(i+2)>=97)&(TheText(i+2)<=122),        % lowercase letters
     xnorm=[0.45 0.65];
     ynorm=[0.90 0.80];
    elseif (TheText(i+2)>=65)&(TheText(i+2)<=90),     % uppercase letters
     xnorm=[0.40 0.60];
     ynorm=[1.20 1.10];
    end
    xcoord2=xnorm*cos(teta)-ynorm*sin(teta);
    ycoord2=xnorm*sin(teta)+ynorm*cos(teta);
    plot(x+xcoord2*scale, y+ycoord2*scale,'LineWidth',thick);
    i=i+2;
    if (TheText(i)==105),                             % special case for the i
     xnorm=[0.3 0.7 0.5 0.5 0.7 0.3]; 
     ynorm=[0.1 0.1 0.1 0.7 0.7 0.7];
     xcoord2=xnorm*cos(teta)-ynorm*sin(teta);
     ycoord2=xnorm*sin(teta)+ynorm*cos(teta);
     plot(x+xcoord2*scale, y+ycoord2*scale,style,'LineWidth',thick);
     i=i+1; x=x+scale*xcoord; y=y+scale*ycoord;
    end
   end;
  end;
  if (i+3<=TextCol),
   if (TheText(i+1)==103)&(TheText(i+2)==101)&...
      (TheText(i+3)==113),                            % \geq
    xnorm=[0.1 0.9 0.1];
    ynorm=[0.3 0.6 0.9];
    xcoord2=xnorm*cos(teta)-ynorm*sin(teta);
    ycoord2=xnorm*sin(teta)+ynorm*cos(teta);
    plot(x+xcoord2*scale, y+ycoord2*scale,'LineWidth',thick);
    xnorm=[0.15 0.9];
    ynorm=[0.1+0.5*(0.15-0.1) 0.4];
    xcoord2=xnorm*cos(teta)-ynorm*sin(teta);
    ycoord2=xnorm*sin(teta)+ynorm*cos(teta);
    plot(x+xcoord2*scale, y+ycoord2*scale,'LineWidth',thick);
    x=x+scale*xcoord; y=y+scale*ycoord; i=i+4;
   elseif (TheText(i+1)==108)&(TheText(i+2)==101)&...
      (TheText(i+3)==113),                            % \leq
    xnorm=1-[0.1 0.9 0.1];
    ynorm=[0.3 0.6 0.9];
    xcoord2=xnorm*cos(teta)-ynorm*sin(teta);
    ycoord2=xnorm*sin(teta)+ynorm*cos(teta);
    plot(x+xcoord2*scale, y+ycoord2*scale,'LineWidth',thick);
    xnorm=1-[0.15 0.9];
    ynorm=[0.1+0.5*(0.15-0.1) 0.4];
    xcoord2=xnorm*cos(teta)-ynorm*sin(teta);
    ycoord2=xnorm*sin(teta)+ynorm*cos(teta);
    plot(x+xcoord2*scale, y+ycoord2*scale,'LineWidth',thick);
    x=x+scale*xcoord; y=y+scale*ycoord; i=i+4;
   end
  end;
  if (i+4<=TextCol),
   if (TheText(i+1)==98)&(TheText(i+2)==101)&...
      (TheText(i+3)==116)&(TheText(i+4)==97),         % \beta
    i=i+5;
   end;
  end;
  if (i+5<=TextCol),
   if (TheText(i+1)==97)&(TheText(i+2)==108)&(TheText(i+3)==112)&...
      (TheText(i+4)==104)&(TheText(i+5)==97), % \alpha
   end;
  end;
 elseif (TheText(i)==94)&(i<TextCol),
  if (TheText(i+1)==123)&(i+2<=TextCol),              % exponent
   i=i+2; dx=0; dy=0.5;
   Deltax=dx*cos(teta)-dy*sin(teta);
   Deltay=dx*sin(teta)+dy*cos(teta);
   x=savedx+scale*Deltax; y=savedy+scale*Deltay; scale=scale*0.6;
   while (TheText(i)~=125),
    draw_char(x,y,TheText(i),scale,angle,style,thick);
    x=x+scale*xcoord; y=y+scale*ycoord; i=i+1;
   end
   dx=0; dy=-0.5; i=i+1;
   Deltax=dx*cos(teta)-dy*sin(teta);
   Deltay=dx*sin(teta)+dy*cos(teta);
   scale=savedscale; x=x+scale*Deltax; y=y+scale*Deltay;
  else
   draw_char(x,y,TheText(i),scale,angle,style,thick);
   x=x+scale*xcoord; y=y+scale*ycoord; i=i+1;
   savedx=x; savedy=y;
  end
 elseif (TheText(i)==95)&(i<TextCol),
  if (TheText(i+1)==123)&(i+2<=TextCol),              % subscript
   i=i+2; dx=0; dy=-0.2;
   Deltax=dx*cos(teta)-dy*sin(teta);
   Deltay=dx*sin(teta)+dy*cos(teta);
   x=savedx+scale*Deltax; y=savedy+scale*Deltay; scale=scale*0.6;
   cont=(TheText(i)~=125);
   while (i<TextCol)&cont,
    draw_char(x,y,TheText(i),scale,angle,style,thick);
    x=x+scale*xcoord; y=y+scale*ycoord; i=i+1;
    cont=(TheText(i)~=125);
   end
   dx=0; dy=+0.2; i=i+1;
   Deltax=dx*cos(teta)-dy*sin(teta);
   Deltay=dx*sin(teta)+dy*cos(teta);
   scale=savedscale; x=x+scale*Deltax; y=y+scale*Deltay;
  else
   draw_char(x,y,TheText(i),scale,angle,style,thick);
   x=x+scale*xcoord; y=y+scale*ycoord; i=i+1;
   savedx=x; savedy=y;
  end
 else
  if (TheText(i)~=32), draw_char(x,y,TheText(i),scale,angle,style,thick); end;
  x=x+scale*xcoord; y=y+scale*ycoord; i=i+1;
  savedx=x; savedy=y;
 end
end
xout=x; yout=y;
return
