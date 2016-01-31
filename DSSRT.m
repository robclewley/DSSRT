function DSSRT(networkObjDirPar,parFile)
% DSSRT.m (loader for DSSRT version 1.xx)
% (c) Robert Clewley, Center for Biodynamics, Boston University, 2003-2004
%
% Usage: DSSRT [ <ODE network object folder> ]

% Last structural update: November 2005, for Matlab R14

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constants, do not change!
UNIX    = 0;
WINDOWS = 1;
MAC     = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==1
    parFile = '';
end
if nargin==0
    parFile = '';
    networkObjDirPar = '';
end
if nargin > 2
    disp('DSSRT:  Too many arguments specified. Maximum of 2')
    return
end

% system details for this version
DSSRTversionStr = '1.32';
DSSRTverDateStr = 'Feb 06';
systemFilenames = {'StepNet','FuncNet','ReadConfigFile','AttEst_epochs','AttEst_regimes','RegimeDet','GlobalConstants', ...
        'linkStateMapPots','emptyStateMap','linkStateMapActs','linkStateMapPots','questdlg','inputdlg','GetColStyleSet'};
globalsFile = 'GlobalConstants'; % script containing common global constant definitions
compileMe = {'ReadConfigFile','FuncNet'}; % the order is not important because all are preprocessed before compilation
% each set in compileFiles corresponds to the compilation files associated with the corresponding entry in compileMe, for that platform
compileFiles.win  = { {'readconfigfile.c', 'readconfigfile_mex.c', 'linkstatemapacts.h', 'linkstatemappots.h', ...
                         'emptystatemap.h', 'fgetl.h', 'num2str.h', 'readconfigfile.h', 'str2num.h'}, ...
                      {'funcnet.c', 'funcnet_mex.c', 'funcnet.h', 'num2str.h', 'readconfigfile.h'}};
compileFiles.unix = { {''}, ...
                      {''} }; % for future use
compileFiles.mac  = { {''}, ...
                      {''} }; % for future use
% the following string is searched for in the compilation object files to find position to insert the definitions
% it should be always be the name of the globals file followed by the identifying comment
compileGlobalsStr = [ globalsFile ' % this line will be replaced by the contents of this script file when compiled (see DSSRT.m)'];

%
% Preamble
%
clc % clear command window
disp(   '  DSSRT: Dominant-Scale System Reduction Tool for ODE analysis')
disp(   '  (c) Robert Clewley, Center for Biodynamics, Boston University')
disp(   '              and Department of Mathematics, Cornell University')
disp(   '  Copyright 2002 - 2006')
fprintf('  This is version %s (%s) running on %s\n\n',DSSRTversionStr,DSSRTverDateStr,datestr(now))
disp('DSSRT:  Initializing system')

%
% Check Matlab version
%
vStr = version;
vVal = str2num(vStr(1:3)); % assumes only a one-digit leading number, i.e. change this at Version 10 of Matlab :-)
if vVal < 6.5
    beep
    disp('DSSRT:  Warning! Old version of Matlab being used. Compatibility is not guaranteed!')
    disp('        Press any key to continue (or CTRL-C to quit) ...')
    pause
end

%
% Get user settings
%
if exist( 'DSSRT_useroptions.m', 'file' ) == 2
    DSSRT_useroptions % get user settings
else
    beep
    disp('DSSRT:  Missing `DSSRT_useroptions.m` initialization file')
    disp('        This file needs to be in the same directory as this script')
    return
end

%
% Check all system files are present (DSSRT_useroptions.m already checked)
%
for i=1:length(systemFilenames)
    sysfname = systemFilenames{i};
    if exist( [sysfname '.m'], 'file' ) ~= 2
        beep
        fprintf('DSSRT:  File %s.m is missing from the DSSRT root directory.\n',sysfname)
        disp(   '        Cannot continue')
        return
    end
end

%
% Check user settings
%
if exist(DSSRT_ROOTPATH,'file')==7
    cd(DSSRT_ROOTPATH)
    proceed = true;
else
    disp('DSSRT: The DSSRT_ROOTPATH specified in DSSRT_useroptions.m is not a valid directory')
    disp('       Cannot continue')
    return
end

if exist( 'COMPILE_ONLY' ) ~= 1 % then parameter not set up ok
    disp('DSSRT:  COMPILE_ONLY flag not set properly by DSSRT_useroptions')
    disp('        (defaulting to OFF)')
    COMPILEONLY = false;
end


if exist( 'START_VERBOSE' ) ~= 1 % then parameter not set up ok
    disp('DSSRT:  START_VERBOSE flag not set properly by DSSRT_useroptions')
    disp('        (defaulting to OFF)')
    START_VERBOSE = false;
end

if exist( 'KEEP_LOG' ) ~= 1 % then parameter not set up ok
    disp('DSSRT:  KEEP_LOG flag not set properly by DSSRT_useroptions')
    disp('        (defaulting to OFF)')
    KEEP_LOG = false;
end

%
% Setup platform-specific stuff if will compile
%
cstr = computer;
if strcmp(cstr,'PCWIN')
    PLATFORM = WINDOWS;
    compileFilesDel = compileFiles.win;
elseif strcmp(cstr,'MAC')
    PLATFORM = MAC;
    compileFilesDel = compileFiles.mac;
else % various possibilities
    PLATFORM = UNIX;
    compileFilesDel = compileFiles.unix;
end

if vVal > 7.
    if FORCE_COMPILE & ~FORCE_NOTCOMPILE
        disp('DSSRT:  Overriding FORCE_COMPILE user option because compilation to MEX')
        disp('         is no longer necessary in Matlab R14')
        disp('         (you may like to edit the user options file for this system)')
    end
    FORCE_NOTCOMPILE = true;
end

if FORCE_NOTCOMPILE
    compileFilesDel = [];
end

%
% Refresh / create compiled versions of .m files for local platform
%
if proceed
    done = false;
    while ~done
        if FORCE_NOTCOMPILE
            disp('DSSRT:  Development option: Forcing no compilation to MEX files')
            done = true;
        else
            updated = false;
            if FORCE_COMPILE
                disp('DSSRT:  Development option: Forcing compilation to MEX files')
                doCompile = ones(1,length(compileMe));
            else
                doCompile = zeros(1,length(compileMe));
                for i=1:length(compileMe)
                    if ~(exist([compileMe{i} '.' mexext],'file')==3)
                        doCompile(i) = 1;
                    end
                end
            end
            if sum(doCompile) > 0
                % preprocess all files before compiling, so all the .h headers are correct
                disp('DSSRT:  Preprocessing compilation object .m files for global constant definitions ...')
                for i=1:length(compileMe)
                    if ~doCompile(i)
                        continue % skip to next file
                    end
                    compFname = compileMe{i};
                    [status mess messId] = copyfile([compFname '.m'],[compFname '.temp'],'f');
                    if status
                        result = PreprocessForComp([compFname '.m'],[globalsFile '.m'],compileGlobalsStr);
                        if ~result % then preprocessing failed
                            disp(['DSSRT:  Preprocessing of ' compFname '.m failed'])
                            disp('DSSRT:  Cannot continue')
                            return
                        end
                    else
                        disp('DSSRT:  Copy failed in preprocessing. Matlab error follows ...')
                        fprintf(' %s\n',mess)
                        disp('        If there`s a big problem getting compilation to work on your system')
                        disp('        try switching `FORCE_NOTCOMPILE = true;` in DSSRTuseroptions.m')
                        return
                    end
                end
                for i=1:length(compileMe)
                    if ~doCompile(i)
                        continue % skip to next file
                    end
                    compFname = compileMe{i};
                    if doCompile(i)
                        fprintf('         Compiling %s.m to a MEX file ... ',compFname)
                        try
                            eval(['mcc -x -i ' compFname '.m'])
                            updated = true;
                            for delIx = 1:length(compileFilesDel{i})
                                if exist( compileFilesDel{i}{delIx} )
                                    delete(compileFilesDel{i}{delIx})
                                else
                                    disp(' ')
                                    beep
                                    disp('DSSRT:  Internal error: compilation temp file not found. Quitting ...')
                                    return
                                end
                            end
                            if isempty(compileFilesDel{1})
                                disp('DSSRT:  You will manually have to remove temporary compilation files')
                                disp('        from the system folder.')
                            end
                            done = true;
                            fprintf('done\n')
                        catch
                            fprintf('failed\n\n')
                            disp('DSSRT: Matlab reported the following error during compilation ...')
                            fprintf('%s\n\n',lasterr)
                            disp('DSSRT:  If a Matlab compiler is not present try `FORCE_NOTCOMPILE = true;` in DSSRTuseroptions.m')
                            disp('        (temporarily changing to this compile option to continue')
                            FORCE_NOTCOMPILE = true;
                            done = false;
                            break % out of this for loop
                        end
                    else
                        done = true;
                    end
                end
                if updated & done % then all files were successfully compiled, so tidy up .temp files
                    for i=1:length(compileMe)
                        if ~doCompile(i)
                            continue % skip to next file
                        end
                        compFname = compileMe{i};
                        delete([compFname '.m'])
                        [status mess messId] = movefile([compFname '.temp'],[compFname '.m'],'f');
                        if ~status
                            beep
                            disp(' ')
                            disp('DSSRT:  Processing failed at final step. Matlab error follows ...')
                            fprintf(' %s\n',mess)
                            return
                        end
                    end
                end
            else % nothing to do
                done = true;
            end % if sum(doCompile)>0
            if ~updated
                disp('         No MEX files updated')
            end
        end % if FORCE_NOTCOMPILE
    end % while ~done
else
    return
end % if proceed

if COMPILE_ONLY
    disp('DSSRT:  COMPILE_ONLY option selected. Returning ...')
    return
end

warning off all
% warning off MATLAB:rmpath:DirNotFound

disp('DSSRT:  System initialization complete')

systemChosen = false; % initial value only
if ~isempty(networkObjDirPar) && ischar(networkObjDirPar)
    networkObjectDir = networkObjDirPar;
else
    disp('DSSRT:  Bad network object path argument')
    return
end

dir_set = nargin;
while ~systemChosen
    if dir_set == 0
        % get user input for `networkObjectDir`
        result = DialogBox('Network name (a directory name), or <cancel> quits');
        if ~isempty(result)
            networkObjectDir = result;
        else % user pressed cancel, so quit
            disp('DSSRT:  No network specified. Quitting...')
            return
        end
    end
    if exist([DSSRT_ROOTPATH networkObjectDir],'file') == 7 % checks existence of directory, exist returns 7 in that case
        fprintf('DSSRT:  Specified directory path `%s` is OK\n', networkObjectDir)
        systemChosen = true;
        % ensure 'Figures' directory is present in network object dir, for
        % auto-saving of figure files during DSSRT
        if exist([DSSRT_ROOTPATH networkObjectDir '/Figures'],'file')~=7
            mkdir([DSSRT_ROOTPATH networkObjectDir],'Figures')
        end
    else
        disp('DSSRT:  Specified directory path for network does not exist.')
        dir_set = 0;
    end
end % while ~systemChosen

%
% Check parFile if provided
%
if ~ischar(parFile)
    disp('DSSRT:  Bad start-up parameter file argument')
    return
end
parFile = EndStripStr(parFile,'.');
if ~isempty(parFile) && exist([DSSRT_ROOTPATH networkObjectDir '/' parFile '.par']) ~= 2
    fprintf('DSSRT:  Start-up parameter file %s.par not found on network object path\n',parFile)
    return
end

%
% Setup history session log
% If used, delete this history file regularly!
%
if KEEP_LOG
    disp('DSSRT:  History file for session being kept')
    disp('DSSRT:  Starting StepNet')
    diary DSSRT_History.txt
else
    disp('DSSRT:  History file for session not being kept')
    disp('DSSRT:  Starting StepNet')
end

disp(' ')

%
% Run StepNet
%
options.rootPath         = DSSRT_ROOTPATH;
options.networkObjectDir = networkObjectDir;
options.parFile          = parFile;
options.startVerbose     = START_VERBOSE;
StepNet(options)

if KEEP_LOG
    diary off
end

%
% Clean up before quitting
%
warning on all % turn warnings back on
% the following warnings are off by default, in Matlab
warning off MATLAB:nonScalarConditional
warning off MATLAB:usinglongnames

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = DialogBox(promptTitle,defVal)
%% This function assumes that negative numbers will not be legal
%% values for any parameters to be got from user input!
result = '';

comLineStr1 = 'DSSRT:  Enter ';
comLineStr2 = ' in dialog box';
disp([ comLineStr1 lower(promptTitle) comLineStr2 ]);

prompt  = [ promptTitle '?' ];
ptitle  = 'Network object folder:';
lines   = 1;
if nargin == 2
    if isnumeric(defVal)
        def = {num2str(defVal)};
    elseif ischar(defVal)
        def = {defVal};
    else
        disp('DialogBox:  Internal error. `defVal` must be a numeric or string')
        return
    end
else
    def = {''};
end

answer  = inputdlg(prompt,ptitle,lines,def,'off');
if ~isempty(answer)
    if ~isempty(answer{1})
        result = answer{1};
    end
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function success = PreprocessForComp(mFile,globFile,testStr)
success = false; % initial value
if nargin ~=3
    disp('Incorrect number of arguments to PreprocessForComp()')
    return
end

fid = fopen(mFile,'r');
file_done = false;
globals_found = false;
lineCount = 0;
outputText = {};

while ~file_done & ~globals_found
    fline = fgetl(fid);
    if ~ischar(fline) % then EOF (could use feof)
        file_done = true; % redundant!
        disp('DSSRT:  Premature end of .m file reached without finding globals script command')
        disp('        Cannot continue')
        fclose(fid);
        return
    end
    foundTestStr = strfind(fline,testStr);
    if isempty(foundTestStr)
        % then a regular line
        lineCount = lineCount + 1;
        outputText = {outputText{:}, fline};
    else % found the point at which to write in globals definitions
        globals_found = true;
    end
end

% now add the global constant definitions from the file
globid = fopen(globFile,'r');
addedGlobsDone = false;
while ~addedGlobsDone
    gline = fgetl(globid);
    if ~ischar(gline) % then EOF (could use feof)
        addedGlobsDone = true;
        fclose(globid);
    else
        outputText = {outputText{:}, gline};
        lineCount = lineCount + 1;
    end
end

% file_done is still not true if we get to this point
while ~file_done % a redundant check
    fline = fgetl(fid);
    if ~ischar(fline) % then EOF (could use feof)
        file_done = true;
    end
    outputText = {outputText{:}, fline}; % finish reading in .m file
    lineCount = lineCount + 1;
end
fclose(fid);

fidw = fopen( mFile ,'w');
if fidw ~= -1
	for line=1:lineCount
        fprintf(fidw, '%s\r\n', outputText{line});
	end
	fclose(fidw);
else
    fprintf('DSSRT:  Problem opening file %s for writing. Cannot continue\n',mFile);
end

success = true; % if we've got this far
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% strip trailing extension in a filename, using findChar to mark the beginning of the extension
function stripStr = EndStripStr(inputStr,findChar)

ix = findstr(findChar,inputStr);
if ~isempty(ix)
    if max(ix)>1
        stripStr = inputStr(1:max(ix)-1); % keep all other occurrences of findChar except the final one
    else
        stripStr = inputStr;
    end
else
    stripStr = inputStr;
end
return
