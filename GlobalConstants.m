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


