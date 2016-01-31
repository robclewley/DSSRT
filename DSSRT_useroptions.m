%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    DSSRTuseroptions.m --- user-editable system options script    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Option for program development diagnostics -- forces recompilation of principal .m files
% (Matlab release 13 and earlier only)
% Usually `false`
FORCE_COMPILE = true;

% Option for program development diagnostics (over-rides FORCE_COMPILE option)
% or if no Matlab compiler is available, etc.
% (Only set false for Matlab release 13 and earlier)
FORCE_NOTCOMPILE = true;

% Option to do / check for any compilations, but do not start DSSRT
% (Matlab release 13 and earlier only)
% Usually `false`
COMPILE_ONLY = false;

% User-specification of DSSRT home directory
% Edit this next line according to your operating system and installation directory
DSSRT_ROOTPATH = 'c:\DSSRT\';

% Use this option to control the default verbosity of console output
% (can be toggled within DSSRT using command `?`)
% Setting this to be `true` is good for newbies
START_VERBOSE = true;

% Use this option to keep a 'diary' log of all console output in `DSSRT_History.txt`
% This is useful in case of an unrecoverable program or Matlab crash.
KEEP_LOG = false;
