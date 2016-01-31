function mapP = linkStateMapPots
global ST_ERR NOT_AP POTENT ACTIVE
% ST_ERR  =-1;
% NOT_AP  = 0;
% POTENT  = 1;
% ACTIVE  = 2;

mapP = {  [ [ NOT_AP POTENT ]; [ POTENT ST_ERR ]; [ ACTIVE ACTIVE ] ]  }; % must return cell array, even if only one entry