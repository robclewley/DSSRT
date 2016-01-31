function mapA = linkStateMapActs
global ST_ERR NOT_AP POTENT ACTIVE
% ST_ERR  =-1;
% NOT_AP  = 0;
% POTENT  = 1;
% ACTIVE  = 2;

mapA = {  [ [ NOT_AP ACTIVE ]; [ POTENT ST_ERR ]; [ ACTIVE ST_ERR ] ]  }; % must return cell array, even if only one entry