function mapA_U = nodeStateMapActs_U

STATE_ERR        =-1;
NODE_NONE_AP     = 0;
NODE_POT_U       = 1;
NODE_ACT_U       = 2;

% maps first column value to second column's value
% Note that each row's first column is just the row index!
% ... You can list these mappings in any order if you
% want to.
mapA_U = [[ NODE_NONE_AP     NODE_ACT_U ]; ...
          [ NODE_POT_U       STATE_ERR ]; ...
          [ NODE_ACT_U       STATE_ERR ];   ];
  