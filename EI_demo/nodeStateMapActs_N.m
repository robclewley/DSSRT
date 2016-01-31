function mapA_N = nodeStateMapActs_N

STATE_ERR        =-1;
NODE_NONE_AP     = 0;
NODE_EITH_P      = 1;
NODE_ONLY_ACT_M  = 2;
NODE_ONLY_ACT_N  = 3;
NODE_BOTH_A      = 4;
NODE_ACT_M_POT_N = 5;
NODE_ACT_N_POT_M = 6;

mapA_N = [[ NODE_NONE_AP     NODE_ONLY_ACT_N ]; ...
          [ NODE_EITH_P      STATE_ERR ]; ...
          [ NODE_ONLY_ACT_M  NODE_BOTH_A ]; ...
          [ NODE_ONLY_ACT_N  STATE_ERR ]; ...
          [ NODE_BOTH_A      NODE_BOTH_A ];
          [ NODE_ACT_M_POT_N STATE_ERR   ]; ...
          [ NODE_ACT_N_POT_M STATE_ERR   ]         ];
