function sinf = sinf_e(V)
sa = 5.0*(1+tanh(V/4));
sb = 0.6;
sinf = sa/(sa+sb);
return