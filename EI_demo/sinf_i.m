function sinf = sinf_i(V)
sa = 2.0*(1+tanh(V/4));
sb = 0.0667;
sinf = sa/(sa+sb);
return