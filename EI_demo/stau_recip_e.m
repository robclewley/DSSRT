function stau_recip = stau_recip_e(V)
sa = 5.0*(1+tanh(V/4));
sb = 0.6;
stau_recip = sa+sb;
return