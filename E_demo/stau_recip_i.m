function stau_recip = stau_recip_i(V)
sa = 2.0*(1+tanh(V/4));
sb = 0.0667;
stau_recip = sa+sb;
return