function mtau_recip = mtau_recip(V)
if V == -54
    ma = 1.2800;
else
    ma=0.32*(V+54)/(1-exp(-(V+54)/4));
end
if V == -27
    mb = 1.400;
else
    mb=0.28*(V+27)/(exp((V+27)/5)-1);
end
mtau_recip = ma+mb;
return
