function minf = minf(V)
if V == -54
    minf = 0.14423680;
    return
else
    ma=0.32*(V+54)/(1-exp(-(V+54)/4));
end
if V == -27
    minf = 0.8606983;
    return
else
    mb=0.28*(V+27)/(exp((V+27)/5)-1);
end
minf = ma/(ma+mb);
return
