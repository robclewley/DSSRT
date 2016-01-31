function hinf = hinf(V)
ha=.128*exp(-(50+V)/18);
hb=4/(1+exp(-(V+27)/5));
hinf = ha/(ha+hb);
return
