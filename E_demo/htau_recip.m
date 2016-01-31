function htau_recip = htau_recip(V)
ha=.128*exp(-(50+V)/18);
hb=4/(1+exp(-(V+27)/5));
htau_recip = ha+hb;
return
