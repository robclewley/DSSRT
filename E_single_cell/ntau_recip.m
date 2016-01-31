function ntau_recip = ntau_recip(V)
if V == -52
    na=0.1600;
else
    na=.032*(V+52)/(1-exp(-(V+52)/5));
end
nb=.5*exp(-(57+V)/40);
ntau_recip = na+nb;
return
