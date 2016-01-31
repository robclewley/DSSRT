function ninf = ninf(V)
if V == -52
    na = 0.16;
else
    na=.032*(V+52)/(1-exp(-(V+52)/5));
end
nb=.5*exp(-(57+V)/40);
ninf = na/(na+nb);
return
