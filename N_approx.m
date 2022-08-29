function Neq = N_approx(L0,N0,Lt,k,b0)
    Neq = N0*Lt/(12*L0)*(14 + exp(-3*k/b0)+3*exp(-k/b0));
end