function [G] = HalfrootNyquistFilter(T,beta,f)

if  abs(f) <= (1-beta)/(2*T)
    G = sqrt(T);
elseif abs(f) <= (1+beta)/(2*T) 
    G = sqrt(T/2*(1+cos(pi*T/beta*(abs(f)-(1-beta)/(2*T)))));

else 
    G = 0;
end
end
