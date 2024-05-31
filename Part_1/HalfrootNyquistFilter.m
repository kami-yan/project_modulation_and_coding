function [G] = HalfrootNyquistFilter(T,beta,f)

if abs(f) >= 0 && abs(f) < (1-beta)/(2*T)
    G = sqrt(T);
end
if abs(f) <= (1+beta)/(2*T) && abs(f) >= (1-beta)/(2*T)
    G = sqrt(T/2*(1+cos(pi*T/beta*(abs(f)-(1-beta)/(2*T)))));
end
if abs(f) > (1+beta)/(2*T)
    G = 0;
end
end
