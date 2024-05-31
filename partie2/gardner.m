function [corrected_signal,epsilon] = gardner(received_signal,gardner_upsample,k_gardner)
n = length(received_signal);
real_size = n/gardner_upsample;
corrected_signal = zeros(1,real_size);
corrected_signal(1) = received_signal(1);
epsilon = zeros(1,real_size);
for n = 1:real_size-1
    points = ((n-1)*gardner_upsample+1:(n+1)*gardner_upsample);
    sequence = received_signal(points);
    OneHalf_n_estimated = (n-1/2)*gardner_upsample+1-epsilon(n);
    actuel_n_estimated = n*gardner_upsample+1-epsilon(n);
    corrected_signal(n+1) = interp1(points,sequence,actuel_n_estimated,'linear');
    corrected_signal_OneHalf = interp1(points,sequence,OneHalf_n_estimated,'linear');
    epsilon(n+1)= epsilon(n)+2*k_gardner*real(corrected_signal_OneHalf*(conj(corrected_signal(n+1))-conj(corrected_signal(n)))) ;
end
    epsilon = epsilon/2; % epsilon is actually the time error because we didnt include T in the expression
end