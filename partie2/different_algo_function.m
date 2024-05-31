function [n_estimate,estimate_delta_f,phi_0_estimate, max_t,total_vector] = different_algo_function (pilot_vector_mapped,received_signal,T,K,nbr_sym_data)    

    N = length(pilot_vector_mapped);
    n_possible =  length(received_signal);
    total_vector = zeros(1,n_possible);
    
    for n = 1:(n_possible-N-1)
        sum_dk = 0;
        for k = 1:K
            term = 0;
            for l = k+1 : N
                term1 = (conj(received_signal(n+l))*pilot_vector_mapped(l))*conj(conj(received_signal(n+l-k))*pilot_vector_mapped(l-k));
                term = term +term1;
            end 
            dk = 1/(N-k)*term;
            sum_dk = sum_dk + abs(dk);
        end
        total_vector(n) = sum_dk;
    end
   
    
    [max_value, n_estimate] = max(total_vector) ;
    n_estimate = n_estimate +1;
    max_t = max_value;
    estimate_delta_f = 0;
    for k = 1:K
            term = 0;
            for l = k+1 : N
                term1 = (conj(received_signal(n_estimate+l))*pilot_vector_mapped(l))*conj(conj(received_signal(n_estimate+l-k))*pilot_vector_mapped(l-k));
                term = term +term1;
            end 
    dk = 1/(N-k)*term;
    estimate_delta_f = estimate_delta_f + angle(dk)/(2*pi*k*T);
    end
    estimate_delta_f = (-1/K)*estimate_delta_f;

    phi_0_estimate = angle(received_signal(n_estimate)*exp(-1j*2*pi*estimate_delta_f*T)/pilot_vector_mapped(1)); % a verifier
   
     % plot(total_vector);
     % figure;
end