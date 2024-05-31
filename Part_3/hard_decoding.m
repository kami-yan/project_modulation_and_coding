function [bits_stream] = decoding(stream_encoded,H,k_limit)


[k,n] = size(H);
bits_stream_temp = zeros(length(stream_encoded)/n,n);

% tanner
for m = 1:length(stream_encoded)/n 
    c = zeros(1,n);
    for w = 1:n
        c(w) = stream_encoded((m-1)*n+w);
    end
   
    while_cond = 1;
    loop = 0;
    while while_cond == 1
        loop = loop +1;
        
        c_bis = zeros(k+2,n);
    
        for v = 1:k
            for w = 1:n
                if H(v,w) ~=0 
                    c_bis(k+2,w) = c_bis(k+2,w) + 1; % we have compute c(w) for the c_bis(k+1,w) times
                    c_bis(c_bis(k+2,w),w) = c * H(v,:)' - c(w); % cancel the participation to himself (due to identity matrix)
                    c_bis(c_bis(k+2,w),w) = mod(c_bis(c_bis(k+2,w),w),2);           

                    
                end
            end
        end
    
        c_bis(k+1,:) = c;

        y = zeros(1,n);
        for w = 1:n
            for v = 1:k+1
                y(w) = y(w) + c_bis(v,w);
            end
            y(w) = y(w) / (c_bis(k+2,w)+1); % (plus 1 because we had the received signal)
        end

        u = zeros(1,n);
        for w = 1:n
            vote = y(w);
            if vote> 0.5
                bits_stream_temp(m,w) = 1; % majority of 1, so the decision is to put a 1
                u(w)=1;
            else
                bits_stream_temp(m,w) = 0; % majority of 0, so the decision is to put a 0
                u(w)=0;
            end
        end
        
        
       
        s = u * H';
        s = mod(s,2);

        if nnz(s) == 0 || loop >= k_limit
            while_cond = 0;     
        end

    c = u; % the information decided become the new information received for the next iteration
      
    end % end of while loop
    for w = 1:n
        u(w) = bits_stream_temp(m,w);
    end

end % end loop for m = 1:length(stream_encoded) - n
bits_stream = reshape(bits_stream_temp',[1 size(bits_stream_temp,1) * size(bits_stream_temp,2)]);


end % end of the decoding function
