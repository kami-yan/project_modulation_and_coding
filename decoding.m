function [bits_stream] = decoding(stream_encoded,P,H)


[k,n] = size(H);
y = zeros(3,n);

bits_stream_temp = zeros(length(stream_encoded)/n,n);

% tanner
for m = 1:length(stream_encoded)/n 
    c = zeros(1,n);
    c_bis = zeros(2,n);
    for w = 1:n
        c(w) = stream_encoded((m-1)*n+w);
    end
   
    while_cond = 1;
    loop = 0;
    while while_cond == 1
        loop = loop +1;
        y(1,:) = c; % what we received
        
    
        for v = 1:k
            for w = 1:n
                if H(v,w) ~=0 && c_bis(2,w) == 0 % To be sure that we went for the first time trhough c(w)
                    c_bis(2,w) = 1; % we have compute c(w) for the first time
                    c_bis(1,w) = c * H(v,:)' - c(w); % cancel the participation to himself
                    c_bis(1,w) = mod(c_bis(1,w),2);           
                end
            end
        end
    
        y(2,:) = c_bis(1,:); 
    
        for v = 1:k
            for w = 1:n
                if H(v,w) ~=0
                    c_bis(1,w) = c * H(v,:)' - c(w); % cancel the participation to himself
                    c_bis(1,w) = mod(c_bis(1,w),2);
                end
            end
        end
    
        y(3,:) = c_bis(1,:); 
        
        for w = 1:n
            bits_stream_temp(m,w) = y(3,w) + y(2,w) + y(1,w);
            if bits_stream_temp(m,w)>=2
                bits_stream_temp(m,w) = 1; % majority of 1, so the decision is to put a 1
            else
                bits_stream_temp(m,w) = 0; % majority of 0, so the decision is to put a 0
            end
        end
        
        u = zeros(1,n);
        for w = 1:n
            u(w) = bits_stream_temp(m,w);
        end
    
        s = u * H';
        s = mod(s,2);

        if nnz(s) == 0 || loop > 10
            while_cond = 0;
            c = u; % the information decided become the new information received for the next iteration
        end


    end % end of while loop

end % end loop for m = 1:length(stream_encoded) - n
bits_stream = reshape(bits_stream_temp',[1 size(bits_stream_temp,1) * size(bits_stream_temp,2)]);


end % end of the decoding function
