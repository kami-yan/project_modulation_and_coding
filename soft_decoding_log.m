function [bits_stream] = soft_decoding(stream_encoded,H,variance,k_limit)

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


    L_q = -2.*c(1,:)./variance;

    while while_cond == 1
        loop = loop +1;
        
        r_bis = zeros(k,n);  
        L_r = zeros(1,k);
        for v = 1:k
       
            coord = find(H(v,:));

            L_r(v) = prod(sign(L_q(coord)));          
            
        
            for l = 1:length(coord) % issue with dimension for L_r
                store = L_q(coord(l));
                L_q(coord(l)) = max(L_q); % so that we do not take the walue by error with the min
                L_q_min = min(L_q);

                r_bis(v,coord(l)) = L_r(v) * sign(store) * L_q_min; % cancel the participation of qij => sign^2 and its reliability is stored apart
                L_q(coord(l)) = store;

            end
        end
        
    
        

        y = L_q;
        for w = 1:n
            for v = 1:k
                y(w) = y(w) + r_bis(v,w); % L(Q)
            end
        end

        u = zeros(1,n);
        for w = 1:n
            if y(w)< 0 % strangely, it is this inverse of what inside the course
                bits_stream_temp(m,w) = 0; 
                u(w)=0;
            else
                bits_stream_temp(m,w) = 1; 
                u(w)=1;
            end
        end
        
        
       
        s = u * H';
        s = mod(s,2);

        if nnz(s) == 0 || loop >= k_limit
            while_cond = 0;     
        end

    L_q = y; % the information decided become the new information received for the next iteration
      
    end % end of while loop


end % end loop for m = 1:length(stream_encoded) - n
bits_stream = reshape(bits_stream_temp',[1 size(bits_stream_temp,1) * size(bits_stream_temp,2)]);


end