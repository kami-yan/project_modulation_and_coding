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
    while while_cond == 1
        loop = loop +1;
        
        c_bis = zeros(k+2,n);
    
        for v = 1:k
            q_0 = 1/2;
            q_0_bis = 1;
            proba = zeros(1,n);
            coord = find(H(v,:));

            for w = 1:length(coord)
            
                proba(coord(w)) = (1 - 2*(1 - 1 / (1 + exp(2*c(coord(w))/variance)))) ; % (1-2q(1))
                q_0_bis = q_0_bis * proba(coord(w)); % compute PI(pour tout i) (1-2q(1))
                c_bis(k+2,coord(w)) = c_bis(k+2,coord(w)) + 1; % we have compute c(w) for the c_bis(k+1,w) times
            end

            for l = 1:length(coord)
                c_bis(k+1,coord(l)) = proba(coord(l));
                c_bis(v,coord(l)) = 1/2 * q_0_bis / proba(coord(l)); % cancel the participation of ci => compute 1/2*PI(pour tout i' excluant i) (1-2q(1))
                c_bis(v,coord(l)) = q_0 + q_0_bis; % +1/2 : r(0)
            end
          
        end
    
        

        y = zeros(3,n);
        for w = 1:n
            y(1,w) = 1;
            y(2,w) = 1;
            for v = 1:k
                if c_bis(v,w)~=0
                    y(1,w) = y(1,w) * c_bis(v,w); % q(0)
                end
            end
            y(2,w) = y(1,w);
            
            y(1,w) = y(1,w) * c_bis(k+1,w); % Q(0) not normalized
            y(2,w) = y(2,w) * (1-c_bis(k+1,w)); % Q(1) not normalized
            y(3,w) = y(1,w) / (y(1,w)+y(2,w)); % Q(0) normalized
          
        end

        
         % the k+1 th is q(1)
         % need to multiply by q(1) (first row of y) then by q(0) (second
         % row) then normalize

        u = zeros(1,n);
        for w = 1:n
            vote = y(3,w);
            if vote> 0.5
                bits_stream_temp(m,w) = 0; % majority of 1, so the decision is to put a 1
                u(w)=0;
            else
                bits_stream_temp(m,w) = 1; % majority of 0, so the decision is to put a 0
                u(w)=1;
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
