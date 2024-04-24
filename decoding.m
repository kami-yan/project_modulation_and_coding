function [bits_stream] = decoding(stream_encoded,P,H)

[k,n] = size(P);
%H = [eye(n) P'];
[k,n] = size(H);
y = zeros(4,n);
bits_stream_temp = zeros(size(stream_encoded));

% tanner
for m = 1:length(stream_encoded) - n
    c = zeros(1,n);
    c_bis = zeros(2,n);
    for w = 1:n
        c(w) = stream_encoded(m+w-1);
    end
    
    % f = c * H';
    % 
    % % Modulo 2
    % for l = 1:length(f)
    %     if rem(f(l),2) == 0
    %         f(l) = 0;
    %     else
    %         f(l) = 1; 
    %     end
    % end
    % % end of modulo 2

    y(1,:) = c; % what we received


    for v = k
        for w = 1:n
            if H(v,w) ~=0 && c_bis(2,w) == 0 % To be sure that we went for the first time trhough c(w)
                c_bis(2,w) = 1; % we have compute c(w) for the first time
                c_bis(1,w) = c * H(v,:)' - c(w); % cancel the participation to himself
                % modulo 2
                if rem(c_bis(1,w),2) == 0
                        c_bis(1,w) = 0;
                    else
                        c_bis(1,w) = 1; 
                end
            end
        end
    end

    y(2,:) = c_bis(1,:); 

    for v = k
        for w = 1:n
           
            c_bis(1,w) = c * H(v,:)' - c(w); % cancel the participation to himself
            % modulo 2
            if rem(c_bis(1,w),2) == 0
                    c_bis(1,w) = 0;
                else
                    c_bis(1,w) = 1; 
            end
            
        end
    end

    y(3,:) = c_bis(1,:); 
    
    for w = 1:n
        bits_stream_temp(m+n-1) = y(3,n) + y(2,n) + y(1,n);
        if bits_stream_temp(m+n-1)>=2
            bits_stream_temp(m+n-1) = 1; % majority of 1, so the decision is to put a 1
        else
            bits_stream_temp(m+n-1) = 0; % majority of 0, so the decision is to put a 0
        end
    end
    

end % end loop for m = 1:length(stream_encoded) - n

% decoding
bits_stream = zeros(length(bits_stream_temp)/n,k);
for w = 1:length(bits_stream_temp)/n
    u = zeros(1,n);
    for m = 1:n
        u(m) = bits_stream_temp((w-1)*n+m);
    end
    s = u * H';
    for m = 1:k
        if rem(s(m),2) == 0
            s(m) = 0;
        else
            s(m) = 1; 
        end
    end
    bits_stream(w,:) = s;
end
bits_stream = bits_stream(:);

end
