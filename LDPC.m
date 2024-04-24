
rate_error = []; % list of error rate depending on BER

f_carrier = 2E9; % Hz
T = 200E-9; % (duration of a symbol)


nb_bits = 2^5; % Number of bits send 2^6=64
bits_p_sym = 4; % number of bits per symbol, so the number of symbol is 512/8 = 64

vect = [randi([0,1],5,1)';
    randi([0,1],5,1)';
    randi([0,1],5,1)';]; % random vector of bits


% % Create initial parity check matrix of size 128 x 256
% H0 = generate_ldpc(128, 256,0,1,3);
% % Compute parity bits (128 x nb of packets) and generate final parity check matrix
% [paritybits, H] = encode_ldpc(vect, H0, 0);

H =[1 1 0 1 1 0 0 1 0 0;
    0 1 1 0 1 1 1 0 0 0;
    0 0 0 1 0 0 0 1 1 1;
    1 1 0 0 0 1 1 0 1 0;
    0 0 1 0 0 1 0 1 0 1];



P =[0 1 1 1 0;
    1 0 1 0 0;
    1 0 1 0 1;
    0 0 1 1 1;
    1 1 0 0 1].'; 

I_5 = eye(5);

Gen = [P I_5];

vect_enc = vect * Gen;

        
vect_enc = vect_enc(:)';
for m = 1:length(vect_enc)
    if rem(vect_enc(m),2) == 0
        vect_enc(m) = 0;
    else
        vect_enc(m) = 1; 
    end
end

decoded = decoding(vect_enc,P,H);
vect_verif = vect(:);






% 
% vect_enc = [vect_enc,0,1];
% 
% 
% 
% 
% for l = 1:length(vect_enc)
%     if rem(vect_enc(l),2) == 0
% 
%         vect_enc(l) = 0;
%     else
% 
%         vect_enc(l) = 1;
% 
%     end
% 
% end
% 
% vect_enc = vect_enc';
% 
% map = mapping(vect_enc,bits_p_sym,'qam');
% ma2= [];
% for i = map'
%     ma2 = [ma2,0,0,0,0,i,0,0,0,0,0];
% end
% % demap = demapping(map,4,'qam');
% % 64 bit per symbol
% % 5 mega symbol per second
% 
% beta = 0.3 ; % (0<=beta<=1)
% f = -25E6:25E4:25E6-25E4;
% G = [];
% for l = f
%     G = [G,HalfrootNyquistFilter(T,beta,l)];
% end
% 
% G = fftshift(G); % for the symetric
% g = ifft(G);
% g = ifftshift(g); % for the signal in positive timelin
% 
% 
% s= conv(ma2,g,'same');
% 
% V = var(s) / 2; % variance of the signal
% 
% 
% Eb = V * T  / bits_p_sym; % Calculation of a bit energy
% 
% 
% nb_point_BER = -1.5:0.5:13; % BER that will be explore to draw curves
% for k = nb_point_BER 
%     % We adjust the noise to make the curve
%     N0 = Eb / (10^(k/10));
% 
%     Y = sqrt(N0/2 * 50E6) .* (randn(1,length(s)) + 1i * randn(1,length(s)));  % Creating white gaussian noise
%     t = s + Y; % Signal with added noise
% 
%     g_invert = conj(flip(g)); % flip function inverse the order of the vecteur g(-t). In fact g(-t) = g(t)
%     r_10 = conv(t,g_invert,'same'); % convolution with g(-t). Caution with rate !
%     r_10 = r_10 / max(conv(g,g)); %normalization
% 
% 
%     % To go back to the intial sampling rate, we take one symbol over 10
%     r=[];
% 
%     count = 6; % The fourth is the first to be sampled
%     for l = r_10
%         if rem(count, 10) == 0
%             r = [r,l];
% 
%         end
%         count = count + 1;
%     end
% 
%     r = r';
% 
% 
%     demap = demapping(r,bits_p_sym,'qam');
% 
% 
%     % We count the number of error
%     error_vct = vect_enc - demap;
% 
% 
%     nb_error = 0;
% 
%     for l = error_vct'
%         if l ~= 0
%             nb_error = nb_error + 1;
%         end
%     end
% 
%     rate_error = [rate_error,nb_error/length(error_vct)];
% 
% end
% 
% 
% %semilogy(nb_point_BER,rate_error)
% 
% 
% 
% 
