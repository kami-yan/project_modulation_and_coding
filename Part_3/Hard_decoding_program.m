K_H = 2^6; % 2^6
N_H = K_H*2;



%_____________________adjustible variables_________________________________

bits_p_sym = 4; % number of bits per symbol. 4
nbr_sym_data = 2^15;


%______________________definition of error related vector___________________
rate_error_uncoded = []; % list of error rate depending on BER
rate_error_1it = [];
rate_error_2it = [];
rate_error_5it = [];
rate_error_10it = [];
rate_error_1it_soft = [];
rate_error_2it_soft  = [];
rate_error_3it_soft  = [];
rate_error_4it_soft  = [];

%______________________declaration of varibales_____________________________



f_carrier = 600E6; % Hz




%______________________generation of signal with pilot________________


vect = randi([0,1],bits_p_sym*nbr_sym_data,1); % random vector of bits




%___________________encoding________________________________

final_vect = [];
infobits = reshape(vect,[K_H,length(vect)/K_H]);


% Create initial parity check matrix of size 128 x 256
H0 = generate_ldpc(K_H, N_H,0,1,3);

% Compute parity bits (128 x nb of packets) and generate final parity check matrix
[paritybits, H] = encode_ldpc(infobits, H0, 0);

infobits = infobits';
paritybits = paritybits';

for l = 1:length(vect)/K_H
    final_vect = [final_vect,paritybits(l,:),infobits(l,:)];
end 

final_vect = final_vect';

%______________________mapping_______________________________



map = mapping(final_vect,bits_p_sym,'qam');
ma2= [];
for l = map'
    ma2 = [ma2,0,0,0,0,l,0,0,0,0,0];
end

map3 = mapping(vect,bits_p_sym,'qam');
ma3= [];
for i_uncoded = map3'
    ma3 = [ma3,0,0,0,0,i_uncoded,0,0,0,0,0];
end





% demap = demapping(map,4,'qam');
% 64 bit per symbol
% 5 mega symbol per second
T = 200E-9; % (duration of a symbol)
beta = 0.3 ; % (0<=beta<=1)
f = -25E6:25E4:25E6-25E4;
G = [];
for l = f
    G = [G,HalfrootNyquistFilter(T,beta,l)];
end

G = fftshift(G); % for the symetric
g = ifft(G);
g = ifftshift(g); % for the signal in positive timelin


s= conv(ma2,g,'same');
s_uncoded = conv(ma3,g,'same');

fcarrier = 600e6;


V = var(s) / 2; % variance of the signal
V_uncoded = var(s_uncoded) / 2; % variance of the signal


Eb = V * T  / bits_p_sym; % Calculation of a bit energy
Eb_uncoded = V_uncoded * T  / bits_p_sym; % Calculation of a bit energy



nb_point_BER = -3:1:14; % BER that will be explore to draw curves
nIterations  = length(nb_point_BER);
for k = 1:nIterations

    % We adjust the noise to make the curve
    N0 = Eb / (10^(nb_point_BER(k)/10));
    N0_uncoded = Eb_uncoded / (10^(nb_point_BER(k)/10));

    Y = sqrt(N0/2 * 50E6) .* (randn(1,length(s))+1i*randn(1,length(s))) ;  % Creating white gaussian noise
    Y_uncoded = sqrt(N0_uncoded/2 * 50E6) .* (randn(1,length(s_uncoded))+1i*randn(1,length(s_uncoded))) ;  % Creating white gaussian noise
    t = s + Y; % Signal with added noise
    t_uncoded = s_uncoded + Y_uncoded; % Signal with added noise
    

    g_invert = flip(g); % flip function inverse the order of the vecteur g(-t). In fact g(-t) = g(t)
    r_10 = conv(t,g_invert,'same'); % convolution with g(-t). Caution with rate !
    r_10_uncoded = conv(t_uncoded,g_invert,'same'); % convolution with g(-t). Caution with rate !
    r_10 = r_10 / max(conv(g,g)); %normalization
    r_10_uncoded = r_10_uncoded / max(conv(g,g)); %normalization

   
    % To go back to the intial sampling rate, we take one symbol over 10
    r=[];

    count = 6; % The fourth is the first to be sampled
    for l = r_10
        if rem(count, 10) == 0
            r = [r,l];

        end
        count = count + 1;
    end

    r = r';
% To go back to the intial sampling rate, we take one symbol over 10
    r_uncoded=[];

    count = 6; % The fourth is the first to be sampled
    for l = r_10_uncoded
        if rem(count, 10) == 0
            r_uncoded = [r_uncoded,l];

        end
        count = count + 1;
    end

    r_uncoded = r_uncoded';

%______________________Soft hard_decoding here_____________________________________________
    % 
    % [bits_stream_1it_soft] = soft_decoding_log(r,H,var(Y),1);
    % [bits_stream_2it_soft] = soft_decoding_log(r,H,var(Y),2);
    % [bits_stream_3it_soft] = soft_decoding_log(r,H,var(Y),3);
    % [bits_stream_4it_soft] = soft_decoding_log(r,H,var(Y),4);

%________________________demapping____________________________________________________

    demap_uncoded = demapping(r_uncoded,bits_p_sym,'qam');
    demap = demapping(r,bits_p_sym,'qam');


%_______________________hard_decoding_____________________________________________________


    [bits_stream_1it] = hard_decoding(demap,H,1);
    [bits_stream_2it] = hard_decoding(demap,H,2);
    [bits_stream_5it] = hard_decoding(demap,H,5);
    [bits_stream_10it] = hard_decoding(demap,H,10);
    


        
%____________________error computation_________________________________________________ 
    
    
    % We count the number of error
    error_vct_uncoded = vect - demap_uncoded;
    error_vct_1it = final_vect' - bits_stream_1it;
    error_vct_2it = final_vect' - bits_stream_2it;
    error_vct_5it = final_vect' - bits_stream_5it;
    error_vct_10it = final_vect' - bits_stream_10it;
    
    % error_vct_1it_soft = final_vect' - bits_stream_1it_soft;
    % error_vct_2it_soft = final_vect' - bits_stream_2it_soft;
    % error_vct_3it_soft = final_vect' - bits_stream_3it_soft;
    % error_vct_4it_soft = final_vect' - bits_stream_4it_soft;

    nb_error = zeros(1,9);
    error_vct_uncoded = paddata(error_vct_uncoded,length(error_vct_1it));
    error_vct = [error_vct_uncoded';error_vct_1it;error_vct_2it;error_vct_5it;error_vct_10it];
    % error_vct = [error_vct_uncoded';error_vct_1it_soft;error_vct_2it_soft;error_vct_3it_soft;error_vct_4it_soft;];

    for p = 1:5
        for l = error_vct(p,:)
            if l ~= 0
                nb_error(p) = nb_error(p) + 1;
            end
        end
    end

    rate_error_uncoded = [rate_error_uncoded,nb_error(1)/length(error_vct_uncoded)];
    rate_error_1it = [rate_error_1it,nb_error(2)/length(error_vct_1it)];
    rate_error_2it = [rate_error_2it,nb_error(3)/length(error_vct_2it)];
    rate_error_5it = [rate_error_5it,nb_error(4)/length(error_vct_5it)];
    rate_error_10it = [rate_error_10it,nb_error(5)/length(error_vct_10it)];

    % rate_error_1it_soft = [rate_error_1it_soft ,nb_error(2)/length(error_vct_1it_soft )];
    % rate_error_2it_soft  = [rate_error_2it_soft ,nb_error(3)/length(error_vct_2it_soft )];
    % rate_error_3it_soft  = [rate_error_3it_soft ,nb_error(4)/length(error_vct_3it_soft )];
    % rate_error_4it_soft  = [rate_error_4it_soft ,nb_error(5)/length(error_vct_4it_soft )];

    




end



semilogy(nb_point_BER,rate_error_uncoded)
hold on 
grid on 
semilogy(nb_point_BER,rate_error_1it)
semilogy(nb_point_BER,rate_error_2it)
semilogy(nb_point_BER,rate_error_5it)
semilogy(nb_point_BER,rate_error_10it)

% semilogy(nb_point_BER,rate_error_1it_soft)
% semilogy(nb_point_BER,rate_error_2it_soft)
% semilogy(nb_point_BER,rate_error_3it_soft)
% semilogy(nb_point_BER,rate_error_4it_soft)

legend('uncoded','soft 1it','soft 2it','soft 3it','soft 4it')
legend('uncoded','1it','2it','5it','10it')
hold off
xlabel('SNR')
ylabel('BER')
