K_H = 2^6;
N_H = K_H*2;

h = waitbar(0, 'Starting...', 'Name', 'Loading Progress');

%_____________________adjustible variables_________________________________

bits_p_sym = 4; % number of bits per symbol
nbr_sym_data = 2^10;


%______________________definition of error related vector___________________
rate_error_uncoded = []; % list of error rate depending on BER
rate_error_1it = [];
rate_error_5it = [];
rate_error_10it = [];
rate_error_20it = [];

%______________________declaration of varibales_____________________________



f_carrier = 600E6; % Hz



delta_w = 2*pi*10^-5; % Here delta_w = 10 ppm (10^-5)





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

for i = 1:length(vect)/K_H
    final_vect = [final_vect,paritybits(i,:),infobits(i,:)];
end 

final_vect = final_vect';

%______________________mapping_______________________________



map = mapping(final_vect,bits_p_sym,'qam');
ma2= [];
for i = map'
    ma2 = [ma2,0,0,0,0,i,0,0,0,0,0];
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
% plot(f,G);
g = ifft(G);
g = ifftshift(g); % for the signal in positive timelin


s= conv(ma2,g,'same');
s_uncoded = conv(ma3,g,'same');

fcarrier = 600e6;


V = var(s) / 2; % variance of the signal
V_uncoded = var(s_uncoded) / 2; % variance of the signal


Eb = V * T  / bits_p_sym; % Calculation of a bit energy
Eb_uncoded = V_uncoded * T  / bits_p_sym; % Calculation of a bit energy

%Eb = trapz(abs(s).^2)*T / bits_p_sym;

nb_point_BER = -1.5:0.5:13; % BER that will be explore to draw curves
nIterations  = length(nb_point_BER);
for k = 1:nIterations

    % We adjust the noise to make the curve
    N0 = Eb / (10^(nb_point_BER(k)/10));
    N0_uncoded = Eb_uncoded / (10^(nb_point_BER(k)/10));

    Y = sqrt(N0/2 * 50E6) .* (randn(1,length(s)) + 1i * randn(1,length(s)));  % Creating white gaussian noise
    Y_uncoded = sqrt(N0_uncoded/2 * 50E6) .* (randn(1,length(s)) + 1i * randn(1,length(s)));  % Creating white gaussian noise
    t = s + Y; % Signal with added noise
    t_uncoded = s + Y_uncoded; % Signal with added noise
    
    progress = k / nIterations;
    waitbar(progress, h, sprintf('Progress: %d%%', round(progress * 100)));
    % We adjust the noise to make the curve

    % Vaut mieux utiliser randn, plus facile à controler
    % La puissance serait juste la variance de l'enveloppe divisé par 2 ?




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


%________________________demapping____________________________________________________

    demap_uncoded = demapping(r_uncoded,bits_p_sym,'qam');
    demap = demapping(r,bits_p_sym,'qam');


%_______________________decoding_____________________________________________________


    [bits_stream_1it] = decoding(demap,H,1);
    [bits_stream_5it] = decoding(demap,H,5);
    [bits_stream_10it] = decoding(demap,H,10);
    [bits_stream_20it] = decoding(demap,H,100);
    % 
    % infob = [];
    % 
    % for i = 1:length(bits_stream)/N_H
    %     infob = [infob,bits_stream(i+K_H:i+N_H-1)];
    % end

        
%____________________error computation_________________________________________________ 
    
    
    % We count the number of error
    error_vct_uncoded = vect - demap_uncoded;
    error_vct_1it = final_vect' - bits_stream_1it;
    error_vct_5it = final_vect' - bits_stream_5it;
    error_vct_10it = final_vect' - bits_stream_10it;
    error_vct_20it = final_vect' - bits_stream_20it;
    

    nb_error = zeros(1,5);
    error_vct = [error_vct_uncoded';error_vct_1it;error_vct_5it;error_vct_10it;error_vct_20it];
    for i = 1:5
    for l = error_vct(i,:)
        if l ~= 0
            nb_error(i) = nb_error(i) + 1;
        end
    end
    end

    rate_error_uncoded = [rate_error_uncoded,nb_error(1)/length(error_vct_uncoded)];
    rate_error_1it = [rate_error_1it,nb_error(2)/length(error_vct_1it)];
    rate_error_5it = [rate_error_5it,nb_error(3)/length(error_vct_5it)];
    rate_error_10it = [rate_error_10it,nb_error(4)/length(error_vct_10it)];
    rate_error_20it = [rate_error_20it,nb_error(5)/length(error_vct_20it)];

    
%partie apres l'algorithme




end


loglog(rate_error_uncoded)
hold on 
loglog(rate_error_1it)
hold on 
loglog(rate_error_5it)
hold on 
loglog(rate_error_10it)
hold on 
loglog(rate_error_20it)
grid on 
legend('uncoded','1it','5it','10it','20it')


