clearvars ;

h = waitbar(0, 'Starting...', 'Name', 'Loading Progress');
waitbar(0, h, sprintf('Progress: %d%%', 0));
%_____________________adjustible variables________________________________
bits_p_sym = 2; % number of bits per symbol
nbr_sym_data = 2^14;
up_sample_factor = 10;
f_cut = 5E6; % Hz
taps =101;
beta = 0.3 ; % (0<=beta<=1)
fcarrier = 600e6;

%______________________declaration of varibales_____________________________
f_sampling = f_cut*up_sample_factor;
T = 1/f_cut;

%______________________generation of signal________________
vect = randi([0,1],bits_p_sym*nbr_sym_data,1); % random vector of bits

%______________________mapping_______________________________
map3 = mapping(vect,bits_p_sym,'qam');
ma3= upsample(map3,up_sample_factor);

f = linspace(-f_cut/2,f_cut/2,taps);
G = zeros(1,taps);
for l = 1:taps
    G(l) = HalfrootNyquistFilter(T,beta,f(l));
end

G = fftshift(G); % for the symetric
g = ifft(G);
g = ifftshift(g); % for the signal in positive timeline

s_uncoded = conv(ma3',g,'same');

V_uncoded = var(s_uncoded) / 2; % variance of the signal
Eb_uncoded = V_uncoded * T  / bits_p_sym; % Calculation of a bit energy

nb_point_BER = -1.5:0.1:16; % BER that will be explore to draw curves
nIterations  = length(nb_point_BER);
%______________________definition of error related vector___________________
rate_error_uncoded = ones(1,nIterations); % list of error rate depending on BER


for k = 1:nIterations
    progress = k / nIterations;
    waitbar(progress, h, sprintf('Progress: %d%%', round(progress * 100)));

    % We adjust the noise to make the curve
    N0_uncoded = Eb_uncoded / (10^(nb_point_BER(k)/10));

    
    Y_uncoded = sqrt(N0_uncoded/2 * 50E6) .* (randn(1,length(s_uncoded)) + 1i * randn(1,length(s_uncoded)));  % Creating white gaussian noise
    t_uncoded = s_uncoded + Y_uncoded; % Signal with added noise
    
    
    % We adjust the noise to make the curve

    % Vaut mieux utiliser randn, plus facile à controler
    % La puissance serait juste la variance de l'enveloppe divisé par 2 ?

    g_invert = flip(g); % flip function inverse the order of the vecteur g(-t). In fact g(-t) = g(t)
    r_10_uncoded = conv(t_uncoded,g_invert,'same'); % convolution with g(-t). Caution with rate !
    r_10_uncoded = r_10_uncoded / max(conv(g,g)); %normalization

    % To go back to the intial sampling rate, we take one symbol over 10
 
% To go back to the intial sampling rate, we take one symbol over 10
    r_uncoded=downsample(r_10_uncoded,up_sample_factor);
    r_uncoded = r_uncoded';

%________________________demapping____________________________________________________
    demap_uncoded = demapping(r_uncoded,bits_p_sym,'qam');

%____________________error computation_________________________________________________ 
    % We count the number of error
    error_vct_uncoded = vect' - demap_uncoded';
    nb_error =0;
    for l = error_vct_uncoded
        if l ~= 0
            nb_error = nb_error + 1;
        end
    end

    rate_error_uncoded(k) = nb_error/length(error_vct_uncoded);
    
end

loglog(rate_error_uncoded)
legend('uncoded')
hold off
close(h);
