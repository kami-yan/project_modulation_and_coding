clearvars ;

h = waitbar(0, 'Starting...', 'Name', 'Loading Progress');
waitbar(0, h, sprintf('Progress: %d%%', 0));

bits_p_sym_list = [1,2,4,6];
technique = ["pam","qam","qam","qam"];

for bits_p_sym_indice = 1:length(technique)
%_____________________adjustible variables________________________________
bits_p_sym = bits_p_sym_list(bits_p_sym_indice); % number of bits per symbol
nbr_sym_data = 2^17;
up_sample_factor = 10;
f_cut = 5E6; % Hz
taps =51;
beta = 0.3 ; % (0<=beta<=1)
fcarrier = 600e6;

%______________________declaration of varibales_____________________________
f_sampling = f_cut*up_sample_factor;
T = 1/f_cut;

%______________________generation of signal________________
vect = randi([0,1],bits_p_sym*nbr_sym_data,1); % random vector of bits

%______________________mapping_______________________________
map3 = mapping(vect,bits_p_sym,technique(bits_p_sym_indice));

ma3= upsample(map3,up_sample_factor);

f = linspace(-f_cut/2,f_cut/2,taps);
G = zeros(1,taps);
for l = 1:taps
    G(l) = HalfrootNyquistFilter(T,beta,f(l));
end

G = ifftshift(G); % for the symetric
g = ifft(G);
g = fftshift(g); % for the signal in positive timeline

s_uncoded = conv(ma3',g,'same');

V_uncoded = var(s_uncoded) / 2; % variance of the signal
Eb_uncoded = V_uncoded * T  / bits_p_sym; % Calculation of a bit energy

point_BER = -1.5:0.1:16; % BER that will be explore to draw curves
nIterations  = length(point_BER);
%______________________definition of error related vector___________________
rate_error_uncoded = ones(1,nIterations); % list of error rate depending on BER


for k = 1:nIterations
    progress = (bits_p_sym_indice-1)/length(bits_p_sym_list) + 0.25*k / nIterations ;
    waitbar(progress, h, sprintf('Progress: %d%%', round(progress * 100)));

    % We adjust the noise to make the curve
    N0_uncoded = Eb_uncoded / (10^(point_BER(k)/10));

    
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
    if bits_p_sym == 1
        r_uncoded = real(r_uncoded);
    end
    demap_uncoded = demapping(r_uncoded,bits_p_sym,technique(bits_p_sym_indice));

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
semilogy(point_BER,rate_error_uncoded)
hold on
end

legend('BPSK','QPSK','16QAM','64QAM')
hold off
close(h);
