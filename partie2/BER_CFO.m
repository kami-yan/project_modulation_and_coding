clearvars ;

h = waitbar(0, 'Starting...', 'Name', 'Loading Progress');
waitbar(0, h, sprintf('Progress: %d%%', 0));

bits_p_sym_list = [4];
bits_p_sym_indice = 1;

angles = [0 1 2 3 4 5 10 20 30];
ppms = [0 2 10];
angle_index = 1;
nb_avg = 10;
point_BER = -1.5:0.5:16; % BER that will be explore to draw curves
nIterations  = length(point_BER);


for ppm_index = 1:length(ppms) 
    rate_avg = zeros(1,nIterations);
for avg = 1:nb_avg
    angle_value = angles(angle_index);
    ppm = ppms(ppm_index);
%_____________________adjustible variables________________________________
bits_p_sym = bits_p_sym_list(bits_p_sym_indice); % number of bits per symbol
nbr_sym_data = 2^15;
up_sample_factor = 20;
f_cut = 5E6; % Hz
taps =11;
beta = 0.3 ; % (0<=beta<=1)
fcarrier = 1e6;
k_gardner = 0.05;

phase_offset = angle_value/180*pi;

%______________________declaration of varibales_____________________________
f_sampling = f_cut*2*up_sample_factor;
T = 1/f_cut;
CFO_phase = 2*pi*ppm*fcarrier*1e-6;

%______________________generation of signal________________
vect = randi([0,1],bits_p_sym*nbr_sym_data,1); % random vector of bits

%______________________mapping_______________________________
if bits_p_sym ==1
    technique = "pam";
else
    technique = "qam";
end


map3 = mapping(vect,bits_p_sym,technique);

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


%______________________definition of error related vector___________________
rate_error_uncoded = ones(1,nIterations); % list of error rate depending on BER
rate_error_uncoded_cfo = ones(1,nIterations);

elt = (0:length(s_uncoded)-1)*1/f_sampling;
CFO = exp(1j*CFO_phase*elt+phase_offset);
CFO_neg = exp(-1j*CFO_phase*elt+phase_offset);

for k = 1:nIterations
    progress = (ppm_index-1)/length(ppms)+(avg-1)/nb_avg*1/length(ppms)+ 1/length(ppms)*1/nb_avg*k / nIterations  ;
    waitbar(progress, h, sprintf('Progress: %d%%', floor(progress * 100)));

    % We adjust the noise to make the curve
    N0_uncoded = Eb_uncoded / (10^(point_BER(k)/10));
    
    Y_uncoded = sqrt(N0_uncoded/2 * 50E6) .* (randn(1,length(s_uncoded)) + 1i * randn(1,length(s_uncoded)));  % Creating white gaussian noise
    t_uncoded = s_uncoded + Y_uncoded; % Signal with added noise
    
    t_uncoded_cfo = t_uncoded.*CFO;
    % We adjust the noise to make the curve

    % Vaut mieux utiliser randn, plus facile à controler
    % La puissance serait juste la variance de l'enveloppe divisé par 2 ?

    g_invert = flip(g); % flip function inverse the order of the vecteur g(-t). In fact g(-t) = g(t)
    % r_10_uncoded = conv(t_uncoded,g_invert,'same'); % convolution with g(-t). Caution with rate !
    % r_10_uncoded = r_10_uncoded / max(conv(g,g)); %normalization
    % r_uncoded = downsample(r_10_uncoded,up_sample_factor);
    % r_uncoded = r_uncoded';

    %cfo inverse filtering
    
    r_10_uncoded_cfo = conv(t_uncoded_cfo,g_invert,'same'); % convolution with g(-t). Caution with rate !
    r_10_uncoded_cfo = r_10_uncoded_cfo / max(conv(g,g)); %normalization
    %________________________garner algorithm correction _______________________
% we downsample the received signal by half of original sampling factor
    
    % half_downsampling_factor = up_sample_factor/2;
    % r_uncoded_half_downsampled=downsample(r_10_uncoded_cfo,half_downsampling_factor);
    % 
    % [r_uncoded_cfo,epsilon] = gardner(r_uncoded_half_downsampled,2,k_gardner);

    r_uncoded_cfo = downsample(r_10_uncoded_cfo,up_sample_factor);

    r_uncoded_cfo = r_uncoded_cfo';

%________________________demapping____________________________________________________
    if bits_p_sym == 1
        % r_uncoded = real(r_uncoded);
        r_uncoded_cfo = real(r_uncoded_cfo);
    end
    % demap_uncoded = demapping(r_uncoded,bits_p_sym,technique);
    demap_uncoded_cfo = demapping(r_uncoded_cfo,bits_p_sym,technique);

%____________________error computation_________________________________________________ 
    % We count the number of error
    % error_vct_uncoded = vect' - demap_uncoded';
    % nb_error =0;
    % for l = error_vct_uncoded
    %     if l ~= 0
    %         nb_error = nb_error + 1;
    %     end
    % end
    % 
    % rate_error_uncoded(k) = nb_error/length(error_vct_uncoded);

 %___________________error computation with garner_________________________
    error_vct_uncoded_cfo = vect' - demap_uncoded_cfo';
    nb_error =0;
    for l = error_vct_uncoded_cfo
        if l ~= 0
            nb_error = nb_error + 1;
        end
    end

    rate_error_uncoded_cfo(k) = nb_error/length(error_vct_uncoded_cfo);

    
end
    rate_avg = rate_avg+rate_error_uncoded_cfo;
end
rate_avg = rate_avg/nb_avg;

semilogy(point_BER,rate_avg)
    hold on

end



legend('ppm = 0','ppm = 2','ppm = 10')
title('BER with a CFO 16QAM')
ylabel('BER')
xlabel('Eb/N0')
hold off
close(h);

