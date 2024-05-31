clearvars; close all;
h = waitbar(0, 'Starting...', 'Name', 'Loading Progress');
waitbar(0, h, sprintf('Progress: %d%%', 0));

bits_p_sym_list = [4];
bits_p_sym_indice = 1;

%_____________________adjustible variables________________________________
bits_p_sym = bits_p_sym_list(bits_p_sym_indice); % number of bits per symbol
nbr_sym_data = 2^9;
up_sample_factor = 20;
f_cut = 5E6; % Hz
taps = 11;
beta = 0.3 ; % (0<=beta<=1)
fcarrier = 1e6;


k_gardners = [ 0.02 ];
angles = [0 1 2 3 4 5 10 20 30];
ppms = [0 2 10 1000];
time_shifts = [ 0 1 2 3 7];

angle_index = 1;
ppm_index = 1;
nb_avg = 100;
point_BER = 5; % BER that will be explore to draw curves
nIterations  = length(point_BER);
time_shift_index = 5;


for error_weight_index = 1:length(k_gardners)
    epsilon_list = zeros(nb_avg,nbr_sym_data);
    rate_avg = zeros(1,nIterations);
for avg = 1:nb_avg
    angle_value = angles(angle_index);
    ppm = ppms(ppm_index);
    time_shift = time_shifts(time_shift_index);
    error_weight = k_gardners(error_weight_index);




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

df = (1/taps)*f_sampling;    
fmax = df*(taps-1)/2;

f = linspace(-fmax,fmax,taps);
G = zeros(1,taps);
for l = 1:taps
    G(l) = HalfrootNyquistFilter(T,beta,f(l));
end

% G = ifftshift(G); % for the symetric
% g = ifft(G);
% g = fftshift(g); % for the signal in positive timeline

G_rc = ifftshift(G.^2);
G = sqrt(G_rc);       
g_rc = ifft(G_rc);
g = ifft(G);
g = fftshift(g/sqrt(g_rc(1)));
 

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
    progress = (error_weight_index-1)/length(k_gardners)+(avg-1)/nb_avg*1/length(k_gardners)+ 1/length(k_gardners)*1/nb_avg*k / nIterations  ;
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
    r_10_uncoded_cfo = circshift(r_10_uncoded_cfo,time_shift);

    %________________________garner algorithm correction _______________________
% we downsample the received signal by half of original sampling factor
    
     if time_shift  == 0
          r_uncoded_cfo = downsample(r_10_uncoded_cfo,up_sample_factor);
    else
    half_downsampling_factor = up_sample_factor/2;
    r_uncoded_half_downsampled=downsample(r_10_uncoded_cfo,half_downsampling_factor);

    [r_uncoded_cfo,epsilon] = gardner(r_uncoded_half_downsampled,2,error_weight);
    epsilon_list(avg,:) = epsilon; 
    end

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

mean_epsilon = mean(epsilon_list)+time_shift/up_sample_factor; % mean of time error plus an offset so that the value is positive
std_epsilon = std(epsilon_list);

color_special = ['r','g','Black'];

plot(mean_epsilon,'Color',color_special(error_weight_index),'LineStyle','-')
hold on
plot(mean_epsilon+std_epsilon,'Color','b','LineStyle','--')
plot(mean_epsilon-std_epsilon,'Color','b','LineStyle','--')
end

h1 = plot([1 1.1], [0.5 0.55], 'r'); % Red line
h2 = plot([1 1.1], [0.5 0.55], 'g'); % Green line
h3 = plot([1 1.1], [0.5 0.55], 'black'); % Black line
% text = {sprintf('error weight = %d ',k_gardners(1)),sprintf('error weight = %d ',k_gardners(2)),sprintf('error weight = %d ',k_gardners(3))};
% legend([h1,h2,h3],text)
legend('ppm = 7 and error weight = 0.02')
title('time error with cfo 16QAM')
axis([0 nbr_sym_data 0 0.4])%legend will generate a unwanted plot so i will limit it
xlabel('symbols')
ylabel('time error mean +/- std')
hold off
close(h);

