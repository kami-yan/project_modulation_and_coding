function  rate_avg = error_calculation(N,K_window,ToA,with_cfo,nbavg,cfo_calcul)

%_____________________adjustible variables_________________________________

bits_p_sym = 4; % number of bits per symbol
nbr_sym_pilot = N;
nb_train = 2;


point_BER = -1.5:0.5:16;
point_BER = [0,2,4,6,8,10,12,14,16];

rate_avg= zeros(nbavg,length(point_BER));
mean = zeros(1,length(point_BER));


for avg = 1:nbavg

%______________________definition of error related vector___________________
rate_error = []; % list of error rate depending on BER
rate_error_cfo = [];
rate_error_phase_shift = [];


%deuxieme partie

rate = [];

%______________________declaration of varibales_____________________________
bits_p_sym = 2; % number of bits per symbol
nbr_sym_data = 2^15;
up_sample_factor = 10;
delta_w = 10^-5;
f_cut = 5E6; % Hz
taps =101;
beta = 0.3 ; % (0<=beta<=1)
f_carrier = 600e6;
phase_error =0;

nbr_sym_data=108;

%______________________declaration of varibales_____________________________
f_sampling = f_cut*up_sample_factor;
T = 1/f_cut;

pilot_size = bits_p_sym*nbr_sym_pilot; % it should only be multiple of bits_sym


nb_bits = (nbr_sym_data*bits_p_sym+pilot_size); 




%_____________________adaptation des grandeurs_______________________

bits_p_sym_with_pilot =   bits_p_sym;





%______________________generation of signal with pilot________________


pilot = randi([0,1],pilot_size,1);

unused = randi([0,1],(ToA-1)*bits_p_sym,1); 
vect = randi([0,1],bits_p_sym*nbr_sym_data*nb_train,1); % random vector of bits


piloted_vect = [unused;pilot;vect];

%______________________mapping_______________________________

pilot_mapped = mapping(pilot,bits_p_sym,'qam');

map = mapping(piloted_vect,bits_p_sym_with_pilot,'qam');
ma2= upsample(map,up_sample_factor);

df = (1/taps)*f_sampling;    
fmax = df*(taps-1)/2;

f = linspace(-fmax,fmax,taps);
G = zeros(1,taps);
for l = 1:taps
    G(l) = HalfrootNyquistFilter(T,beta,f(l));
end

G_rc = ifftshift(G.^2);
G = sqrt(G_rc);       
g_rc = ifft(G_rc);
g = ifft(G);
g = fftshift(g/sqrt(g_rc(1)));


s= conv(ma2',g,'same');

V = var(s) / 2; % variance of the signal
Eb = V * T / bits_p_sym; % Calculation of a bit energy

 % BER that will be explore to draw curves
nIterations  = length(point_BER);


for k = 1:nIterations
    % We adjust the noise to make the curve

    % Vaut mieux utiliser randn, plus facile à controler
    % La puissance serait juste la variance de l'enveloppe divisé par 2 ?

    N0 = Eb / (10^(point_BER(k)/10));
    Y = sqrt(N0/2 * 50E6) .* (randn(1,length(s)) + 1i * randn(1,length(s)));  % Creating white gaussian noise
    t = s + Y; % Signal with added noise

    
    t_cfo = [];
    t_phase_err = [];
    count = 7; 
     for m = t
         t_cfo = [t_cfo,m*exp(1i*f_carrier*delta_w*T*count/10)];   
         t_phase_err = [t_phase_err,m*exp(1i*phase_error)];
    end

    
    g_invert = flip(g); % flip function inverse the order of the vecteur g(-t). In fact g(-t) = g(t)
    r_10 = conv(t,g_invert,'same'); % convolution with g(-t). Caution with rate !
    r_10 = r_10 / max(conv(g,g)); %normalization

    r_10_cfo = conv(t_cfo,g_invert,'same'); % convolution with g(-t). Caution with rate !
    r_10_cfo = r_10_cfo / max(conv(g,g)); %normalization

    r_10_phase = conv(t_phase_err,g_invert,'same'); % convolution with g(-t). Caution with rate !
    r_10_phase = r_10_phase / max(conv(g,g)); %normalization

    % To go back to the intial sampling rate, we take one symbol over 10
    r=downsample(r_10,up_sample_factor);
    r_cfo =downsample(r_10_cfo,up_sample_factor);
    r_phase_err =downsample(r_10_phase,up_sample_factor);
   
    r = r';
    r_cfo = r_cfo';
    r_phase_err = r_phase_err';

    % scatter(real(map),imag(map))
    % hold on
    % scatter(real(r_cfo),imag(r_cfo))
    % hold off
%______________________error correction (differential algo)_____________________________________________
if with_cfo == true
    received_signal = r_cfo';
else
    received_signal = r';
end


[n_estimate,estimate_delta_f,phi_0_estimate,max_value,total_vector] = different_algo_function (pilot_mapped',received_signal,T,K_window,nbr_sym_data); %K ne doit pas depasser la taille du pilot

r_without_cfo = zeros(1,length(received_signal));
for k1 = 1:length(r)
   r_without_cfo(k1) = received_signal(k1)*exp(-1j*2*pi*estimate_delta_f*k1*T);
end

r_correct = r_without_cfo(n_estimate:length(r));

% final signal r_correct still with pilots

%real  delta_f

real_delta_f = angle(received_signal(1)-received_signal(2))/(2*pi*T); 

%________________________demapping____________________________________________________

    demap = demapping(r,bits_p_sym_with_pilot,'qam');
    demap_cfo = demapping(r_cfo,bits_p_sym_with_pilot,'qam');
    demap_phase_shift = demapping(r_phase_err,bits_p_sym_with_pilot,'qam');



%____________________error computation_________________________________________________ 
    
    if cfo_calcul 
        error = abs((10^-6)*(estimate_delta_f));
    else
        error = abs(n_estimate-ToA);
    end

    % We count the number of error
    error_vct = piloted_vect - demap;
    error_vct_cfo = piloted_vect - demap_cfo;
    error_vct_phase_shift = piloted_vect - demap_phase_shift;


    nb_error = 0;

    for l = error_vct'
        if l ~= 0
            nb_error = nb_error + 1;
        end
    end

    nb_error_cfo = 0;
    for l = error_vct_cfo'
        if l ~= 0
            nb_error_cfo = nb_error_cfo + 1;
        end
    end

    nb_error_phase_shift = 0;
    for l = error_vct_phase_shift'
        if l ~= 0
            nb_error_phase_shift = nb_error_phase_shift + 1;
        end
    end


    rate = [rate, error];
    rate_error = [rate_error,nb_error/length(error_vct)];
    rate_error_cfo = [rate_error_cfo,nb_error_cfo/length(error_vct_cfo)];
    rate_error_phase_shift = [rate_error_phase_shift,nb_error_phase_shift/length(error_vct_phase_shift)];

%partie apres l'algorithme




end
mean = mean + rate;

rate_avg(avg,:) = rate;
end

mean = mean/nbavg;

for i = 1:nbavg
    rate_avg(i,:) = rate_avg(i,:)-mean;
end


rate_avg = sqrt(1/(nbavg-1)*sum((rate_avg)).^2);
end