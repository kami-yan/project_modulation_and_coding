
rate_error = []; % list of error rate depending on BER
rate_error_cfo = [];
rate_error_phase_shift = [];

f_carrier = 2E9; % Hz
T = 200E-9; % (duration of a symbol)
phase_error = 2;

delta_f = 10*10^-6; % Here delta_f = 10 ppm (10^-5)




nb_bits = 2^12; % Number of bits send 2^6=64
bits_p_sym = 4; % number of bits per symbol, so the number of symbol is 512/8 = 64
vect = randi([0,1],nb_bits,1); % random vector of bits

map = mapping(vect,bits_p_sym,'qam');
ma2= [];
for i = map'
    ma2 = [ma2,0,0,0,0,i,0,0,0,0,0];
end
% demap = demapping(map,4,'qam');
% 64 bit per symbol
% 5 mega symbol per second

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

V = var(s) / 2; % variance of the signal


Eb = V * T  / bits_p_sym; % Calculation of a bit energy
%Eb = trapz(abs(s).^2)*T / bits_p_sym;

nb_point_BER = -1.5:0.5:13; % BER that will be explore to draw curves
for k = nb_point_BER 
    % We adjust the noise to make the curve
    N0 = Eb / (10^(k/10));

    Y = sqrt(N0/2 * 50E6) .* (randn(1,length(s)) + 1i * randn(1,length(s)));  % Creating white gaussian noise
    t = s + Y; % Signal with added noise
    
    
    t_cfo = [];
    t_phase_err = [];
    count = 0; 
    for m = t
        t_cfo = [t_cfo,m*exp(1i*2*pi*f_carrier*delta_f*T*count)];   % What changed : count/10 => count
        t_phase_err = [t_phase_err,m*exp(1i*phase_error)];
        count = count +1;
    end




    g_invert = conj(flip(g)); % flip function inverse the order of the vecteur g(-t). In fact g(-t) = g(t)
    r_10 = conv(t,g_invert,'same'); % convolution with g(-t). Caution with rate !
    r_10 = r_10 / max(conv(g,g)); %normalization

    r_10_cfo = conv(t_cfo,g_invert,'same'); % convolution with g(-t). Caution with rate !
    r_10_cfo = r_10_cfo / max(conv(g,g)); %normalization

    
    count = 0; 
    for m = 1:length(r_10_cfo)
        r_10_cfo(m) = r_10_cfo(m)*exp(-1i*2*pi*f_carrier*delta_f*T*count);    % We want only the ISI, not the phase shift
        count = count +1;
    end

    r_10_phase = conv(t_phase_err,g_invert,'same'); % convolution with g(-t). Caution with rate !
    r_10_phase = r_10_phase / max(conv(g,g)); %normalization

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

    % To go back to the intial sampling rate, we take one symbol over_cfo 10
    r_cfo=[];

    count = 6; % The fourth is the first to be sampled
    for l = r_10_cfo
        if rem(count, 10) == 0   
            r_cfo = [r_cfo,l];

        end
        count = count + 1;
    end

    r_cfo = r_cfo';

    % To go back to the intial sampling rate, we take one symbol over_cfo 10
    r_phase_err = [];

    count = 6; % The fourth is the first to be sampled
    for l = r_10_phase
        if rem(count, 10) == 0
            r_phase_err = [r_phase_err,l];

        end
        count = count + 1;
    end

    r_phase_err = r_phase_err';
    
    

    demap = demapping(r,bits_p_sym,'qam');
    demap_cfo = demapping(r_cfo,bits_p_sym,'qam');
    demap_phase_shift = demapping(r_phase_err,bits_p_sym,'qam');
    

    % We count the number of error
    error_vct = vect - demap;
    error_vct_cfo = vect - demap_cfo;
    error_vct_phase_shift = vect - demap_phase_shift;


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

    

    rate_error = [rate_error,nb_error/length(error_vct)];
    rate_error_cfo = [rate_error_cfo,nb_error_cfo/length(error_vct_cfo)];
    rate_error_phase_shift = [rate_error_phase_shift,nb_error_phase_shift/length(error_vct_phase_shift)];
 

end

semilogy(nb_point_BER,rate_error)

grid on
hold on 
semilogy(nb_point_BER,rate_error_cfo)
semilogy(nb_point_BER,rate_error_phase_shift)
legend('no CFO', 'CFO', 'Phase shift')
hold off


