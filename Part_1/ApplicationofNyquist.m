vect = randi([0,1],64,1);
map = mapping(vect,4,'qam');

ma2 = [];
for l = map'
    ma2 = [ma2,0,0,0,0,l,0,0,0,0,0];
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
g = ifftshift(g); % for the signal in positive timeline

% plot(g);

s= conv(ma2,g,'same');
% plot(real(s));


t = s;

r_10 = conv(t,g,'same'); % convolution with g(-t)
r_10 = r_10 / max(conv(g,g));

r=[];
count = 6;
for i = r_10
    if rem(count, 10) == 0
        r = [r,i];
    end
    count = count + 1;
end
r = r';

% demapping
demap = demapping(r,4,'qam');

tiledlayout(5,1)
nexttile
plot(conv(g,g))
title("Nyquist filter")
nexttile
plot(real(map))
title("Original signal after mapping")

nexttile
plot(real(r))
title("Signal received before demapping")
nexttile
plot(vect)
title("Original signal before mapping")
nexttile
plot(demap)
title("Signal received after demapping")

figure
G = fftshift(G);
plot(f, abs(G));
title('Limited communication bandwidth of the nyquiest filter');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;

