clearvars; close all;

h = waitbar(0, 'Starting...', 'Name', 'Loading Progress');
waitbar(0, h, sprintf('Progress: %d%%', 0));

ToA = 34;
nbavg = 10;
nb_point_BER_log = [0,2,4,6,8,10,12,14,16];
point_BER = -1.5:0.5:16;
count =12;

%calcul of ToA without cfo
cfo_calcul = false;
with_cfo = false;

rate_N10_K8  = error_calculation(10,8,ToA,with_cfo,nbavg,cfo_calcul);
waitbar( 1/count, h, sprintf('Progress: %d%%', round(1/count * 100)));

rate_N20_K8  = error_calculation(20,8,ToA,with_cfo,nbavg,cfo_calcul);
waitbar( 2/count, h, sprintf('Progress: %d%%', round(2/count * 100)));
rate_N40_K8  = error_calculation(40,8,ToA,with_cfo,nbavg,cfo_calcul);
waitbar( 3/count, h, sprintf('Progress: %d%%', round(3/count * 100)));
rate_N20_K1  = error_calculation(20,1,ToA,with_cfo,nbavg,cfo_calcul);
waitbar( 4/count, h, sprintf('Progress: %d%%', round(4/count * 100)));
rate_N20_K16 = error_calculation(20,16,ToA,with_cfo,nbavg,cfo_calcul);
waitbar( 5/count, h, sprintf('Progress: %d%%', round(5/count * 100)));

with_cfo = true;

rate_cfo_ToA  = error_calculation(40,8,ToA,with_cfo,nbavg,cfo_calcul);
waitbar( 6/count, h, sprintf('Progress: %d%%', round(6/count * 100)));

%calcul of cfo estimation error
cfo_calcul = true;
with_cfo = false;

rate_N10_K8_ppm  = error_calculation(10,8,ToA,with_cfo,nbavg,cfo_calcul);
waitbar( 7/count, h, sprintf('Progress: %d%%', round(7/count * 100)));
rate_N20_K8_ppm  = error_calculation(20,8,ToA,with_cfo,nbavg,cfo_calcul);
waitbar( 8/count, h, sprintf('Progress: %d%%', round(8/count * 100)));
rate_N40_K8_ppm  = error_calculation(40,8,ToA,with_cfo,nbavg,cfo_calcul);
waitbar( 9/count, h, sprintf('Progress: %d%%', round(9/count * 100)));
rate_N20_K1_ppm  = error_calculation(20,1,ToA,with_cfo,nbavg,cfo_calcul);
waitbar( 10/count, h, sprintf('Progress: %d%%', round(10/count * 100)));
rate_N20_K16_ppm = error_calculation(20,16,ToA,with_cfo,nbavg,cfo_calcul);
waitbar( 11/count, h, sprintf('Progress: %d%%', round(11/count * 100)));


with_cfo = true;

rate_cfo_ppm  = error_calculation(40,8,ToA,with_cfo,nbavg,cfo_calcul);
waitbar( 12/count, h, sprintf('Progress: %d%%', round(12/count * 100)));
close(h)


tiledlayout(2,3)

nexttile
plot(nb_point_BER_log,rate_N10_K8)
hold on
plot(nb_point_BER_log,rate_N20_K8)
hold on
plot(nb_point_BER_log,rate_N40_K8)
hold off
grid on
legend('N=10 K=8','N=20 K=8','N=40 K=8')

nexttile
plot(nb_point_BER_log,rate_N20_K1)
hold on
plot(nb_point_BER_log,rate_N20_K8)
hold on
plot(nb_point_BER_log,rate_N20_K16)
hold off
grid on
legend('N=20 K=1','N=20 K=8','N=20 K=16')

nexttile
plot(nb_point_BER_log,rate_cfo_ToA)
hold on
plot(nb_point_BER_log,rate_N40_K8)
hold off
grid on
legend('CFO ToA','without CFO')

nexttile
plot(nb_point_BER_log,rate_N10_K8_ppm)
hold on
plot(nb_point_BER_log,rate_N20_K8_ppm)
hold on
plot(nb_point_BER_log,rate_N40_K8_ppm)
hold off
grid on
legend('CFO N=10 K=8 ppm','CFO N=20 K=8 ppm','CFO N=40 K=8 ppm')

nexttile
plot(nb_point_BER_log,rate_N20_K1_ppm)
hold on
plot(nb_point_BER_log,rate_N20_K8_ppm)
hold on
plot(nb_point_BER_log,rate_N20_K16_ppm)
hold off
grid on
legend('CFO N=20 K=1 ppm','CFO N=20 K=8 ppm','CFO N=20 K=16 ppm')


nexttile
plot(nb_point_BER_log,rate_cfo_ppm)
hold on
plot(nb_point_BER_log,rate_N40_K8_ppm)

hold off
grid on
legend('with cfo ppm','without cfo ppm')



