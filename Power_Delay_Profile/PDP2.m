close all;
clear all;

fid =fopen(['Waveform.txt'],'r');
x = fscanf(fid, '%f\n');
fclose(fid);

% DTFT
f = [-0.5:0.01:0.5];
n = 0:length(x)-1;
W = -j*2*n*pi;
DTFTx = zeros(1,length(f));
for i = 1:length(x)
   DTFTx = DTFTx + x(i)*exp(W(i)*f);
end
figure;
plot(f,10*log10(abs(DTFTx)));
grid on;

% In-phase
% down conversion
IF = 70;
SAMPLING_RATE = 200;
n = [0:length(x)-1];
x_I = x'.*cos(2*pi*IF/SAMPLING_RATE*n);
[h,err,res] = firgr(80,[0 0.125 0.2 1],[1 1 0 0],[1,5]);
x_I_LPF = conv(x_I, h);

% Quadrature
IF = 70;
SAMPLING_RATE = 200;
n = [0:length(x)-1];
x_Q = x'.*sin(2*pi*IF/SAMPLING_RATE*n);
[h,err,res] = firgr(80,[0 0.125 0.2 1],[1 1 0 0],[1,5]);
x_Q_LPF = conv(x_Q, h);

% DTFT
f = [-0.5:0.01:0.5];
n = 0:length(x_I_LPF)-1;
W = -j*2*n*pi;
DTFTx_I_LPF = zeros(1,length(f));
for i = 1:length(x_I_LPF)
   DTFTx_I_LPF = DTFTx_I_LPF + x_I_LPF(i)*exp(W(i)*f);
end
figure;
plot(f,10*log10(abs(DTFTx_I_LPF)));
grid on;

% DTFT
f = [-0.5:0.01:0.5];
n = 0:length(x_Q_LPF)-1;
W = -j*2*n*pi;
DTFTx_Q_LPF = zeros(1,length(f));
for i = 1:length(x_Q_LPF)
   DTFTx_Q_LPF = DTFTx_Q_LPF + x_Q_LPF(i)*exp(W(i)*f);
end
figure;
plot(f,10*log10(abs(DTFTx_Q_LPF)));
grid on;

% correlatioon
MSEQ = mseq(2, 8);
NSAMP = 40;
MSEQ_up = zeros(1, NSAMP*length(MSEQ));
for i=1:length(MSEQ)
    MSEQ_up(NSAMP*(i-1)+1:NSAMP*i) = [MSEQ(i) zeros(1, NSAMP-1)];
end
Tx_filter = rcosfir(0, [-length(MSEQ) length(MSEQ)], NSAMP, 1,'sqrt')*sqrt(NSAMP);
x_cr = conv(MSEQ_up, Tx_filter);

z_I = conv(x_I_LPF, x_cr);
z_Q = conv(x_Q_LPF, x_cr);

figure;
plot(z_I);
figure;
plot(z_Q);

figure;
plot(10*log10((sqrt(z_I.^2 + z_Q.^2))));
grid on;


