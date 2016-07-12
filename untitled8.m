
N = 83;
B = 1000;
T = [0:2*pi/5:2*pi-0.1];

Z = exp(1j*T(unidrnd(numel(T),[N,B])));
mag = abs(mean(Z));
ang = angle(mean(Z));

%%

