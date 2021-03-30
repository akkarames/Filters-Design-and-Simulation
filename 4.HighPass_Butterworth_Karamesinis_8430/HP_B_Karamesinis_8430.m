%% High Pass - Butterworth
%%Karamesinis Antonios-Rafail

%% AEM
a_1 = 8;
a_2 = 4;
a_3 = 3;
a_4 = 0;

%% Prodiagrafes
m = 0;
fp = (4 + m) * 10^3;
fs = fp / 2.6;

a_min = 24  + a_3*(6/9);
a_max = 0.5 + a_4/36;

w_p = 2*pi*fp;
w_s = 2*pi*fs;


Wp = 1;
Ws = w_p/w_s;

%% sintelestes a, e, taksi filtrou
e = sqrt(10^(a_max/10)-1);
n = (log10( (10^(a_min/10)-1) / (10^(a_max/10)-1) ) / (2*log10(Ws/Wp)) );

n = ceil(n);

%% sixnotita hmiseias isxios

Wo = Wp / (10^(a_max/10) - 1)^(1/(2*n));
wo = w_p / Wo;
fo = wo / (2*pi);

%% gwnies butterworth
y_k1 = 0;
y_k2 = +36;
y_k3 = -36;
y_k4 = +72;
y_k5 = -72;

%% poloi
p_1 = -1;
p_2 = -0.809 + j*0.5877;
p_3 = -0.809 - j*0.5877;
p_4 = -0.309 + j*0.951;
p_5 = -0.309 - j*0.951;


Q_1 = 0.5;

Q_2 = 0.618;

Q_3 = 1.618;

%--------------------------
% Sallen Key Strategy (1)
% Capacitance 1uF
% 10 dB

%% Monada 1

C11 = 1; 
R11 = 1;
kf1 = wo;
p1 = 1;
C11 = 1 * 10^(-6);
km1 = 1 / (kf1 * C11);
R11 = R11 * km1;

%% Monada 2

C21 = 1;
C22 = 1;
R21 = 1;
R22 = 1;
C21 = 1 * 10^(-6);
C22 = C21;
k2 = 3 - (1/Q_2);
kf2 = wo;
km2 = 1 / (kf2 * C22);
R21 = R21 * km2;
R22 = R22 * km2;
r21 =1;
r22 = 2 - (1/Q_2);
r21= r21*km2;
r22= r22*km2;

%% Monada 3

C31 = 1;
C32 = 1;
R31 = 1;
R32 = 1;
k3 = 3-(1/Q_3);
kf3 = wo;
C31 = 1 * 10^(-6);
C32 = C31;
km3 = 1 / (kf3 * C32);
R31 = km3 * R31;
R32 = km3 * R32;
r31 =1;
r32 = 2 - (1/Q_3);
r31= r31*km3;
r32= r32*km3;

%% Ktot
Ktot = k2 * k3 ;


%%Rythmisi Kerdous
a= 10^(1/2)/Ktot;
Z2=1;
Z3= (a*Z2)/(1-a);

%%Logika tha anevoyme klimaka sta kOhm gia na exoyme diaforetiki taksi ston
%%diaireti tasis

%% Transfer functions

T1 = tf ([1 0], [1 wo])
T2 = tf ([ 1 0 0], [1 wo/Q_2 wo^2])
T3 = tf ([ 1 0 0], [1 wo/Q_3 wo^2])

%T_test=T1*T2
T_total = T1*T2*T3 ;

%plot_transfer_function(T1, [fp, fs])

%plot_transfer_function(T2, [fp, fs])

%plot_transfer_function(T3, [fp, fs])

%plot_transfer_function(T_test, [fp, fs])

%plot_transfer_function(T_total, [fp, fs])

ltiview({'bodemag'}, T1, T2, T3, T_total)
 
InvSys_new = inv (T_total)
%plot_transfer_function(InvSys_new, [fp fs])

%% fourier analysis

f11= (0.2*w_s) / (2*pi);
f12= (0.7*w_s) / (2*pi);
f13= (1.6*w_p) / (2*pi);
f14= (2.4*w_p) / (2*pi);
f15= (3.5*w_p) / (2*pi) ;
fs= 200*10^3;
T=0.002;
dt=1/fs;
t=0:dt:(T);

u1= cos(2*pi*f11*t)+0.6*cos(2*pi*f12*t)+1.5*cos(2*pi*f13*t)+0.7*cos(2*pi*f14*t)+0.4*cos(2*pi*f15*t);
figure
plot(u1)

N=T/dt;
figure
lsim(T_total,u1,t)
xt=lsim(T_total,u1,t);
figure
plot(t,xt)
n=2^nextpow2(N);
xfourier= fft(xt,n);
p2=abs(xfourier/n);
p1=p2(1:n/2+1);
p1(2:end-1)=2*p1(2:end-1);
f=200*10^3*(0:(n/2))/n;
figure
plot(f,p1)

nfft=n;
y=fft(u1,nfft);
y = y(1:nfft/2); 
y_mag=abs(y);
f = (0:nfft/2-1)*fs/nfft; 
figure
plot(f,y_mag)

