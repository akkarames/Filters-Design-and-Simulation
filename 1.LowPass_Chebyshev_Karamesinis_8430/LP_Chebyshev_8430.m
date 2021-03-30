%% Low Pass - Chebyshev
%%Karamesinis Antonios-Rafail

clear

%% AEM
a_1 = 8;
a_2 = 4;
a_3 = 3;
a_4 = 0;

%% Prodiagrafes
m = 2;
fp = 1.1*(3 + m) * 10^3;
fs = 1.6*fp ;

a_min = 21.5 + [max(1,a_3)-5]*0.5;
a_max = 0.55 + [max(1,a_4)-5]/10;

w_p = 2*pi*fp;
w_s = 2*pi*fs;

%%Kanonikopoihsh
Wp = 1;
Ws = w_s/w_p;

%% Taksi n kai e 
e = sqrt(10^(a_max/10)-1);

n = acosh(((10^(a_min/10)-1)/(10^(a_max/10)-1))^(1/2))/(acosh(Ws));
n= ceil(n);

a = (1/n)*asinh(1/e);
w_hp = cosh((1/n)*acosh(1/e));
f_hp = w_hp * fp;

%%butterworth gwnies
y_k1 = 0;
y_k2 = 36;
y_k3 = -36;
y_k4 = 72;
y_k5 = -72;

%%poloi
p_k1 = -sinh(a);
W_01 = (sqrt(real(p_k1)^2 + imag(p_k1)^2));
Q1 = 0.5;

p_k2 = -sinh(a)*cosd(y_k2) + 1i*cosh(a)*sind(y_k2);
p_k3 = -sinh(a)*cosd(y_k2) - 1i*cosh(a)*sind(y_k2);
W_023 = (sqrt(real(p_k2)^2 + imag(p_k2)^2));
Q23 = 1.177;

p_k4 = -sinh(a)*cosd(y_k4) + 1i*cosh(a)*sind(y_k4);
p_k5 = -sinh(a)*cosd(y_k4) - 1i*cosh(a)*sind(y_k4);
W_045 = (sqrt(real(p_k4)^2 + imag(p_k4)^2));
Q45 = 4.54;


%the real value of the poles is:
w_01 = W_01 * 2 * pi * fp;
w_023 = W_023 * 2 * pi * fp;
w_045 = W_045 * 2 * pi * fp;

%% Stratigiki 1 - ison piknoton ison antistaseon

% Prwti Monada
R11 = 1;
C11 = 1; 
kf1 = w_01;
C11 = 1 * 10^(-6);
km1 = 1 / ( C11 * kf1 );
R11 = R11 * km1;

% Deuteri Monada
C21 = 1;
C22 = 1;
R21 = 1;
R22 = 1;
C21 = 1 * 10^(-6);
C22 = C21;
k2 = 3 - (1/Q23);
r21 =1;
r22 = 2 - (1/Q23);
C21 = 1 * 10^(-6);
C22 = C21;
kf2 = w_023;
km2 = 1 / (kf2 * C22);
R21 = R21 * km2;
R22 = R22 * km2;
r21= r21*km2;
r22= r22*km2;

% Triti Monada
C31 = 1;
C32 = 1;
R31 = 1;
R32 = 1;
C31 = 1 * 10^(-6);
C32 = C21;
k3 = 3 - (1/Q45);
r31 =1;
r32 = 2 - (1/Q45);
C31 = 1 * 10^(-6);
C32 = C31;
kf3 = w_045;
km3 = 1 / (kf3 * C32);
R31 = R31 * km3;
R32 = R32 * km3;
r31= r31*km3;
r32= r32*km3;

%% Rythmisi kerdous
% Theloume gain 5 db ara 20logaK = 5
Ktot= k2 * k3 ;
a_gain = 10^(1/4) / Ktot;

Z2 = R11/a_gain;
Z3 = R11/(1-a_gain);

%% Transfer functions

T1 = tf ([0 w_01], [1 w_01])
T2 = tf ([ 0 0 k2*(w_023^2)], [1 w_023/Q23 w_023^2])
T3 = tf ([ 0 0 k3*(w_045^2)], [1 w_045/Q45 w_045^2])

T_total = a_gain * T1*T2*T3

%plot_transfer_function(T1, [fp, fs])

%plot_transfer_function(T2, [fp, fs])

%plot_transfer_function(T3, [fp, fs])

plot_transfer_function(T_total, [fp, fs])

%ltiview({'bodemag'}, T1, T2, T3, T_total)
 
InvSys_new = inv (T_total)
%plot_transfer_function(InvSys_new, [fp fs])


%% Fourier
% eisodos palmoi 2kHz

x = [ones(1, 20), zeros(1, 60), ones(1,20)];
y=x;
m = 10;
for i=2:m
    y = [y x];
end

t = linspace(0, m*0.5*10^(-3), m*100);

figure, plot(t, y)

z = lsim(T_total,y,t);
hold on, plot(t,z)

L = length(y);
fs = 200000;
Y = fft(y);
Y1 = abs(Y/L);
p = Y1(1:L/2+1);
p(2:end-1)=2*p(2:end-1);
f = fs*(0:(L/2))/L;
figure, plot(f,p)


Z = fft(z);
Z1 = abs(Z/L);
p = Z1(1:L/2+1);
p(2:end-1)=2*p(2:end-1);
f = fs*(0:(L/2))/L;
figure, plot(f,p)
