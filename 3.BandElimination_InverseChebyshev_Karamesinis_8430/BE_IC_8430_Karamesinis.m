%% Band Elimination Inverse Chebyshev
% Karamesinis Antonios-Rafail
%AEM: 8430

clear
clc
%% set parameters
a_1 = 8;
a_2 = 4;
a_3 = 3;
a_4 = 0;

f_0 = 1.8*10^3;
f_1 = 1200 +25*(9-a_4);
f_2 = (f_0^2)/f_1;
D = (f_0^2 - f_1^2)/(1.8*f_1);
f_3 = (-D+sqrt(D^2+4*f_0^2))/2;
f_4 = (f_0^2)/f_3;

w_0 = 2*pi*f_0;
w_1 = 2*pi*f_1;
w_2 = 2*pi*f_2;
w_3 = 2*pi*f_3;
w_4 = 2*pi*f_4;

a_min = 30 - a_3;
a_max = 0.5 +a_4/18;

bw = w_2 - w_1;
q_c = w_0 / bw;

W_p = 1;
W_s = (w_2 - w_1)/(w_4-w_3);

%W_p = 1/W_s;
%W_s = 1;

n = acosh( sqrt( ( 10^( a_min / 10 ) -1 ) / ( 10^( a_max / 10 ) - 1 ) ) ) / acosh( W_s );
n = ceil(n);

e = 1 / sqrt( 10^( a_min / 10 ) -1 );
a = ( 1 / n ) * ( asinh( 1 / e ) );

whp = 1 / (cosh((1/n)*acosh(1/e)));

%n=5
y1 = 0;
y2 = 36;
y3 = -36;
y4 = 72;
y5 = -72;

%Poles
p1 = -sinh(a) * cosd(y1)+ (1i)*cosh(a) * sind(y1);
p2 = -sinh(a) * cosd(y2)+ (1i)*cosh(a) * sind(y2);
p3 = -sinh(a) * cosd(y3) + (1i)*cosh(a) * sind(y3);
p4 = -sinh(a) * cosd(y4) + (cosh(a) * sind (y4))*1i;
p5 = -sinh(a) * cosd(y5) + (cosh(a) * sind (y5))*1i;

W1  =   abs(p1);
W23 = abs(p2);
W45 = abs(p4);

Q1 = W1/(-2*real(p1));
Q23= W23 / ( -2 * real(p2) );
Q45 = W45/(-2*real(p4));

%antistrofi polon
InvW1 = 1 / W1;
InvW23 = 1 / W23;
InvW45 = 1 / W45;

%klimakopoihsh sixnotitas
InvW1= InvW1 * W_s ;
InvW23= InvW23 * W_s;
InvW45= InvW45 * W_s;

%% midenika gia k=1 kai k=3 k=5
Z1 = sec(pi/10);
Z2 = sec(3*pi/10);
Z3 = sec(5*pi/10);

% klimakopoihsh mhdenikwn
Z1 = Z1 * W_s;
Z2 = Z2 * W_s;
Z3 = 0;

%Antistrofi polwn
S1=1/InvW1;
InvW23= 1/ InvW23;
InvW45= 1/InvW45;
S23 = -InvW23 / (2* Q23);
W23 = sqrt( InvW23^2 - S23^2 );
S45 = -InvW45 / (2* Q45);
W45 = sqrt( InvW45^2 - S45^2 );
S1=-S1;

%Antistrofi midenikon
InvZ1=1/Z1;
InvZ2=1/Z2;
InvZ3=0;

%% Algorithm Geffe

%Metasxhmatismos pragmatikou polou
w_01=w_0;
Q1= q_c /(-S1);

%Metasxhmatismos 2ou migadikou polou
p23Inv = S23 + ( W23 * 1i );
C2= S23^2 + W23^2;
D2= -2* S23 / q_c;
E2= 4 + C2/ q_c^2;
G2= sqrt (E2^2 - 4* D2^2);
Q2_3= 1/D2 * sqrt ( 1/2* ( E2+ G2) );
k2= -S23 * Q2_3 /q_c;
W2= k2 + sqrt( k2^2 -1);
w_02 = 1/W2 * w_0;
w_03 = W2 * w_0;

%Metasxhmatismos 3ou migadikou polou
p45Inv = S45 + ( W45 * 1i );
C3= S45^2 + W45^2;
D3= -2* S45 / q_c;
E3= 4 + C3/ q_c^2;
G3= sqrt (E3^2 - 4* D3^2);
Q4_5= 1/D3 * sqrt ( 1/2* ( E3+ G3) );
k2= -S45 * Q4_5 /q_c;
W3= k2 + sqrt( k2^2 -1);
w_04 = 1/W3 * w_0;
w_05 = W3 * w_0;

%Metasxhmatismos 1ou midenikou
K_zero1 = 2 + (InvZ1^2) / (q_c^2);
x1 = ( K_zero1 + sqrt( K_zero1^2 - 4 ) ) / 2;
wz2 = w_0 * ( sqrt(x1) );
wz3 = w_0 / ( sqrt(x1) );

%Metasxhmatismos 2ou midenikou
K_zero2 = 2 + (InvZ2^2) / (q_c^2);
x2 = ( K_zero2 + sqrt( K_zero2^2 - 4 ) ) / 2;
wz4 = w_0 * ( sqrt(x2) );
wz5 = w_0 / ( sqrt(x2) );

wz1 = w_0;

%% Monades

% Monada1-Notch
k11 = 0;
R11=1;
R13=1;
R12= Q1^2* (k11+2)^2;
R14= (k11 +2 )* Q1^2 ;
C= 1/ (Q1 * ( 2 + k11) );
k12= ((k11 +2) * Q1^2)/((k11+2) *Q1^2 +1);

%Klimakopoihsh 1uF
kf1=w_01;
km1= C/ (kf1 * 10^(-6));
C11= 10^(-6);
C12=C11;
R11=R11* km1;
R12= R12 * km1;
R13= R13 * km1;
R14= R14 * km1;
k_unit1=k12;

% Monada 2 - LPN
wz02= wz2/ w_02;
R21=1;
R22= 4 * (Q2_3)^2;
R23= wz02^2 / ( 2 * (Q2_3)^2) ;
R24= 1;
R25= 4 * (Q2_3)^2 / ( wz02^2 -1);
C21= 1/ (2 * Q2_3);
k_unit2 = 1 / ( 1 +  wz02^2 / ( 2 * Q2_3^2) );
C22=C21;

%Klimakopoihsh
kf2= w_02;
km2= C21/ (kf2 * 10^ (-6));
R21=R21* km2;
R22= R22 * km2;
R23= R23 * km2;
R24= R24 * km2;
R25 = R25 * km2;
C21=10^(-6);
C22=C21;

% Monada 3 - HPN
wz03 = (wz3)/(w_03);
k3_1 = (w_03/wz3)^2 - 1;
k3_2=((Q2_3^2)*(k3_1+2))/((Q2_3^2)*(k3_1+2)+1);
k_unit3=k3_2*(w_03/wz3)^2; 
R31 = 1;
R32 = (Q2_3^2)*(k3_1+2)^2;
R33 = 1;
R34 = (Q2_3^2)*(k3_1+2);
C31 = 1/(Q2_3*(k3_1+2));
C32= k3_1*C31;
C33 = C31;

%Klimakopoihsh
kf3=w_03;
km3= C31/ (kf3 * 10^ (-6));
R31=R31*km3;
R32=R32*km3;
R33=R33*km3;
R34=R34*km3;
C31=10^ (-6);
C32=C32/(km3*kf3);
C33=C31;

% Monada 4 - LPN
wz04= wz4/ w_04;
R41=1;
R42= 4 * (Q4_5)^2;
R43= wz04^2 / ( 2 * (Q4_5)^2) ;
R44= 1;
R45= 4 * (Q4_5)^2 / ( wz04^2 -1);
C41= 1/ (2 * Q4_5);
k_unit4 = 1 / ( 1 +  wz04^2 / ( 2 * Q4_5^2) );
C42=C41;

%Klimakopoihsh
kf4= w_04;
km4= C41/ (kf4 * 10^ (-6));
R41=R41* km4;
R42= R42 * km4;
R43= R43 * km4;
R44= R44 * km4;
R45 = R45 * km4;
C41=10^(-6);
C42=C41;

% Monada 5 - HPN
wz05 = (wz5)/(w_05);
k5_1 = (w_05/wz5)^2 - 1;
k5_2=((Q4_5^2)*(k5_1+2))/((Q4_5^2)*(k5_1+2)+1);
k_unit5=k5_2*(w_05/wz5)^2; 
R51 = 1;
R52 = (Q4_5^2)*(k5_1+2)^2;
R53 = 1;
R54 = (Q4_5^2)*(k5_1+2);
C51 = 1/(Q4_5*(k5_1+2));
C52= k5_1*C51;
C53 = C51;

%Klimakopoihsh
kf5=w_05;
km5= C51/ (kf5 * 10^ (-6));
R51=R51*km5;
R52=R52*km5;
R53=R53*km5;
R54=R54*km5;
C51=10^ (-6);
C52=C52/(km5*kf5);
C53=C51;

%% Transfer functions
T_1 = tf( [k_unit1 0 ( k_unit1 * w_01^2 ) ], [ 1 ( w_01 / Q1 ) w_01^2 ] );
T_2 = tf( [k_unit2 0 ( k_unit2 * wz2^2 ) ], [ 1 ( w_02 / Q2_3 ) w_02^2 ] );
T_3 = tf( [k_unit3 0 ( k_unit3 * wz3^2 ) ], [ 1 ( w_03 / Q2_3 ) w_03^2 ] );
T_4 = tf( [k_unit4 0 ( k_unit4 * wz4^2 ) ], [ 1 ( w_04 / Q4_5 ) w_04^2 ] );
T_5 = tf( [k_unit5 0 ( k_unit5 * wz5^2 ) ], [ 1 ( w_05 / Q4_5 ) w_05^2 ] );

%T_total
T_total = T_1*T_2*T_3*T_4*T_5;

%Rythmisi Kerdous
K_total = k_unit1*k_unit2*k_unit3*k_unit4*k_unit5;
a = 10^(10/20) / K_total;

Z1=10^4;
Z2=a*Z1;

T_total = a * T_total;


%plot_transfer_function(T_1, [f_1 f_2 f_3 f_4])
%plot_transfer_function(T_2, [f_1 f_2 f_3 f_4]) 
%plot_transfer_function(T_3, [f_1 f_2 f_3 f_4]) 
%plot_transfer_function(T_4, [f_1 f_2 f_3 f_4]) 
%plot_transfer_function(T_5, [f_1 f_2 f_3 f_4]) 
%plot_transfer_function(T_total, [f_1 f_2 f_3 f_4])  
%plot_transfer_function(inv(T_total), [f_1 f_2 f_3 f_4]) 
%ltiview({'bodemag'}, T_1,T_2,T_3,T_4,T_5,T_total) 

%% Signals 

f11 = (w_0 - (w_0-w_3)/2) / (2*pi);
f12 = (w_0 + (w_0+w_3)/2) / (2*pi);
f13 = 0.5*w_1 / (2*pi);
f14 = (2.4*w_2) / (2*pi);
f15 = (3.5*w_2) / (2*pi) ;

fs= 200*10^3;
T=0.002;
dt=1/fs;
t=0:dt:(T);

f_in = 0.8*cos(2*pi*f11*t) + 1*cos(2*pi*f12*t) + cos(2*pi*f13*t) + 0.8*cos(2*pi*f14*t) + 0.4*cos(2*pi*f15*t);

figure
plot(t,f_in)
title('Input') 
xlabel('t (sec)') 
ylabel('Amplitude') 

f_out=lsim(T_total,f_in,t);
figure
plot(t,f_out)
title('Output') 

%% Fourier Analysis

%Input Signal
input_fft=fft(f_in);
L=length(f_in);
P2 = abs(input_fft/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;

figure;
plot(f,P1);
axis([0.01 10000 0 inf]);
title('FFT of Input signal');

%Output Signal
output_fft=fft(f_out);
L=length(f_out);

P2 = abs(output_fft/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;
figure;
plot(f,P1);
axis([0.01 10000 0 inf]);
title('FFT of Output signal');