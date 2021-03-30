%Band Pass filter -- Inverse Chebyshev
%Karamesinis Antonios-Rafail -- AEM 8430

clear
%% parameters
a1=8;
a2=4;
a3=3;
a4=0;

f0=1*1000;
f1=650 + 25*a4;
D= 2.1*(f0^2-f1^2)/f1;
f2= f0^2/f1;
f3= (-D+ sqrt(D^2+4*f0^2))/2;
f4=f0^2/f3;

amin= 35 - a3;
amax= 0.4 + a4/36;

wo = 2*pi*f0;
w1 = 2 * pi * f1;
w2= 2 * pi * f2;
w3= 2* pi * f3;
w4 = 2 * pi * f4;


Wp = 1;
Ws = ( w4 - w3 ) / ( w2 - w1 );

BW = 2 * pi * ( f2 - f1 );
qc = wo / BW;

n = acosh( sqrt( ( 10^( amin / 10 ) -1 ) / ( 10^( amax / 10 ) - 1 ) ) ) / acosh( Ws );
n = ceil(n);

e = 1 / sqrt( 10^( amin / 10 ) -1 );
a = ( 1 / n ) * ( asinh( 1 / e ) );

whp = 1 / cosh( acosh( 1 / e ) / n);

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
%%

%antistrofi polon
InvW1 = 1 / W1;
InvW23 = 1 / W23;
InvW45 = 1 / W45;

%klimakopoihsh sixnotitas
InvW1= InvW1 * Ws ;
InvW23= InvW23 * Ws;
InvW45= InvW45 * Ws;

S1 = InvW1 / ( 2 * Q1);
S23= InvW23 / ( 2 * Q23);
S45= InvW45 / ( 2 * Q45);
%%
W1 = sqrt( InvW1^2 - S1^2 );
W23 = sqrt( InvW23^2 - S23^2 );
W45 = sqrt( InvW45^2 - S45^2 );

%% midenika gia k=1 kai k=3 k=5
Z1 = sec(pi/10);
Z2 = sec(3*pi/10);
Z3 = sec(5*pi/10);

% klimakopoihsh mhdenikwn
Z1 = Z1 * Ws;
Z2 = Z2 * Ws;
Z3 = 0;

%% Transform 1st pole
Q1_Geffe = qc/S1;
y_1 = acosd(1/(2*Q1_Geffe));
p1_1 = wo*(-cosd(y_1)+(1i)*sin(y_1));
p1_2 = wo*(-cosd(y_1)-(1i)*sin(y_1));
w_01 = wo;

%% Transform 2nd  pole
p2_Inv = -S23 + ( W23 * 1i );
C2= S23^2 + W23^2;
D2= 2* S23 / qc;
E2= 4 + C2/ qc^2;
G2= sqrt (E2^2 - 4* D2^2);
Q2_Geffe= 1/D2 * sqrt ( 1/2* ( E2+ G2) );
k2= S23 * Q2_Geffe /qc;
W2= k2 + sqrt( k2^2 -1);
w_02 = 1/W2 * wo;
w_03 = W2 * wo;

%% Transform 3rd  pole
p3_Inv = -S45 + ( W45 * 1i );
C3= S45^2 + W45^2;
D3= 2* S45 / qc;
E3= 4 + C3/ qc^2;
G3= sqrt (E3^2 - 4* D3^2);
Q3_Geffe= 1/D3 * sqrt ( 1/2* ( E3+ G3) );
k3= S45 * Q3_Geffe /qc;
W3= k3 + sqrt( k3^2 -1);
w_04 = 1/W3 * wo;
w_05 = W3 * wo;

%% Transform Z1-Z2-Z3
K_Z1 = 2 + (Z1)^2/ (qc)^2;
x1 = ( K_Z1 + sqrt( K_Z1^2 - 4 ) ) / 2;
w_z2 = wo * ( sqrt(x1) );
w_z3 = wo / ( sqrt(x1) );

K_Z2 = 2 + (Z2)^2/ (qc)^2;
x2 = ( K_Z2 + sqrt( K_Z2^2 - 4 ) ) / 2;
w_z4 = wo * ( sqrt(x2) );
w_z5 = wo / ( sqrt(x2) );

w_z1 = wo;

%% Unit-1 Delyiannis-Fried Strategy 2 
Q1=Q1_Geffe;
C11=1;
C12=1;
R11=1/(2*Q1);
R12=2*(Q1);
H11=2*(Q1)^2;
%%
% Klimakopoihsh
kf1= wo;
Cnew = 0.01*10^(-6);
km1 = 1 / (Cnew*kf1);
C11 = Cnew;
C12 = Cnew;
R11 = R11*km1;
R12 = R12*km1;
k_unit1 = H11;

w0_unit1 = wo;
wz_unit1 = wo;


%% Unit-2 HPN
wz02 = (w_z3)/(w_04);
k2_1 = (w_04/w_z3)^2 - 1;
Q2 = Q3_Geffe;
k2_2=((Q2^2)*(k2_1+2))/((Q2^2)*(k2_1+2)+1);
k_unit2=k2_2*(w_04/w_z3)^2; 

R21 = 1;
R22 = (Q2^2)*(k2_1+2)^2;
R23 = 1;
R24 = (Q2^2)*(k2_1+2);
C2 = 1/(Q2*(k2_1+2));
C21= k2_1*C2;
C22 = 1/(Q2*(k2_1+2));

w0_unit2=w_04;
wz_unit2=w_z3;
%%
%Klimakopoihsh
Cnew = 0.01*10^(-6);
kf2 = w0_unit2;
km2 = (C2) / (Cnew*kf2);
R21 = R21*km2;
R22 = R22*km2;
R23 = R23*km2;
R24 = R24*km2;
C21 = C21/(km2*kf2);
C2 = Cnew;
C22 = C2;

%% Unit-3 LPN
wz03 = (w_z4)/(w_05);
Q3=Q3_Geffe;
R31=1;
R32=4*Q3^2;
R33 = (wz03)^2 / (2*Q3^2);
R34=1;
R35 = (4*Q3^2)/(wz03^2-1);
C3=1/(2*Q3);
k_unit3 = 1/(1+(wz03^2)/(2*Q3^2));

w0_unit3 = w_05;
wz_unit3= w_z4;

%Klimakopoihsh
Cnew = 0.01*10^(-6);
kf3 = w0_unit3;
km3 = C3/(Cnew*kf3);
R31 = R31*km3;
R32 = R32*km3;
R33 = R33*km3;
R34 = R34*km3;
R35 = R35*km3;
C3 = Cnew;

%% Unit-4 HPN
wz04 = (w_z5)/(w_02);
Q4 = Q2_Geffe;
k4_1 = (w_02/w_z5)^2 - 1;
k4_2=((Q4^2)*(k4_1+2))/((Q4^2)*(k4_1+2)+1);
k_unit4=k4_2*(w_02^2)/(w_z5^2); 
R41 = 1;
R42 = (Q4^2)*(k4_1+2)^2;
R43 = 1;
R44 = (Q4^2)*(k4_1+2);
C4 = 1/(Q4*(k4_1+2));
C41=k4_1*C4;
C42 = C4;

w0_unit4 = w_02;
wz_unit4 = w_z5;

%Klimakopoihsh
Cnew = 0.01*10^(-6);
kf4 = w0_unit4;
km4 = (C4) / (Cnew*kf4);
R41=1*km4;
R42 = R42*km4;
R43 = R43*km4;
R44 = R44*km4;
C41 = C41/(km4*kf4);
C4 = Cnew;
C42 = C4;

%% Unit-5 LPN
wz05 = (w_z2)/(w_03);
Q5=Q2_Geffe;
R51 = 1;
R52 = 4*Q5^2;
R53 = (wz05)^2 / (2*Q5^2);
R54 = 1;
R55 = (4*Q5^2)/((wz05^2)-1);
C5=1/(2*Q5);
k_unit5 = 1/(1+(wz05^2)/(2*Q5^2)); 

w0_unit5 = w_03;
wz_unit5 = w_z2;

%Klimakopoihsh
Cnew = 0.01*10^(-6);
kf5 = w0_unit5;
km5 = C5/(Cnew*kf5);
R51 = R51*km5;
R52 = R52*km5;
R53 = R53*km5;
R54 = R54*km5;
R55 = R55*km5;
C5 = Cnew;

%% Transfer functions 
T_1 = tf( [ -k_unit1*(w0_unit1)/(Q1) 0 ],[ 1 (w0_unit1/ Q1) w0_unit1^2 ] );
T_2 = tf( [ k_unit2 0 (k_unit2*wz_unit2^2) ],[ 1 (w0_unit2 / Q2) w0_unit2^2 ] );
T_3 = tf( [ k_unit3 0 (k_unit3*wz_unit3^2) ],[ 1 (w0_unit3 / Q3) w0_unit3^2 ] );
T_4 = tf( [ k_unit4 0 (k_unit4*wz_unit4^2) ],[ 1 (w0_unit4 / Q4) w0_unit4^2 ] );
T_5 = tf( [ k_unit5 0 (k_unit5*wz_unit5^2) ],[ 1 (w0_unit5 / Q5) w0_unit5^2 ] );

% Kerdos
K_total = k_unit1*k_unit2*k_unit3*k_unit4*k_unit5;

%T_total
T_total = T_1*T_2*T_3*T_4*T_5;

gain = abs(evalfr(T_total, wo * 1i));
a = 10^(0/20) / gain;

%a<1
Z1 = 10^4;
Z2 = Z1*a;

T_total = a * T_total;

%plot_transfer_function(T_1, [f1 f2 f3 f4])
%plot_transfer_function(T_2, [f1 f2 f3 f4]) 
%plot_transfer_function(T_3, [f1 f2 f3 f4]) 
%plot_transfer_function(T_4, [f1 f2 f3 f4]) 
%plot_transfer_function(T_5, [f1 f2 f3 f4]) 
%plot_transfer_function(T_total, [f1 f2 f3 f4])  
%plot_transfer_function(inv(T_total), [f1 f2 f3 f4]) 
%ltiview({'bodemag'}, T_1,T_2,T_3,T_4,T_5,T_total) 


%% Signals 

f11 = (wo - (wo-w1)/2) / (2*pi);
f12 = (wo + (wo+w1)/3) / (2*pi);
f13 = 0.4*w3 / (2*pi);
f14 = (2.5*w4) / (2*pi);
f15 = (3*w4) / (2*pi) ;

fs= 200*10^3;
T=0.01;
dt=1/fs;
t=0:dt:T;

f_in = cos(2*pi*f11*t) + 0.8*cos(2*pi*f12*t) + 0.8*cos(2*pi*f13*t) + 0.6*cos(2*pi*f14*t) + 0.5*cos(2*pi*f15*t);

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
% Input
Input_fft=fft(f_in);
L=length(f_in);
P2 = abs(Input_fft/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;
% Plot fft input
figure;
plot(f,P1);
axis([0.01 30000 0 inf]);
title('FFT of Input signal');

% Output
Output_fft=fft(f_out);
L=length(f_out);
P2 = abs(Output_fft/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;
% Plot fft output
figure;
plot(f,P1);
axis([0.01 30000 0 inf]);
title('FFT of Output signal');