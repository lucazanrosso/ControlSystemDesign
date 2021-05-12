%%  Nominal model definitions

%   mechanical parameters

%   mechanical parameters
Jm  = 0.000026;         %   motor rotor inertia
Jd1 = 1.1e-5;          %   small disk inertia
Jd2 = 0.00015;         %   large disk inertia
J   = Jm+Jd1+Jd2;           %   tot load inertia
b   = 3e-4;      %   viscous damping coeff
%   electrical parameters
Ki  = 2;                % current loop transconductance [A/V]
Kt  = 0.071;            % motor torque const [Nm/A]
A=[0 1 0
    0 -b/J -Kt/J
    0 0 0]
B=[0 Kt/J 0]'
C=[1 0 0]
D=0


% sys1 = con B nominale
% sys2 = con B nullo

%     Tra mettere B nullo o no, non cambia nulla sul
%     test di bartlett. Guardando i disturbo stimato 
%     si nota che la corrente di disturbo Id è 
%     maggiore con B nullo poiché devo compensare di
%     più rispetto al caso nominale.

sys1=ss(A,B,C,D)
b=0;
A=[0 1 0
    0 -b/J -Kt/J
    0 0 0];
sys2=ss(A,B,C,D)  


% sys1 = con B nominale
% sys2 = con B nullo

sys1d=c2d(sys2,0.001,'zoh')
[F,G,H,DD]=ssdata(sys1d)
alpha=-b/J;
beta= Kt/J;
T=0.001;


% Calcolo i graminani per entrambi i rumori

WT1=[T^3/3*beta^2 T^2/2*beta^2+T^3/3*alpha*beta^2 0
    T^2/2*beta^2+T^3/3*alpha*beta^2 T*beta^2+alpha^2*beta^2*T^3/3+2*alpha*beta^2*T^2/2 0
    0 0 0]
WT2=[0 0 0
    0 beta^2*T^3/3 -beta*T^2/2
    0 -beta*T^2/2 T]



%     Inizialmente simulo il modello preoccupandomi di WT1 (Lo faccio 
%     utilizando KF_playwith), taro quel peso li affinchè sia una retta,
%     poi passo a WT2 usando lo schema 3x3 (KF3_playwithreal), ma non
%     uso lo schema LQG, perché inizialmente mi preoccupo di stimare
%     i rumori interni del sistema. Purtroppo bisogna portare il segnale
%     ai livelli più alti di intensità compatibilmente con la saturazione
%     (sotto la corrente max) e poi rendere l'innovazione bianca. Così 
%     stimato tutto. Ci teniamo Q e R ben stretti, ora li usiamo nella
%     progettazione del controllore in LGQ_toplaywith

Q=(10/2^12)^2/12*WT1*1000+WT2*0.01

R=(2*pi/2000)^2/12

L=(dlqr(F',H',Q,R))'
Fes=F-L*H;
disp('run your KF model')

pause

yykf=(cumsum(abs(fft(detrend(innovation(1000:8000,2))))));

%%
figure(1)
% f = 1:3500;
yy = abs(fft(detrend(innovation(1000:8000,2))));
f = (0:length(yy)-1)*1000/length(yy);
% plot(f(1:3500), 0.1*ones(1,3500))
hold on
plot(f(1:3500), yy(1:3500))
hold off

%%
figure(2)
% f = 1:3500;
f = (0:length(yykf)-1)*1000/length(yykf);
plot(f(1:3500),f(1:3500)/500)
hold on
plot(f(1:3500), yykf(1:3500)/yykf(3500));
% hold off
% A=[0 1 0
%     0 -b/J -Kt*Ki/J
%     0 0 0]
% B=[0 Kt*Ki/J 0]'
% C=[1 0 0]
% D=0
% 
% sys1=ss(A,B,C,D)
% sys1d=c2d(sys1,0.001,'zoh')
% [F,G,H,DD]=ssdata(sys1d)
% Q=eye(3)
% R=1
% L=(dlqr(F',H',Q,R))'