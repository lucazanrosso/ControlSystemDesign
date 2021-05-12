%%  Nominal model definitions

%   mechanical parameters
Jm  = 0.000026;         %   motor rotor inertia
Jd1 = 1.1e-5;          %   small disk inertia
Jd2 = 0.00015;         %   large disk inertia
J   = Jm+Jd1+Jd2;           %   tot load inertia
b   = 3e-4;      %   viscous damping coeff


%   electrical parameters
Ki  = 2;                % current loop transconductance [A/V]
Kt  = 0.071;            % motor torque const [Nm/A]

% Ipotizzo comando in corrente - devo dividerlo per 2 prima di applicarlo
% al modello del motore

A=[0 1
    0 -b/J]
B=[0 Kt/J]'
C=[1 0]
D=0

sys1=ss(A,B,C,D)
sys1d=c2d(sys1,0.001,'zoh')
[F,G,H,DD]=ssdata(sys1d)
alpha=-b/J;
beta=Kt/J;
T=0.001;
WT=[T^3/3*beta^2 T^2/2*beta^2+T^3/3*alpha*beta^2
    T^2/2*beta^2+T^3/3*alpha*beta^2 T*beta^2+alpha^2*beta^2*T^3/3+2*alpha*beta^2*T^2/2]

A=[0 1
    0 -b/J]
Q=(10/2^12)^2/12*WT*1000




R=(2*pi/2000)^2/12


L=(dlqr(F',H',Q,R))'

Fes=F-L*H;

disp('run your KF model')

pause

yykf=(cumsum(abs(fft(dtrend(innovation(1000:8000,2))))));
plot(yykf(1:3500)/yykf(3500))
hold on
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