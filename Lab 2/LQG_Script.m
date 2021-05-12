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
A=[0 1 0
    0 -b/J -Kt/J
    0 0 0]
An=[0 1 0
    0 -0/J -Kt/J
    0 0 0]
B=[0 Kt/J 0]'
C=[1 0 0]
D=0

sys1=ss(A,B,C,D)
sys11=ss(An,B,C,D)
sys1d=c2d(sys11,0.001,'zoh')
[F,G,H,DD]=ssdata(sys1d)
alpha=-0/J;
beta=Kt/J;
T=0.001;
WT1=[T^3/3*beta^2 T^2/2*beta^2+T^3/3*alpha*beta^2 0
    T^2/2*beta^2+T^3/3*alpha*beta^2 T*beta^2+alpha^2*beta^2*T^3/3+2*alpha*beta^2*T^2/2 0
    0 0 0]
WT2=[0 0 0
    0 beta^2*T^3/3 -beta*T^2/2
    0 -beta*T^2/2 T]

Q=(10/2^12)^2/12*WT1*1000+WT2*100

R=(2*pi/2000)^2/12

% Computing the Estimator gains
L=(dlqr(F',H',Q,R))'

%Trick for real-time implementation
Fes=F-L*H;

%computing the feedback gains for the reachable subsystem
K=dlqr(F(1:2,1:2),G(1:2),[1 0;0 .1],0.01)

%%
% yykf=(cumsum(abs(fft(detrend(innovation(1000:8000,2))))));
% f = (0:length(yykf)-1)*1000/length(yykf);
% figure(1)
% plot(f(1:3500),f(1:3500)/500)
% hold on
% plot(f(1:3500), yykf(1:3500)/yykf(3500));



%K=acker(F(1:2,1:2),G(1:2),[.98 .95])
% z=tf('z',0.001)
% WW=H(1:2)*inv(z*eye(2)-F(1:2,1:2)+G(1:2)*K)*G(1:2)*K(1);
% WW=minreal(WW,.1);
% [n,d]=tfdata(WW,'verbose')
%     n=n(2:3)
    






















