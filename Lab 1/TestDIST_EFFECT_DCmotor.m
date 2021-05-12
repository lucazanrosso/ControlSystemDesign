Jm=0.000026;
J1=1.1e-5;

J2=0.00015;
Jt=0.00019456;
Bm=0.000256 ;
Kt=0.071;
A=[0 1
    0 -Bm/Jt];
B=[0; Kt/Jt];
C=[1 0];

Ts=0.001
sysC=ss(A,B,C,0);
sysD=c2d(sysC,Ts,'zoh');
[F,G,H,D]=ssdata(sysD);

%Primo caso: Stimatore costruito su Motore ideale senza attrito e parametri di coppia e inerzia
%reali
As=[0 1
    0 0];
%C=[1 0];
Cs=eye(2);
Ts=0.001
sysCs=ss(A,B,Cs,0);

% Costruisco il modello discreto
sysDs=c2d(sysCs,Ts,'zoh');

[Fs,Gs,Hs,Ds]=ssdata(sysDs);

% Progetto il controllore per il sistema ideale
Q=[1 0;0 0];
R=0.01;
KK=dlqr(Fs,Gs,Q,R)

lambda = eig(Fs-Gs*KK)
% lambda = lambda(1)^4

% Progetto lo stimatore per il sistema 2x2

L=(acker(Fs',[1 0]',[0.84 0.84]))'
Festar=Fs-L*[1 0]
Gestar=[Gs L]
% Lanciare la simulazione con dlqr_est e verificare gli errori
% dello stimatore 


%pause % premere enter per continuare


% Valuto con le equazioni viste a lezione i valori delle stime a regime
DIST=.1/Kt
%XX2=L(1)/L(2)*DIST-Ts/2*DIST
%XX1=(DIST+KK(2)*XX2)/(-KK(1))
XX2=DIST*(L(1)/L(2)*Gs(2) - Gs(1))/Fs(1,2)
X1diff=-DIST*Gs(2)/L(2) 

% Creo il sistema aumentato per lo stimatore con modello
% di ordine zero per il disturbo, che agisce all'ingresso
% Metto attrito nullo
% Usare il secondo schema simulink dlqr_est_2

% Bm=0;

A1=[0 1 0
    0 -Bm/Jt -Kt/Jt
    0 0 0]
B1= [0 Kt/Jt 0]';
C1=eye(3);
sysC1=ss(A1,B1,C1,0);
sysD1=c2d(sysC1,Ts,'zoh');
[F1,G1,H1,D1]=ssdata(sysD1)

L1=(acker(F1',[1 0 0]',[0.84 0.84 0.84]))'
Festar1=F1-L1*[1 0 0]
Gestar1=[G1 L1]
% Adesso posso lanciare la simulazione... e verifico che 
% le stime sono senza errori

%%%%%%