clear all
clc
%% Original FCS
% \V = \{ u_{\min},...,u_{\max} \} \subset \Z
u_min=0;    %Minimum value for input
u_max=1;    %Maximum value for input

%% Multistep FCS
% u \in \V^n_u \subset \Z^n_u
% \U_k = \V^n_u Np \subset \Z^p; p=n_u Np

n_u = 3;             % Number of inputs u=[Sa Sb Sc]'
N1 = 1;
N = 2;              % Prediction Horizon length
N2 = N-N1+1;
Fs = 10000;          % Sampling Frequency
Ts = 1/Fs;           % Sample Time
p = n_u*N;          % Length of switching sequences U_k
Niter_max = 32000;   % Maximum number of explored nodes for SDA
freqHzSin = 50; %Frequência da senoide de referência em Hz


%% Basic Matrices
I_2 = eye(2);          % Identity matrices
I_3 = eye(3);
O_2 = zeros(2);        % Zero matrices
J = [0 -1;
   1  0];
Tr = sqrt(2/3)*[1     -0.5         -0.5;      % Clarke transform 
              0   sqrt(3)/2   -sqrt(3)/2];
                 

%% Continuous-Time alpha-beta Model -> 2 lvl VSI with output LC filter as example

n_y = 2;             % Number of outputs for control y= [vo_alpha vo_beta]
Lf = 2e-3;           % LC filter inductor
Cf = 50e-6;          % LC filter capacitor
Rl = 60;             % Resistive load

P = 0; %ERRO DE MODELAGEM NOS COMPONENTES
% Vnoise = 2; %RUÍDO NA MEDIDA DA SAÍDA
Vnoise = 0; %RUÍDO NA MEDIDA DA SAÍDA

Lsim = Lf*(1-P/100); % LC filter inductor in simulation
Csim = Cf*(1-P/100); % LC filter capacitor in simulation
Rsim = Rl*(1-P/100); % Resistive load in simulation
w = 2*pi*50;         % Angular frequency (@50 Hz)
Vdc = 400;                  % Dc link voltage
Vdc_sim = Vdc*1;            % Dc link voltage in simulation
Vo_ref_amp = 120*sqrt(2);   % Output voltage amplitude reference
lambda_u = 60.0;            % Weighting factor: Switching effort penalty

Rr=[1  -w*Ts;        % Rotation matrix
     w*Ts  1];

Ac_xy = [O_2        -(1/Lf)*I_2         O_2     ;
      (1/Cf)*I_2      O_2       -(1/Cf)*I_2   ;
       O_2            O_2             w*J     ;];

Bc_xy = [(Vdc/Lf)*I_2 O_2 O_2]'*Tr ;
        
Cxy = [I_2 O_2 O_2;
      O_2 I_2 O_2];   % C for observer
 
ct_sys = ss(Ac_xy,Bc_xy,Cxy,[]);

%% Discrete-Time alpha-beta Model
dt_sys = c2d(ct_sys,Ts);
A = dt_sys.a;
B = dt_sys.b;
C = dt_sys.c(3:4,:);


%% Multistep Matrices (H is lower triangular)
[UpsilonT,Gamma,lSTE,W_inv,H] = MPC_Matrices_l(A,B,C,N,lambda_u); %Multistep matrices

