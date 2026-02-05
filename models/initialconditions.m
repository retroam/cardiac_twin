KL =0.2;%15.2%.271;
KR =  10;%10 [uM] 
KG = 0.9041;%.75[uM] 
alpha_L = 1/32;% 1/32
gamma_L = 1.3; %1
KA = 500e-6;%500e-6 2e-3;
alpha_A = 1; %1 1.5574
gamma_A =  1; %1 
PKAItot         = 0.59;           % (uM) total type 1 PKA
PKAIItot        = 0.025;          % (uM) total type 2 PKA
load KLcalc.mat; load fitalpha.mat;
KG =  2.4131;gamma_L =  0.3762;

%  p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);p(1) = 0;
%  y0 = zeros(29,1);y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
%  [~,y] = ode15s(@daeODE,[0;60*60*1000],y0,[],p);y0Sig = y(end,:);
 options = odeset('RelTol',1e-5,'MaxStep',5e-3,'Stats','on');
 p = daePARAMSc(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);p(1) = 0;
%  load('y0_1Hz.mat');y0(1:29) = y0Sig;
load y0full;
 [t,y] = ode15s(@daeODEc,[0;10*60],y0full,options,p);plot(t,y(:,47));
%  y0full = y(end,:);save y0full y0full;