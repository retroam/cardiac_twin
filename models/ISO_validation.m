% Run file for rat signaling model
% 
%
% Copyright 2011, Cardiac Systems Biology Lab, University of Virginia
%   JS: Jeff Saucerman  <jsaucerman@virginia.edu>
%   JY: Jason Yang      <jhyang@virginia.edu>
%   RA: Robert Amanfu   <rka2p@virginia.edu>
%
% Robert Amanfu
% 11/08/11
%% Receptor parameters to optimize
 clc; clear all;%close all;
           
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




%% cAMP validation with Ca2+ experiments
% load KLcalc.mat; load fitalpha.mat;
% 
% 
% KL =KLcalc(1);
% alpha_L = fitalpha(1);% 0.3037
% KA = KLcalc(20);%2e-3;
% alpha_A =  fitalpha(20);
%  p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);p(1) = 0;
%  y0 = zeros(29,1);y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
%  [t,y] = ode15s(@daeODE,[0;60*60*1000],y0,[],p);y0 = y(end,:);
% tspan = [0;10*60*1000];  % 10 minutes
% 
% 
% 
% PROrange = 10.^[-5:.1:2.5];
% ISO = 0.1;
% p(1:2) = 0;
% for i=1:length(PROrange)
%     disp(i);
%     p(2) = PROrange(i);  % Atot 
%     p(1) = ISO;
%     [t,y] = ode15s(@daeODE,tspan,y0,[],p);
%     cAMP(i) = (y(end,9));
%          
% end
% 
% figure;semilogx(PROrange*1e-6,cAMP);hold all;
% save ISOPRO_cAMP cAMP;
%% Ca2+ validation 
load KLcalc.mat; load fitalpha.mat;
KG =  2.4131;gamma_L =  0.3762;
KL =[KLcalc(1) 0.5 0.8];%KLcalc(1);
alpha_L = fitalpha(1);% 0.3037
KA = KLcalc(20);%2e-3;
alpha_A =  fitalpha(20);
 p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);p(1) = 0;
 y0 = zeros(29,1);y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
 [~,y] = ode15s(@daeODE,[0;60*60*1000],y0,[],p);y0Sig = y(end,:);
tspan = [0;4*60];  % 10 minutes

options = odeset('RelTol',1e-5,'MaxStep',5e-3,'Stats','on');
p = daePARAMSc(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
load y0full;y0 = y0full;


ISOrange = 10.^[-5:.2:1.0];

p(1:2) = 0;

for h = 1
for i=1:length(ISOrange)
    y0 = y0full;p(1:2) = 0;
    disp(i);
  
    p = daePARAMSc(KR,KL(h),KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
   [t1,y1] = ode15s(@daeODEc,[0;2*60],y0,options,p);y0 = y1(end,:);
    
    p(2) = 0;  % Atot 
    p(1) = ISOrange(i);
   [t2,y2] = ode15s(@daeODEc,[2*60;6*60],y0,options,p);
   y = [y1;y2];t = [t1;t2];
   cAMP(h,i) = (y(end,9));
   ISO_Ca{h,i} = y(:,47);
       
end
end

% figure;semilogx(PROrange*1e-6,cAMP);hold all;
% figure;semilogx(PROrange*1e-6,calcium);
save ISO_cAMPII cAMP;
save ISO_CaII ISO_Ca;
% toc
