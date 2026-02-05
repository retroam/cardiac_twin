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
 clc; clear all;%close all ;
           
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
load KLcalc.mat; load fitalpha.mat;load Ki;
KL = KLcalc(1);
alpha_L = fitalpha(1);
KA = KLcalc(17);
alpha_A = fitalpha(17);

options = odeset('RelTol',1e-5,'MaxStep',5e-3,'Stats','on');
p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
 y0 = zeros(29,1); y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
t = [0;1000*60*1000]; p(1:2) = 0;[~,y] = ode15s(@daeODE,t,y0,[],p);y0Sig = y(end,:);
 p = daePARAMSc(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
 load('y0_1Hz.mat');y0(1:29) = y0Sig;
% [~,y] = ode15s(@daeODEc,[0;30*60],y0,options,p); y0 = y(end,:);
p(17) = 0.59;p(23) = 0.025;

ISO1 =0.1 ; ISO2 =10;BB =0.1 ;
tspan = [0;4*60];
p(1) =ISO1; p(2) = BB;  % agonist conc and beta blocker conc
[t,y] = ode15s(@daeODEc,tspan,y0,options,p);
yfinal = y(end,:); 
tspan2 = [4*60 8*60]; 
p(1) = ISO2; p(2) = BB;     % agonist conc and beta blocker conc
[t2,y2]=ode15s(@daeODEc,tspan2,yfinal,options,p);
yfinal2 = y2(end,:);
t = [t;t2]./60;y = [y;y2];
yCell=mat2cell(y,size(y,1),ones(size(y,2),1));
[Ri,G, b1AR_S464,b1AR_S301,...
    GsaGTPtot, GsaGDP, Gsby, AC_GsaGTP ,cAMPtot, PDEp, ...
  RC_I, RCcAMP_I, RCcAMPcAMP_I, RcAMPcAMP_I, PKACI ,PKACI_PKI, ...
  RC_II, RCcAMP_II, RCcAMPcAMP_II, RcAMPcAMP_II, PKACII, PKACII_PKI, ...
  I1p_PP1, I1ptot, LCCap ,LCCbp, PLBp, PLMp, TnIp, ...
    m, h, jo,v, w, x, yo, z,rto, sto ,ssto,rss ,...
    sss,Ca_nsr ,Ca_jsr,Nai ,Ki ,Cai,Vm,trelo]=yCell{:};

subplot(3,3,1);plot(t,rss);title('rss');axis tight;hold all;
subplot(3,3,2);plot(t,sss);title('sss');axis tight;hold all;
subplot(3,3,3);plot(t,Ca_nsr);title('Ca_{nsr}');axis tight;hold all;
subplot(3,3,4);plot(t,Ca_jsr);title('Ca_{jsr}');axis tight;hold all;
subplot(3,3,5);plot(t,Nai);title('Nai');axis tight;hold all;
subplot(3,3,6);plot(t,Ki);title('Ki');axis tight;hold all;
subplot(3,3,7);plot(t,Cai);title('[Ca]');axis tight;hold all;
subplot(3,3,8);plot(t,Vm);title('Vm');axis tight;hold all;
subplot(3,3,9);plot(t, trelo);title('trel');axis tight;hold all;


