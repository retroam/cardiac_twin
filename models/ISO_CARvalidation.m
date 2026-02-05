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
% CARrange = 10.^[-5:.1:2.5];
% ISO = 0.1;
% p(1:2) = 0;
% for i=1:length(CARrange)
%     disp(i);
%     p(2) = CARrange(i);  % Atot 
%     p(1) = ISO;
%     [t,y] = ode15s(@daeODE,tspan,y0,[],p);
%     cAMP(i) = (y(end,9));
%          
% end
% 
% figure;semilogx(CARrange*1e-6,cAMP);hold all;
% save ISOPRO_cAMP cAMP;
% %% ISO KL vs CAR
% load KLcalc.mat; load fitalpha.mat;load Ki;
% KG =  2.4131;gamma_L =  0.3762;
% p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
% RelTol = 1e-13;MaxStep = 1e-1;
% options = odeset('MaxStep',1e3,'NonNegative',[1:29],'RelTol',RelTol);
%  y0 = zeros(29,1); y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
% t = [0;1000*60*1000]; p(1:2) = 0;[~,y] = ode15s(@daeODE,t,y0,options,p);y0 = y(end,:);
% 
% dose = [0.1 1 1];drugs = [20 18 17];
% KLrange = 0.1:.03:1;
% for i=1:length(drugs)
%       for j = 1:length(KLrange)
%   KL = KLrange(j);%KLcalc(1);
%    alpha_L = fitalpha(1);
%    KA = KLcalc(drugs(i)); alpha_A = fitalpha(drugs(i));
%    p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
%     
%     tspan = [0;2*60*1000];
% p(1) = 0; p(2) = 0;  % agonist conc and beta blocker conc
% [~,y] = ode15s(@daeODE,tspan,y0,options,p);
% yfinal = y(end,:); 
% 
% tspan = [2*60*1000;6*60*1000];
% p(1) = 0.1; p(2) = dose(i);  % agonist conc and beta blocker conc
% [~,y2] = ode15s(@daeODE,tspan,yfinal,options,p); 
% yfinal = y2(end,:); 
% 
% tspan2 = [6*60*1000 10*60*1000]; %changed the time
% p(1) = 10; p(2) = dose(i);     % agonist conc and beta blocker conc
% [~,y3]=ode15s(@daeODE,tspan2,yfinal,options,p);
% response(i,j) = max(y3(:,9));
%       end
% 
% end
% 
% plot(KLrange,response);


%% CAR alpha vs CAR
load KLcalc.mat; load fitalpha.mat;load Ki;
% KG =  2.4131;gamma_L =  0.3762;
% p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
% RelTol = 1e-13;MaxStep = 1e-1;
% options = odeset('MaxStep',1e3,'NonNegative',[1:29],'RelTol',RelTol);
%  y0 = zeros(29,1); y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
% t = [0;1000*60*1000]; p(1:2) = 0;[~,y] = ode15s(@daeODE,t,y0,options,p);y0 = y(end,:);
% 
% dose = [0.1 1 1];drugs = [20 18 17];
% alpharange = 0.1:.03:1;
% for i=1:length(drugs)
%       for j = 1:length(KLrange)
%   KL = 0.3;%KLcalc(1);
%    alpha_L = fitalpha(1);
%    KA = KLcalc(drugs(i)); alpha_A = fitalpha(drugs(i));
%    p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
%     
%     tspan = [0;2*60*1000];
% p(1) = 0; p(2) = 0;  % agonist conc and beta blocker conc
% [~,y] = ode15s(@daeODE,tspan,y0,options,p);
% yfinal = y(end,:); 
% 
% tspan = [2*60*1000;6*60*1000];
% p(1) = 0.1; p(2) = dose(i);  % agonist conc and beta blocker conc
% [~,y2] = ode15s(@daeODE,tspan,yfinal,options,p); 
% yfinal = y2(end,:); 
% 
% tspan2 = [6*60*1000 10*60*1000]; %changed the time
% p(1) = 10; p(2) = dose(i);     % agonist conc and beta blocker conc
% [~,y3]=ode15s(@daeODE,tspan2,yfinal,options,p);
% response(i,j) = max(y3(:,9));
%       end
% 
% end
% 
% plot(KLrange,response);
%% Ca2+ validation 
load KLcalc.mat; load fitalpha.mat;
KG =  2.4131;gamma_L =  0.3762;
KL =0.5;%KLcalc(1);
alpha_L = fitalpha(1);% 0.3037
KA = KLcalc(17);%2e-3;
alpha_A =  fitalpha(17);
 p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);p(1) = 0;
 y0 = zeros(29,1);y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
 [~,y] = ode15s(@daeODE,[0;60*60*1000],y0,[],p);y0Sig = y(end,:);
tspan = [0;4*60];  % 10 minutes

options = odeset('RelTol',1e-5,'MaxStep',5e-3,'Stats','on');
p = daePARAMSc(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
load y0full;y0 = y0full;

CARrange = 10.^[-5:.25:2.5];
ISO = 0.1;
p(1:2) = 0;
for i=1:length(CARrange)
    disp(i);
    y0 = y0full;p(1:2) = 0;
    [t1,y1] = ode15s(@daeODEc,[0;2*60],y0,options,p);y0 = y1(end,:);
    p(2) = CARrange(i);  % Atot 
    p(1) = ISO;
   [t2,y2] = ode15s(@daeODEc,[2*60;6*60],y0,options,p);
   y = [y1;y2];t = [t1;t2];
    cAMP(i) = (y(end,9));
    ISOCAR_Ca{i} = y(:,47);
         
end

% figure;semilogx(CARrange*1e-6,cAMP);hold all;
% figure;semilogx(CARrange*1e-6,calcium);
save ISOCAR_cAMP cAMP;
save ISOCAR_Ca ISOCAR_Ca;

