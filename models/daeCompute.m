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


%% Checking receptor species
%  params = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
%  options = odeset('RelTol',1e-7);
% % % params(34) = 6.26e-04;params(33) =1.064e-08;
% % %  params(49) = 0; 
% % %   params(31:34) =  1/8*params(31:34);     % No desensitization or PDE's
% % %   params(35:37) = 0; 
% % 
%  tspan = [0;200*60*1000];
%  y0 = zeros(29,1);%p(17) = 0.59;p(23) = 0.025;
% [~,y] = ode15s(@daeODE,tspan,y0,options,params);
% y0 = y(end,:);
% % %
% % % save ICs y0;
% % % plot(t,y);
% % %  load ICs;
% % 
% params(1) = 0;params(2) = 0;
% tspan = [0;10*60*1000];
% [t,y] = ode15s(@daeODE,tspan,y0,options,params);
% yfinal = y(end,:); tspan2 = [10*60*1000 20*60*1000]; 
% params(1) = 1000;params(2) = 0;
% [t2,y2]=ode15s(@daeODE,tspan2,yfinal,options,params);
% % yfinal2 = y2(end,:);tspan3 = [20*60*1000 30*60*1000]; 
% % params(1) = 0;params(2) = 0;
% % [t3,y3]=ode15s(@daeODE,tspan3,yfinal2,options,params);
% t = [t;t2;]./60e3;
% y = [y;y2;];
% 
% % Plot Outputs
% yCell=mat2cell(y,size(y,1),ones(size(y,2),1));
% [Ra, LRa ,LRi, RaG ,LRaG ,ARi ,ARa, ARaG, b1AR_S464 ,b1AR_S301,...
%     GsaGTPtot, GsaGDP, Gsby, AC_GsaGTP, cAMPtot, PDEp, ...
%   RC_I, RCcAMP_I, RCcAMPcAMP_I, RcAMPcAMP_I, PKACI, PKACI_PKI, ...
%   RC_II ,RCcAMP_II, RCcAMPcAMP_II, RcAMPcAMP_II, PKACII, PKACII_PKI, ...
%   I1p_PP1, I1ptot, LCCap, LCCbp, PLBp, PLMp, TnIp] = yCell{:};
% Ri = 0.0132- Ra - LRi - LRa - RaG - LRaG - ARi - ARa - ARaG;
% Rtot = Ri + Ra + LRi + LRa + RaG + LRaG + ARi + ARa + ARaG;
% figure(2)
% subplot(3,3,1);plot(t,Ra);title('Ra');axis tight;hold all;
% subplot(3,3,2);plot(t,Rtot);title('Rtot');axis tight;hold all;
% subplot(3,3,3);plot(t,LRa);title('LRa');axis tight;hold all;
% subplot(3,3,4);plot(t,LRi);title('LRi');axis tight;hold all;
% subplot(3,3,5);plot(t,RaG);title('RaG');axis tight;hold all;
% subplot(3,3,6);plot(t,LRaG);title('LRaG');axis tight;hold all;
% subplot(3,3,7);plot(t,GsaGTPtot);title('GsaGTPtot');axis tight;hold all;
% subplot(3,3,8);plot(t, b1AR_S464);title('b1AR_{S464}');axis tight;hold all;
% subplot(3,3,9);plot(t,cAMPtot);title('cAMPtot');axis tight;hold all;

%% sensitivity analysis

% tspan = [0;10*60*1000];
% ISOrange = 10.^(-5:.1:3);
% deltaParams = [0.1,1,10];
% Params = {'KR','KL','KA','KG','alpha_L','alpha_A','gamma_L','gamma_A'};
% 
% for i= 1:length(Params)
%     KLp =KL;KRp =  KR;KGp =  KG;alpha_Lp = alpha_L;gamma_Lp = gamma_L; 
%    KAp = KA;alpha_Ap = alpha_A; gamma_Ap = gamma_A; 
%     for j = 1:length(deltaParams)
%          eval([Params{i} 'p = deltaParams(j)*' Params{i}]);
%        [KRp,KLp,KAp,KGp,alpha_Lp,alpha_Ap,gamma_Lp,gamma_Ap]
%         p = daePARAMS(KRp,KLp,KAp,KGp,alpha_Lp,alpha_Ap,gamma_Lp,gamma_Ap);p(1) = 0;%p(31:34) = 0;
%         y0 = zeros(29,1);[~,y] = ode15s(@daeODE,[0;60*60*1000],y0,[],p);y0 = y(end,:);
%         for k = 1:length(ISOrange)
%            
%             p(1) = ISOrange(k);  % Ltot 
%             [t,y] = ode15s(@daeODE,tspan,y0,[],p);
%             cAMP(k) = y(end,9);
%         end
%      subplot(4,2,i);semilogx(ISOrange*1e-6,cAMP/.168);hold all;
%    
%     end
%     title(Params{i});axis tight;
% %       load('Petroff_ISOcAMPdr.dat');xlabel('ISO(M)');ylabel('cAMP(pmol/mg)');
% %      semilogx(10.^(-Petroff_ISOcAMPdr(:,1)),Petroff_ISOcAMPdr(:,2),'o');
% end


%% ?ARK vs. cAMP
%  p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
%  options = odeset('RelTol',1e-7);
%  tspan = [0;200*60*1000];y0 = zeros(29,1);p(1) = 0;p(34) = 1;p(33) =2.2e-6;
%  [~,y] = ode15s(@daeODE,tspan,y0,options,p);y0 = y(end,:);
% 
% 
% 
%       tspan = [0;30*60*1000];p(1) = 1000;
%      [t,y] = ode15s(@daeODE,tspan,y0,options,p); bARK_resens = p(33)*y(:,9)./(p(34) + y(:,9));
%      plot(y(:,9),bARK_resens);
%        title([ 'k2_{bARK} = ' num2str(p(33)) '  Km_{bARK} = '  num2str(p(34))]);
%     
%     ylabel('bARK_{resens}'); xlabel('b1AR_{S464}');
     
%  ISOrange = 10.^(-5:.1:3);
%  for i = 1:length(ISOrange)
%      tspan = [0;20*60*1000];p(1) = ISOrange(i);
%      [t,y] = ode15s(@daeODE,tspan,y0,options,p);
%      cAMP(i) = y(end,9);b1AR_S464(i) = y(end,9);
%  end
%  plot(cAMP,b1AR_S464);title([ 'k2_{bARK} = ' num2str(p(33)) '  Km_{bARK} = '  num2str(p(34))]);
%  xlabel('cAMP'); ylabel('b1AR_{S464}');
%% Fitting KG and KR to match ISO cAMP dose response from 
%  load KG.mat;load gamma_L;%load KR.mat;
%  KL =15.2;%KL should be calculated by Ki(alpha*Kr+ 1)/alpha/(Kr+1)
%  KG = 0.7; gamma_L =  0.0313;
% KR =   5;
% alpha_L = 1/32;% 0.3037
% % gamma_L = 1; %1
% KA =  500e-6;%2e-3;
% alpha_A = 1; %1.5574
% gamma_A = 1; %1 
% p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
% p(1) = 0;
% y0 = zeros(29,1);[~,y] = ode15s(@daeODE,[0;60*60*1000],y0,[],p);
% y0 = y(end,:);RaGo = y(end,4);LRaGo = y(end,5);
% tspan = [0;20*60*1000]; 
% p(1) = 0.01;  
% [t,y] = ode15s(@daeODE,tspan,y0,[],p);
% figure(1);plot(t./60e3,y(:,9)/.168);%conversion from Bers Table 9...for Ca2+ fluxes
% load('Petroff_ISOcAMPtl.dat');hold all;
% plot(Petroff_ISOcAMPtl(:,1),Petroff_ISOcAMPtl(:,2),'o');
% xlabel('time(mins)');ylabel('cAMP(pmol/mg)');
% 
% ISOrange = 10.^[-5:.1:3];
% 
% 
% for i=1:length(ISOrange)
%     p(1) = ISOrange(i);  % Ltot 
%     [t,y] = ode15s(@daeODE,tspan,y0,[],p);
%     cAMP(i) = y(end,9);
% end
% RaG = y(end,4);LRaG = y(end,5);fold_change = (RaG + LRaG)/(RaGo + LRaGo);
% analytical = (KG*(KR + 1) + p(30))/(gamma_L*KG*(alpha_L*KR + 1) + p(30));
% figure(2);semilogx(ISOrange*1e-6,cAMP/.168);hold all;
% load('Petroff_ISOcAMPdr.dat');xlabel('ISO(M)');ylabel('cAMP(pmol/mg)');
% semilogx(10.^(Petroff_ISOcAMPdr(:,1)),Petroff_ISOcAMPdr(:,2),'o');

%% Adenylyl Cyclase Assay- replicating membranes w. no desensitization (skip to next section)
% Check to make sure parameters are correct, add changes in the antagonist
% and generally, get a clue.
% 
%  KG =  2.4131;gamma_L =  1;
%  p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
% tspan = [0;20*60*1000];  % 20 minutes
%  p(51) = 0; p(31:36) = 0;     % No desensitization or PDE's  
%  y0 = zeros(29,1);%p(17) = 0.59;p(23) = 0.025;
%  [t,y] = ode15s(@daeODE,[0;60*60*1000],y0,[],p);plot(t,y);
%  y0 = y(end,:); y0(9) = .0001;plot(t,y);
% 
% [t,y] = ode15s(@daeODE,tspan,y0,[],p);
% p(1)=1000;p(2) = 0;%80 pm ICYP was used in Hoffman paper
% [t2,y2] = ode15s(@daeODE,tspan,y0,[],p);
% 
% figure;plot(t./60e3,y(:,9),t2./60e3,y2(:,9))
% ratio = y2(end,9)/y(end,9)
% 
% % max cAMP;min cAMP; ratio; [b1ARi]; [RaRaG]; [Ra]; [RaG]
% simEpi = y2(end,9);simContr = y(end,9);
% 
% cAMPHoff = [65.4,30.9,76.785,73.335,53.67,50.91,31.59,32.97,...
%             27.105,19.515,29.175,24.345,22.62,22.275,24.345,22.965,...
%             21.24,22.62,19.515,18.825,22.275,23.31]';
% cAMPperc = (cAMPHoff-cAMPHoff(2))./(cAMPHoff(1)-cAMPHoff(2));%changing  to nepi
% cAMPdata = cAMPperc*(simEpi-simContr)+ simContr;
% save cAMPdata cAMPdata;




 %% Peforming parameter estimation on alpha using AC assay
%    y0 = zeros(29,1);
% 
%  KG =  2.4131;gamma_L =  1;
% 
%   params = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
% 
% params(1) = 0;params(51) = 0; params(31:36) = 0;
%   [t,y] = ode15s(@daeODE,[0;60*60*1000],y0,[],params);plot(t,y(:,9));
%     y0 = y(end,:); y0(9) = .0001;
% % 
% tspan = [0;60*20*1000];
%  simOptions = [];%odeset('NonNegative',[1:30]);
% 
% lowerbnd = 0.00001;
% upperbnd = 1000000;
% paramsEst = 1/32;   % use reasonable initial guesse for alpha
% params(1) = 1000;
% load cAMPdata;
% %Isoproterenol Basal Epinephrine Norepinephrine Fenoterol Formoterol
% %Salbutamol Terbutaline Broxaterol Salmeterol BRL-37344 CGP-12177
% %Alprenolol Pindolol SR 59230A Atenolol Carvedilol Metoprolol Bisoprolol
% %Propranolol CGP-20712 ICI-118551
% 
% for i=1:length(cAMPdata)
%     disp(['Fitting data set ',num2str(i),'\n']);
%     paramsEst = (upperbnd-lowerbnd).*rand(size(paramsEst))+lowerbnd;    % use random initial guesses (overwrites above line)
%     data = cAMPdata(i);
% 
%     % gradient based optimization
%     
%     gaOptions = gaoptimset('PlotFcns',@gaplotbestf);
%   numParams = length(paramsEst);
%     fitnessFcn = @(paramsEst) norm(objectiveFcn(KL,KR,gamma_L,alpha_A,gamma_A,paramsEst,params,simOptions,tspan,y0,data));   % uses anonymous function so it can utilize the same objective function as lsqnonlin
%     [paramsEst,resnorm]=ga(fitnessFcn,numParams,[],[],[],[],lowerbnd,upperbnd,[],gaOptions);
%     
%     optimOptions = optimset('Display','iter','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',[],'LargeScale','on','PlotFcns',@optimplotresnorm);
%     objFcn = @(paramsEst) objectiveFcn(KL,KR,gamma_L,alpha_A,gamma_A,paramsEst,params,simOptions,tspan,y0,data);
%     [paramsEst1,resnorm,residual,exitflag,output,lambda,Jacobian]=lsqnonlin(objFcn,paramsEst,lowerbnd,upperbnd,optimOptions);
%     disp(['Algorithm: ',output.algorithm]);
%     disp(['Exit flag: ',num2str(exitflag),', resnorm = ',num2str(resnorm)]);
%     paramsVariance = resnorm*inv(Jacobian'*Jacobian)/length(data(:,1));  % estimates the parameter estimate covariance using the Cramer-Rao inequality
%     paramsStd = sqrt(diag(paramsVariance));   % standard deviations for parameter estimates: 2*std is 95% conf interval
% 
% 
%     % redo gradient based using prior est.
%     [paramsEst2,resnorm,residual,exitflag,output,lambda,Jacobian]=lsqnonlin(objFcn,paramsEst1,lowerbnd,upperbnd,optimOptions);
%     paramsVariance3 = resnorm*inv(Jacobian'*Jacobian)/length(data(:,1));  % estimates the parameter estimate covariance using the Cramer-Rao inequality
%     paramsStd3 = sqrt(diag(paramsVariance3));   % standard deviations for parameter estimates: 2*std is 95% conf interval
%     
%     % save results from each
%     paramFits(i,:)=[paramsEst1,paramsEst2,paramsStd3,resnorm];
%     fitalpha(i) = paramsEst2;
%     disp(paramFits);
% 
%     % run simulation with the final param value
%  
% 
%     alpha_L=paramsEst2;
%     params = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
%     params(1) = 1000;params(51) = 0; params(31:36) = 0;
%     [t,y] = ode15s(@daeODE,tspan,y0,simOptions,params);
%     cAMPsim(i,1)=y(end,9);
%     
% end
% %  save fitting_data paramFits cAMPsim;
% % Then run each of the cases to verify the fit
% close all;bar([cAMPdata,cAMPsim]);
% Ki = [0.224,NaN,3.97,3.57,13.6,1.71,2.44,31.3,1.31,1.6,37.9,4.70e-03,...
%       5.80e-03,2.60e-03,1.64e-02,3.88e-01,1.70e-03,4.70e-02,2.24e-02,...
%       1.80e-03,4.50e-03,4.95e-02];
%  KLcalc = (Ki.*(fitalpha*KR + 1)./(fitalpha*(KR + 1)));
% save KLcalc KLcalc; save fitalpha fitalpha;
%% Receptor binding assays
% load KLcalc.mat; load fitalpha.mat;
% KL =  0.8553;
% alpha_L = fitalpha(1);gamma_L = 1;
% KG =  0.7;gamma_L =  0.3762;
% KA = 95e-6;
% alpha_A = 1;gamma_A = 1;
% p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
% tspan = [0;20*60*1000];  % 20 minutes
% 
% ISOrange = 10.^[-5:.1:2];
% y0 = zeros(29,1); y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
% RelTol = 1e-13;
% MaxStep = 1e-1;
% options = odeset('MaxStep',1e3,'NonNegative',[1:29],'RelTol',RelTol);
% t = [0;1000*60*1000]; p(1:2) = 0;[~,y] = ode15s(@daeODE,t,y0,options,p);y0 = y(end,:);
% p(2) = 10e-6;   % CYP concentration 40 pM is Mason/JBC1999
% p(51) = 0; p(31:36) = 0;
% for i=1:length(ISOrange)
%     p(1) = ISOrange(i);  % Ltot
%     p(30) = 3.83;    % Gtot = 3.83  
%     [t,y] = ode15s(@daeODE,tspan,y0,[],p);
%        for tstep=1:length(t),
%         [~,algvars(tstep,:)]=daeODE(t(tstep),y(tstep,:),p);
%        end
%     algvarsCell=mat2cell(algvars,size(algvars,1),ones(size(algvars,2),1));
%     [Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};
%     CYPboundwoGPP(i)= ARi(end)+ARa(end)+ARaG(end);
%     clear algvars algvarscell;
%       p(30) = 0;    % Gtot = 0  
%     [t,y] = ode15s(@daeODE,tspan,y0,[],p);
%     
%     for tstep=1:length(t),
%         [~,algvars(tstep,:)]=daeODE(t(tstep),y(tstep,:),p);
%     end
%     algvarsCell=mat2cell(algvars,size(algvars,1),ones(size(algvars,2),1));
%     [Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};
%     CYPboundwGPP(i)= ARi(end)+ARa(end)+ARaG(end);
%     clear algvars algvarscell;
%     
% end
% 
% figure(1);
% load('masonJBC1999_fig4.mat');
% semilogx(ISOrange*1e-6,CYPboundwoGPP/max(CYPboundwoGPP)*100,...
%     ISOrange*1e-6,CYPboundwGPP/max(CYPboundwGPP)*100,...
%     GlyIso,GlywoGPP,'bo',GlyIso,GlywGPP,'ro');
% xlabel('Isoproterenol (M)');
% ylabel('CYP bound (%)'); legend('-GPP','+GPP'); title('\beta_1-AR Gly389');
% figure;semilogx(ISOrange*1e-6,CYPboundwoGPP/max(CYPboundwoGPP)*100,...
%     ISOrange*1e-6,CYPboundwGPP/max(CYPboundwGPP)*100,...
%     ArgIso,ArgwoGPP,'bo',ArgIso2,ArgwGPP,'ro');
% xlabel('Isoproterenol (M)');
% ylabel('CYP bound (%)'); legend('-GPP','+GPP'); title('\beta_1-AR Arg389');


%% Adenylyl Cyclase Assay- Iso dose response- compared with Rathz JBC 2003
% load KLcalc.mat; load fitalpha.mat;
% KL =  0.8553;
% alpha_L = fitalpha(1);gamma_L =  0.3762;
% KA = 95e-6;
% alpha_A = 1;gamma_A = 1;
% tspan = [0;20*60*1000];  % 20 minutes
% y0 = zeros(29,1); y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
% RelTol = 1e-13;MaxStep = 1e-1;
% options = odeset('MaxStep',1e3,'NonNegative',[1:29],'RelTol',RelTol);
% pArg = daePARAMS(KR,KL,KA,0.7,alpha_L,alpha_A,gamma_L,gamma_A);
% pGly = daePARAMS(KR,KL,KA, 2.4131,alpha_L,alpha_A,gamma_L,gamma_A);
% pArg(1:2) = 0;pGly(1:2) = 0;
% t = [0;1000*60*1000]; [~,y] = ode15s(@daeODE,t,y0,options,pArg);y0 = y(end,:);
% ISOrange = 10.^[-5:.3:2];
% for i=1:length(ISOrange)
%     pArg(1) = ISOrange(i);  % Ltot
%      pGly(1) = ISOrange(i);  % Ltot
%     [t,y] = ode15s(@daeODE,tspan,y0,[],pArg);
%     cAMParg(i)=y(end,9);
%     [t,y] = ode15s(@daeODE,tspan,y0,[],pGly);
%     cAMPgly(i)=y(end,9);
% end
%  semilogx(ISOrange,cAMPgly/max(cAMPgly)*100,'color','b');xlabel('Isoproterenol (\muM)');
% ylabel('Adenylyl Cyclase Activity (% max)');
% hold on;
% semilogx(ISOrange,cAMParg/max(cAMPgly)*100,'color','r');
% 
% 
% % Plot experimental data from Rathz,Ligget JBC 2003
% load('rathz_ACactivity_fig2a.mat');
% 
%  errorbar(iso,gly_cAMP/max(gly_cAMP)*100,(arg_error-arg_cAMP)/max(arg_cAMP)*100,'o','color','b');
%   errorbar(iso,arg_cAMP/max(gly_cAMP)*100,(arg_error-arg_cAMP)/max(gly_cAMP)*100,'o','color','r');
%   errorbarlogx;
% hold off;

%% cAMP validation
% close all;
% load KLcalc.mat; load fitalpha.mat;
% KG =  2.4131;gamma_L =  0.3762;
% 
% KL = 0.5;%KLcalc(1);
% alpha_L = fitalpha(1);% 0.3037
% KA = KLcalc(20);%2e-3;
% alpha_A =  fitalpha(20);
%  p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);p(1) = 0;
%  y0 = zeros(29,1);y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
%  [t,y] = ode15s(@daeODE,[0;60*60*1000],y0,[],p);
%  y0 = y(end,:);plot(t,y);
% tspan = [0;20*60*1000];  % 20 minutes
% 
%  p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);p(1) = 0;
%  p(1) = .01;   
%  [t,y] = ode15s(@daeODE,tspan,y0,[],p);
% plot(t./60e3,y(:,9)/.168);%conversion from Bers Table 9...for Ca2+ fluxes
% load('Petroff_ISOcAMPtl.dat');hold all;
% plot(Petroff_ISOcAMPtl(:,1),Petroff_ISOcAMPtl(:,2),'o');
% xlabel('time(mins)');ylabel('cAMP(pmol/mg)');
% 
% % save ICscAMP y0;
% 
% ISOrange = 10.^[-5:.3:2];
% 
% p(1) = 0;
% for i=1:length(ISOrange)
%     p(1) = ISOrange(i);  % Ltot 
%     [t,y] = ode15s(@daeODE,tspan,y0,[],p);
%     cAMP(i) = (y(end,9));
%          
% end
% 
% figure;semilogx(ISOrange*1e-6,cAMP/.168);hold all;
% load('Petroff_ISOcAMPdr.dat');xlabel('ISO (M)');ylabel('cAMP (pmol/mg)');
% semilogx(10.^(Petroff_ISOcAMPdr(:,1)),Petroff_ISOcAMPdr(:,2),'o');
% plot(1e-8,Petroff_ISOcAMPtl(end,2),'o');

%% PLB validation
% close all;
% load KLcalc.mat; load fitalpha.mat;
% KG =  2.4131;gamma_L =  0.3762;
% 
% KL = 0.5;%KLcalc(1);
% alpha_L = fitalpha(1);% 0.3037
% KA = KLcalc(20);%2e-3;
% alpha_A =  fitalpha(20);
%  p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);p(1) = 0;
%  y0 = zeros(29,1);y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
% [~,y] = ode15s(@daeODE,[0;60*60*1000],y0,[],p);y0 = y(end,:);
% tspan = [0;20*60*1000];  % 20 minutes
% 
% ISOrange = 10.^[-5:.1:2];
% for i=1:length(ISOrange)
%     p(1) = ISOrange(i);  % Ltot 
%     p(2) = 1e-20;
%     [t,y] = ode15s(@daeODE,tspan,y0,[],p);
%     PLBp(i) = y(end,27);
% end
% 
% semilogx(ISOrange*1e-6,PLBp./max(PLBp));hold all;
% load Vittone.mat;
% semilogx(Vittone(:,1)*1e-9,Vittone(:,2)./100,'o');
% xlabel('ISO (M)');ylabel('PLBp (% of max)');

%% calcium ISO validation
% load KLcalc.mat; load fitalpha.mat;
% RelTol = 1e-6;MaxStep = 0.5;
% options = odeset('MaxStep',MaxStep,'RelTol',RelTol);
%  y0 = zeros(49,1);y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
% p = daePARAMSc(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
% [t,y] = ode15s(@daeODEc,[0;1*60*1000],y0,options,p);y0 = y(end,:);plot(t,y(:,47));
%% desensitization validation
% 
% close all;
% load KLcalc.mat; load fitalpha.mat; 
% y0 = zeros(29,1);
% % % % load ICs.mat;p = daePARAMS;close all;
% p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
% desens_factor =1;
% disp(p(31))
% disp(p(32))
% p(32) = desens_factor*p(32);
% disp(p(31))
% disp(p(32))
% % KL = KLcalc(1);
% % alpha_L = fitalpha(1);% 0.3037
% % KA =  KLcalc(20);
% % alpha_A = fitalpha(20);
% [~,y] = ode15s(@daeODE,[0;60*60*1000],y0,[],p);
% %control
%  yinitial = y(end,:);p(51) = 0; p(31:36) = 0;
% p(1) = 10; [t,y] = ode15s(@daeODE,[0;20*60*1000],yinitial,[],p);
% cAMP_max = y(end,9); figure(1);subplot(2,1,1); plot(t./60e3,y(:,15));
% subplot(2,1,2);plot(t./60e3,y(:,9));
% Rp_max = .0132;
% p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);p(32) =  desens_factor*p(32);
% %5 minute desensitization
% p(1) = 1;tspan = [0;10*60*1000];
%  [t,y] = ode15s(@daeODE,tspan,yinitial,[],p); b(1,1) =  y(end,9);c(1,1) =  y(end,9);
% figure;plot(t./60e3,y(:,15)); hold all
% p(1) = 10;tspan = [0*60*1000;20*60*1000]; y0 = y(end,:); p(51) = 0; p(31:36) = 0;
% [t,y] = ode15s(@daeODE,tspan,y0,[],p); b(1,2) =  y(end,9);c(1,2) =  y(end,9);
% plot(t./60e3,y(:,15));
% 
% p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);p(32) =desens_factor*p(32);
% %30 minute desensitization
% p(1) = 1;tspan = [0;30*60*1000]; 
%  [t,y] = ode15s(@daeODE,tspan,yinitial,[],p); b(2,1) =  y(end,9);
% figure;plot(t./60e3,y(:,15));hold all;
% p(1) = 10;tspan = [0*60*1000;20*60*1000]; y0 = y(end,:); p(51) = 0; p(31:36) = 0;
% [t,y] = ode15s(@daeODE,tspan,y0,[],p); b(2,2) =  y(end,9);c(2,2) =  y(end,9);
% plot(t./60e3,y(:,15));
% b = b./cAMP_max;c = c./Rp_max;c = 1-c;b(1,1) = 0.72;b(2,1) = 0.56;
% figure(5);bar(b);
% legend('expt','model'); xlabel('minutes'); ylabel('cAMP(% control of max)');
% c(1,1) = 0.72;c(2,1) = 0.56;figure;bar(c)
%% time step 
% % tic
% load KLcalc.mat; load fitalpha.mat;load Ki;
% KG =  2.4131;gamma_L =  0.3762;
% KL = KLcalc(1);
% alpha_L = fitalpha(1);
% KA = KLcalc(17);
% alpha_A = fitalpha(17);
% p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
% RelTol = 1e-13;
% MaxStep = 1e-1;
% options = odeset('MaxStep',1e3,'NonNegative',[1:29],'RelTol',RelTol);
%  y0 = zeros(29,1); y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
% t = [0;1000*60*1000]; p(1:2) = 0;[~,y] = ode15s(@daeODE,t,y0,options,p);y0 = y(end,:);
% tspan1 = [0;2*60*1000]; 
% p(1) = 0;p(2) = 0;
% [t1,y1] = ode15s(@daeODE,tspan1,y0,options,p);
% for tstep=1:length(t1),
%     [~,algvars1(tstep,:)]=daeODE(t1(tstep),y1(tstep,:),p);
% end
% options = odeset('MaxStep',MaxStep,'NonNegative',[1:29],'RelTol',RelTol);
% tspan = [2*60*1000; 2.1*60*1000]; y0 = y1(end,:);
% p(1) = 0.1;p(2) = 1;
% [t2,y2] = ode15s(@daeODE,tspan,y0,[],p);
% for tstep=1:length(t2),
%     [~,algvars2(tstep,:)]=daeODE(t2(tstep),y2(tstep,:),p);
% end
% options = odeset('MaxStep',MaxStep,'NonNegative',[1:29],'RelTol',RelTol);
% tspan = [2.1*60*1000; 6*60*1000]; y0 = y2(end,:);
% p(1) = 0.1;p(2) = 1;
% [t3,y3] = ode15s(@daeODE,tspan,y0,[],p);
% for tstep=1:length(t3),
%     [~,algvars3(tstep,:)]=daeODE(t3(tstep),y3(tstep,:),p);
% end
% 
% tspan = [6*60*1000; 10*60*1000]; y0 = y3(end,:);
% p(1) = 10;p(2) = 1;
% [t4,y4] = ode15s(@daeODE,tspan,y0,[],p);
% for tstep=1:length(t4),
%     [~,algvars4(tstep,:)]=daeODE(t4(tstep),y4(tstep,:),p);
% end
% 
% t = [t1;t2;t3;t4]./60e3;
%  y = [y1;y2;y3;y4];
%  algvars = [algvars1;algvars2;algvars3;algvars4];
% yCell=mat2cell(y,size(y,1),ones(size(y,2),1));
% algvarsCell=mat2cell(algvars,size(algvars,1),ones(size(algvars,2),1));
% [Ri,G, b1AR_S464,b1AR_S301,...
%     GsaGTPtot, GsaGDP, Gsby, AC_GsaGTP ,cAMPtot, PDEp, ...
%   RC_I, RCcAMP_I, RCcAMPcAMP_I, RcAMPcAMP_I, PKACI ,PKACI_PKI, ...
%   RC_II, RCcAMP_II, RCcAMPcAMP_II, RcAMPcAMP_II, PKACII, PKACII_PKI, ...
%   I1p_PP1, I1ptot, LCCap ,LCCbp, PLBp, PLMp, TnIp] = yCell{:};
% 
% [Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};
% 
% Rtot = Ri + Ra + LRi + LRa + RaG + LRaG + ARi + ARa + ARaG + b1AR_S464 + b1AR_S301 ;
% figure(1)
% subplot(5,3,1);plot(t,Ra);title('Ra');axis tight;hold all;
% subplot(5,3,2);plot(t,Ri);title('Ri');axis tight;hold all;
% subplot(5,3,3);plot(t,LRa);title('LRa');axis tight;hold all;
% subplot(5,3,4);plot(t,LRi);title('LRi');axis tight;hold all;
% subplot(5,3,5);plot(t,RaG);title('RaG');axis tight;hold all;
% subplot(5,3,6);plot(t,LRaG);title('LRaG');axis tight;hold all;
% subplot(5,3,7);plot(t,ARa);title('ARa');axis tight;hold all;
% subplot(5,3,8);plot(t,ARi);title('ARi');axis tight;hold all;
% subplot(5,3,9);plot(t,ARaG);title('ARaG');axis tight;hold all;
% subplot(5,3,10);plot(t,GsaGTPtot);title('GsaGTPtot');axis tight;hold all;
% subplot(5,3,11);plot(t,GsaGDP);title('GsaGDP');axis tight;hold all;
% subplot(5,3,12);plot(t,Gsby);title('Gsby');axis tight;hold all;
% subplot(5,3,13);plot(t, b1AR_S464);title('b1AR_{S464}');axis tight;hold all;
% subplot(5,3,14);plot(t, b1AR_S301);title('b1AR_{S301}');axis tight;hold all;
% subplot(5,3,15);plot(t,cAMPtot);title('cAMPtot');axis tight;hold all;
% plot(t,cAMPtot);title('cAMPtot');axis tight;hold all;
% toc

%% ISO & PRO
% load KLcalc.mat; load fitalpha.mat;%load ICs.mat;
% options = odeset('MaxStep',1e10,'NonNegative',[1:35],'RelTol',1e-3);
% KL = KLcalc(1);
% alpha_L = fitalpha(1);% 0.3037
% KA =  KLcalc(20);
% alpha_A = fitalpha(20);
% p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
% % % p(1) = 0;   
% % % p(6) = KR;  p(8)  = (alpha_L*KL); 
% % % p(10) = KL;  p(12) = (alpha_L*KR);
% % % p(14) = KG;  p(16) = (alpha_L*gamma_L*KL);
% % % p(18) = (gamma_L*KG); p(20)  = (alpha_L*KL);
%  y0 = zeros(29,1);
%  t = [0;1000*60*1000]; p(1:2) = 0;[~,y] = ode15s(@daeODE,t,y0,options,p);y0 = y(end,:);
% 
% 
% options = odeset('MaxStep',1e10,'NonNegative',[1:35],'RelTol',1e-10);
% tspan = [0;2*60*1000]; 
% p(1) = 0;p(2) = 0;
% [t1,y1] = ode15s(@daeODE,tspan,y0,[],p);
% yfinal = y1(end,:);
% tspan = [2*60*1000;6*60*1000]; 
% p(1) = 0.1;p(2) = 0;
% [t2,y2] = ode15s(@daeODE,tspan,yfinal,[],p);
% yfinal = y2(end,:);
% tspan = [6*60*1000;10*60*1000]; 
% p(1) = 10;p(2) = 0;
% [t3,y3] = ode15s(@daeODE,tspan,yfinal,[],p);
% t = [t1;t2;t3]./60e3;
% y = [y1;y2;y3];
% % plot(t,y(:,15)); hold all;
% yCell=mat2cell(y,size(y,1),ones(size(y,2),1));
% [Ra, LRa ,LRi, RaG ,LRaG ,ARi ,ARa, ARaG, b1AR_S464 ,b1AR_S301,...
%     GsaGTPtot, GsaGDP, Gsby, AC_GsaGTP, cAMPtot, PDEp, ...
%   RC_I, RCcAMP_I, RCcAMPcAMP_I, RcAMPcAMP_I, PKACI, PKACI_PKI, ...
%   RC_II ,RCcAMP_II, RCcAMPcAMP_II, RcAMPcAMP_II, PKACII, PKACII_PKI, ...
%   I1p_PP1, I1ptot, LCCap, LCCbp, PLBp, PLMp, TnIp] = yCell{:};
% 
% 
% Ri = 0.0132- Ra - LRi - LRa - RaG - LRaG - ARi - ARa - ARaG;
% Rtot = Ri + Ra + LRi + LRa + RaG + LRaG + ARi + ARa + ARaG;
% figure(1)
% subplot(5,3,1);plot(t,Ra);title('Ra');axis tight;hold all;
% subplot(5,3,2);plot(t,Rtot);title('Rtot');axis tight;hold all;
% subplot(5,3,3);plot(t,LRa);title('LRa');axis tight;hold all;
% subplot(5,3,4);plot(t,LRi);title('LRi');axis tight;hold all;
% subplot(5,3,5);plot(t,RaG);title('RaG');axis tight;hold all;
% subplot(5,3,6);plot(t,LRaG);title('LRaG');axis tight;hold all;
% subplot(5,3,7);plot(t,ARa);title('ARa');axis tight;hold all;
% subplot(5,3,8);plot(t,ARi);title('ARi');axis tight;hold all;
% subplot(5,3,9);plot(t,ARaG);title('ARaG');axis tight;hold all;
% subplot(5,3,10);plot(t,GsaGTPtot);title('GsaGTPtot');axis tight;hold all;
% subplot(5,3,11);plot(t,GsaGDP);title('GsaGDP');axis tight;hold all;
% subplot(5,3,12);plot(t,Gsby);title('Gsby');axis tight;hold all;
% subplot(5,3,13);plot(t, b1AR_S464);title('b1AR_{S464}');axis tight;hold all;
% subplot(5,3,14);plot(t, b1AR_S301);title('b1AR_{S301}');axis tight;hold all;
% subplot(5,3,15);plot(t,cAMPtot);title('cAMPtot');axis tight;hold all;
% 
% tspan = [0;2*60*1000]; 
% p(1) = 0;p(2) = 0;
% [t1,y1] = ode15s(@daeODE,tspan,y0,options,p);
% yfinal = y1(end,:);
% tspan = [2*60*1000;6*60*1000]; 
% p(1) = 0;p(2) =0.1;
% [t2,y2] = ode15s(@daeODE,tspan,yfinal,options,p);
% yfinal = y2(end,:);
% tspan = [6*60*1000;10*60*1000]; 
% p(1) = 0;p(2) = 0.1;
% [t3,y3] = ode15s(@daeODE,tspan,yfinal,options,p);
% t = [t1;t2;t3]./60e3;
% y = [y1;y2;y3];
% % plot(t,y(:,15)); hold all;xlabel('time(minutes)');ylabel('cAMP(\muM)');
% yCell=mat2cell(y,size(y,1),ones(size(y,2),1));
% [Ra, LRa ,LRi, RaG ,LRaG ,ARi ,ARa, ARaG, b1AR_S464 ,b1AR_S301,...
%     GsaGTPtot, GsaGDP, Gsby, AC_GsaGTP, cAMPtot, PDEp, ...
%   RC_I, RCcAMP_I, RCcAMPcAMP_I, RcAMPcAMP_I, PKACI, PKACI_PKI, ...
%   RC_II ,RCcAMP_II, RCcAMPcAMP_II, RcAMPcAMP_II, PKACII, PKACII_PKI, ...
%   I1p_PP1, I1ptot, LCCap, LCCbp, PLBp, PLMp, TnIp] = yCell{:};
% Ri = 0.0132- Ra - LRi - LRa - RaG - LRaG - ARi - ARa - ARaG;
% Rtot = Ri + Ra + LRi + LRa + RaG + LRaG + ARi + ARa + ARaG;
% subplot(5,3,1);plot(t,Ra);title('Ra');axis tight;hold all;
% subplot(5,3,2);plot(t,Rtot);title('Rtot');axis tight;hold all;
% subplot(5,3,3);plot(t,LRa);title('LRa');axis tight;hold all;
% subplot(5,3,4);plot(t,LRi);title('LRi');axis tight;hold all;
% subplot(5,3,5);plot(t,RaG);title('RaG');axis tight;hold all;
% subplot(5,3,6);plot(t,LRaG);title('LRaG');axis tight;hold all;
% subplot(5,3,7);plot(t,ARa);title('ARa');axis tight;hold all;
% subplot(5,3,8);plot(t,ARi);title('ARi');axis tight;hold all;
% subplot(5,3,9);plot(t,ARaG);title('ARaG');axis tight;hold all;
% subplot(5,3,10);plot(t,GsaGTPtot);title('GsaGTPtot');axis tight;hold all;
% subplot(5,3,11);plot(t,GsaGDP);title('GsaGDP');axis tight;hold all;
% subplot(5,3,12);plot(t,Gsby);title('Gsby');axis tight;hold all;
% subplot(5,3,13);plot(t, b1AR_S464);title('b1AR_{S464}');axis tight;hold all;
% subplot(5,3,14);plot(t, b1AR_S301);title('b1AR_{S301}');axis tight;hold all;
% subplot(5,3,15);plot(t,cAMPtot);title('cAMPtot');axis tight;hold all;
%  

%% ISO and CAR
% tic
% load KLcalc.mat; load fitalpha.mat;load Ki;
% 
% KL = KLcalc(1);
% alpha_L = fitalpha(1);
% KA = KLcalc(17);
% alpha_A = fitalpha(17);
% p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
% RelTol = 1e-13;
% MaxStep = 1e-1;
% options = odeset('MaxStep',1e3,'NonNegative',[1:29],'RelTol',RelTol);
%  y0 = zeros(29,1); y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
% t = [0;1000*60*1000]; p(1:2) = 0;[~,y] = ode15s(@daeODE,t,y0,options,p);y0 = y(end,:);
% tspan1 = [0;2*60*1000]; 
% p(1) = 0;p(2) = 0;
% [t1,y1] = ode15s(@daeODE,tspan1,y0,options,p);
% BBdose = 0;
% for tstep=1:length(t1),
%     [~,algvars1(tstep,:)]=daeODE(t1(tstep),y1(tstep,:),p);
% end
% options = odeset('MaxStep',MaxStep,'NonNegative',[1:29],'RelTol',RelTol);
% tspan = [2*60*1000; 2.0167*60*1000]; y0 = y1(end,:);
% p(1) = 0.1;p(2) = BBdose;
% [t2,y2] = ode15s(@daeODE,tspan,y0,[],p);
% for tstep=1:length(t2),
%     [~,algvars2(tstep,:)]=daeODE(t2(tstep),y2(tstep,:),p);
% end
% options = odeset('MaxStep',MaxStep,'NonNegative',[1:29],'RelTol',RelTol);
% tspan = [2.0167*60*1000; 6*60*1000]; y0 = y2(end,:);
% p(1) = 0.1;p(2) =BBdose;
% [t3,y3] = ode15s(@daeODE,tspan,y0,[],p);
% for tstep=1:length(t3),
%     [~,algvars3(tstep,:)]=daeODE(t3(tstep),y3(tstep,:),p);
% end
% 
% tspan = [6*60*1000; 10*60*1000]; y0 = y3(end,:);
% p(1) = 10;p(2) = BBdose;
% [t4,y4] = ode15s(@daeODE,tspan,y0,[],p);
% for tstep=1:length(t4),
%     [~,algvars4(tstep,:)]=daeODE(t4(tstep),y4(tstep,:),p);
% end
% 
% t = [t1;t2;t3;t4]./60e3;
%  y = [y1;y2;y3;y4];
%  algvars = [algvars1;algvars2;algvars3;algvars4];
% yCell=mat2cell(y,size(y,1),ones(size(y,2),1));
% algvarsCell=mat2cell(algvars,size(algvars,1),ones(size(algvars,2),1));
% [Ri,G, b1AR_S464,b1AR_S301,...
%     GsaGTPtot, GsaGDP, Gsby, AC_GsaGTP ,cAMPtot, PDEp, ...
%   RC_I, RCcAMP_I, RCcAMPcAMP_I, RcAMPcAMP_I, PKACI ,PKACI_PKI, ...
%   RC_II, RCcAMP_II, RCcAMPcAMP_II, RcAMPcAMP_II, PKACII, PKACII_PKI, ...
%   I1p_PP1, I1ptot, LCCap ,LCCbp, PLBp, PLMp, TnIp] = yCell{:};
% 
% [Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};
% 
% Rtot = Ri + Ra + LRi + LRa + RaG + LRaG + ARi + ARa + ARaG + b1AR_S464 + b1AR_S301 ;
% figure(1)
% % subplot(5,3,1);plot(t,Rtot);title('Rtot');axis tight;hold all;
% % subplot(5,3,2);plot(t,Ri);title('Ri');axis tight;hold all;
% % subplot(5,3,3);plot(t,LRa);title('LRa');axis tight;hold all;
% % subplot(5,3,4);plot(t,LRi);title('LRi');axis tight;hold all;
% % subplot(5,3,5);plot(t,RaG);title('RaG');axis tight;hold all;
% % subplot(5,3,6);plot(t,LRaG);title('LRaG');axis tight;hold all;
% % subplot(5,3,7);plot(t,ARa);title('ARa');axis tight;hold all;
% % subplot(5,3,8);plot(t,ARi);title('ARi');axis tight;hold all;
% % subplot(5,3,9);plot(t,ARaG);title('ARaG');axis tight;hold all;
% % subplot(5,3,10);plot(t,GsaGTPtot);title('GsaGTPtot');axis tight;hold all;
% % subplot(5,3,11);plot(t,GsaGDP);title('GsaGDP');axis tight;hold all;
% % subplot(5,3,12);plot(t,Gsby);title('Gsby');axis tight;hold all;
% % subplot(5,3,13);plot(t, b1AR_S464);title('b1AR_{S464}');axis tight;hold all;
% % subplot(5,3,14);plot(t, b1AR_S301);title('b1AR_{S301}');axis tight;hold all;
% % subplot(5,3,15);plot(t,cAMPtot);title('cAMPtot');axis tight;hold all;
% plot(t,cAMPtot);title('cAMPtot');axis tight;hold all;
% toc


%% ISO and MET
% tic
% load KLcalc.mat; load fitalpha.mat;load Ki;
% 
% KL = KLcalc(1);
% alpha_L = fitalpha(1);
% KA = KLcalc(18);
% alpha_A = fitalpha(18);
% p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
% RelTol = 1e-13;
% MaxStep = 1e-1;
% options = odeset('MaxStep',1e3,'NonNegative',[1:29],'RelTol',RelTol);
%  y0 = zeros(29,1); y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
% t = [0;1000*60*1000]; p(1:2) = 0;[~,y] = ode15s(@daeODE,t,y0,options,p);y0 = y(end,:);
% tspan1 = [0;2*60*1000]; 
% p(1) = 0;p(2) = 0;
% [t1,y1] = ode15s(@daeODE,tspan1,y0,options,p);
% for tstep=1:length(t1),
%     [~,algvars1(tstep,:)]=daeODE(t1(tstep),y1(tstep,:),p);
% end
% options = odeset('MaxStep',MaxStep,'NonNegative',[1:29],'RelTol',RelTol);
% tspan = [2*60*1000; 2.0167*60*1000]; y0 = y1(end,:);
% p(1) = 0.1;p(2) = 1;
% [t2,y2] = ode15s(@daeODE,tspan,y0,[],p);
% for tstep=1:length(t2),
%     [~,algvars2(tstep,:)]=daeODE(t2(tstep),y2(tstep,:),p);
% end
% options = odeset('MaxStep',MaxStep,'NonNegative',[1:29],'RelTol',RelTol);
% tspan = [2.0167*60*1000; 6*60*1000]; y0 = y2(end,:);
% p(1) = 0.1;p(2) =1;
% [t3,y3] = ode15s(@daeODE,tspan,y0,[],p);
% for tstep=1:length(t3),
%     [~,algvars3(tstep,:)]=daeODE(t3(tstep),y3(tstep,:),p);
% end
% 
% tspan = [6*60*1000; 10*60*1000]; y0 = y3(end,:);
% p(1) = 10;p(2) = 1;
% [t4,y4] = ode15s(@daeODE,tspan,y0,[],p);
% for tstep=1:length(t4),
%     [~,algvars4(tstep,:)]=daeODE(t4(tstep),y4(tstep,:),p);
% end
% 
% t = [t1;t2;t3;t4]./60e3;
%  y = [y1;y2;y3;y4];
%  algvars = [algvars1;algvars2;algvars3;algvars4];
% yCell=mat2cell(y,size(y,1),ones(size(y,2),1));
% algvarsCell=mat2cell(algvars,size(algvars,1),ones(size(algvars,2),1));
% [Ri,G, b1AR_S464,b1AR_S301,...
%     GsaGTPtot, GsaGDP, Gsby, AC_GsaGTP ,cAMPtot, PDEp, ...
%   RC_I, RCcAMP_I, RCcAMPcAMP_I, RcAMPcAMP_I, PKACI ,PKACI_PKI, ...
%   RC_II, RCcAMP_II, RCcAMPcAMP_II, RcAMPcAMP_II, PKACII, PKACII_PKI, ...
%   I1p_PP1, I1ptot, LCCap ,LCCbp, PLBp, PLMp, TnIp] = yCell{:};
% 
% [Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};
% 
% Rtot = Ri + Ra + LRi + LRa + RaG + LRaG + ARi + ARa + ARaG + b1AR_S464 + b1AR_S301 ;
% figure(1)
% % subplot(5,3,1);plot(t,Rtot);title('Rtot');axis tight;hold all;
% % subplot(5,3,2);plot(t,Ri);title('Ri');axis tight;hold all;
% % subplot(5,3,3);plot(t,LRa);title('LRa');axis tight;hold all;
% % subplot(5,3,4);plot(t,LRi);title('LRi');axis tight;hold all;
% % subplot(5,3,5);plot(t,RaG);title('RaG');axis tight;hold all;
% % subplot(5,3,6);plot(t,LRaG);title('LRaG');axis tight;hold all;
% % subplot(5,3,7);plot(t,ARa);title('ARa');axis tight;hold all;
% % subplot(5,3,8);plot(t,ARi);title('ARi');axis tight;hold all;
% % subplot(5,3,9);plot(t,ARaG);title('ARaG');axis tight;hold all;
% % subplot(5,3,10);plot(t,GsaGTPtot);title('GsaGTPtot');axis tight;hold all;
% % subplot(5,3,11);plot(t,GsaGDP);title('GsaGDP');axis tight;hold all;
% % subplot(5,3,12);plot(t,Gsby);title('Gsby');axis tight;hold all;
% % subplot(5,3,13);plot(t, b1AR_S464);title('b1AR_{S464}');axis tight;hold all;
% % subplot(5,3,14);plot(t, b1AR_S301);title('b1AR_{S301}');axis tight;hold all;
% % subplot(5,3,15);plot(t,cAMPtot);title('cAMPtot');axis tight;hold all;
% plot(t,cAMPtot);title('cAMPtot');axis tight;hold all;
% toc
%% Beta Blocker comparisons
% tic
% load KLcalc.mat; load fitalpha.mat;load Ki;
% KG =  2.4131;gamma_L =  0.3762;
% p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
% RelTol = 1e-13;MaxStep = 1e-1;
% options = odeset('MaxStep',1e3,'NonNegative',[1:29],'RelTol',RelTol);
%  y0 = zeros(29,1); y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
% t = [0;1000*60*1000]; p(1:2) = 0;[~,y] = ode15s(@daeODE,t,y0,options,p);y0 = y(end,:);
% tspan = [0;2*60*1000]; 
% 
% [~,y] = ode15s(@daeODE,tspan,y0,options,p);
% bg(1,1) = max(y(end,9));
% yfinal = y(end,:); 
% 
% tspan = [2*60*1000;6*60*1000];
% p(1) = 0.1; p(2) = 0;  % agonist conc and beta blocker conc
% [~,y2] = ode15s(@daeODE,tspan,yfinal,options,p); 
% bg(1,2) = max(y2(end,9));
% yfinal = y2(end,:); 
% 
% tspan2 = [6*60*1000 10*60*1000]; %changed the time
% p(1) = 10; p(2) = 0;     % agonist conc and beta blocker conc
% [~,y3]=ode15s(@daeODE,tspan2,yfinal,options,p);
% bg(1,3) = max(y3(end,9));
% 
% dose = [0.1 1 1];drugs = [20 18 17];
% 
% for i=1:length(drugs)
%    
%   KL =  KLcalc(1);
%    alpha_L = fitalpha(1);
%    KA = KLcalc(drugs(i)); alpha_A = fitalpha(drugs(i));
%    p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
%     
%     tspan = [0;2*60*1000];
% p(1) = 0; p(2) = 0;  % agonist conc and beta blocker conc
% [~,y] = ode15s(@daeODE,tspan,y0,options,p);
% bg(i+1,1) = max(y(end,9));
% yfinal = y(end,:); 
% 
% tspan = [2*60*1000;6*60*1000];
% p(1) = 0.1; p(2) = dose(i);  % agonist conc and beta blocker conc
% [~,y2] = ode15s(@daeODE,tspan,yfinal,options,p); 
% bg(i+1,2) = max(y2(end,9));
% yfinal = y2(end,:); 
% 
% tspan2 = [6*60*1000 10*60*1000]; %changed the time
% p(1) = 10; p(2) = dose(i);     % agonist conc and beta blocker conc
% [~,y3]=ode15s(@daeODE,tspan2,yfinal,options,p);
% bg(i+1,3) = max(y3(end,9));
% 
% end
% 
% bar(bg);

%% Beta Blocker comparisons calcium
% tic
% load KLcalc.mat; load fitalpha.mat;
% KG =  2.4131;gamma_L =  0.3762;
%  p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);p(1) = 0;
%  y0 = zeros(29,1);y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
%  [~,y] = ode15s(@daeODE,[0;60*60*1000],y0,[],p);y0Sig = y(end,:);
% 
% 
% options = odeset('RelTol',1e-5,'MaxStep',5e-3,'Stats','on');
% p = daePARAMSc(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
% load y0full;y0 = y0full;
% 
% 
% tspan = [0;2*60];
% [t,y] = ode15s(@daeODEc,tspan,y0,options,p);
% bgcAMP{1,1} = y(:,9);
% bgCa{1,1} = y(:,47);
% time{1,1} = t;
% yfinal = y(end,:); 
% 
% tspan = [2*60;6*60];
% p(1) = 0.1; p(2) = 0;  % agonist conc and beta blocker conc
% [t2,y2] = ode15s(@daeODEc,tspan,yfinal,options,p); 
% bgcAMP{1,2} = y2(:,9);
% bgCa{1,2} = y2(:,47);
% time{1,2} = t2;
% yfinal = y2(end,:); 
% 
% tspan2 = [6*60 10*60]; %changed the time
% p(1) = 10; p(2) = 0;     % agonist conc and beta blocker conc
% [t3,y3]=ode15s(@daeODEc,tspan2,yfinal,options,p);
% bgcAMP{1,3} = y3(:,9);
% bgCa{1,3} = y3(:,47);
% time{1,3} = t3;
% 
% dose = [0.1 1 1];drugs = [20 18 17];
% 
% for i=1:length(drugs)
%    
%   KL = KLcalc(1);
%    alpha_L = fitalpha(1);
%    KA = KLcalc(drugs(i)); alpha_A = fitalpha(drugs(i));
%    p = daePARAMSc(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
%     
%     tspan = [0;2*60];
% p(1) = 0; p(2) = 0;  % agonist conc and beta blocker conc
% [t,y] = ode15s(@daeODEc,tspan,y0,options,p);
% 
% bgcAMP{i+1,1} = y(:,9);
% bgCa{i+1,1} = y(:,47);
% time{i+1,1} = t;
% yfinal = y(end,:); 
% 
% tspan = [2*60;6*60];
% p(1) = 0.1; p(2) = dose(i);  % agonist conc and beta blocker conc
% [t2,y2] = ode15s(@daeODEc,tspan,yfinal,options,p); 
% bgcAMP{i+1,2} = y2(:,9);
% bgCa{i+1,2} = y2(:,47);
% time{i+1,2} = t2;
% yfinal = y2(end,:); 
% tspan2 = [6*60 10*60]; %changed the time
% p(1) = 10; p(2) = dose(i);     % agonist conc and beta blocker conc
% [t3,y3]=ode15s(@daeODEc,tspan2,yfinal,options,p);
% bgcAMP{i+1,3} = y3(:,9);
% bgCa{i+1,3} = y3(:,47);
% time{i+1,3} = t3;
% end
% cd('A:\Robert\Dissertation\Beta Blocker modeling\final_model validation\data');
% save bgcAMP bgcAMP; save bgCa bgCa; save time time;
%% Polymorphisms simulation


% load KLcalc.mat; load fitalpha.mat;load Ki;
% KG =  2.4131;gamma_L =  0.3762;
% KL = KLcalc(1);
% alpha_L = fitalpha(1);
% KA = KLcalc(20);
% alpha_A = fitalpha(20);
% p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
% RelTol = 1e-13;
% MaxStep = 1e-1;
% options = odeset('MaxStep',1e3,'NonNegative',[1:29],'RelTol',RelTol);
%  y0 = zeros(29,1); y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
% t = [0;1000*60*1000]; p(1:2) = 0;[~,y] = ode15s(@daeODE,t,y0,options,p);y0 = y(end,:);
% tspan1 = [0;2*60*1000]; 
% p(1) = 0;p(2) = 0;
% [t1,y1] = ode15s(@daeODE,tspan1,y0,options,p);
% 
% options = odeset('MaxStep',MaxStep,'NonNegative',[1:29],'RelTol',RelTol);
% tspan = [2*60*1000; 2.1*60*1000]; y0 = y1(end,:);
% p(1) = 0.1;p(2) = 0.1;
% [t2,y2] = ode15s(@daeODE,tspan,y0,[],p);
% options = odeset('MaxStep',MaxStep,'NonNegative',[1:29],'RelTol',RelTol);
% tspan = [2.1*60*1000; 6*60*1000]; y0 = y2(end,:);
% p(1) = 0.1;p(2) = 0.1;
% [t3,y3] = ode15s(@daeODE,tspan,y0,[],p);
% tspan = [6*60*1000; 10*60*1000]; y0 = y3(end,:);
% p(1) = 10;p(2) = 0.1;
% [t4,y4] = ode15s(@daeODE,tspan,y0,[],p);
% 
% t = [t1;t2;t3;t4]./60e3;
%  y = [y1;y2;y3;y4];
% 
% plot(t,y(:,9));title('cAMPtot');axis tight;hold all;
% 
% 
% KG =  0.7;gamma_L =  0.3762;
% KL = KLcalc(1);
% alpha_L = fitalpha(1);
% KA = KLcalc(20);
% alpha_A = fitalpha(20);
% p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
% RelTol = 1e-13;
% MaxStep = 1e-1;
% options = odeset('MaxStep',1e3,'NonNegative',[1:29],'RelTol',RelTol);
%  y0 = zeros(29,1); y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
% t = [0;1000*60*1000]; p(1:2) = 0;[~,y] = ode15s(@daeODE,t,y0,options,p);y0 = y(end,:);
% tspan1 = [0;2*60*1000]; 
% p(1) = 0;p(2) = 0;
% [t1,y1] = ode15s(@daeODE,tspan1,y0,options,p);
% 
% options = odeset('MaxStep',MaxStep,'NonNegative',[1:29],'RelTol',RelTol);
% tspan = [2*60*1000; 2.1*60*1000]; y0 = y1(end,:);
% p(1) = 0.1;p(2) = 0.1;
% [t2,y2] = ode15s(@daeODE,tspan,y0,[],p);
% options = odeset('MaxStep',MaxStep,'NonNegative',[1:29],'RelTol',RelTol);
% tspan = [2.1*60*1000; 6*60*1000]; y0 = y2(end,:);
% p(1) = 0.1;p(2) = 0.1;
% [t3,y3] = ode15s(@daeODE,tspan,y0,[],p);
% tspan = [6*60*1000; 10*60*1000]; y0 = y3(end,:);
% p(1) = 10;p(2) = 0.1;
% [t4,y4] = ode15s(@daeODE,tspan,y0,[],p);
% 
% t = [t1;t2;t3;t4]./60e3;
%  y = [y1;y2;y3;y4];
% 
% plot(t,y(:,9));title('cAMPtot');axis tight;


%%

%  load fitalpha.mat;load Ki;
%  drugList = {'Isoproterenol'; 'Basal' ;'Epinephrine'; 'Norepinephrine'; 'Fenoterol' ;'Formoterol';
% 'Salbutamol' ;'Terbutaline' ;'Broxaterol' ;'Salmeterol'; 'BRL-37344'; 'CGP-12177';
% 'Alprenolol' ;'Pindolol'; 'SR 59230A'; 'Atenolol' ;'Carvedilol' ;'Metoprolol'; 'Bisoprolol';
% 'Propranolol' ;'CGP-20712' ;'ICI-118551'};
% 
%  Ki = [0.224,NaN,3.97,3.57,13.6,1.71,2.44,31.3,1.31,1.6,37.9,4.70e-03,...
%       5.80e-03,2.60e-03,1.64e-02,3.88e-01,1.70e-03,4.70e-02,2.24e-02,...
%       1.80e-03,4.50e-03,4.95e-02];
%   KLcalc =Ki.*(fitalpha*KR + 1)./(fitalpha*(KR + 1));
%    KLcalc = (0.2/KLcalc(1))*KLcalc; 
% KG = 2.4131;gamma_L =  0.3762;
% p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
% 
% RelTol = 1e-13;MaxStep = 1e-1;
% options = odeset('MaxStep',1e3,'NonNegative',[1:29],'RelTol',RelTol);
%  y0 = zeros(29,1); y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
% t = [0;1000*60*1000]; p(1:2) = 0;[~,y] = ode15s(@daeODE,t,y0,options,p);y0 = y(end,:);
% tspan = [0;2*60*1000]; 
% drugs = 1:22;
% dose = ones(size(drugs))*1e-0;%dose([17 18 20])  = [ 1 1 0.1];
% % dose = Ki*(0.1/.01 + 1)*10;
% constraint = 1e2;
% for i=1:length(drugs)
%    
%   KL =  KLcalc(1);
%    alpha_L = fitalpha(1);
%    KA = KLcalc(drugs(i)); alpha_A = fitalpha(drugs(i));
%   KAdrugs(i) = KLcalc(drugs(i)); alphadrugs(i) = fitalpha(drugs(i));
%   p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
%     tspan = [0;2*60*1000];
% p(1) = 0; p(2) = 0;  % agonist conc and beta blocker conc
% [t1,y] = ode15s(@daeODE,tspan,y0,options,p);
% bg(i,1) = max(y(end,9));
% yfinal = y(end,:); 
% tspan = [2*60*1000;6*60*1000];
% p(1) = 0.1; p(2) = dose(i);  % agonist conc and beta blocker conc
% [t2,y2] = ode15s(@daeODE,tspan,yfinal,options,p); 
% bg(i,2) = max(y2(end,9));
% yfinal = y2(end,:); 
% tspan2 = [6*60*1000 10*60*1000]; %changed the time
% p(1) = 10; p(2) = dose(i);     % agonist conc and beta blocker conc
% [t3,y3]=ode15s(@daeODE,tspan2,yfinal,options,p);
% bg(i,3) = max(y3(end,9));
% 
% end
% y = [y;y2;y3];t = [t1;t2;t3];
% for tstep=1:length(t),
%     [~,algvars(tstep,:)]=daeODE(t(tstep),y(tstep,:),p);
% end
% algvarsCell=mat2cell(algvars,size(algvars,1),ones(size(algvars,2),1));
% [Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};
% 
% figure(3);bar(dose);ylabel('dose(\muM)');
% set(gca,'XTick',1:length(drugList),'XTickLabel',drugList);xticklabel_rotate([],45);
% 
% sensitivity = (bg(:,3)- bg(:,2));figure(1); bar(bg);set(gca,'XTick',1:length(drugList),'XTickLabel',drugList);xticklabel_rotate([],45);
% sensitivity2 = (bg(:,3)- bg(:,2))./bg(:,3);sensitivity3 = (bg(:,3)- bg(:,2))./bg(1,3);
% sensitivity(2:4) = [];KAdrugs(2:4) = [];alphadrugs(2:4) = [];drugList(2:4) = [];
% figure(2);subplot(2,2,1:2); bar(sensitivity);xlabel('ligand');ylabel('sens.(cAMP \muM)');box off;axis tight;
% set(gca,'XTick',1:length(drugList),'XTickLabel',drugList);xticklabel_rotate([],45);
% subplot(2,2,3);scatter(KAdrugs*1e-6,sensitivity);ylabel('sens.(cAMP \muM)');xlabel('K_L (M)');
% subplot(2,2,4);scatter(alphadrugs,sensitivity);ylabel('sens.(cAMP \muM)');xlabel('\alpha');
% sensitivity2(2:4) = [];sensitivity3(2:4) = [];

% dx =1e-12; dy = 0.01;text(KAdrugs*1e-6+dx,sensitivity+dx,drugList)
%%
% 
% load KLcalc; load fitalpha;
% KL = KLcalc(1);
% alpha_L = fitalpha(1);
% KA = KLcalc(20);
% alpha_A = fitalpha(20);
% KG = 2.4131;gamma_L =  0.3762;
% p = daePARAMSc(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
% 
% options = odeset('RelTol',1e-5,'MaxStep',5e-3,'Stats','on');
% load y0full;y0 = y0full;
% 
% tspan1 = [0;2*60]; 
% p(1) = 0;p(2) = 0;
% [t1,y1] = ode15s(@daeODEc,tspan1,y0,options,p);
% for tstep=1:length(t1),
%     [~,algvars1(tstep,:)]=daeODEc(t1(tstep),y1(tstep,:),p);
% end
% 
% tspan = [2*60; 2.1*60]; y0 = y1(end,:);
% p(1) = 0.1;p(2) = 0.1;
% [t2,y2] = ode15s(@daeODEc,tspan,y0,options,p);
% for tstep=1:length(t2),
%     [~,algvars2(tstep,:)]=daeODEc(t2(tstep),y2(tstep,:),p);
% end
% 
% 
% tspan = [2.1*60; 6*60]; y0 = y2(end,:);
% p(1) = 0.1;p(2) = 0.1;
% [t3,y3] = ode15s(@daeODEc,tspan,y0,options,p);
% 
% for tstep=1:length(t3),
%     [~,algvars3(tstep,:)]=daeODEc(t3(tstep),y3(tstep,:),p);
% end
% tspan = [6*60; 10*60]; y0 = y3(end,:);
% p(1) = 10;p(2) = 0.1;
% [t4,y4] = ode15s(@daeODEc,tspan,y0,options,p);
% for tstep=1:length(t4),
%     [~,algvars4(tstep,:)]=daeODEc(t4(tstep),y4(tstep,:),p);
% end
% 
% tGly = [t1;t2;t3;t4]./60;
%  yGly = [y1;y2;y3;y4];
% algvarsGly = [algvars1;algvars2;algvars3;algvars4];
% figure(12);
% 
% yCell=mat2cell(yGly,size(yGly,1),ones(size(yGly,2),1));
% 
% [Ri,G, b1AR_S464,b1AR_S301,...
%     GsaGTPtot, GsaGDP, Gsby, AC_GsaGTP ,cAMPtot, PDEp, ...
%   RC_I, RCcAMP_I, RCcAMPcAMP_I, RcAMPcAMP_I, PKACI ,PKACI_PKI, ...
%   RC_II, RCcAMP_II, RCcAMPcAMP_II, RcAMPcAMP_II, PKACII, PKACII_PKI, ...
%   I1p_PP1, I1ptot, LCCap ,LCCbp, PLBp, PLMp, TnIp, ...
%     m, h, jo,v, w, x, yo, z,rto, sto ,ssto,rss ,...
%     sss,Ca_nsr ,Ca_jsr,Nai ,Ki ,Cai,Vm,trelo] = yCell{:};
%  Rtot = sum(algvarsGly,2) + b1AR_S464 + b1AR_S301 ;
%  index = Rtot> .0132+1e-4;algvarsGly(index,:) = NaN;talgvars = tGly; talgvars(index) = NaN;
%  algvarsCell=mat2cell(algvarsGly,size(algvarsGly,1),ones(size(algvarsGly,2),1));
% [Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};
% 
%  subplot(3,3,1);plot(talgvars(~isnan(talgvars)),Ra(~isnan(Ra)));title('Ra');axis tight;hold all;
% subplot(3,3,2);plot(talgvars(~isnan(talgvars)),LRa(~isnan(LRa)));title('LRa');axis tight;hold all;
% subplot(3,3,3);plot(talgvars(~isnan(talgvars)),RaG(~isnan(RaG)));title('RaG');axis tight;hold all;
% subplot(3,3,4);plot(talgvars(~isnan(talgvars)),LRaG(~isnan(LRaG)));title('LRaG');axis tight;hold all;
% subplot(3,3,5);plot(talgvars(~isnan(talgvars)),ARi(~isnan(LRaG)));title('ARi');axis tight;hold all;
% subplot(3,3,6);plot(tGly,100*(b1AR_S464+b1AR_S301)/0.0132);box off;axis tight;
% subplot(3,3,7);plot(tGly,cAMPtot);title('cAMPtot');axis tight;hold all;
% subplot(3,3,8);plot(tGly,PLBp);title('PLBp');axis tight;hold all;
% subplot(3,3,9);plot(tGly,Cai*1000);box off;axis tight;
% xlabel('time (mins)');ylabel('Ca^2+ (\muM)');hold all;