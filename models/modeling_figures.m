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
tic
 clc; clear all;close all;
load fitalpha.mat;

KR =  10;
Ki = [0.224,NaN,3.97,3.57,13.6,1.71,2.44,31.3,1.31,1.6,37.9,4.70e-03,...
      5.80e-03,2.60e-03,1.64e-02,3.88e-01,1.70e-03,4.70e-02,2.24e-02,...
      1.80e-03,4.50e-03,4.95e-02];
  KLcalc =Ki.*(fitalpha*KR + 1)./(fitalpha*(KR + 1));

fileName = 'finalfigures';
KL = KLcalc(1);

KG =  2.4131;
alpha_L = fitalpha(1);
gamma_L =  0.3762;
KA = 500e-6;
alpha_A = 1; 
gamma_A =  1;
PKAItot         = 0.59;           % (uM) total type 1 PKA
PKAIItot        = 0.025;          % (uM) total type 2 PKA
%% Receptor binding assays

KA = 95e-6;
alpha_A = 1;gamma_A = 1;
p = daePARAMS(KR,KL,KA,0.7,alpha_L,alpha_A,gamma_L,gamma_A);
tspan = [0;20*60*1000];  % 20 minutes

ISOrange = 10.^[-5:.1:2];
y0 = zeros(29,1); y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
RelTol = 1e-13;
MaxStep = 1e-1;
options = odeset('MaxStep',1e3,'NonNegative',[1:29],'RelTol',RelTol);
t = [0;1000*60*1000]; p(1:2) = 0;[~,y] = ode15s(@daeODE,t,y0,options,p);y0 = y(end,:);
p(2) = 40e-6;   % CYP concentration 40 pM is Mason/JBC1999
p(51) = 0; p(31:36) = 0;

for i=1:length(ISOrange)
    p(1) = ISOrange(i);  % Ltot
    p(30) = 3.83;    % Gtot = 3.83  
    [t,y] = ode15s(@daeODE,tspan,y0,[],p);
       for tstep=1:length(t),
        [~,algvars(tstep,:)]=daeODE(t(tstep),y(tstep,:),p);
       end
    algvarsCell=mat2cell(algvars,size(algvars,1),ones(size(algvars,2),1));
    [Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};
    CYPboundwoGPP(i)= ARi(end)+ARa(end)+ARaG(end);
    clear algvars algvarscell;
      p(30) = 0;    % Gtot = 0  
    [t,y] = ode15s(@daeODE,tspan,y0,[],p);
    
    for tstep=1:length(t),
        [~,algvars(tstep,:)]=daeODE(t(tstep),y(tstep,:),p);
    end
    algvarsCell=mat2cell(algvars,size(algvars,1),ones(size(algvars,2),1));
    [Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};
    CYPboundwGPP(i)= ARi(end)+ARa(end)+ARaG(end);
    clear algvars algvarscell;
    
end

figure(1);
load('masonJBC1999.mat');
subplot(1,2,1);
semilogx(ISOrange*1e-6,CYPboundwoGPP/max(CYPboundwoGPP)*100,...
    ISOrange*1e-6,CYPboundwGPP/max(CYPboundwGPP)*100,'LineWidth',2);hold all;
errorbar(10.^masonArgnoGPP(:,1),masonArgnoGPP(:,2),masonArgnoGPP(:,3),'bo','MarkerSize',5,'LineWidth',2);
errorbar(10.^masonArgGPP(:,1),masonArgGPP(:,2),masonArgGPP(:,3),'ro','MarkerSize',5,'LineWidth',2); 
xlabel('isoproterenol (M)');
ylabel('CYP bound (%)'); title('\beta_1-AR Arg389');
axis tight; set(gca,'Xtick',[1e-010 1e-08,1e-06,0.0001,0.01;]);box off;

p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
tspan = [0;20*60*1000];  % 20 minutes
ISOrange = 10.^[-5:.1:2];
y0 = zeros(29,1); y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
options = odeset('MaxStep',1e3,'NonNegative',[1:29],'RelTol',RelTol);
t = [0;1000*60*1000]; p(1:2) = 0;[~,y] = ode15s(@daeODE,t,y0,options,p);y0 = y(end,:);
p(2) = 40e-6;   % CYP concentration 40 pM is Mason/JBC1999
p(51) = 0; p(31:36) = 0;

for i=1:length(ISOrange)
    p(1) = ISOrange(i);  % Ltot
    p(30) = 3.83;    % Gtot = 3.83  
    [t,y] = ode15s(@daeODE,tspan,y0,[],p);
       for tstep=1:length(t),
        [~,algvars(tstep,:)]=daeODE(t(tstep),y(tstep,:),p);
       end
    algvarsCell=mat2cell(algvars,size(algvars,1),ones(size(algvars,2),1));
    [Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};
    CYPboundwoGPP(i)= ARi(end)+ARa(end)+ARaG(end);
    clear algvars algvarscell;
      p(30) = 0;    % Gtot = 0  
    [t,y] = ode15s(@daeODE,tspan,y0,[],p);
    
    for tstep=1:length(t),
        [~,algvars(tstep,:)]=daeODE(t(tstep),y(tstep,:),p);
    end
    algvarsCell=mat2cell(algvars,size(algvars,1),ones(size(algvars,2),1));
    [Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};
    CYPboundwGPP(i)= ARi(end)+ARa(end)+ARaG(end);
    clear algvars algvarscell;
    
end
subplot(1,2,2);
semilogx(ISOrange*1e-6,CYPboundwoGPP/max(CYPboundwoGPP)*100,...
    ISOrange*1e-6,CYPboundwGPP/max(CYPboundwGPP)*100,'LineWidth',2);hold all;
errorbar(10.^masonGlynoGPP(:,1),masonGlynoGPP(:,2),masonGlynoGPP(:,3),'bo','MarkerSize',5,'LineWidth',2);
errorbar(10.^masonGlyGPP(:,1),masonGlyGPP(:,2),masonGlyGPP(:,3),'ro','MarkerSize',5,'LineWidth',2); 
xlabel('isoproterenol (M)');
ylabel('CYP bound (%)'); title('\beta_1-AR Gly389');axis tight;
 set(gca,'Xtick',[1e-010 1e-08,1e-06,0.0001,0.01;]);box off

 set(gcf, 'PaperPositionMode', 'manual');
 set(gcf, 'PaperUnits', 'inches');
 set(gcf, 'PaperPosition', [2.0 4.5 5.5  2.0]);
cd('A:\Robert\Dissertation\Beta Blocker modeling\final_model validation\figures');
print -dpdf receptorva


%% Adenylyl Cyclase Assay- Iso dose response- compared with Rathz JBC 2003
cd('A:\Robert\Dissertation\Beta Blocker modeling\final_model validation');
KA = 95e-6;
alpha_A = 1;gamma_A = 1;
tspan = [0;20*60*1000];  % 20 minutes
y0 = zeros(29,1); y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
RelTol = 1e-13;MaxStep = 1e-1;
options = odeset('MaxStep',1e3,'NonNegative',[1:29],'RelTol',RelTol);
pArg = daePARAMS(KR,KL,KA,0.7,alpha_L,alpha_A,gamma_L,gamma_A);
pGly = daePARAMS(KR,KL,KA, 2.4131,alpha_L,alpha_A,gamma_L,gamma_A);
pArg(1:2) = 0;pGly(1:2) = 0;
t = [0;1000*60*1000]; [~,y] = ode15s(@daeODE,t,y0,options,pArg);y0 = y(end,:);
ISOrange = 10.^[-5:.3:2];
for i=1:length(ISOrange)
    pArg(1) = ISOrange(i);  % Ltot
     pGly(1) = ISOrange(i);  % Ltot
    [t,y] = ode15s(@daeODE,tspan,y0,[],pArg);
    cAMParg(i)=y(end,9);
    [t,y] = ode15s(@daeODE,tspan,y0,[],pGly);
    cAMPgly(i)=y(end,9);
end
figure(2);
semilogx(ISOrange*1e-6,cAMParg/max(cAMPgly)*100,'color','r','LineWidth',2); 
xlabel('isoproterenol (\muM)');ylabel('AC activity (% max of Gly)');hold on;
semilogx(ISOrange*1e-6,cAMPgly/max(cAMPgly)*100,'color','b','LineWidth',2);axis tight;
% Plot experimental data from Rathz,Ligget JBC 2003
load('rathz_ACactivity_fig2a.mat');
errorbar(iso*1e-6,gly_cAMP/max(gly_cAMP)*100,(arg_error-arg_cAMP)/max(gly_cAMP)*100,'o','color','b','LineWidth',2);
  errorbar(iso*1e-6,arg_cAMP/max(gly_cAMP)*100,(gly_error-arg_cAMP)/max(gly_cAMP)*100,'o','color','r','LineWidth',2);
  errorbarlogx;
ylabel('AC activity (% max of Gly)');
xlabel('isoproterenol (M)');axis tight;

 set(gca,'Xtick',[1e-010 1e-08,1e-06,0.0001,0.01;]);box off
  set(gcf, 'PaperPositionMode', 'manual');
 set(gcf, 'PaperUnits', 'inches');
 set(gcf, 'PaperPosition', [2.0 4.5 3.1  2.7]);
cd('A:\Robert\Dissertation\Beta Blocker modeling\final_model validation\figures');
print -dpdf ACvalid
%% cAMP validation
cd('A:\Robert\Dissertation\Beta Blocker modeling\final_model validation');
KLcalc = (0.2/KLcalc(1))*KLcalc; 
KL = KLcalc(1);
alpha_L = fitalpha(1);% 0.3037
KA = KLcalc(20);%2e-3;
alpha_A =  fitalpha(20);
 p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);p(1) = 0;
 y0 = zeros(29,1);y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
 [t,y] = ode15s(@daeODE,[0;60*60*1000],y0,[],p);
 y0 = y(end,:);plot(t,y);
tspan = [0;20*60*1000];  % 20 minutes

 p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);p(1) = 0;
 p(1) = .01;   
 [t,y] = ode15s(@daeODE,tspan,y0,[],p);
 cAMPlow = y(:,9);tlow = t;
 
figure(3);plot(tlow./60e3,cAMPlow/.168,'LineWidth',2);%conversion from Bers Table 9...for Ca2+ fluxes
load VillaPetrofftl.mat;hold all;
errorbar(VillaPetrofftl(:,1),VillaPetrofftl(:,2),VillaPetrofftl(:,3),'o','LineWidth',2);xlim([0 20]);ylim([4 20]);
xlabel('time (mins)');ylabel('cAMP (pmol/mg)');box off;
 set(gcf, 'PaperPositionMode', 'manual');
 set(gcf, 'PaperUnits', 'inches');
 set(gcf, 'PaperPosition', [2.0 4.5 3.8  2.0]);
cd('A:\Robert\Dissertation\Beta Blocker modeling\final_model validation\figures');
print -dpdf cAMPtlvalid
cd('A:\Robert\Dissertation\Beta Blocker modeling\final_model validation');

% save ICscAMP y0;

ISOrange = 10.^[-5:.3:2];

p(1) = 0;
for i=1:length(ISOrange)
    p(1) = ISOrange(i);  % Ltot 
    [t,y] = ode15s(@daeODE,tspan,y0,[],p);
    cAMP(i) = (y(end,9));
         
end

figure(4);semilogx(ISOrange*1e-6,cAMP/.168,'LineWidth',2);hold all;
load('VillaPetroffdr');xlabel('isoproterenol (M)');ylabel('cAMP (pmol/mg)');
errorbar(10.^VillaPetroffdr(:,1),VillaPetroffdr(:,2),VillaPetroffdr(:,3),'o','LineWidth',2);
 set(gca,'Xtick',[1e-010 1e-08,1e-06,0.0001,0.01;]);box off;axis tight
 
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2.0 4.5 2.8  2.3]);
cd('A:\Robert\Dissertation\Beta Blocker modeling\final_model validation\figures');
print -dpdf cAMPdrvalid

%% PLB validation
cd('A:\Robert\Dissertation\Beta Blocker modeling\final_model validation');

KL = KLcalc(1);
alpha_L = fitalpha(1);% 0.3037
KA = KLcalc(20);%2e-3;
alpha_A =  fitalpha(20);
 p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);p(1) = 0;
 y0 = zeros(29,1);y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
[~,y] = ode15s(@daeODE,[0;60*60*1000],y0,[],p);y0 = y(end,:);
tspan = [0;20*60*1000];  % 20 minutes

ISOrange = 10.^[-4:.1:1];
for i=1:length(ISOrange)
    p(1) = ISOrange(i);  % Ltot 
    p(2) = 1e-20;
    [t,y] = ode15s(@daeODE,tspan,y0,[],p);
    PLBp(i) = y(end,27);
end

figure(5);semilogx(ISOrange*1e-6,100*PLBp./max(PLBp),'LineWidth',2);hold all;
load Vittone.mat;
errorbar(10.^Vittone(:,1),Vittone(:,2),Vittone(:,3),'o','LineWidth',2);
xlabel('isoproterenol (M)');ylabel('PLBp (% of max)');
 set(gca,'Xtick',[1e-010 1e-08,1e-06,0.0001,0.01;]);box off; axis tight;
  set(gcf, 'PaperPositionMode', 'manual');
 set(gcf, 'PaperUnits', 'inches');
 set(gcf, 'PaperPosition', [2.0 4.5 2.8  2.3]);
cd('A:\Robert\Dissertation\Beta Blocker modeling\final_model validation\figures');
print -dpdf PLBvalid

%% Beta Blocker comparisons calcium

 p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);p(1) = 0;
 y0 = zeros(29,1);y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
 [~,y] = ode15s(@daeODE,[0;60*60*1000],y0,[],p);y0Sig = y(end,:);


options = odeset('RelTol',1e-5,'MaxStep',5e-3,'Stats','on');
p = daePARAMSc(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
load y0full;y0 = y0full;


tspan = [0;2*60];
[t,y] = ode15s(@daeODEc,tspan,y0,options,p);
bgcAMP{1,1} = y(:,9);
bgCa{1,1} = y(:,47);
time{1,1} = t;
yfinal = y(end,:); 

tspan = [2*60;6*60];
p(1) = 0.1; p(2) = 0;  % agonist conc and beta blocker conc
[t2,y2] = ode15s(@daeODEc,tspan,yfinal,options,p); 
bgcAMP{1,2} = y2(:,9);
bgCa{1,2} = y2(:,47);
time{1,2} = t2;
yfinal = y2(end,:); 

tspan2 = [6*60 10*60]; %changed the time
p(1) = 10; p(2) = 0;     % agonist conc and beta blocker conc
[t3,y3]=ode15s(@daeODEc,tspan2,yfinal,options,p);
bgcAMP{1,3} = y3(:,9);
bgCa{1,3} = y3(:,47);
time{1,3} = t3;

dose = [0.1 1 1];drugs = [20 18 17];

for i=1:length(drugs)
   
  KL = KLcalc(1);
   alpha_L = fitalpha(1);
   KA = KLcalc(drugs(i)); alpha_A = fitalpha(drugs(i));
   p = daePARAMSc(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
    
    tspan = [0;2*60];
p(1) = 0; p(2) = 0;  % agonist conc and beta blocker conc
[t,y] = ode15s(@daeODEc,tspan,y0,options,p);

bgcAMP{i+1,1} = y(:,9);
bgCa{i+1,1} = y(:,47);
time{i+1,1} = t;
yfinal = y(end,:); 

tspan = [2*60;6*60];
p(1) = 0.1; p(2) = dose(i);  % agonist conc and beta blocker conc
[t2,y2] = ode15s(@daeODEc,tspan,yfinal,options,p); 
bgcAMP{i+1,2} = y2(:,9);
bgCa{i+1,2} = y2(:,47);
time{i+1,2} = t2;
yfinal = y2(end,:); 
tspan2 = [6*60 10*60]; %changed the time
p(1) = 10; p(2) = dose(i);     % agonist conc and beta blocker conc
[t3,y3]=ode15s(@daeODEc,tspan2,yfinal,options,p);
bgcAMP{i+1,3} = y3(:,9);
bgCa{i+1,3} = y3(:,47);
time{i+1,3} = t3;
end

ISO = [bgCa{1,1}; bgCa{1,2} ;bgCa{1,3}];ISOt = [time{1,1}; time{1,2} ;time{1,3}]./60;
PRO= [bgCa{2,1}; bgCa{2,2} ;bgCa{2,3}];PROt = [time{2,1}; time{2,2} ;time{2,3}]./60;
MET = [bgCa{3,1}; bgCa{3,2} ;bgCa{3,3}];METt = [time{3,1}; time{3,2} ;time{3,3}]./60;
CAR = [bgCa{4,1}; bgCa{4,2} ;bgCa{4,3}];CARt = [time{4,1}; time{4,2} ;time{4,3}]./60;

figure(6);plot(ISOt,ISO,'LineWidth',2);hold all;plot(METt,MET,'LineWidth',2);
xlabel('time (mins)');ylabel('Ca^{2+} (\muM)');box off;
set(gcf, 'PaperPositionMode', 'manual');
 set(gcf, 'PaperUnits', 'inches');
 set(gcf, 'PaperPosition', [2.0 4.5 3.8  2.0]);
cd('A:\Robert\Dissertation\Beta Blocker modeling\final_model validation\figures');
print -dpdf METsens
figure(7);
plot(ISOt(ISOt<=0.0184),ISO(ISOt<=0.0184),'LineWidth',2);hold all;
plot(METt(METt>=9.9816)-9.9816,MET(METt>=9.9816),'LineWidth',2);axis tight;box off;
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2.0 4.5 2.8  2.3]);
print -dpdf METsensexample
% figure(6);plot(ISOt,ISO);hold all;plot(PROt,PRO);plot(METt,MET);plot(CARt,CAR);
% ISO = (ISO - min(ISO))./min(ISO);[ISO,locsISO] = findpeaks(ISO);ISO = ISO./ISO(1);
% PRO = (PRO - min(PRO))./min(PRO);[PRO,locsPRO] = findpeaks(PRO);PRO = PRO./PRO(1);
% MET = (MET - min(MET))./min(MET);[MET, locsMET]= findpeaks(MET);MET = MET./MET(1);
% CAR  = (CAR  - min(CAR))./min(CAR);[CAR,locsCAR] = findpeaks(CAR);CAR = CAR./CAR(1);
% ISOt =ISOt(locsISO);PROt =PROt(locsPRO);METt =METt(locsMET);CARt =CARt(locsCAR);
% xlabel('time (mins)');ylabel('Ca^{2+} (\muM)');legend('ISO','ISO+PRO','ISO+MET','ISO+CAR');
% 
% figure(7);plot(ISOt,ISO);hold all;plot(PROt,PRO);plot(METt,MET);plot(CARt,CAR);legend('ISO','ISO+PRO','ISO+MET','ISO+CAR');
% xlabel('time (mins)');ylabel('Ca^{2+} amplitude (fold change)');export_fig(fileName,'-pdf','-append');
% bg(1,1) = max(ISO(ISOt <2));bg(1,2) = max(ISO(ISOt >5.5 & ISOt <6 ));bg(1,3) = max(ISO(ISOt >6));
% bg(2,1) = max(PRO(PROt <2));bg(2,2) = max(PRO(PROt >5.5 & PROt <6 ));bg(2,3) = max(PRO(PROt >6));
% bg(3,1) = max(MET(METt <2));bg(3,2) = max(MET(METt >5.5 & METt <6 ));bg(3,3) = max(MET(METt >6));
% bg(4,1) = max(CAR(CARt <2));bg(4,2) = max(CAR(CARt >5.5 & CARt <6 ));bg(4,3) = max(CAR(CARt >6));figure(8);bar(bg);
% set(gca,'XTickLabel',{'ISO','ISO+PRO','ISO+MET','ISO+CAR'}, 'XTick',[1 2 3 4]);legend('ctrl','0.1 \muM ISO','10 \muM ISO');

% cd('A:\Robert\Dissertation\Beta Blocker modeling\final_model validation\figures');
%% polymorphisms simulation
KL = KLcalc(1);
alpha_L = fitalpha(1);
KA = KLcalc(20);
alpha_A = fitalpha(20);
p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
RelTol = 1e-13;
MaxStep = 1e-1;
options = odeset('MaxStep',1e3,'NonNegative',[1:29],'RelTol',RelTol);
 y0 = zeros(29,1); y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
t = [0;1000*60*1000]; p(1:2) = 0;[~,y] = ode15s(@daeODE,t,y0,options,p);y0 = y(end,:);
tspan1 = [0;2*60*1000]; 
p(1) = 0;p(2) = 0;
[t1,y1] = ode15s(@daeODE,tspan1,y0,options,p);

options = odeset('MaxStep',MaxStep,'NonNegative',[1:29],'RelTol',RelTol);
tspan = [2*60*1000; 2.1*60*1000]; y0 = y1(end,:);
p(1) = 0.1;p(2) = 0.1;
[t2,y2] = ode15s(@daeODE,tspan,y0,[],p);
options = odeset('MaxStep',MaxStep,'NonNegative',[1:29],'RelTol',RelTol);
tspan = [2.1*60*1000; 6*60*1000]; y0 = y2(end,:);
p(1) = 0.1;p(2) = 0.1;
[t3,y3] = ode15s(@daeODE,tspan,y0,[],p);
tspan = [6*60*1000; 10*60*1000]; y0 = y3(end,:);
p(1) = 10;p(2) = 0.1;
[t4,y4] = ode15s(@daeODE,tspan,y0,[],p);

t = [t1;t2;t3;t4]./60e3;
 y = [y1;y2;y3;y4];

figure(9);plot(t,y(:,9),'LineWidth',2);hold all;


KG =  0.7;gamma_L =  0.3762;
KL = KLcalc(1);
alpha_L = fitalpha(1);
KA = KLcalc(20);
alpha_A = fitalpha(20);
p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
RelTol = 1e-13;
MaxStep = 1e-1;
options = odeset('MaxStep',1e3,'NonNegative',[1:29],'RelTol',RelTol);
 y0 = zeros(29,1); y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
t = [0;1000*60*1000]; p(1:2) = 0;[~,y] = ode15s(@daeODE,t,y0,options,p);y0 = y(end,:);
tspan1 = [0;2*60*1000]; 
p(1) = 0;p(2) = 0;
[t1,y1] = ode15s(@daeODE,tspan1,y0,options,p);

options = odeset('MaxStep',MaxStep,'NonNegative',[1:29],'RelTol',RelTol);
tspan = [2*60*1000; 2.1*60*1000]; y0 = y1(end,:);
p(1) = 0.1;p(2) = 0.1;
[t2,y2] = ode15s(@daeODE,tspan,y0,[],p);
options = odeset('MaxStep',MaxStep,'NonNegative',[1:29],'RelTol',RelTol);
tspan = [2.1*60*1000; 6*60*1000]; y0 = y2(end,:);
p(1) = 0.1;p(2) = 0.1;
[t3,y3] = ode15s(@daeODE,tspan,y0,[],p);
tspan = [6*60*1000; 10*60*1000]; y0 = y3(end,:);
p(1) = 10;p(2) = 0.1;
[t4,y4] = ode15s(@daeODE,tspan,y0,[],p);

t = [t1;t2;t3;t4]./60e3;
 y = [y1;y2;y3;y4];

plot(t,y(:,9),'LineWidth',2);box off;axis tight;ylim([0.3 4.2]);
xlabel('time (mins)');ylabel('cAMP (\muM)');
set(gcf, 'PaperPositionMode', 'manual');
 set(gcf, 'PaperUnits', 'inches');
 set(gcf, 'PaperPosition', [2.0 4.5 5.5  2]);
cd('A:\Robert\Dissertation\Beta Blocker modeling\final_model validation\figures');
print -dpdf polymorphismsSim


%%
load fitalpha.mat;load Ki;
 drugList = {'Isoproterenol'; 'Basal' ;'Epinephrine'; 'Norepinephrine'; 'Fenoterol' ;'Formoterol';
'Salbutamol' ;'Terbutaline' ;'Broxaterol' ;'Salmeterol'; 'BRL-37344'; 'CGP-12177';
'Alprenolol' ;'Pindolol'; 'SR 59230A'; 'Atenolol' ;'Carvedilol' ;'Metoprolol'; 'Bisoprolol';
'Propranolol' ;'CGP-20712' ;'ICI-118551'};

p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);

RelTol = 1e-13;MaxStep = 1e-1;
options = odeset('MaxStep',1e3,'NonNegative',[1:29],'RelTol',RelTol);
 y0 = zeros(29,1); y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
t = [0;1000*60*1000]; p(1:2) = 0;[~,y] = ode15s(@daeODE,t,y0,options,p);y0 = y(end,:);
tspan = [0;2*60*1000]; 
drugs = 1:22;
dose = ones(size(drugs))*1e-0;
for i=1:length(drugs)
   
  KL =  KLcalc(1);
   alpha_L = fitalpha(1);
   KA = KLcalc(drugs(i)); alpha_A = fitalpha(drugs(i));
  KAdrugs(i) = KLcalc(drugs(i)); alphadrugs(i) = fitalpha(drugs(i));
  p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
    tspan = [0;2*60*1000];
p(1) = 0; p(2) = 0;  % agonist conc and beta blocker conc
[t1,y] = ode15s(@daeODE,tspan,y0,options,p);
bg(i,1) = max(y(end,9));
yfinal = y(end,:); 
tspan = [2*60*1000;6*60*1000];
p(1) = 0.1; p(2) = dose(i);  % agonist conc and beta blocker conc
[t2,y2] = ode15s(@daeODE,tspan,yfinal,options,p); 
bg(i,2) = max(y2(end,9));
yfinal = y2(end,:); 

tspan2 = [6*60*1000 10*60*1000]; %changed the time
p(1) = 10; p(2) = dose(i);     % agonist conc and beta blocker conc
[t3,y3]=ode15s(@daeODE,tspan2,yfinal,options,p);
bg(i,3) = max(y3(end,9));

end
y = [y;y2;y3];t = [t1;t2;t3];
for tstep=1:length(t),
    [~,algvars(tstep,:)]=daeODE(t(tstep),y(tstep,:),p);
end
algvarsCell=mat2cell(algvars,size(algvars,1),ones(size(algvars,2),1));
[Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};

figure(12);bar(dose);ylabel('dose(\muM)');
set(gca,'XTick',1:length(drugList),'XTickLabel',drugList);xticklabel_rotate([],45);

sensitivity = (bg(:,3)- bg(:,2));figure(13); bar(bg);set(gca,'XTick',1:length(drugList),'XTickLabel',drugList);xticklabel_rotate([],45);
toc

sensitivity(2:4) = [];KAdrugs(2:4) = [];alphadrugs(2:4) = [];drugList(2:4) = [];
figure(14);subplot(2,2,1:2); bar(sensitivity);xlabel('ligand');ylabel('sens.(cAMP \muM)');box off;axis tight;
set(gca,'XTick',1:length(drugList),'XTickLabel',drugList);xticklabel_rotate([],45);
subplot(2,2,3);scatter(KAdrugs*1e-6,sensitivity);ylabel('sens.(cAMP \muM)');xlabel('K_L (M)');
subplot(2,2,4);scatter(alphadrugs,sensitivity);ylabel('sens.(cAMP \muM)');xlabel('\alpha');
