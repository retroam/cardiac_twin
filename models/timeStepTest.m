close all; clear all;
tic;
load KLcalc.mat; load fitalpha.mat;
KR =  10;
Ki = [0.224,NaN,3.97,3.57,13.6,1.71,2.44,31.3,1.31,1.6,37.9,4.70e-03,...
      5.80e-03,2.60e-03,1.64e-02,3.88e-01,1.70e-03,4.70e-02,2.24e-02,...
      1.80e-03,4.50e-03,4.95e-02];
  KLcalc =Ki.*(fitalpha*KR + 1)./(fitalpha*(KR + 1));

fileName = 'KL0p2scaled';
KL = KLcalc(1);

KG =  2.4131;
alpha_L = fitalpha(1);
gamma_L =  0.3762;
KA =KLcalc(20);
alpha_A = fitalpha(20); 
gamma_A =  1;
PKAItot         = 0.59;           % (uM) total type 1 PKA
PKAIItot        = 0.025;          % (uM) total type 2 PKA

 p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);p(1) = 0;
 y0 = zeros(29,1);y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
 [~,y] = ode15s(@daeODE,[0;60*60*1000],y0,[],p);y0Sig = y(end,:);


options = odeset('RelTol',1e-5,'MaxStep',5e-3,'Stats','on');
p = daePARAMSc(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
load y0full;y0 = y0full;


tspan = [0;2*60];
[t1,y1] = ode15s(@daeODEc,tspan,y0,options,p);
yfinal = y1(end,:); 
for tstep=1:length(t1),
    [~,algvars1(tstep,:)]=daeODE(t1(tstep),y1(tstep,:),p);
end
options = odeset('RelTol',1e-6,'MaxStep',5e-4,'Stats','on');
tspan = [2*60;6*60];
p(1) = 0.1; p(2) = 0.1;  % agonist conc and beta blocker conc
[t2,y2] = ode15s(@daeODEc,tspan,yfinal,options,p); 
yfinal = y2(end,:); 
for tstep=1:length(t2),
    [~,algvars2(tstep,:)]=daeODE(t2(tstep),y2(tstep,:),p);
end
options = odeset('RelTol',1e-5,'MaxStep',5e-3,'Stats','on');
tspan2 = [6*60 10*60]; %changed the time
p(1) = 10; p(2) = 0.1;     % agonist conc and beta blocker conc
[t3,y3]=ode15s(@daeODEc,tspan2,yfinal,options,p);
for tstep=1:length(t3),
    [~,algvars3(tstep,:)]=daeODE(t3(tstep),y3(tstep,:),p);
end

t = [t1;t2;t3];y = [y1;y2;y3];

 algvars = [algvars1;algvars2;algvars3;];
yCell=mat2cell(y,size(y,1),ones(size(y,2),1));
algvarsCell=mat2cell(algvars,size(algvars,1),ones(size(algvars,2),1));
[Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};

[Ri,G, b1AR_S464,b1AR_S301,...
    GsaGTPtot, GsaGDP, Gsby, AC_GsaGTP ,cAMPtot, PDEp, ...
  RC_I, RCcAMP_I, RCcAMPcAMP_I, RcAMPcAMP_I, PKACI ,PKACI_PKI, ...
  RC_II, RCcAMP_II, RCcAMPcAMP_II, RcAMPcAMP_II, PKACII, PKACII_PKI, ...
  I1p_PP1, I1ptot, LCCap ,LCCbp, PLBp, PLMp, TnIp, ...
    m, h, jo,v, w, x, yo, z,rto, sto ,ssto,rss ,...
    sss,Ca_nsr ,Ca_jsr,Nai ,Ki ,Cai,Vm,trelo]=yCell{:};
toc

Rtot = Ri + Ra + LRi + LRa + RaG + LRaG + ARi + ARa + ARaG + b1AR_S464 + b1AR_S301 ;
figure(1)
subplot(5,3,1);plot(t,Rtot);title('Rtot');axis tight;hold all;
subplot(5,3,2);plot(t,Ri);title('Ri');axis tight;hold all;
subplot(5,3,3);plot(t,LRa);title('LRa');axis tight;hold all;
subplot(5,3,4);plot(t,LRi);title('LRi');axis tight;hold all;
subplot(5,3,5);plot(t,RaG);title('RaG');axis tight;hold all;
subplot(5,3,6);plot(t,LRaG);title('LRaG');axis tight;hold all;
subplot(5,3,7);plot(t,ARa);title('ARa');axis tight;hold all;
subplot(5,3,8);plot(t,ARi);title('ARi');axis tight;hold all;
subplot(5,3,9);plot(t,ARaG);title('ARaG');axis tight;hold all;
subplot(5,3,10);plot(t,GsaGTPtot);title('GsaGTPtot');axis tight;hold all;
subplot(5,3,11);plot(t,GsaGDP);title('GsaGDP');axis tight;hold all;
subplot(5,3,12);plot(t,Gsby);title('Gsby');axis tight;hold all;
subplot(5,3,13);plot(t, b1AR_S464);title('b1AR_{S464}');axis tight;hold all;
subplot(5,3,14);plot(t, b1AR_S301);title('b1AR_{S301}');axis tight;hold all;
subplot(5,3,15);plot(t,cAMPtot);title('cAMPtot');axis tight;hold all;