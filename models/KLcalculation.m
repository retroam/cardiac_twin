 clear all;clc;
load fitalpha;
KR =  10;
gamma = 1;
KG = 0.9041;
Gstot = 3.83;
Ki = [0.224,NaN,3.97,3.57,13.6,1.71,2.44,31.3,1.31,1.6,37.9,4.70e-03,...
      5.80e-03,2.60e-03,1.64e-02,3.88e-01,1.70e-03,4.70e-02,2.24e-02,...
      1.80e-03,4.50e-03,4.95e-02];
% KLcalc(:,1) = Ki.*(fitalpha*KR + 1)./(fitalpha*(KR + 1));
% KLcalc(:,2) = Ki.*(fitalpha*KR + 1)./fitalpha/(KR + 1);
% KLcalc(:,3) = (gamma*KG.*(fitalpha*KR + Ki) + Ki*Gstot)./((gamma*KG)*(KR + Ki).*fitalpha);
% KLcalc(:,4) = (gamma*KG.*(fitalpha*KR + Ki) + Ki*Gstot)./fitalpha./((gamma*KG)*(KR + Ki));
% KLcalc(:,5) = (gamma*KG*Ki.*(fitalpha*KR + 1) + Ki*Gstot)./(gamma*(fitalpha.*KG*KR + Gstot));
factor1 =1;factor2 =1;
 fitalpha(1) =  1/32;
KLcalc = Ki./factor1.*(fitalpha./factor2*KR + 1)./(fitalpha./factor2*(KR + 1));
% KLcalc(1) = 0.2;
%  KLcalc = Ki*(10/11);
% ISO = KLcalc(1,2)
% PRO = KLcalc(20,2)
% KLcalc = (gamma*KG*Ki.*(fitalpha*KR + 1) + Ki*Gstot)./(fitalpha./gamma*(KG*KR + KG + Gstot));
%  KLcalc = (0.2/KLcalc(1))*KLcalc; 
% save KLcalc KLcalc; 
 %%
 L = 1e-3:1e-2:1;
 alpha = fitalpha(1);
 Ri = 0.0132/2;
 KL = 0.2;
 Gstot = 3.83;
 analyticala = (L*Ri/alpha*KL*KR + L*Ri/KL + L*Ri*Gstot/alpha*gamma*KL*KG*KR)./...
    (Ri + Ri/KR + Ri*Gstot/KG*KR + L*Ri/alpha*KL*KR + L*Ri/KL + L*Ri*Gstot/alpha*gamma*KL*KG*KR);
 analyticalb = (L*Ri/alpha*KL*KR + L*Ri/KL)./...
    (Ri + Ri/KR  + L*Ri/alpha*KL*KR + L*Ri/KL);
analyticalc = (gamma*KG*L*(alpha*KR + 1) + L*Gstot)/(alpha/gamma*(KG*KR + KG + Gstot));
analyticald = L.*(alpha*KR + 1)./alpha/(KR + 1);
subplot(2,2,1);semilogx(L,[analyticala; analyticalb]);axis tight;box off;
subplot(2,2,2);semilogx(L,[analyticalc; analyticald]);axis tight; box off;
Gstot = 0;
 analyticala = (L*Ri/alpha*KL*KR + L*Ri/KL + L*Ri*Gstot/alpha*gamma*KL*KG*KR)./...
    (Ri + Ri/KR + Ri*Gstot/KG*KR + L*Ri/alpha*KL*KR + L*Ri/KL + L*Ri*Gstot/alpha*gamma*KL*KG*KR);
 analyticalb = (L*Ri/alpha*KL*KR + L*Ri/KL)./...
    (Ri + Ri/KR  + L*Ri/alpha*KL*KR + L*Ri/KL);
analyticalc = (gamma*KG*L*(alpha*KR + 1) + L*Gstot)/(alpha/gamma*(KG*KR + KG + Gstot));
subplot(2,2,3);semilogx(L,[analyticala; analyticalb]);axis tight;box off;
subplot(2,2,4);semilogx(L,[analyticalc; analyticald]);axis tight; box off;

%% Using model equations
close all;
KL =0.2;KR =  10; KG = 0.9041;alpha_L = 1/32;gamma_L = 1.3; 
KA = 500e-6;alpha_A = 1; gamma_A =  1; 
p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
y0 = zeros(2,1);p(30)= 0;
[t,y] = ode15s(@dae_receptorODE,[0;60*60*1000],y0,[],p);plot(t,y);

for tstep=1:length(t),
    [~,algvars(tstep,:)]=dae_receptorODE(t(tstep),y(tstep,:),p);
end

yCell=mat2cell(y,size(y,1),ones(size(y,2),1));
algvarsCell=mat2cell(algvars,size(algvars,1),ones(size(algvars,2),1));
[Ri,G] = yCell{:};
[Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};

Rtot = Ri + Ra + LRi + LRa + RaG + LRaG + ARi + ARa + ARaG ;
figure(1)
subplot(2,2,1);plot(t,Rtot);title('Rtot');axis tight;hold all;
subplot(2,2,2);plot(t,Ri);title('Ri');axis tight;hold all;
subplot(2,2,3);plot(t,Ra);title('Ra');axis tight;hold all;
subplot(2,2,4);plot(t,RaG);title('RaG');axis tight;hold all;


Ri = Ri(end);
 L = 1e-3:1e-2:1;
 figure(2);
 analyticala = (L*Ri/alpha*KL*KR + L*Ri/KL + L*Ri*Gstot/alpha*gamma*KL*KG*KR)./...
    (Ri + Ri/KR + Ri*Gstot/KG*KR + L*Ri/alpha*KL*KR + L*Ri/KL + L*Ri*Gstot/alpha*gamma*KL*KG*KR);
 analyticalb = (L*Ri/alpha*KL*KR + L*Ri/KL)./...
    (Ri + Ri/KR  + L*Ri/alpha*KL*KR + L*Ri/KL);
analyticalc = (gamma*KG*L*(alpha*KR + 1) + L*Gstot)/(alpha/gamma*(KG*KR + KG + Gstot));
analyticald = L.*(alpha*KR + 1)./alpha/(KR + 1);
subplot(2,2,1);semilogx(L,[analyticala; analyticalb]);axis tight;box off;
subplot(2,2,2);semilogx(L,[analyticalc; analyticald]);axis tight; box off;
Gstot = 0;
 analyticala = (L*Ri/alpha*KL*KR + L*Ri/KL + L*Ri*Gstot/alpha*gamma*KL*KG*KR)./...
    (Ri + Ri/KR + Ri*Gstot/KG*KR + L*Ri/alpha*KL*KR + L*Ri/KL + L*Ri*Gstot/alpha*gamma*KL*KG*KR);
 analyticalb = (L*Ri/alpha*KL*KR + L*Ri/KL)./...
    (Ri + Ri/KR  + L*Ri/alpha*KL*KR + L*Ri/KL);
analyticalc = (gamma*KG*L*(alpha*KR + 1) + L*Gstot)/(alpha/gamma*(KG*KR + KG + Gstot));
subplot(2,2,3);semilogx(L,[analyticala; analyticalb]);axis tight;box off;
subplot(2,2,4);semilogx(L,[analyticalc; analyticald]);axis tight; box off;

%% dose response

clear all;
KL =0.2;KR =  10; KG = 0.9041;alpha_L = 1/32;gamma_L = 1.3; 
KA = 500e-6;alpha_A = 1; gamma_A =  1; 
 L = 1e-3:1e-3:1;Gstot = [3.83 0];
 p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
for i =1:2 
    for j = 1:length(L)
        y0 = zeros(2,1);p(1) = L(j);p(30)= Gstot(i);
        [t,y] = ode15s(@dae_receptorODE,[0;60*60*1000],y0,[],p);
       for tstep=1:length(t),
           [~,algvars(tstep,:)]=dae_receptorODE(t(tstep),y(tstep,:),p);
       end
       yCell=mat2cell(y,size(y,1),ones(size(y,2),1));
       algvarsCell=mat2cell(algvars,size(algvars,1),ones(size(algvars,2),1));
       [Ri,G] = yCell{:};
       [Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};
       A = Ri + Ra + LRi + LRa + RaG + LRaG + ARi + ARa + ARaG ;Rtot(i,j) = A(end);  
       B = LRi + LRa + LRaG + ARi + ARa + ARaG ;Lbound(i,j) = B(end); 
       C = ARi + ARa + ARaG ;Ctot(i,j) = C(end);
        clear algvars;clear y;
    end
end
%check if the dual receptor module changes anything
semilogx(L,Lbound./Rtot); hold all; 
Ri = 0.00866;
 analyticala = (L*Ri/alpha_L*KL*KR + L*Ri/KL + L*Ri*Gstot(1)/alpha_L*gamma_L*KL*KG*KR)./...
    (Ri + Ri/KR + Ri*Gstot(1)/KG*KR + L*Ri/alpha_L*KL*KR + L*Ri/KL + L*Ri*Gstot(1)/alpha_L*gamma_L*KL*KG*KR);
Ri = 0.012;
analyticalb = (L*Ri/alpha_L*KL*KR + L*Ri/KL)./...
    (Ri + Ri/KR  + L*Ri/alpha_L*KL*KR + L*Ri/KL);   
  semilogx(L,[analyticala; analyticalb],'--');  
  
  %%
  clear all;
KL =0.2;KR =  10; KG = 0.9041;alpha_L = 1/32;gamma_L = 1.3; 
KA = 500e-6;alpha_A = 1; gamma_A =  1; 
 L = 1e-3:1e-3:1;Gstot = [0 0];
 p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
for i =1:1 
    for j = 1:length(L)
        y0 = zeros(2,1);p(1) = L(j);p(30)= Gstot(i);
        [t,y] = ode15s(@dae_receptorODE,[0;60*60*1000],y0,[],p);
       for tstep=1:length(t),
           [~,algvars(tstep,:)]=dae_receptorODE(t(tstep),y(tstep,:),p);
       end
       yCell=mat2cell(y,size(y,1),ones(size(y,2),1));
       algvarsCell=mat2cell(algvars,size(algvars,1),ones(size(algvars,2),1));
       [Ri,G] = yCell{:};
       [Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};
       A = Ri + Ra + LRi + LRa + RaG + LRaG + ARi + ARa + ARaG ;Rtot(i,j) = A(end);  
       B = LRi + LRa + LRaG + ARi + ARa + ARaG ;Lbound(i,j) = B(end); 
       C = ARi + ARa + ARaG ;Ctot(i,j) = C(end);
        clear algvars;clear y;
    end
end
%check if the dual receptor module changes anything
 Ri = 0.00866;
semilogx(L,Lbound./Rtot); hold all;
analyticalb = (L*Ri/alpha_L*KL*KR)./(Ri + Ri/KR  + L*Ri/alpha_L*KL*KR); 
% analyticalb = (L*Ri/KL)./(Ri + Ri/KR  + L*Ri/KL);
  semilogx(L, analyticalb,'--'); 
%% Hoffmann Fig1
  clear all;
KL =0.2;KR =  10; KG = 0.9041;alpha_L = 1/32;gamma_L = 1.3; 
KA = 68e-6;alpha_A = 1; gamma_A =  1; 
L = 0:1e-6:400e-6;
p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
p(29) = 4e-4;
tspan = [0;20*60*1000];  % 20 minutes
y0 = zeros(29,1);      
for i=1:length(L)
    p(2) = L(i);  % Ltot
    p(30) = 3.83;    % Gtot = 3.83  
    [t,y] = ode15s(@daeODE,tspan,y0,[],p);
       for tstep=1:length(t),
        [~,algvars(tstep,:)]= daeODE(t(tstep),y(tstep,:),p);
       end
    algvarsCell=mat2cell(algvars,size(algvars,1),ones(size(algvars,2),1));
    [Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};
    CYPboundwoGPP(i)= ARi(end)+ARa(end)+ARaG(end);
    clear algvars algvarscell;
      p(30) = 0;    % Gtot = 0  
    [t,y] = ode15s(@daeODE,tspan,y0,[],p);
    
    for tstep=1:length(t),
        [~,algvars(tstep,:)] = daeODE(t(tstep),y(tstep,:),p);
    end
    algvarsCell=mat2cell(algvars,size(algvars,1),ones(size(algvars,2),1));
    [Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};
    CYPboundwGPP(i)= ARi(end)+ARa(end)+ARaG(end);
    clear algvars algvarscell;
    
end
% 
plot(L,CYPboundwGPP);hold all;
 plot(L,CYPboundwoGPP);hold all;
load('HoffmannFig1.dat'); plot(HoffmannFig1(:,1)*1e-6,HoffmannFig1(:,2)*1e-6,'o');
%% Hoffmann Fig3
clear all; clc;
load KLcalc.mat; load fitalpha.mat;
KR = 10;KG = .09041;
KL = 0.2;
alpha_L = fitalpha(1);gamma_L = 1.3;
KA = 95e-6;
alpha_A = 1;gamma_A = 1;
PKAItot = 0.59;PKAIItot = 0.025;  
tspan = [0;20*60*1000];  % 20 minutes
y0 = zeros(29,1); y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
RelTol = 1e-13;MaxStep = 1e-1;
options = odeset('MaxStep',1e3,'NonNegative',[1:29],'RelTol',RelTol);
p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
t = [0;1000*60*1000]; p(1:2) = 0;[~,y] = ode15s(@daeODE,t,y0,options,p);y0 = y(end,:);
p(29) = 8e-4;p(51) = 0; p(31:36) = 0;     % No desensitization or PDE's
ISOrange = 10.^[-3:.3:1];
for i=1:length(ISOrange)
    p(1) = ISOrange(i);  % Ltot
    [t,y] = ode15s(@daeODE,tspan,y0,[],p);
    cAMParg(i)=y(end,9);
end
semilogx(ISOrange,cAMParg/max(cAMParg)*100); xlabel('Isoproterenol (\muM)');
ylabel('Adenylyl Cyclase Activity (% max)');hold all;
load('HoffmannFig2.dat');semilogx(HoffmannFig2(:,1),HoffmannFig2(:,2),'o');
%% Effect of I-CYP on L binding
clear all;
ICYP = 10.^[-1:.5:6].*1e-6;
for k = 1:length(ICYP)
KL =0.2;KR =  10; KG = 0.9041;alpha_L = 1/32;gamma_L = 1.3; 
KA = 68e-6;alpha_A = 1; gamma_A =  1; 
L =  10.^[-4:.1:7];
p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
p(29) = 4e-4;
tspan = [0;20*60*1000];  % 20 minutes
y0 = zeros(29,1); %[t,y] = ode15s(@daeODE,[0;60*60*1000],y0,[],p);y0 = y(end,:);
p(2) = ICYP(k);
for i=1:length(L)
    p(1) = L(i);  % Ltot
    p(30) = 3.83;    % Gtot = 3.83  
    [t,y] = ode15s(@daeODE,tspan,y0,[],p);
       for tstep=1:length(t),
        [~,algvars(tstep,:)]= daeODE(t(tstep),y(tstep,:),p);
       end
    algvarsCell=mat2cell(algvars,size(algvars,1),ones(size(algvars,2),1));
    [Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};
    LboundwoGPP(i)= LRi(end)+LRa(end)+LRaG(end);
    clear algvars algvarscell;
      p(30) = 0;    % Gtot = 0  
    [t,y] = ode15s(@daeODE,tspan,y0,[],p);
    
    for tstep=1:length(t),
        [~,algvars(tstep,:)] = daeODE(t(tstep),y(tstep,:),p);
    end
    algvarsCell=mat2cell(algvars,size(algvars,1),ones(size(algvars,2),1));
    [Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};
    LboundwGPP(i)= LRi(end)+LRa(end)+LRaG(end);
    clear algvars algvarscell;
    
end
% 
%   semilogx(L,LboundwoGPP);hold all;
figure(1);semilogx(L,LboundwGPP);hold all; axis tight;

ic50i = find(LboundwGPP >max(LboundwGPP)/2);
maxgpp = LboundwGPP(ic50i(1)-1);ic50gpp(k) = L(ic50i(1)-1);diffgpp = maxgpp-max(LboundwGPP)/2;
ic50i = find(LboundwoGPP >max(LboundwoGPP)/2);
maxnogpp = LboundwoGPP(ic50i(1)-1);ic50nogpp(k) = L(ic50i(1)-1);diffnogpp = maxnogpp-max(LboundwoGPP)/2;

end
table = [ICYP.*1e6 ;ic50gpp ;ic50nogpp]';figure(2);semilogx(table(:,1),table(:,3));hold all;
clearvars -except ic50nogpp ic50gpp ICYP table

%% %% Effect of I-CYP on L binding
clear all;
ICYP = 10.^[2:.01:4].*1e-6;
for k = 1:length(ICYP)
KL =0.2;KR =  10; KG = 0.9041;alpha_L = 1/32;gamma_L = 1.3; 
KA = 68e-6;alpha_A = 1; gamma_A =  1; 
L =  10.^[-5:.03:6];
p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
p(29) = 4e-4;
tspan = [0;20*60*1000];  % 20 minutes
y0 = zeros(29,1); %[t,y] = ode15s(@daeODE,[0;60*60*1000],y0,[],p);y0 = y(end,:);
p(2) = ICYP(k);
for i=1:length(L)
    p(1) = L(i);  % Ltot
    p(30) = 3.83;    % Gtot = 3.83  
    [t,y] = ode15s(@daeODE,tspan,y0,[],p);
       for tstep=1:length(t),
        [~,algvars(tstep,:)]= daeODE(t(tstep),y(tstep,:),p);
       end
    algvarsCell=mat2cell(algvars,size(algvars,1),ones(size(algvars,2),1));
    [Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};
    LboundwoGPP(i)= ARi(end)+ARa(end)+ARaG(end);
    clear algvars algvarscell;
      p(30) = 0;    % Gtot = 0  y
      
    [t,y] = ode15s(@daeODE,tspan,y0,[],p);
    
    for tstep=1:length(t),
        [~,algvars(tstep,:)] = daeODE(t(tstep),y(tstep,:),p);
    end
    algvarsCell=mat2cell(algvars,size(algvars,1),ones(size(algvars,2),1));
    [Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};
    LboundwGPP1(i)= LRi(end)+LRa(end)+LRaG(end);
    LboundwGPP2(i)= ARi(end)+ARa(end)+ARaG(end);
    clear algvars algvarscell;
    
end
% 
%   semilogx(L,LboundwoGPP);hold all;
figure(1);semilogx(L,LboundwGPP1);hold all; axis tight;
ic502(k) = IC50(L,LboundwGPP2);

ic50i = find(LboundwGPP1 >max(LboundwGPP1)/2);
maxgpp = LboundwGPP1(ic50i(1)-1);ic501(k) = L(ic50i(1)-1);
% ic50i = find(LboundwGPP >max(LboundwGPP)/2);
% maxgpp = LboundwGPP(ic50i(end));ic50gpp(k) = L(ic50i(end));diffgpp = maxgpp-max(LboundwGPP)/2;
% ic50i = find(LboundwoGPP >max(LboundwoGPP)/2);
% maxnogpp = LboundwoGPP(ic50i(end));ic50nogpp(k) = L(ic50i(end));diffnogpp = maxnogpp-max(LboundwoGPP)/2;
end
table = [ICYP.*1e6 ;ic501;ic502]';figure(2);semilogx(table(:,1),table(:,2));hold all;semilogx(table(:,1),table(:,3));
% clearvars -except ic50nogpp ic50gpp ICYP table