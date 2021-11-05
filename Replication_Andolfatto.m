clear all;close all;clc;


%% Calibration and steady state

% Targets
y = 1;
r = 1.01;
n = 0.57;
l = 0.33;

% Parameters exogenously set
delta = 0.025;
theta = 0.36;
sigma = 0.15;
alpha = 0.6;
eta=2;
rhoa = 0.95;
sdepsilona = 0.007;
e = l*(l/2);
q = 0.9;
psi = 0.01;

% Parameters calibrated
beta = 1/r;
k = theta/(r + delta - 1);
c = 1 - delta*k - psi;
m = sigma*n;
v = m/q;
kappa = psi/v;
mu = psi/(c*alpha*m);
p = 1 - q;
phi_1 = ((1-l)^eta)*(1-theta)*y/(n*l*c);
phi_2 = ((1-eta)/((1-e)^(1-eta)))*(-mu/beta + phi_1*((1-l)^(1-eta))/(1-eta) + (1-theta)*y/(n*c) + mu*(1 - sigma -(1-alpha)*p));
w = (1-theta - psi*(1 -(1-sigma)*beta)/(beta*sigma))/(n*l); 

param=[beta  k  c  m  v  kappa  mu  p  phi_1  phi_2  w];
disp('     beta        k         c       m          v        kappa      mu        p       phi_1       phi_2       w');
disp(param);




%% Linearized Model in Matrix Form

M1=zeros(6,6);
M2=zeros(6,5);
M3I=zeros(5,5);
M3L=zeros(5,5);
M4I=zeros(5,6);
M4L=zeros(5,6);
M5=zeros(5,1);

M1(1,1) = 1;
M1(1,3) = -(1 + eta*l/(1-l));
M1(2,2) = 1;
M1(2,5) = -1;
M1(3,1) = y*(1-alpha)*(1-theta)/(w*n*l);
M1(3,3) = (alpha*c*(1-l)^(-eta))/w -1;
M1(3,4) = (1-alpha)*kappa*v/(w*l*(1-n));
M1(3,6) = -1;
M1(4,1) = 1;
M1(4,3) = -(1-theta);
M1(5,2) = -alpha;
M1(5,5) = 1;
M1(6,4) = 1;
M1(6,5) = -1;

M2(1,1) = 1;
M2(1,4) = 1;
M2(2,1) = 1;
M2(2,3) = 1;
M2(3,1) = -1 + y*(1-alpha)*(1-theta)/(w*n*l);
M2(3,3) = -(1-alpha)*kappa*v/(w*l*(1-n));
M2(3,4) =y*(1-alpha)*(1-theta)/(w*n*l);
M2(4,2) = theta;
M2(4,4) = 1 - theta;
M2(4,5) = 1;
M2(5,4) = -(1-alpha)*n/(1-n);
M2(6,4) = n/(1-n);

M3I(1,1) = 1;
M3I(1,2) = beta*theta*y/k;
M3I(2,1) = phi_1*l/(mu*(1-l)^eta);
M3I(2,3) = -(1-sigma-(1-alpha)*p);
M3I(2,4) = phi_1*l/(mu*(1-l)^eta);
M3I(3,2) = 1;
M3I(4,4) = 1;
M3I(5,5) = 1;

M3L(1,1) = -1;
M3L(2,3) = 1/beta;
M3L(3,1) = c/k;
M3L(3,2) = -(1-delta);
M3L(4,4) = -(1-sigma);
M3L(5,5) = -rhoa;

M4I(1,1) = beta*theta*y/k;
M4I(2,1) = phi_1*l/(mu*(1-l)^eta);
M4I(2,3) = -phi_1*l/(mu*(1-l)^eta);
M4I(2,4) = -(1-alpha)*p;

M4L(3,1) = y/k;
M4L(3,2) = -kappa*v/k;
M4L(4,5) = m/n;

M5(5,1) = 1;

%% Solve the Model
% Reduced Form

L0=M3I-M4I*inv(M1)*M2;
L1=M4L*inv(M1)*M2-M3L;
L2=M5;

W=L0\L1;
Q=L0\L2;


%Eigenvalues


[PP,MU]=eig(W);


[llambda,kk]=sort(diag(MU));  

P1=PP(:,kk);
P1I=inv(P1);
MU1=P1I*W*P1;

[abf abf] = size(MU1);

ab=0;
for i = 1:abf;
if abs(MU1(i,i)) < 1
ab=ab+1;
end
end

af = abf- ab;

for i = 1:abf;
if abs(MU1(i,i)) == 1,
disp('Unit root')
else;
end;
end;

disp('backward    ')
disp(ab)
disp('forward')
disp(af)

disp('Eigenvalues')
disp(' ')
disp(diag(MU1));


% Saddle Path Condition 

% Forward looking variables are consumption and mu, ordered in first and
% third position respectively

% Extract saddle path condtion parameters which are on the two last rows
P1IF=[P1I(4,1),P1I(4,3);P1I(5,1),P1I(5,3)];
P1IB=[P1I(4,2),P1I(4,4:5);P1I(5,2),P1I(5,4:5)];

SP = -inv(P1IF)*P1IB;
SPC=SP(1,:);
SPMU=SP(2,:);


% Decision rules 


% Dynamics of the state variable S(t+1) / S(t)
SPK = [W(2,2), W(2,4), W(2,5)] + [W(2,1), W(2,3)]*SP;
SPN = [W(4,2), W(4,4), W(4,5)] + [W(4,1), W(4,3)]*SP;
SPZ = [W(5,2), W(5,4), W(5,5)] + [W(5,1), W(5,3)]*SP;
PIB = [SPK;SPN;SPZ];

% Decision rules for control variables d(t) / S(t)


M2B = [M2(1:6,2),M2(1:6,4:5)];
PIC = M1\(M2B + M2(1:6,1)*SPC + M2(1:6,3)*SPMU);
disp('State variables Dynamics') 
disp(' ')
disp(PIB)


disp('(c,v,l,p,m,w) control variables decision rules ')
disp(' ')
disp(PIC)

%% IRF

nrep=60;

CHOC=[0; 0 ; 1];

for j=1:nrep
RC=(PIB^(j-1))*CHOC;
RCK(j)=RC(1);
RCN(j)=RC(2);
RCZ(j) = RC(3);
RCC(j) = SPC*[RC(1);RC(2);RC(3)];
end

for j=1:nrep;
RC=(PIC*PIB^(j-1))*CHOC;
 RCY(j)=RC(1);
 RCV(j)=RC(2);
 RCL(j)=RC(3);
 RCP(j)=RC(4);
 RCM(j)=RC(5);
 RCW(j)=RC(6);
end;

% Investment, total hours and capital rental rate
for j=2:nrep
    RCI(j) = RCK(j)/delta - (1-delta)*RCK(j-1)/delta;
end

RCR=RCY-RCK;
RTH= RCL+RCN;

figure
subplot(221),plot(RCY(1:nrep))
title('Output')
xlabel('quarters')
ylabel('% Dev.   ')

subplot(222),plot(RCZ(1:nrep))
title('TPF')
xlabel('quarters')
ylabel('% Dev.   ')

subplot(223),plot(RCK(1:nrep))
title('Capital')
xlabel('quarters')
ylabel('% Dev.   ')

subplot(224),plot(RCI(2:nrep))
title('Investment')
xlabel('quarters')
ylabel('% Dev.   ')


figure
subplot(221),plot(RCC(1:nrep))
title('Consumption')
xlabel('quarters')
ylabel('% Dev.   ')

subplot(222),plot(RTH(1:nrep))
title('Total hours')
xlabel('quarters')
ylabel('% Dev.   ')


subplot(223),plot(RCW(1:nrep))
title('Wages')
xlabel('quarters')
ylabel('% Dev.   ')

subplot(224),plot(RCR(1:nrep))
title('Rental rate of capital')
xlabel('quarters')
ylabel('% Dev.   ')


figure
subplot(221),plot(RCN(1:nrep))
title('Extensive margin')
xlabel('quarters')
ylabel('% Dev.   ')

subplot(222),plot(RCL(1:nrep))
title('Intensive margin')
xlabel('quarters')
ylabel('% Dev.   ')


subplot(223),plot(RCV(1:nrep))
title('Vacancies')
xlabel('quarters')
ylabel('% Dev.   ')

subplot(224),plot(RCM(1:nrep))
title('Matches')
xlabel('quarters')
ylabel('% Dev.   ')



%% Stochastic simulation

nsimul=1000;
nlong=500;


for j=1:nsimul;

disp('simulation')
disp(j)


aleaa(1:nlong,j)=randn(nlong,1);


for i=1:nlong;
epsa(i)= aleaa(i,j) * (sdepsilona);
end;



CHT(1)=epsa(1);

for i=2:nlong;
CHT(i)=rhoa*CHT(i-1)+epsa(i);
end;


KC(1)=0;
NC(1)=0;
% parties cycliques
for i=2:nlong;
    KC(i) = PIB(1,:)*[KC(i-1); NC(i-1); CHT(i-1)];
    NC(i) = PIB(2,:)*[KC(i-1); NC(i-1); CHT(i-1)];
end;

for i=1:nlong;
    
CC(i)=SPC*[KC(i); NC(i); CHT(i)]; 

YC(i)=PIC(1,:)*[KC(i); NC(i); CHT(i)]; 

VC(i)=PIC(2,:)*[KC(i); NC(i); CHT(i)]; 

LC(i)=PIC(3,:)*[KC(i); NC(i); CHT(i)]; 

PC(i)=PIC(4,:)*[KC(i); NC(i); CHT(i)]; 

MC(i)=PIC(5,:)*[KC(i); NC(i); CHT(i)]; 

WC(i)=PIC(6,:)*[KC(i); NC(i); CHT(i)]; 
end;

% Investment 
for i=1:nlong-1;
    IC(i) = KC(i+1)/delta - (1-delta)*KC(i)/delta;
end
Kterminal = PIB(1,:)*[KC(nlong); NC(nlong); CHT(nlong)];
IC(nlong) = Kterminal/delta - (1-delta)*KC(nlong)/delta;

% Total hours
THC = LC+NC;
% Wage bill
WBC = WC+THC;
% Labour share
LSC = WBC-YC;
% Productivity
PRC = YC-THC;
% Unemployment
UC=-(n/(1-n)).*NC;
% Labour market tightness
LMTC=VC-UC;


% Table 1: Business Cycles Properties
[YHPT,YHP]=hpfilter(YC,1600);
[CHPT,CHP]=hpfilter(CC,1600);
[IHPT,IHP]=hpfilter(IC,1600);
[THHPT,THHP]=hpfilter(THC,1600);
[NHPT,NHP]=hpfilter(NC,1600);
[LHPT,LHP]=hpfilter(LC,1600);
[WBCHPT,WBHP]=hpfilter(WBC,1600);
[LSCHPT,LSHP]=hpfilter(LSC,1600);
[PRCHPT,PRHP]=hpfilter(PRC,1600);
[WHPT,WHP]=hpfilter(WC,1600);

ETYHP(j)=std(YHP(1:nlong));
ETCHP(j)=std(CHP(1:nlong)); 
ETIHP(j)=std(IHP(1:nlong));
ETTHHP(j)=std(THHP(1:nlong));
ETNHP(j)=std(NHP(1:nlong));
ETLHP(j)=std(LHP(1:nlong));
ETWBHP(j)=std(WBHP(1:nlong));
ETLSHP(j)=std(LSHP(1:nlong));
ETPRHP(j)=std(PRHP(1:nlong));
ETWHP(j)=std(WHP(1:nlong));


RHO=corrcoef(YHP,CHP);
RHOCHP(j)=RHO(1,2);
RHO=corrcoef(YHP,IHP);
RHOIHP=RHO(1,2);
RHO=corrcoef(YHP,THHP);
RHOTHHP=RHO(1,2);
RHO=corrcoef(YHP,NHP);
RHONHP(j)=RHO(1,2);
RHO=corrcoef(YHP,LHP);
RHOLHP(j)=RHO(1,2);
RHO=corrcoef(YHP,WBHP);
RHOWBHP(j)=RHO(1,2);
RHO=corrcoef(YHP,LSHP);
RHOLSHP(j)=RHO(1,2);
RHO=corrcoef(YHP,PRHP);
RHOPRHP(j)=RHO(1,2);
RHO=corrcoef(YHP,WHP);
RHOWHP=RHO(1,2);



xxY=YHP(2:nlong);
yyY=YHP(1:nlong-1);
RHO=corrcoef(xxY,yyY);
RHOYYHP(j)=RHO(1,2);

xxC=CHP(2:nlong);
yyC=CHP(1:nlong-1);
RHO=corrcoef(xxC,yyC);
RHOCCHP(j)=RHO(1,2);

xxI=IHP(2:nlong);
yyI=IHP(1:nlong-1);
RHO=corrcoef(xxI,yyI);
RHOIIHP(j)=RHO(1,2);

xxTH=THHP(2:nlong);
yyTH=THHP(1:nlong-1);
RHO=corrcoef(xxTH,yyTH);
RHOTHTHHP(j)=RHO(1,2);

xxN=NHP(2:nlong);
yyN=NHP(1:nlong-1);
RHO=corrcoef(xxN,yyN);
RHONNHP(j)=RHO(1,2);

xxL=LHP(2:nlong);
yyL=LHP(1:nlong-1);
RHO=corrcoef(xxL,yyL);
RHOLLHP(j)=RHO(1,2);

xxWB=WBHP(2:nlong);
yyWB=WBHP(1:nlong-1);
RHO=corrcoef(xxWB,yyWB);
RHOWBWBHP(j)=RHO(1,2);

xxLS=LSHP(2:nlong);
yyLS=LSHP(1:nlong-1);
RHO=corrcoef(xxLS,yyLS);
RHOLSLSHP(j)=RHO(1,2);

xxPR=PRHP(2:nlong);
yyPR=PRHP(1:nlong-1);
RHO=corrcoef(xxPR,yyPR);
RHOPRPRHP(j)=RHO(1,2);

xxW=WHP(2:nlong);
yyW=WHP(1:nlong-1);
RHO=corrcoef(xxW,yyW);
RHOWWHP(j)=RHO(1,2);



% table 2: total hours, productivity and real wages
CORPR(j,:)=xcorr(PRHP, THHP, 4, 'normalized')';

CORW(j,:)=xcorr(WHP, THHP, 4, 'normalized')';

% table 3: unemployment and vacancies
[VHPT,VHP]=hpfilter(VC,1600);
[UHPT,UHP]=hpfilter(UC,1600);

CORUU(j,:)=xcorr(UHP, 4, 'normalized');
CORUV(j,:)=xcorr(VHP, UHP, 4, 'normalized');


% Shimer puzzle
[LMTHPT,LMTHP]=hpfilter(LMTC,1600);
ETVHP(j)=std(VHP(1:nlong));
ETUHP(j)=std(UHP(1:nlong));
ETLMTHP(j)=std(LMTHP(1:nlong));
end;

%table 1
mETYHP=mean(ETYHP);
mETCHP=mean(ETCHP);
mETIHP=mean(ETIHP);
mETTHHP=mean(ETTHHP);
mETNHP=mean(ETNHP);
mETLHP=mean(ETLHP);
mETWBHP=mean(ETWBHP);
mETLSHP=mean(ETLSHP);
mETPRHP=mean(ETPRHP);
mETWHP=mean(ETWHP);

mRHOCHP=mean(RHOCHP);
mRHOIHP=mean(RHOIHP);
mRHOTHHP=mean(RHOTHHP);
mRHONHP=mean(RHONHP);
mRHOLHP=mean(RHOLHP);
mRHOWBHP=mean(RHOWBHP);
mRHOLSHP=mean(RHOLSHP);
mRHOPRHP=mean(RHOPRHP);
mRHOWHP=mean(RHOWHP);

mRHOYYHP=mean(RHOYYHP);
mRHOCCHP=mean(RHOCCHP);
mRHOIIHP=mean(RHOIIHP);
mRHOTHTHHP=mean(RHOTHTHHP);
mRHONNHP=mean(RHONNHP);
mRHOLLHP=mean(RHOLLHP);
mRHOWBWBHP=mean(RHOWBWBHP);
mRHOLSLSHP=mean(RHOLSLSHP);
mRHOPRPRHP=mean(RHOPRPRHP);
mRHOWWHP=mean(RHOWWHP);


mET1=[mETYHP mETCHP mETIHP mETTHHP mETNHP mETLHP mETWBHP mETLSHP mETPRHP mETWHP]./mETYHP;
mRHO1=[1 mRHOCHP mRHOIHP mRHOTHHP mRHONHP mRHOLHP mRHOWBHP mRHOLSHP mRHOPRHP mRHOWHP];
mARHO1=[mRHOYYHP mRHOCCHP mRHOIIHP mRHOTHTHHP mRHONNHP mRHOLLHP mRHOWBWBHP mRHOLSLSHP mRHOPRPRHP mRHOWWHP];

disp('variable order: Y - C - I - TH - N - L - WB - LS - PR - W ')
disp(' ')

disp('relative standard deviation')
disp(' ')
disp(mET1)


disp('correlation')
disp(' ')
disp(mRHO1)


disp('serial correlation')
disp(' ')
disp(mARHO1)


%table 2
mCORPR = mean(CORPR);
mCORW = mean(CORW);

disp('Cross correlation hours productivity')
disp('')
disp(mCORPR)

disp('Cross correlation hours wages')
disp('')
disp(mCORW)

%table 3
mCORUU=mean(CORUU);
mCORUV=mean(CORUV);

disp('Cross correlation unemployment')
disp('')
disp(mCORUU)

disp('Cross correlation unemployment vacancies')
disp('')
disp(mCORUV)


% ACF output TFP
Ygrowth=(YC(2:end)-YC(1:end-1))';
TPFgrowth=(CHT(2:end)-CHT(1:end-1))';

ACF=autocorr(Ygrowth,24);
ACFY=[ACF(2:end)];

ACF=autocorr(TPFgrowth,24);
ACFTPF=[ACF(2:end)];

figure
plot([ACFY, ACFTPF])
title('ACF')
xlabel('Lags')
ylabel('Autocorrelations ')
legend = legend('Output growth', 'TPF growth');
set(legend,'FontSize',12,'Interpreter','Latex')


% Table for Shimer puzzle
mETVHP=mean(ETVHP);  
mETUHP=mean(ETUHP);
mETLMTHP=mean(ETLMTHP); 


% latex table 
a = [mET1; mRHO1; mARHO1];
text = ["(1)","(2)","(3)"];
table11 = fopen('file.tex', 'w');
fprintf(table11, '\\begin{tabular}{|cccccc|}\\hline \n');
fprintf(table11, ' & Output & Consumption & Investment & Total hours & Employment \\\\ \\hline \n');
for i=1:3
    fprintf(table11, '%s & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ ', text(i), a(i,1), a(i,2), a(i,3), a(i,4), a(i,5));
    if i==3
        fprintf(table11, '\\hline ');
    end
    fprintf(table11, '\n');
end
fprintf(table11, '  & Hours per worker & Wage Bill & Labour share & Productivity & Real wages \\\\ \\hline \n');
for i=1:3
    fprintf(table11, '%s & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f   \\\\ ', text(i), a(i,6), a(i,7), a(i,8), a(i,9), a(i,10));
    if i==3
        fprintf(table11, '\\hline ');
    end
    fprintf(table11, '\n');
end
fprintf(table11, '\\end{tabular}\n');
fclose(table11);


table2 = fopen('file.tex', 'w');
fprintf(table2, '\\begin{tabular}{|cccccccccc|}\\hline \n');
fprintf(table2, 'Lags & -4 & -3 & -2 & -1 & 0 & 1 & 2 & 3 & 4 \\\\ \\hline \n');
fprintf(table2, '%s & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ ', "Productivity", mCORPR(1), mCORPR(2), mCORPR(3), mCORPR(4), mCORPR(5), mCORPR(6), mCORPR(7), mCORPR(8), mCORPR(9));
fprintf(table2, '\n');
fprintf(table2, '%s & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ ', "Real wages", mCORW(1), mCORW(2), mCORW(3), mCORW(4), mCORW(5), mCORW(6), mCORW(7), mCORW(8), mCORW(9));
fprintf(table2, '\\hline ');
fprintf(table2, '\n');
fprintf(table2, '\\end{tabular}\n');
fclose(table2);


table3 = fopen('file.tex', 'w');
fprintf(table3, '\\begin{tabular}{|cccccccccc|}\\hline \n');
fprintf(table3, 'Lags & -4 & -3 & -2 & -1 & 0 & 1 & 2 & 3 & 4 \\\\ \\hline \n');
fprintf(table3, '%s & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ ', "Unemployment", mCORUU(1), mCORUU(2), mCORUU(3), mCORUU(4), mCORUU(5), mCORUU(6), mCORUU(7), mCORUU(8), mCORUU(9));
fprintf(table3, '\n');
fprintf(table3, '%s & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f  \\\\ ', "Vacancies", mCORUV(1), mCORUV(2), mCORUV(3), mCORUV(4), mCORUV(5), mCORUV(6), mCORUV(7), mCORUV(8), mCORUV(9));
fprintf(table3, '\\hline ');
fprintf(table3, '\n');
fprintf(table3, '\\end{tabular}\n');
fclose(table3);


table4 = fopen('file.tex', 'w');
fprintf(table4, '\\begin{tabular}{|ccc|}\\hline \n');
fprintf(table4, '  Vacancies & Unemployment & Tightness \\\\ \\hline \n');
fprintf(table4, ' %8.2f & %8.2f & %8.2f   \\\\ ', mETVHP, mETUHP, mETLMTHP);
fprintf(table4, '\\hline ');
fprintf(table4, '\n');
fprintf(table4, '\\end{tabular}\n');
fclose(table4);