function [f,g]=finalSizeMulti2sub(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,isdual,numseed,eps,del,cross,randic,tauend)
%This is an upgrade on "MA2subSpace(Debug)"
%Parameters:
%solvetype: ODE only
%eps=.3;
%cross=.6;
%randic=1;
ICdepAge=1;
if isdual==0
    beta=betaS;
elseif isdual==1
    beta=betaD;
else
    beta=betaI;
end
demog=1;%Ageing: 1=on
%tauend=200;
plotTau=0;%Plot curve for this year: plotTau<1 for no plot
time=(1:tauend);
lt=length(time);
t0=0; tend=1800;
mu=0;
%%
%Initial condition (independent of space):
%Default - total immunity:
S0=[1;0;0];
%IC with immunity propto population density (1 year IC/1YIC):
%{
Srand=rand(3,1);
Srand=Srand/sum(Srand);
S0=.2*Srand;%(1:3);
%}
%%
phi1=1; phi2=0;
%numseed=10^(-8); 
seed=numseed;
%
%Theta:
Prow1=[0,0,0,0,.5-del,.5+del];%zeros(1,21); Prow1(end)=1;
lp1=length(Prow1);
Prow1=Prow1/sum(Prow1);
Prow2=[0,0,0,1-eps,eps];%Prow1=Prow2; lp1=length(Prow1);
lp2=length(Prow2);
Prow2=Prow2/sum(Prow2);
%%
vover=[1/5,1/14,1/45,0];
V=kron(vover',ones(n,1)); V1=repmat(V,1,lp1); V2=repmat(V,1,lp2);
VV=kron(vover',ones(n,1)); Stake123=repmat(VV,3,1);
wover=[0,0,0,1/16];
WW=kron(wover',ones(n,1)); Stake4=repmat(WW,3,1);
%%
Ni=repmat(NNrep,1,nbar); Nj=Ni';
Niover=1./Ni; Niover(Ni==0)=1; Njover=Niover';
Mj=(Kbar')*NNbar;
Mj(Mj==0)=1;
Mjover=1./Mj;
Mjover=repmat(Mjover',nbar,1);
if isdual==0
    D=Kbar.*Mjover.*Cbar.*Nj;
elseif isdual==1
    matrix1=K1.*Mjover;
    matrix2=Kbar'.*Cbar;
    D=matrix1*matrix2;
    D=D.*Nj;
elseif isdual==2
    D=Kbar'.*Niover.*Cbar.*Nj;
elseif isdual==3
    D=Kbar.*Cbar;
else
    error('Invalid isdual identifier')
end
%%
A4=[1,1,1,1;1,0,1,0;0,1,0,1]; B4=[1,0,1,0;0,1,0,1;1,0,1,0;0,1,0,1];
A=kron(A4,D);
B=kron(B4,D);
%%
znbar=zeros(nbar,1);
if randic==1
    %From 1D model:
    if ICdepAge==1
        r=1;
        maxYears1=lp1-1;
        Prows=zeros(nbar,maxYears1);
        numNonZero=floor((maxYears1+1)*rand(nbar,1));%0-6 years - pick whow many are non zero
        for i=1:nbar
            rsi=randsample(maxYears1,numNonZero(i));
            Prows(i,rsi)=1;
        end
        findNZ1=find(Prows);
        lNZ=length(findNZ1);
        randNZ=rand(lNZ,1);
        Prows(findNZ1)=randNZ;
        rand2=r*rand(nbar,1); repRand=repmat(rand2,1,maxYears1);%rand2=r*rand()
        sumP=sum(Prows,2); repP=repmat(sumP,1,maxYears1); repP(repP==0)=1;
        X=Prows./repP.*repRand;
        %
        maxYears2=lp2-1;
        Prows=zeros(nbar,maxYears2);
        numNonZero=floor((maxYears2+1)*rand(nbar,1));%0-6 years - pick whow many are non zero
        for i=1:nbar
            rsi=randsample(maxYears2,numNonZero(i));
            Prows(i,rsi)=1;
        end
        findNZ2=find(Prows);
        lNZ=length(findNZ2);
        randNZ=rand(lNZ,1);
        Prows(findNZ2)=randNZ;
        rand2=r*rand(nbar,1); repRand=repmat(rand2,1,maxYears2);%rand2=r*rand()
        sumP=sum(Prows,2); repP=repmat(sumP,1,maxYears2); repP(repP==0)=1;
        Y=Prows./repP.*repRand;
        %b1=[znbar,X.*repmat(NNbar,1,maxYears1)]; b2=[znbar,Y.*repmat(NNbar,1,maxYears2)];
        b1=[X.*repmat(NNbar,1,maxYears1),znbar]; b2=[Y.*repmat(NNbar,1,maxYears2),znbar];
        %
        %Remove overlap:
        s01=sum(X,2); s10=sum(Y,2);%s11/s00 surrently not specified
        Z=s01+s10;
        s00=Z-1; s00(s00<0)=0;
        findDef=find(s00);
        Xfrac=s01(findDef)./Z(findDef); Yfrac=s10(findDef)./Z(findDef);
        %b1frac=ones(nbar,1); b1frac(findDef)=Xfrac;
        %b2frac=ones(nbar,1); b2frac(findDef)=Yfrac;
        %b1(:,1:end-1)=b1(:,1:end-1).*repmat(b1frac,1,maxYears1); b2(:,1:end-1)=b2(:,1:end-1).*repmat(b2frac,1,maxYears2); 
        %b1(:,2:end)=b1(:,2:end).*repmat(b1frac,1,maxYears1); b2(:,2:end)=b2(:,2:end).*repmat(b2frac,1,maxYears2); 
        Xtake=Xfrac.*s00(findDef); Ytake=Yfrac.*s00(findDef);
        s01(findDef)=s01(findDef)-Xtake; s10(findDef)=s10(findDef)-Ytake;
        %
        %Define s00:
        %s11=1-s10-s01; 
        Zmin=min([s01,s10],[],2);
        s00add=rand(nbar,1).*Zmin;
        s00=s00+s00add;%/3?
        s10=s10-s00add; s01=s01-s00add;
        s11=1-s10-s01-s00;
        %removed2sub
        S0star=[s11,s10,s01,s00];
        S0star(S0star<0)=0;
        S0sum=sum(S0star,2); S0sum(S0sum==0)=1;
        S0star=S0star./repmat(S0sum,1,4);
        %S0star=S0star./sum(S0star);
        S0=reshape(S0star(:,1:3),3*nbar,1).*repmat(NNbar,3,1);
    else
        error('This bit of code isnt written yet')
    end
    clear Ni Nj Mj Niover Njover Mjover K1 Kbar Cbar
    clear X X1 Y Y1 Z Z1 Nratio1 Nratio2 Nratio3 rsi randNZ rand2 randNZ2 Sdiff sub1 sub2 %s00 s01 s10 s11
    clear maxYears2 maxYears2 numNonZero Prows repRand sumP repP lNZ lpmax
else
    %Default - no spatial/age dependence:
    s00=1-sum(S0); s00=s00*NNbar;%Immune to both
    S0=kron(S0,NNbar);%/sum(S0);%Number sus to each distinct set of subtypes
    %
    b1=zeros(nbar,lp1); b2=zeros(nbar,lp2);
    %1YIC:
    %{
    b1=repmat(s00+S0(2*nbar+1:3*nbar),1,lp1).*[repmat(Prow1(2:end),nbar,1),znbar];
    b2=repmat(s00+S0(nbar+1:2*nbar),1,lp2).*[repmat(Prow2(2:end),nbar,1),znbar];
    %}
end
A1=zeros(3*nbar,lt);
A2=zeros(2*nbar,lt); A2(:,1)=zeros(2*nbar,1);%?? Why re-define 1st column?
%%
Shat=[S0(1:nbar);S0];%This is Shat(0)
zn=zeros(length(Shat),1);
options=odeset('refine',1);
NNrep(NNrep==0)=1;%********26/9/18
NNrep3=repmat(NNrep,3,1); NNrep4=[NNrep3;NNrep];
seeddot=1;
for tau=1:lt
%%
y0=[S0;zn;zn];
%{
if tau==1
    cross1=0;
else
    cross1=cross;
end
%}
cross1=cross;
[tout,yout]=ode45(@(t,y)integr8all(t,y,beta,gamma,nbar,A,B,NN,NNrep3,NNrep4,seed,phi1,phi2,tau,seeddot,cross1),[t0,tend],y0,options);
if tau==plotTau
    fs=12; lw=2;
    figure
    Y=yout(:,3*nbar+1:7*nbar);
    Yages=Y(:,1:nbar)+Y(:,nbar+1:2*nbar)+Y(:,2*nbar+1:3*nbar)+Y(:,3*nbar+1:end);
    Ysum=sum(Y,2)/sum(NN);
    Y1=Y(:,1:nbar)+Y(:,2*nbar+1:3*nbar);
    Y1=Y1(:,1:n)+Y1(:,n+1:2*n)+Y1(:,2*n+1:3*n)+Y1(:,3*n+1:end);
    Y2=Y(:,nbar+1:2*nbar)+Y(:,3*nbar+1:end);
    Y2=Y2(:,1:n)+Y2(:,n+1:2*n)+Y2(:,2*n+1:3*n)+Y2(:,3*n+1:end);
    Y1sum=sum(Y1,2); Y2sum=sum(Y2,2);
    hold on
    plot(tout,[Y1sum,Y2sum],'-','linewidth',lw);
    xlabel('Time (days)','FontSize',fs);
    ylabel('Infectives','FontSize',fs);
    set(gca,'FontSize',fs);
    maxY=max(max(max([Y1sum,Y2sum])))+.01;%Y
    axis ([0,360,0,maxY])%tend
    legend('H1','H2','location','NE')
    grid on
    grid minor
    hold off
end
    Rt=yout(end,7*nbar+1:end)+yout(end,3*nbar+1:7*nbar);
    Rt=Rt';
    Rt1=Rt(1:nbar); Rt2=Rt(nbar+1:2*nbar); Rt3=Rt(2*nbar+1:3*nbar); Rt4=Rt(3*nbar+1:end);
    A1tx=S0-[Rt1+Rt2;Rt3;Rt4]+[zeros(nbar,1);Rt2;Rt1]; 
    toSum=reshape(A1tx,nbar,3); sumA1tx=sum(toSum,2);
    %
    A1t=[A1tx;NNbar-sumA1tx];
    A2t=[Rt1+Rt3;Rt2+Rt4];
    A2(:,tau)=A2t;
    %
    minAttack=min(A2t); maxAttack=max(A2t./repmat(NNrep,2,1));
    if minAttack<0 || maxAttack>1
        error('Attack rate not in [0,1]')
    end
    %
    b1=b1+repmat((Rt1+Rt3),1,lp1).*repmat(Prow1,nbar,1);
    b2=b2+repmat((Rt2+Rt4),1,lp2).*repmat(Prow2,nbar,1);
    w1=b1(:,1); b1=[b1(:,2:end),znbar];
    w2=b2(:,1); b2=[b2(:,2:end),znbar];
    %%
    A1t2=A1t(nbar+1:2*nbar); A1t3=A1t(2*nbar+1:3*nbar); A1t4=A1t(3*nbar+1:end);
    A1t34div=A1t3+A1t4; A1t34div(A1t34div==0)=1;
    A1t24div=A1t2+A1t4; A1t24div(A1t24div==0)=1;
    %Waning - new method - get negative S_00:
    from01=A1t3./A1t34div.*w1;
    from00s1=A1t4./A1t34div.*w1;
    from10=A1t2./A1t24div.*w2;
    from00s2=A1t4./A1t24div.*w2;
    A1t4div=A1t4; A1t4div(A1t4==0)=1;
    frac1=from00s1./A1t4div; frac2=from00s2./A1t4div;
    from00both=frac1.*frac2.*A1t4;
    from00both(isnan(from00both)==1)=0;
    from00s1=from00s1-from00both;
    from00s2=from00s2-from00both;
    A1t=A1t+[from10+from01+from00both;
             from00s1-from10;
             from00s2-from01;
             -from00both-from00s1-from00s2];
    %%
    %Update S0:
    S0=A1t(1:3*nbar);%*NN;%Assume duration of immunity unrelated to cross infection
    %
    %Optional different output - susceptibility at start of season:
    %Need to comment out "A1(:,t)=..." above
    %%
    %Age the populations:
    %Only works if age dist assumed uniform!
    if demog==1
        b1take=b1.*V1; b1add=circshift(b1take,n,1);
        b1(end-n+1:end,:)=b1(end-n+1:end,:)*15/16;
        b1=b1-b1take+b1add;
        b2take=b2.*V2; b2add=circshift(b2take,n,1);
        b2(end-n+1:end,:)=b2(end-n+1:end,:)*15/16;
        b2=b2-b2take+b2add;
        %Prep S0 separately - because of "immune to both" class:
        S0rest=NNbar-sum(reshape(S0,nbar,3),2);
        SageTake=S0.*Stake123; SageAdd=circshift(SageTake,n);
        Sdeath=S0.*Stake4; Sbirth=reshape(Sdeath,n,3*na); Sbirth=sum(Sbirth,2);
        S0=S0-SageTake+SageAdd-Sdeath; S0(1:n)=S0(1:n)+Sbirth+S0rest(end-n+1:end)/16;
    end
    A1(:,tau)=S0;
end

f=A1;
g=A2;
end

%%

function f=integr8all(t,y,beta,gamma,nbar,A,B,NN,NNrep3,NNrep4,seed,phi1,phi2,tau,seeddot,cross)
mu=1/1800;%5*360=1800
phi=1;%phi1-phi2*cos(pi*t/180);
S1=y(1:nbar); S2=y(nbar+1:2*nbar); S3=y(2*nbar+1:3*nbar);
S=[S1;(1-cross)*S2;(1-cross)*S3];%XXXX
Shat=[S(1:nbar);S];
I=y(3*nbar+1:7*nbar);
%seedf=seed*S*exp(-t).*seeddot;%.* if seeddot is a vector (non-trivial case)
%seedfhat=[seedf(1:nbar);seedf];
if t<30
    seedf=seed.*S;%*exp(-t);
    seedfhat=repmat(seed(nbar+1:end),2,1).*Shat;%*exp(-t);
else
    seedf=zeros(3*nbar,1);
    seedfhat=zeros(4*nbar,1);
end
Iscaled=I./NNrep4;
Sdot=-beta*S.*(A*Iscaled)*phi-seedf;%+mu*R;
Rdot=gamma*I;%-mu*R;
Idot=beta*Shat.*(B*Iscaled)*phi+seedfhat-Rdot;
f=[Sdot;Idot;Rdot];
end

%%

function [value,isterminal,direction]=evZero(t,y,nbar)
value=y(3*nbar+1:7*nbar);
isterminal=ones(4*nbar,1);
direction=-ones(4*nbar,1);
end