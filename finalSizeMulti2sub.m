function [f,g]=finalSizeMulti2sub(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,isdual)
%This is an upgrade on "MA2subSpace(Debug)"
%Parameters:
%R0=1.8;
eps=.2;
cross=0;
if isdual==0
    beta=betaS;
elseif isdual==1
    beta=betaD;
else
    beta=betaI;
end
%
%solvetype: ODE only
demog=1;%Ageing: 1=on
amax=80;
tauend=50;
plotTau=0;%Plot curve for this year: plotTau<1 for no plot
time=(1:tauend);
lt=length(time);
t0=0; tend=360;
mu=0;
S0=[0;1;0];
%
phi1=1; phi2=0;
numseed=10^(-8); seed=numseed;
%
%Theta:
Prow1=[0,0,0,0,.5,.5];
lp1=length(Prow1);
Prow1=Prow1/sum(Prow1);
Prow2=[0,0,0,1-eps,eps];
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
end

%%
A4=[1,1,1,1;1,0,1,0;0,1,0,1]; B4=[1,0,1,0;0,1,0,1;1,0,1,0;0,1,0,1];
A=kron(A4,D);%.*repD; is repD necessary here? Also, got dimension mismatch error - odd
B=kron(B4,D);%.*repD;
%
Rzero=[1-sum(S0);flipud(S0)];
Rzero=kron(Rzero,NNbar);
%
S0=kron(S0,NNbar)/sum(S0);%Number sus to each distinct set of subtypes
Shat=[S0(1:nbar);S0];%This is Shat(0)
zn=zeros(length(Shat),1);
znbar=zeros(nbar,1);
%
A1=zeros(3*nbar,lt);
A2=zeros(2*nbar,lt); A2(:,1)=zeros(2*nbar,1);%?? Why re-define 1st column?
%b1=zeros(1,lp1); b2=zeros(1,lp2);
b1=repmat(S0(nbar+1:2*nbar),1,lp1).*[repmat(Prow1(2:end),nbar,1),znbar];
b2=repmat(S0(2*nbar+1:3*nbar),1,lp2).*[repmat(Prow2(2:end),nbar,1),znbar];

%Check method/S0 definition above - *Prow if IMMUNE
%b1=zeros(nbar,lp1);
%b2=zeros(nbar,lp2); z1nbar=ones(nbar,1);
%b2(:,1)=[z1nbar];%,.3*z1nbar,.1*z1nbar];

amod=1/amax;%(amax-1)/amax;
options=odeset('refine',1);
%xoverFcn=@(t,y)evZero(t,y,nbar); 
%options=odeset(options,'Events',xoverFcn);
%thresh=1-1/R0;
NNrep(NNrep==0)=1;%********26/9/18
NNrep3=repmat(NNrep,3,1); NNrep4=[NNrep3;NNrep];
%Rzero=zn;
seeddot=1;
for tau=1:lt
%%
y0=[S0;zn;zn];
[tout,yout]=ode45(@(t,y)integr8all(t,y,beta,gamma,nbar,A,B,NN,NNrep3,NNrep4,seed,phi1,phi2,tau,seeddot,cross),[t0,tend],y0,options);
if tau==plotTau
    fs=15; lw=2;
    figure
    Y=yout(:,3*nbar+1:7*nbar);
    Yages=Y(:,1:nbar)+Y(:,nbar+1:2*nbar)+Y(:,2*nbar+1:3*nbar)+Y(:,3*nbar+1:end);
    Ysum=sum(Y,2)/sum(NN);
    hold on
    plot(tout,Yages,'-','linewidth',lw);
    xlabel('Time (days)','FontSize',fs);
    ylabel('Infectives','FontSize',fs);
    set(gca,'FontSize',fs);
    maxY=max(max(max(Y)))+.01;
    axis tight%([0,tend,0,maxY])
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
    Rzero=zn;%Vector of Rs
    %Rzero=[S0(3*nbar+1:end);S0(2*nbar+1:3*n);S0(nbar+1:2*n);S0(1:nbar)];
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
%S=[y(1:nbar);cross*y(nbar+1:2*nbar);cross*y(2*nbar+1:3*nbar)];%XXXX possibility of cross12 and cross21
S1=y(1:nbar); S2=y(nbar+1:2*nbar); S3=y(2*nbar+1:3*nbar);
S=[S1;(1-cross)*S2;(1-cross)*S3];%XXXX
Shat=[S(1:nbar);S];
I=y(3*nbar+1:7*nbar);
seedf=seed*S*exp(-t).*seeddot;%.* if seeddot is a vector (non-trivial case)
seedfhat=[seedf(1:nbar);seedf];
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