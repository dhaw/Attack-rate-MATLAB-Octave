function [f,g,other]=MA1D(solvetype,start,eps)
%isdual: 0=SM, 1=DM, 2=IM
%solvetype: 1=FSC, 2=ODE, 3=SCM
%Parameters:
R0=1.8; gamma=1/2.6;
NN=1000; n=length(NN);
demog=1;%Ageing: 1=on
agemix=1;%Total homog: =0; else specify v - twice
tauend=20100;
time=(1:tauend);
lt=length(time);
t0=0; tend=3600;
mu=0;
phi1=1; phi2=0;
%NNprob=NNbar/sum(NN); %NNprob=ones(nbar,1)/sum(NNbar);
numseed=10^(-3); seed=numseed;%*NNprob;
addbit=0;%seed;%********SEED IN FSC********
randic=1;
%
%Theta:
%eps=.15;
%Prow=exppdf((1:10),5);
%Prow=poisspdf((0:10),5);
%Prow=binopdf((0:10),4500,1/1000);
%Prow=[0,0,0,.1,.5,.3,.1];%Chaotic?
%Prow=ones(1,100);
%Prow=[0,0,.5*eps,eps,1-3*eps,eps,.5*eps];
Prow=[0,0,0,eps,1-eps];
%Prow=[eps,eps,eps,eps,eps,1-8*eps,eps,eps,eps,eps,eps];
lp=length(Prow);
Prow=Prow/sum(Prow);
%%
if agemix==0
    C=1; NNbar=NN;
    amax=80;
else
    %C=[6,.5;1,.5];%Steven's paper
    Cnum=[6.92,.25,.77,.45;.19,3.51,.57,.2;.42,.38,1.4,.17;.36,.44,1.03,1.83];
    Cdur=[3.88,.28,1.04,.49;.53,2.51,.75,.5;1.31,.8,1.14,.47;1,.85,.88,1.73];
    C=Cnum.*Cdur;
    %Age distribution (nbar=4):
    v=[5,14,45,16]/80;
    NNbar=[NN*v(1);NN*v(2);NN*v(3);NN*v(4)]; %NNbar=round(NNbar);%For discrete populations only (SCM/ABM)
    if solvetype==3
        NNbar=round(NNbar);%For discrete populations only (SCM/ABM)
        NN=sum(NNbar);
    end
end
nbar=length(C);%Number of age groups
NNrep=repmat(NN,nbar,1);
NN0=NNrep; NN0(NNrep==0)=1;
Nages=NNbar./NN0;
%
Sstart=repmat(NNbar,1,nbar);
%ABM needs matrix "ages" - omitted here
%R0 calculation:
D1=Sstart.*C/NN; G=1/gamma*D1;
d=eigs(G,1); R0a=max(d); beta=R0/R0a;
%%
D=C; Pmat=repmat(Prow,nbar,1);
%if ageoff==1
    %V=repmat(1/NN,1,lp); W=V;
%else
if agemix~=0
    v=v*80;%v=[5,14,45,16];
    vover=1./v; V=vover'; V=repmat(V,1,lp);
    w=[0,v(1:end-1)];%Only kill off from oldest
    woverv=w.*vover; W=woverv'; W=repmat(W,1,lp);
end
A1=zeros(1,lt);
A2=A1; other=A1;
%%
%Brand=rand(nbar,1); %Brand=Brand./repmat(sum(Brand,2),1,5);
B=zeros(nbar,lp);
%B(:,2)=.1*rand(nbar,1);%Non-trivial IC %B(:,2)=.2*Brand;
%
%Z0=sum(B(:,2:end),2)*(1-mu);
Z0=start*ones(nbar,1);%Initial immunity
%Z0=zeros(nbar,1)+start*rand(nbar,1);
if randic==1
    maxYears=lp-1;
    numNonZero=ceil(maxYears*rand);%1-6 years - pick whow many are non zero
    here=randsample(maxYears,numNonZero); pend=max(here);%Pick which ones are non-zero
    vals=rand(1,numNonZero);%Random number for each non-zero
    vals=rand*vals/sum(vals);%Total proportion immune=rand
    ProwRand=zeros(1,pend);
    ProwRand(here)=vals;
    Brand=rand(nbar,1).*NNbar./NN0;
    pdiff=maxYears-pend;
    if pdiff>0
        ProwRand(pend+1:maxYears)=0;
    end
    B=[zeros(nbar,1),repmat(ProwRand,nbar,1).*repmat(Brand,1,maxYears)];
    Z0=sum(B(:,2:end),2)*(1-mu);
end
%%
zn=zeros(nbar,1);
for tau=1:lt
    if solvetype==1%FSC
    %IC=.5+.1*rand(nbar,1);
    IC=XODEsolveAllMulti(gamma,NN,n,nbar,NNbar,NNrep,NN0,D,Z0,beta,t0,tend,zn,phi1,phi2,seed,tau,2); IC=IC./NN0;
    if agemix~=1
        IC=IC+.2*(1+rand(nbar,1)).*(1-IC);
    end
    other(tau)=IC(1);
    options=optimset('Display','off');
    funt=@(Zi)solveZi(Zi,Z0,beta,gamma,D,Nages,addbit);
    Zsol=fsolve(funt,IC,options);
    %
    else%solvetype=2/3
    Zsol=XODEsolveAllMulti(gamma,NN,n,nbar,NNbar,NNrep,NN0,D,Z0,beta,t0,tend,zn,phi1,phi2,seed,tau,solvetype);
    Zsol=Zsol./NN0;
    end
    nu=Zsol-Z0;
    A1(:,tau)=sum(reshape(Zsol,1,nbar),2);%Prop immune for spatial cell (before antigenic drift)
    A2(:,tau)=sum(reshape(nu,1,nbar),2);%AR for spatial cell
    %Assign immunity:
    Bhat=[B(:,2:end),zeros(nbar,1)];
    B=Bhat+repmat(nu,1,lp).*Pmat;
    %Age the populations:
    if demog==1
        if agemix==0
            B=B*(amax-1)/amax;
        else
            BVtake=B.*V; Bsus=sum(BVtake(1,2:end),2);
            BVadd=[zeros(1,lp);BVtake(1:nbar-1,:)].*W;
            B=B-BVtake+BVadd; B(1,1)=B(1,1)+Bsus;
        end
    end
    Z0=sum(B(:,2:end),2)*(1-mu);%SCM: need to round NNbar?
end
%%
f=A1;
g=A2;
end

function f=solveZi(Zi,Z0,beta,gamma,D,Nages,addbit)
mu2=0;
lZi=log(Nages-Zi); lZi(Zi>1)=NaN;
f=lZi-log(Nages-Z0*(1-mu2))+beta/gamma*D*(Zi-Z0*(1-mu2))+addbit;
end
%%
function f=XODEsolveAllMulti(gamma,NN,n,nbar,NNbar,NNrep,NN0,D,ic,beta,t0,tend,zn,phi1,phi2,seed,tau,solvetype)
icR=ic.*NNrep;
y0=[NNbar-icR;zn;icR];
%
%Deterministic:
if solvetype==2
[tout,yout]=ode45(@(t,y)integr8all(t,y,beta,gamma,nbar,NN,NN0,D,seed,phi1,phi2),[t0,tend],y0);
%
%Incidence curve in here (code at bottom)
%{
if tau==71
figure
n=length(NN);
fs=15;
hold on
plot(tout,yout(:,1),'-','linewidth',2,'color',[0,.7,0])
plot(tout,yout(:,2),'-','linewidth',2,'color',[.7,0,0;])
plot(tout,yout(:,3),'k-','linewidth',2)
xlabel('Time (days)','FontSize',fs);
ylabel('Population','FontSize',fs);
set(gca,'FontSize',fs);
axis([0,tend,0,1])
legend('S','I','R')
end
%}
%
rec0=yout(end,2*nbar+1:end);
f=rec0';
%}
else
%Stochastic (don't use with FSC!):
%ic1=round(NNbar/1000); neg=y0(1:nbar)-ic1; neg(neg>0)=0; ic1=ic1+neg;
%ic1=round(y0(1:nbar)/500);
y1=round(y0);%+[-ic1;ic1;zn];
f=stochSim(y1,beta,gamma,n,nbar,NN,NN0,D,seed,phi1,phi2,tau);
%}
end
end
%%
function f=integr8all(t,y,beta,gamma,nbar,NN,NN0,D,seed,phi1,phi2)
mu=1/1800;%5*360=1800
phi=1;%phi1-phi2*cos(pi*t/180);
S=y(1:nbar);
I=y(nbar+1:2*nbar);
R=y(2*nbar+1:end);
Sfoi=beta*S.*(D*(I./NN0))*phi;
seed=seed*S./NN0*exp(-t);
Sdot=-Sfoi-seed;%+mu*R;
Idot=Sfoi+seed-gamma*I;
Rdot=gamma*I;%-mu*R;
f=[Sdot;Idot;Rdot];
end
%%
%%
function f=stochSim(y,beta,gamma,n,nbar,NN,N0,D,seed,phi1,phi2,tau)
mu=1/1800;%5*360=1800
factor=1;
tend=360*factor; beta=beta/factor; gamma=gamma/factor;
%
S=y(1:nbar);
I=y(nbar+1:2*nbar);
R=y(2*nbar+1:end);
%Different from here:
i=1;
while i<tend && (i<30 || sum(I)>0)%At least 30 time-steps
phi=1;%phi1-phi2*cos(pi*i*f1/180);%Just seed here
Sout=1-exp(-phi*(beta*(D*(I./N0))+seed));%+mu*R;
Sout=binornd(S,Sout); Sout(S==0)=0;
S=S-Sout; I=I+Sout;
%
Iout=1-exp(-gamma);
Iout=binornd(I,Iout);
I=I-Iout; R=R+Iout;
i=i+1;
end
f=R;
if sum(isnan(R))>0
    fuck=1;
end
end
%%
%Incidence curve:
%{
if tau==10
figure
fs=15;
Y=yout(:,nbar+1:2*nbar);
Y=Y(:,1:n)+Y(:,n+1:2*n)+Y(:,2*n+1:3*n)+Y(:,3*n+1:end);
Nover=1./NN; Nover(NN==0)=0; Y=Y.*repmat(Nover',length(tout),1);
X=[NN,Y']; X=sortrows(X,1); Nsort=X(:,1); X=X(:,2:end); X=X';
hold on
k=100; nk=floor(n/k);
cmap=parula(nk);
colormap(cmap);
for i=1:nk
    plot(tout,X(:,i*100),'-','linewidth',2,'color',cmap(i,:));
end
xlabel('Time (days)','FontSize',fs);
ylabel('Infectives','FontSize',fs);
set(gca,'FontSize',fs);
axis([0,tend,0,max(max(X))])
end
%}