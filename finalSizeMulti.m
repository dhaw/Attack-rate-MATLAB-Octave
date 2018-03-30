<<<<<<< HEAD
function [f,g]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,beta,isdual,solvetype,numseed)%,NNbar)
%ZfinalSizeAllMulti2
%isdual: 0=SM, 1=DM, 2=IM
%solvetype: 1=FSC, 2=ODE, 3=SCM
tauend=1;
time=(1:tauend);
lt=length(time);
mu=0;%In ODE code
NN0=NNrep; NN0(NNrep==0)=1;
Nages=NNbar./NN0;
t0=0; tend=360;

R0=1.8;

%
%Theta:
eps=.13;
%Prow=poisspdf((0:10),5);
%Prow=binopdf((0:10),45,1/10);
Prow=[0,0,0,0,eps,1-eps];%Chaotic?
%Prow=[0,0,0,eps,1-2*eps,eps];
%Prow=[0,0,0,0,0,eps,1-2*eps,eps];
%Prow=1/11*ones(1,11);
%
lp=length(Prow);
Prow=Prow/sum(Prow);
Pmat=repmat(Prow,nbar,1);
v=[5,14,45,16]; vover=1./v; V=kron(vover',ones(n,1)); V=repmat(V,1,lp);
w=[0,5,14,45]; woverv=w.*vover; W=kron(woverv',ones(n,1)); W=repmat(W,1,lp);
%
cen=0; ages0=0; %maxN=0;
A1=zeros(n,lt);
A2=A1;
N0=NN; N0(NN==0)=1;
%
Ni=repmat(NNrep,1,nbar); Nj=Ni';
Niover=1./Ni; Niover(Ni==0)=1; Njover=Niover';
Mj=(Kbar')*NNbar;%C in denom?? .*Cbar %(Kbar')*
Mj(Mj==0)=1;
Mjover=1./Mj;
Mjover=repmat(Mjover',nbar,1);
if isdual==0
    D=Kbar.*Mjover.*Cbar.*Nj;
elseif isdual==1
    D=(K1.*Mjover)*(Kbar'.*Cbar);
    D=D.*Nj;
elseif isdual==2
    D=Kbar'.*Niover.*Cbar.*Nj;
end
%%
%Brand=rand(nbar,1); %Brand=Brand./repmat(sum(Brand,2),1,5);
B=zeros(nbar,lp);
%B(:,2)=.1*rand(nbar,1);%Non-trivial IC
%B(:,2)=.2*Brand;
%%
Z0=sum(B(:,2:end),2)*(1-mu);%=0;
%For final sizes (if IC not in loop):
%IC=XODEsolveAllMulti(gamma,NN,n,nbar,NNbar,NNrep,NN0,maxN,cen,D,Z0,beta,ages0,1); IC=IC./NN0;
%IC=.5*.2*rand(nbar,1);
%%
zn=zeros(nbar,1);
phi1=1; phi2=0;
%NNprob=NNbar/sum(NN);
%NNprob=ones(nbar,1)/sum(NNbar);
seed=numseed;%*NNprob;

thresh=1-1/R0;

for t=1:lt    
    if solvetype==1
    %Final size:
    addbit=0;%seed;%tend*
    %If seed=0, no epidemic happens - because use ODE solver for initial condition
    IC=XODEsolveAllMulti(gamma,NN,n,nbar,NNbar,NNrep,NN0,minNind,maxNind,D,Z0,beta,ages0,t,t0,tend,zn,phi1,phi2,seed,2,thresh); IC=IC./NN0;%New IC here (if ever necessary)
    options = optimset('Display','off');
    funt=@(Zi)solveZi(Zi,Z0,beta,gamma,D,Nages,addbit);
    Zsol=fsolve(funt,IC,options);
    %
    else%solvetype=2/3
    Zsol=XODEsolveAllMulti(gamma,NN,n,nbar,NNbar,NNrep,NN0,minNind,maxNind,D,Z0,beta,ages0,t,t0,tend,zn,phi1,phi2,seed,solvetype,thresh);
    Zsol=Zsol./NN0;
    end
    nu=Zsol-Z0;
    A1(:,t)=sum(reshape(Zsol,n,na),2);%Prop immune for spatial cell (before antigenic drift)
    A2(:,t)=sum(reshape(nu,n,na),2);%AR for spatial cell
    %
    %Assign immunity:
    Bhat=[B(:,2:end),zeros(nbar,1)];
    B=Bhat+repmat(nu,1,lp).*Pmat;
    %Age the populations:
    %
    BVtake=B.*V; Bsus=sum(BVtake(1:n,2:end),2);
    BVadd=[zeros(n,lp);BVtake(1:3*n,:)].*W;
    B=B-BVtake+BVadd; B(1:n,1)=B(1:n,1)+Bsus;
    %}
    Z0=sum(B(:,2:end),2)*(1-mu);
end
%%
%}
f=A1;
g=A2;
end

function f=solveZi(Zi,Z0,beta,gamma,D,Nages,addbit)
lZi=log(Nages-Zi); lZi(Zi>1)=NaN;
f=lZi-log(Nages-Z0)+beta/gamma*D*(Zi-Z0)+addbit;
end
%%
function f=XODEsolveAllMulti(gamma,NN,n,nbar,NNbar,NNrep,NN0,minNind,maxNind,D,ic,beta,ages0,tau,t0,tend,zn,phi1,phi2,seed,solvetype,thresh)
icR=ic.*NNrep;
y0=[NNbar-icR;zn;icR];
cond=sum(icR<thresh);
if cond==0
    seed=0;
end
%
if solvetype==2
[tout,yout]=ode45(@(t,y)integr8all(t,y,beta,gamma,n,nbar,NN,NN0,D,seed,phi1,phi2),[t0,tend],y0);
%Incidence curve in here:
%{
figure
%fs=15; lw=2
fs=12; lw=2;
Y=yout(:,nbar+1:2*nbar); %Z=yout(:,2*nbar+1:3*nbar);
%Y=sum(Y,2);
Y1=Y(:,minNind)+Y(:,minNind+n)+Y(:,minNind+2*n)+Y(:,minNind+3*n); Y1=Y1/NN(minNind);
Y2=Y(:,maxNind)+Y(:,maxNind+n)+Y(:,maxNind+2*n)+Y(:,maxNind+3*n); Y2=Y2/NN(maxNind);
hold on
plot(tout,Y1,'--','linewidth',lw,'color',[.165,.31,.431]);%[.165,.31,.431][.447,.553,.647]
plot(tout,Y2,'-','linewidth',lw,'color',[.165,.31,.431]);
xlabel('Time (days)','FontSize',fs);
ylabel('Infectives','FontSize',fs);
set(gca,'FontSize',fs);
maxY=max(max([Y1;Y2]));
axis([0,tend,0,maxY])
legend('Min','Max','location','NE')
hold off
%}
rec0=yout(end,2*nbar+1:end);%+yout(end,nbar+1:2*nbar);%******** Truncation - add infectious to rercovered?
f=rec0';

elseif solvetype==3
%ic1=round(NNbar/1000); neg=y0(1:nbar)-ic1; neg(neg>0)=0; ic1=ic1+neg;
%ic1=round(y0(1:nbar)/500);
y1=round(y0);%+[-ic1;ic1;zn];
f=stochSim(y1,beta,gamma,n,nbar,NN,NN0,D,seed,phi1,phi2,tau);
%}
end
end
%%
function f=integr8all(t,y,beta,gamma,n,nbar,NN,NN0,D,seed,phi1,phi2)
mu=1/1800;%5*360=1800
alpha=1;
%phi=phi1-phi2*cos(pi*t/180);
phi=1;%+sin(pi*t/180);
S=y(1:nbar);
I=y(nbar+1:2*nbar);
R=y(2*nbar+1:end);
seed1=seed*S./NN0;%*exp(-t);
Sfoi=phi*(beta*S.*(D*(I./NN0)).^alpha+seed1);
Sdot=-Sfoi;%+mu*R;
Idot=Sfoi-gamma*I;
Rdot=gamma*I;%-mu*R;
f=[Sdot;Idot;Rdot];
end
%%
function f=stochSim(y,beta,gamma,n,nbar,NN,N0,D,seed,phi1,phi2,tau)
%Feed in mu if required
factor=24;
tend=360*factor; beta=beta/factor; gamma=gamma/factor;
Vec=zeros(nbar,tend);
%
S=y(1:nbar);
I=y(nbar+1:2*nbar);
R=y(2*nbar+1:end);

%Test seed:
%[Nmax,maxInd]=max(NN);
%ii=4*(maxInd-1)+3;
%seed=0; I(ii)=10; S(ii)=S(ii)-10;
%
%I=floor(N0/10000); S=S-I;

%Different from here:
i=1;
threshold=30;%Number of time-steps for necessary simulation/seed
while i<tend && (i<threshold || sum(I)>0)%At least 30 time-steps
phi=1;%phi1-phi2*cos(pi*i*f1/180);%Just seed here
Sout=1-exp(-phi*(beta*(D*(I./N0))+seed*heaviside(threshold-i)));%+mu*R;
Sout=binornd(S,Sout); Sout(S==0)=0;
S=S-Sout; I=I+Sout;
%
Iout=1-exp(-gamma);
Iout=binornd(I,Iout);
I=I-Iout; R=R+Iout;
Vec(:,i)=I;
i=i+1;
end

Vsum=sum(Vec,1);%/sum(NN);
maxV=max(Vsum);
fs=15; lw=2;
figure
plot(1:tend/factor,Vsum(factor:factor:end),'-','linewidth',2,'color',[.447,.553,.647]);
xlabel('Time (days)','FontSize',fs);
ylabel('Infectives','FontSize',fs);
axis([0,tend/factor,0,maxV]); xlabel('Time (years)'); ylabel('Incidence'); set(gca,'fontsize',15)

f=R;
if sum(isnan(R))>0
    fuck=1;
end
end

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

figure
fs=15;
Y=yout(:,nbar+1:2*nbar);
Y=sum(Y,2);
plot(tout,Y,'-','linewidth',2,'color',[.447,.553,.647]);
xlabel('Time (days)','FontSize',fs);
ylabel('Infectives','FontSize',fs);
set(gca,'FontSize',fs);
axis([0,tend,0,max(Y)])
=======
function [f,g]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,beta,isdual,solvetype,numseed)%,NNbar)
%ZfinalSizeAllMulti2
%isdual: 0=SM, 1=DM, 2=IM
%solvetype: 1=FSC, 2=ODE, 3=SCM
tauend=1;
time=(1:tauend);
lt=length(time);
mu=0;%In ODE code
NN0=NNrep; NN0(NNrep==0)=1;
Nages=NNbar./NN0;
t0=0; tend=360;

R0=1.8;

%
%Theta:
eps=.13;
%Prow=poisspdf((0:10),5);
%Prow=binopdf((0:10),45,1/10);
Prow=[0,0,0,0,eps,1-eps];%Chaotic?
%Prow=[0,0,0,eps,1-2*eps,eps];
%Prow=[0,0,0,0,0,eps,1-2*eps,eps];
%Prow=1/11*ones(1,11);
%
lp=length(Prow);
Prow=Prow/sum(Prow);
Pmat=repmat(Prow,nbar,1);
v=[5,14,45,16]; vover=1./v; V=kron(vover',ones(n,1)); V=repmat(V,1,lp);
w=[0,5,14,45]; woverv=w.*vover; W=kron(woverv',ones(n,1)); W=repmat(W,1,lp);
%
cen=0; ages0=0; %maxN=0;
A1=zeros(n,lt);
A2=A1;
N0=NN; N0(NN==0)=1;
%
Ni=repmat(NNrep,1,nbar); Nj=Ni';
Niover=1./Ni; Niover(Ni==0)=1; Njover=Niover';
Mj=(Kbar')*NNbar;%C in denom?? .*Cbar %(Kbar')*
Mj(Mj==0)=1;
Mjover=1./Mj;
Mjover=repmat(Mjover',nbar,1);
if isdual==0
    D=Kbar.*Mjover.*Cbar.*Nj;
elseif isdual==1
    D=(K1.*Mjover)*(Kbar'.*Cbar);
    D=D.*Nj;
elseif isdual==2
    D=Kbar'.*Niover.*Cbar.*Nj;
end
%%
%Brand=rand(nbar,1); %Brand=Brand./repmat(sum(Brand,2),1,5);
B=zeros(nbar,lp);
%B(:,2)=.1*rand(nbar,1);%Non-trivial IC
%B(:,2)=.2*Brand;
%%
Z0=sum(B(:,2:end),2)*(1-mu);%=0;
%For final sizes (if IC not in loop):
%IC=XODEsolveAllMulti(gamma,NN,n,nbar,NNbar,NNrep,NN0,maxN,cen,D,Z0,beta,ages0,1); IC=IC./NN0;
%IC=.5*.2*rand(nbar,1);
%%
zn=zeros(nbar,1);
phi1=1; phi2=0;
%NNprob=NNbar/sum(NN);
%NNprob=ones(nbar,1)/sum(NNbar);
seed=numseed;%*NNprob;

thresh=1-1/R0;

for t=1:lt    
    if solvetype==1
    %Final size:
    addbit=0;%seed;%tend*
    %If seed=0, no epidemic happens - because use ODE solver for initial condition
    IC=XODEsolveAllMulti(gamma,NN,n,nbar,NNbar,NNrep,NN0,minNind,maxNind,D,Z0,beta,ages0,t,t0,tend,zn,phi1,phi2,seed,2,thresh); IC=IC./NN0;%New IC here (if ever necessary)
    options = optimset('Display','off');
    funt=@(Zi)solveZi(Zi,Z0,beta,gamma,D,Nages,addbit);
    Zsol=fsolve(funt,IC,options);
    %
    else%solvetype=2/3
    Zsol=XODEsolveAllMulti(gamma,NN,n,nbar,NNbar,NNrep,NN0,minNind,maxNind,D,Z0,beta,ages0,t,t0,tend,zn,phi1,phi2,seed,solvetype,thresh);
    Zsol=Zsol./NN0;
    end
    nu=Zsol-Z0;
    A1(:,t)=sum(reshape(Zsol,n,na),2);%Prop immune for spatial cell (before antigenic drift)
    A2(:,t)=sum(reshape(nu,n,na),2);%AR for spatial cell
    %
    %Assign immunity:
    Bhat=[B(:,2:end),zeros(nbar,1)];
    B=Bhat+repmat(nu,1,lp).*Pmat;
    %Age the populations:
    %
    BVtake=B.*V; Bsus=sum(BVtake(1:n,2:end),2);
    BVadd=[zeros(n,lp);BVtake(1:3*n,:)].*W;
    B=B-BVtake+BVadd; B(1:n,1)=B(1:n,1)+Bsus;
    %}
    Z0=sum(B(:,2:end),2)*(1-mu);
end
%%
%}
f=A1;
g=A2;
end

function f=solveZi(Zi,Z0,beta,gamma,D,Nages,addbit)
lZi=log(Nages-Zi); lZi(Zi>1)=NaN;
f=lZi-log(Nages-Z0)+beta/gamma*D*(Zi-Z0)+addbit;
end
%%
function f=XODEsolveAllMulti(gamma,NN,n,nbar,NNbar,NNrep,NN0,minNind,maxNind,D,ic,beta,ages0,tau,t0,tend,zn,phi1,phi2,seed,solvetype,thresh)
icR=ic.*NNrep;
y0=[NNbar-icR;zn;icR];
cond=sum(icR<thresh);
if cond==0
    seed=0;
end
%
if solvetype==2
[tout,yout]=ode45(@(t,y)integr8all(t,y,beta,gamma,n,nbar,NN,NN0,D,seed,phi1,phi2),[t0,tend],y0);
%Incidence curve in here:
%{
figure
%fs=15; lw=2
fs=12; lw=2;
Y=yout(:,nbar+1:2*nbar); %Z=yout(:,2*nbar+1:3*nbar);
%Y=sum(Y,2);
Y1=Y(:,minNind)+Y(:,minNind+n)+Y(:,minNind+2*n)+Y(:,minNind+3*n); Y1=Y1/NN(minNind);
Y2=Y(:,maxNind)+Y(:,maxNind+n)+Y(:,maxNind+2*n)+Y(:,maxNind+3*n); Y2=Y2/NN(maxNind);
hold on
plot(tout,Y1,'--','linewidth',lw,'color',[.165,.31,.431]);%[.165,.31,.431][.447,.553,.647]
plot(tout,Y2,'-','linewidth',lw,'color',[.165,.31,.431]);
xlabel('Time (days)','FontSize',fs);
ylabel('Infectives','FontSize',fs);
set(gca,'FontSize',fs);
maxY=max(max([Y1;Y2]));
axis([0,tend,0,maxY])
legend('Min','Max','location','NE')
hold off
%}
rec0=yout(end,2*nbar+1:end);%+yout(end,nbar+1:2*nbar);%******** Truncation - add infectious to rercovered?
f=rec0';

elseif solvetype==3
%ic1=round(NNbar/1000); neg=y0(1:nbar)-ic1; neg(neg>0)=0; ic1=ic1+neg;
%ic1=round(y0(1:nbar)/500);
y1=round(y0);%+[-ic1;ic1;zn];
f=stochSim(y1,beta,gamma,n,nbar,NN,NN0,D,seed,phi1,phi2,tau);
%}
end
end
%%
function f=integr8all(t,y,beta,gamma,n,nbar,NN,NN0,D,seed,phi1,phi2)
mu=1/1800;%5*360=1800
alpha=1;
%phi=phi1-phi2*cos(pi*t/180);
phi=1;%+sin(pi*t/180);
S=y(1:nbar);
I=y(nbar+1:2*nbar);
R=y(2*nbar+1:end);
seed1=seed*S./NN0;%*exp(-t);
Sfoi=phi*(beta*S.*(D*(I./NN0)).^alpha+seed1);
Sdot=-Sfoi;%+mu*R;
Idot=Sfoi-gamma*I;
Rdot=gamma*I;%-mu*R;
f=[Sdot;Idot;Rdot];
end
%%
function f=stochSim(y,beta,gamma,n,nbar,NN,N0,D,seed,phi1,phi2,tau)
%Feed in mu if required
factor=24;
tend=360*factor; beta=beta/factor; gamma=gamma/factor;
Vec=zeros(nbar,tend);
%
S=y(1:nbar);
I=y(nbar+1:2*nbar);
R=y(2*nbar+1:end);

%Test seed:
%[Nmax,maxInd]=max(NN);
%ii=4*(maxInd-1)+3;
%seed=0; I(ii)=10; S(ii)=S(ii)-10;
%
%I=floor(N0/10000); S=S-I;

%Different from here:
i=1;
threshold=30;%Number of time-steps for necessary simulation/seed
while i<tend && (i<threshold || sum(I)>0)%At least 30 time-steps
phi=1;%phi1-phi2*cos(pi*i*f1/180);%Just seed here
Sout=1-exp(-phi*(beta*(D*(I./N0))+seed*heaviside(threshold-i)));%+mu*R;
Sout=binornd(S,Sout); Sout(S==0)=0;
S=S-Sout; I=I+Sout;
%
Iout=1-exp(-gamma);
Iout=binornd(I,Iout);
I=I-Iout; R=R+Iout;
Vec(:,i)=I;
i=i+1;
end

Vsum=sum(Vec,1);%/sum(NN);
maxV=max(Vsum);
fs=15; lw=2;
figure
plot(1:tend/factor,Vsum(factor:factor:end),'-','linewidth',2,'color',[.447,.553,.647]);
xlabel('Time (days)','FontSize',fs);
ylabel('Infectives','FontSize',fs);
axis([0,tend/factor,0,maxV]); xlabel('Time (years)'); ylabel('Incidence'); set(gca,'fontsize',15)

f=R;
if sum(isnan(R))>0
    fuck=1;
end
end

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

figure
fs=15;
Y=yout(:,nbar+1:2*nbar);
Y=sum(Y,2);
plot(tout,Y,'-','linewidth',2,'color',[.447,.553,.647]);
xlabel('Time (days)','FontSize',fs);
ylabel('Infectives','FontSize',fs);
set(gca,'FontSize',fs);
axis([0,tend,0,max(Y)])
>>>>>>> acfad3e4da850e525286e08138beef9f06a8c2cc
%}