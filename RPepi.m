function [f,g]=RPepi(gamma,NN,n,K,L,Vi,Viover,beta,kangle,h,w,rhohat,isflat,numseed,sample,boxLat,boxLong)%,loc)
%mu=0;%In ODE code
t0=0; tend=360;
plotto=tend;
logplotmin=-4;
plotmin=10^(logplotmin);
plotmax=10^4;
logplot=1;
incidence=1;%1 for incidence, 0 for prev.
plotfig=0;
outputGlobal=1;%1 for global incidence, 0 (or ither) for local incidence
%
alpha=1;
a=boxLat;
b=boxLong;
fa=floor((a-1)/2); fb=floor((b-1)/2);
%loc=[a*fb+fa,a*fb+fa+1,a*(1+fb)+fa,a*(1+fb)+fa+1];%[5050,5051,5150,5151];%ceil(n/2)-49;%Seed location(s)
loc=1;
Z0=zeros(n,1);
%%
if isflat==1
    Leff=1/rhohat/kangle*L;
else
    Leff=w/kangle*Viover.*L;
end
L=Leff;
%%
zn=zeros(n,1);
phi1=1; phi2=0;
seed=numseed;%*NNprob;
seedvec=zeros(n,1); seedvec(loc)=seed;
[f,g]=XODEsolveAllMulti(gamma,NN,n,L,beta,Z0,t0,tend,zn,phi1,phi2,seed,alpha,seedvec,plotto,sample,logplotmin,plotmin,plotmax,logplot,incidence,plotfig,outputGlobal);
%{
[tout,Zsol]=XODEsolveAllMulti(gamma,NN,n,L,beta,Z0,t0,tend,zn,phi1,phi2,seed,alpha,seedvec,plotto,sample,logplotmin,plotmin,plotmax,logplot,incidence);
Zsol=Zsol./NN;%Prop immune for spatial cell (before antigenic drift)
Zsol(NN==0)=0;
f=Zsol;
%}
end
%%
function [f,g]=XODEsolveAllMulti(gamma,NN,n,L,beta,Z0,t0,tend,zn,phi1,phi2,seed,alpha,seedvec,plotto,sample,logplotmin,plotmin,plotmax,logplot,incidence,plotfig)
icR=Z0.*NN;
y0=[NN-icR;zn;icR];
%
[tout,yout]=ode45(@(t,y)integr8all(t,y,gamma,NN,n,L,beta,seed,phi1,phi2,alpha,seedvec),[t0,tend],y0);
if incidence==1
    Y=yout(:,1:n);
    Y=movmean(Y,3,1);
    Y=-diff(Y).*repmat(diff(tout),1,n);
    tout(1)=[];
else
    Y=yout(:,n+1:2*n); %Z=yout(:,2*nbar+1:3*nbar);
    Y=movmean(Y,5,1);
end
%
Ysum=sum(Y,2);
tend=find(tout>=plotto);
tend=tend(1);
tvec=tout(1:tend);
%
if plotfig==1
%int=floor(n/50);
%sample=randsample(n,100);%(int:int:end);
figure
fs=18; lw=2;
%Y=sum(Y,2);
%Y1=Y(:,minNind)+Y(:,minNind+n)+Y(:,minNind+2*n)+Y(:,minNind+3*n); Y1=Y1/NN(minNind);
%Y2=Y(:,maxNind)+Y(:,maxNind+n)+Y(:,maxNind+2*n)+Y(:,maxNind+3*n); Y2=Y2/NN(maxNind);
%
%Logged plots:
if logplot==1
    semilogy(tvec,Ysum(1:tend,:),'k','linewidth',lw);
    hold on
    semilogy(tvec,Y(1:tend,sample));
    hold off
    maxY=max(Ysum);
    logMaxY=ceil(log10(maxY));
    tickvec=logspace(logplotmin,logMaxY,logMaxY-logplotmin+1);
    maxY=10^logMaxY;
    axis([0,min(tvec(end),plotto),plotmin,maxY])
    set(gca,'ytick',tickvec)
else
    plot(tout,Ysum,'k','linewidth',lw);
    hold on
    plot(tout,Y);
    hold off
    maxY=max(Ysum);
    axis ([0,min(tvec(end),plotto),plotmin,maxY])
end
%
xlabel('Time (days)','FontSize',fs);
if incidence==1
    ylabel('Incidence','FontSize',fs);
else
    ylabel('Prevalence','FontSize',fs);
end
set(gca,'FontSize',fs);
grid on
grid minor
hold off
end
%
rec0=yout(end,2*n+1:end);%+yout(end,nbar+1:2*nbar);%******** Truncation - add infectious to rercovered?
%f=rec0';
f=tout;
if outputGlobal==1
    g=Ysum;
else
    g=Y;
end
%%
function f=integr8all(t,y,gamma,NN,n,L,beta,seed,phi1,phi2,alpha,seedvec)
mu=1/1800;%5*360=1800
%phi=phi1-phi2*cos(pi*t/180);
phi=1;%+sin(pi*t/180);
S=y(1:n);
I=y(n+1:2*n);
R=y(2*n+1:end);
%seed1=seed*S./NN0;%*exp(-t);
seed1=seedvec.*S./NN;%*exp(-t);
Sfoi=phi*(beta*S.*(L*I.^alpha).^alpha+seed1);
Sdot=-Sfoi;%+mu*R;
Idot=Sfoi-gamma*I;
Rdot=gamma*I;%-mu*R;
f=[Sdot;Idot;Rdot];
end