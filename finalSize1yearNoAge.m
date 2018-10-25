function [f,g,z]=finalSize1yearNoAge(gamma,n,NN,K,betaS,betaI,betaD,isdual,solvetype,sample,glob)%,numseed)
%This code is set up to output incidence - modification required for
%prevalence
%isdual: 0=SM, 1=DM, 2=IM
%glob=1 to output/plot global incidence
%solvetype: 1=FSC NOT AN OPTION, 2=ODE, 3=SCM
plotInc=0;%=1 to plot incidence
%
seed=10^(-5);%Scalar - proportion of local population
%Vector of seed locations:
seedhere=ones(n,1);%Seed everywhere
%
%seedhere=zeros(n,1);
%[maxn,maxind]=max(NN);
%seedhere(sample)=1;%Seed in most populated cell
%
%seedhere=sample;%Input seed locations
%
%Mobility assumption:
if isdual==0
    beta=betaS;
elseif isdual==1
    beta=betaD;
else
    beta=betaI;
end
t0=0; tend=360;
mu=0;%In ODE code
phi1=1; phi2=0;
NN0=NN; NN0(NN==0)=1;
%
Ni=repmat(NN,1,n); Nj=Ni';
Niover=1./Ni; Niover(Ni==0)=1;
Mj=K'*NN;%C in denom?? .*Cbar %(Kbar')*
Mj(Mj==0)=1;
Mjover=1./Mj;
Mjover=repmat(Mjover',n,1);
if isdual==0
    D=K.*Mjover.*Nj;
elseif isdual==1
    D=(K.*Mjover)*K';
    D=D.*Nj;
elseif isdual==2
    D=K'.*Niover.*Nj;
end
%NNprob=NNbar/sum(NN); NNprob=ones(nbar,1)/sum(NNbar);
seedvec=seed*seedhere;

zn=zeros(n,1);
Z0=zn;%Proportion immune at start
[f,g,z]=XODEsolveAllMulti(gamma,NN,n,NN0,D,Z0,beta,t0,tend,zn,phi1,phi2,seed,solvetype,seedvec,glob);
if plotInc==1
figure
fs=12; lw=2;
plot(f,g,'-','linewidth',lw)%,'color',[.447,.553,.647]);
maxY=max(max(g));
axis([0,f(end),0,maxY]);
xlabel('Time (days)'); ylabel('Incidence'); set(gca,'fontsize',fs)
grid on
grid minor
box on
end
end
%%
function [f,g,z]=XODEsolveAllMulti(gamma,NN,n,NN0,D,ic,beta,t0,tend,zn,phi1,phi2,seed,solvetype,seedvec,glob)
icR=ic.*NN;
y0=[NN-icR;zn;icR];
%
if solvetype==2
[tout,yout]=ode45(@(t,y)integr8all(t,y,beta,gamma,n,NN,NN0,D,seed,phi1,phi2,seedvec),[t0,tend],y0);
%
%Y=yout(:,n+1:2*n); Y=sum(Y,2);
Y=yout(:,1:n);
Y=-diff(Y,1,1); tdiff=diff(tout); Y=Y./repmat(tdiff,1,n); tout(1)=[]; Y=movmean(Y,3,1);
%dt=mean(diff(tout)); dy=-gradient(Y,dt); tout=dt*(1:length(dy)); Y=dy;
if glob==1
    Y=sum(Y,2);
end
f=tout;
g=Y;
z=sum(yout(end,n+1:end))/sum(NN);
%
elseif solvetype==3
y1=round(y0);
[f1,g1,z1]=stochSim(y1,beta,tend,gamma,n,NN,NN0,D,seed,phi1,phi2);
inc=movmean(g1,3,1);
inc=diff(inc,1,1)*(f1(2)-f1(1));
if glob==1
    inc=sum(inc,2);
end
f=f1(2:end);
g=inc;
z=z1;
end
end
%
function f=integr8all(t,y,beta,gamma,n,NN,NN0,D,seed,phi1,phi2,seedvec)
%mu=1/1800;%5*360=1800
%phi=phi1-phi2*cos(pi*t/180);
phi=1;
S=y(1:n);
I=y(n+1:2*n);
R=y(2*n+1:end);
%seed1=seed*S./NN0;%*exp(-t);
seed1=seedvec.*S./NN0;%*exp(-t);
Sfoi=phi*(beta*S.*(D*(I./NN0))+seed1);
Igam=gamma*I;
Sdot=-Sfoi;%+mu*R;
Idot=Sfoi-Igam;
Rdot=Igam;%-mu*R;
f=[Sdot;Idot;Rdot];
end
%%
function [f,g,z]=stochSim(y,beta,tendin,gamma,n,NN,N0,D,seed,phi1,phi2)
%Feed in mu if required
factor=6;
tend=tendin*factor; beta=beta/factor; gamma=gamma/factor;
Vec=zeros(n,tend);
%
S=y(1:n);
I=y(n+1:2*n);
R=y(2*n+1:end);
%Different from here:
i=1;
threshold=30;%Number of time-steps for necessary simulation/seed
while i<tend && (i<threshold || sum(I)>0)%At least 30 time-steps
phi=1;%phi1-phi2*cos(pi*i*f1/180);%Just seed here
Sout=1-exp(-phi*(beta*(D*(I./N0))+seed*heaviside(threshold-i)));%+mu*R;.^alpha
Sout(Sout>1)=1;
Sout=binornd(S,Sout); Sout(S==0)=0;
S=S-Sout; I=I+Sout;
%
Iout=1-exp(-gamma);
Iout=binornd(I,Iout);
I=I-Iout; R=R+Iout;
Vec(:,i)=-S;%For incidence
i=i+1;
end
tend=i-1;
%
%Vsum=sum(Vec(:,1:tend),1);%/sum(NN);
%maxV=max(Vsum);
tvec=(1:tend)/factor;
f=tvec';
v=Vec(:,1:tend); g=v';
z=(sum(I)+sum(R))/sum(NN);
if sum(isnan(R))>0
    fuck=1;
    print('Somenting is NaN')
end
end