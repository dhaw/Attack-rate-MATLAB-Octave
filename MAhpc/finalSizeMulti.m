function [f,g,D]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,isdual,solvetype,numseed,eps,randic,tauend)
%ZfinalSizeAllMulti2
%isdual: 0=SM, 1=DM, 2=IM
%solvetype: 1=FSC, 2=ODE, 3=SCM
%eps=.3;
ICdepAge=1;
if isdual==0
    beta=betaS;
elseif isdual==1
    beta=betaD;
elseif isdual==2
    beta=betaI;
else
    beta=beta3;
end
demog=1;
%tauend=1;
plotTau=1;
time=(1:tauend);
lt=length(time);
t0=0; tend=360;
mu=1/80;%In ODE code
phi1=1; phi2=0;
NN0=NNrep; NN0(NNrep==0)=1;
Nages=NNbar./NN0;
alpha=1;%TSIR/sub-exp parameter
%randic=1;
thetaGeom=0;
rate=eps;
%
%Theta:
%Prow=[0,0,0,1-eps,eps];
Prow=[0,0,0,0,1-eps,eps];
lp=length(Prow);
Prow=Prow/sum(Prow);
Pmat=repmat(Prow,nbar,1);
%v=[5,14,45,16]; vover=1./v; V=kron(vover',ones(n,1)); V=repmat(V,1,lp);
%w=[0,5,14,45]; woverv=w.*vover; W=kron(woverv',ones(n,1)); W=repmat(W,1,lp);
vover=[1/5,1/14,1/45,0];
V=kron(vover',ones(n,1)); V=repmat(V,1,lp);
%
%cen=0; 
ages0=0; %maxN=0;
%
Ni=repmat(NNrep,1,nbar); Nj=Ni';
Niover=1./Ni; Niover(Ni==0)=1;
Mj=(Kbar')*NNbar;%C in denom?? .*Cbar %(Kbar')*
Mj(Mj==0)=1;
Mjover=1./Mj;
Mjover=repmat(Mjover',nbar,1);
if isdual==0
    D=Kbar.*Mjover.*Cbar.*Nj;
elseif isdual==1
    D=(K1.*Mjover)*((Kbar').*Cbar);
    D=D.*Nj;
elseif isdual==2
    D=Kbar'.*Niover.*Cbar.*Nj;
elseif isdual==3
    D=Kbar.*Cbar;
end
%%
%Initial condition:
B=zeros(nbar,lp);
if randic==1
    %{
    maxYears=lp-1;
    R=rand(nbar,maxYears);
    R=round(R); rsum=sum(sum(R)); S=rand(rsum,1);
    R(R==1)=S;
    Rsum=sum(R,2); Rsum(Rsum==0)=1; R=R./repmat(Rsum,1,maxYears);
    T=rand(nbar,1).*NNbar./NN0; Trep=repmat(T,1,maxYears);
    B=[zeros(nbar,1),R.*Trep];
    %}
    %{
    ProwRand=rand(1,lp+1); ProwRand(1)=[]; ProwRand=ProwRand/sum(ProwRand);
    Brand=rand(nbar,1).*NNbar./NN0; %Brand=rand(1)*Brand./repmat(sum(Brand,2),1,lp);
    B=repmat(ProwRand,nbar,1).*repmat(Brand,1,lp);
    %}
    %
    if ICdepAge==1
        maxYears=lp-1;
        Prows=zeros(nbar,maxYears);
        numNonZero=floor((maxYears+1)*rand(nbar,1));%0-6 years - pick whow many are non zero
        for i=1:nbar
            rsi=randsample(maxYears,numNonZero(i));
            Prows(i,rsi)=1;
        end
        findNZ=find(Prows);
        lNZ=length(findNZ);
        randNZ=rand(lNZ,1);
        Prows(findNZ)=randNZ;
        rand2=rand(nbar,1); repRand=repmat(rand2,1,maxYears);
        sumP=sum(Prows,2); repP=repmat(sumP,1,maxYears);
        Nratio=repmat(NNbar./NNrep,1,maxYears);
        B=Prows./repP.*repRand.*Nratio;
        B(isinf(B)==1)=0; B(isnan(B)==1)=0;
        B=[zeros(nbar,1),B];
    else
        maxYears=lp-1;
        Prows=zeros(n,maxYears);
        numNonZero=floor((maxYears+1)*rand(n,1));%0-6 years - pick whow many are non zero
        for i=1:n
            rsi=randsample(maxYears,numNonZero(i));
            Prows(i,rsi)=1;
        end
        findNZ=find(Prows);
        lNZ=length(findNZ);
        randNZ=rand(lNZ,1);
        Prows(findNZ)=randNZ;
        rand2=rand(n,1); repRand=repmat(rand2,1,maxYears);
        sumP=sum(Prows,2); repP=repmat(sumP,1,maxYears);
        B=Prows./repP.*repRand;
        B(isinf(B)==1)=0; B(isnan(B)==1)=0;
        NNratio=repmat(NNbar./NNrep,1,maxYears);
        B=[zeros(nbar,1),repmat(B,4,1).*NNratio];
    end
end
%}
%
Z0=sum(B(:,2:end),2)*(1-mu);

%Z0(3*n+1:end)=eps*NNbar(3*n+1:end)./NN0(3*n+1:end);
%Z0(2*n+1:3*n)=12/45*eps*NNbar(2*n+1:3*n)./NN0(2*n+1:3*n);

%%
zn=zeros(nbar,1);
A1=zeros(n,lt);
A2=A1;
N0=NN; N0(NN==0)=1;
%NNprob=NNbar/sum(NN); NNprob=ones(nbar,1)/sum(NNbar);
seed=10^(-numseed);%*NNprob;
seedvec=zeros(nbar,1); seedvec(2*n+1:3*n)=seed*ones(n,1);
%seedvec=seed*ones(nbar,1);
thresh=0;%Remove from ODE solver
for t=1:lt    
    if solvetype==1
    %Final size:
    addbit=0;%seed;%tend*
    %If seed=0, no epidemic happens - because use ODE solver for initial condition
    IC=XODEsolveAllMulti(gamma,NN,n,nbar,NNbar,NNrep,NN0,minNind,maxNind,D,Z0,beta,ages0,t,t0,tend,zn,phi1,phi2,seed,2,thresh,alpha,seedvec); IC=IC./NN0;%New IC here (if ever necessary)
    options = optimset('Display','off');
    funt=@(Zi)solveZi(Zi,Z0,beta,gamma,D,Nages,addbit);
    Zsol=fsolve(funt,IC,options);
    %
    else%solvetype=2/3
    Zsol=XODEsolveAllMulti(gamma,NN,n,nbar,NNbar,NNrep,NN0,minNind,maxNind,D,Z0,beta,ages0,t,t0,tend,zn,phi1,phi2,seed,solvetype,thresh,alpha,seedvec,plotTau);
    minAttack=min(Zsol); maxAttack=max(Zsol);
%
    if minAttack<0 
        error('Attack rate not in [0,1]')
    end
%}
    Zsol=Zsol./NN0;
    end
    nu=Zsol-Z0;
    A1(:,t)=sum(reshape(Zsol,n,na),2);%Prop immune for spatial cell (before antigenic drift)
    A2(:,t)=sum(reshape(nu,n,na),2);%AR for spatial cell
    %
    if thetaGeom==1
        if demog==1
            Z0=(1-rate)*(1-mu)*Zsol;
        else
            Z0=(1-rate)*Zsol;
        end
    else
    %Assign immunity:
    Bhat=[B(:,2:end),zeros(nbar,1)];
    B=Bhat+repmat(nu,1,lp).*Pmat;
    %Age the populations:
    if demog==1
        %{
        BVtake=B.*V; Bsus=sum(BVtake(1:n,2:end),2);
        BVadd=[zeros(n,lp);BVtake(1:3*n,:)].*W;
        B=B-BVtake+BVadd; B(1:n,1)=B(1:n,1)+Bsus;
        %}
        Btake=B.*V; Badd=circshift(Btake,n,1);
        B(end-n+1:end,:)=B(end-n+1:end,:)*15/16;
        B=B-Btake+Badd;
    end
    Z0=sum(B(:,2:end),2);%*(1-mu);
    end
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
function f=XODEsolveAllMulti(gamma,NN,n,nbar,NNbar,NNrep,NN0,minNind,maxNind,D,ic,beta,ages0,tau,t0,tend,zn,phi1,phi2,seed,solvetype,thresh,alpha,seedvec,plotTau)
icR=ic.*NNrep;
y0=[NNbar-icR;zn;icR];
%{
cond=sum(icR<thresh);
if cond==0
    seed=0;
end
%}
if solvetype==2
    [tout,yout]=ode45(@(t,y)integr8all(t,y,beta,gamma,n,nbar,NN,NN0,D,seed,phi1,phi2,alpha,seedvec),[t0,tend],y0);
    %Incidence curve in here:
    %
    if tau==plotTau
        figure
        fs=12; lw=2;
        Y=yout(:,nbar+1:2*nbar); %Z=yout(:,2*nbar+1:3*nbar);
        %Y=sum(Y,2);
        %Y1=Y(:,minNind)+Y(:,minNind+n)+Y(:,minNind+2*n)+Y(:,minNind+3*n); Y1=Y1/NN(minNind);
        %Y2=Y(:,maxNind)+Y(:,maxNind+n)+Y(:,maxNind+2*n)+Y(:,maxNind+3*n); Y2=Y2/NN(maxNind);
        Ysum=sum(Y,2);
        Yall=Y(:,1:n)+Y(:,n+1:2*n)+Y(:,2*n+1:3*n)+Y(:,3*n+1:end);
        Yall=abs(Yall);%****cheat ;)
        %
        %Unlogged plots:
        hold on
        %plot(tout,Y1,'--','linewidth',lw,'color',[.165,.31,.431]);%[.165,.31,.431][.447,.553,.647]
        %plot(tout,Y2,'-','linewidth',lw,'color',[.165,.31,.431]);
        plot(tout,Ysum,'k','linewidth',lw);
        plot(tout,Yall);
        %}
        %{
        %Logged plots:
        semilogy(tout,Ysum,'k','linewidth',lw);
        hold on
        semilogy(tout,Yall);
        %}
        xlabel('Time (days)','FontSize',fs);
        ylabel('Prevalence','FontSize',fs);
        set(gca,'FontSize',fs);
        maxY=max(Ysum);
        %axis([0,tend,0,maxY])
        axis ([0,tend,0,maxY])
        %legend('Min','Max','location','NE')
        grid on
        grid minor
        hold off
    end
%}
rec0=yout(end,2*nbar+1:end)+yout(end,nbar+1:2*nbar);%******** Truncation - add infectious to rercovered?
f=rec0';

elseif solvetype==3
%ic1=round(NNbar/1000); neg=y0(1:nbar)-ic1; neg(neg>0)=0; ic1=ic1+neg;
%ic1=round(y0(1:nbar)/500);
y1=round(y0);%+[-ic1;ic1;zn];
f=stochSim(y1,beta,gamma,n,nbar,NN,NN0,D,seed,phi1,phi2,tau,alpha,plotTau);
%}
end
end
%%
function f=integr8all(t,y,beta,gamma,n,nbar,NN,NN0,D,seed,phi1,phi2,alpha,seedvec)
mu=1/1800;%5*360=1800
%phi=phi1-phi2*cos(pi*t/180);
phi=1;%+sin(pi*t/180);
S=y(1:nbar);
I=y(nbar+1:2*nbar);
R=y(2*nbar+1:end);
%seed1=seed*S./NN0;%*exp(-t);
seed1=seedvec.*S./NN0;%*exp(-t);
Sfoi=phi*(beta*S.*(D*(I./NN0))+seed1);%seed1);
Sdot=-Sfoi;%+mu*R;
Idot=Sfoi-gamma*I;
Rdot=gamma*I;%-mu*R;
f=[Sdot;Idot;Rdot];
end
%%
function f=stochSim(y,beta,gamma,n,nbar,NN,N0,D,seed,phi1,phi2,tau,alpha,plotTau)
%Feed in mu if required
factor=6;
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
threshold=30*factor;%Number of time-steps for necessary simulation/seed
while i<tend && (i<threshold || sum(I)>0)%At least 30 time-steps
    phi=1;%phi1-phi2*cos(pi*i*f1/180);%Just seed here
    Sout=1-exp(-phi*(beta*(D*(I./N0).^alpha)+seed*heaviside(threshold-i)));%+mu*R;.^alpha
    Sout(Sout>1)=1;
    Sout=binornd(S,Sout); Sout(S==0)=0;
    S=S-Sout; I=I+Sout;
    %
    Iout=1-exp(-gamma);
    Iout=binornd(I,Iout);
    I=I-Iout; R=R+Iout;
    Vec(:,i)=I;
    i=i+1;
end
if tau==plotTau
    figure
    fs=12; lw=2;
    Y=Vec'; %Z=yout(:,2*nbar+1:3*nbar);
    tout=(1:size(Y,1))/factor;
    %Y=sum(Y,2);
    %Y1=Y(:,minNind)+Y(:,minNind+n)+Y(:,minNind+2*n)+Y(:,minNind+3*n); Y1=Y1/NN(minNind);
    %Y2=Y(:,maxNind)+Y(:,maxNind+n)+Y(:,maxNind+2*n)+Y(:,maxNind+3*n); Y2=Y2/NN(maxNind);
    Ysum=sum(Y,2);
    Yall=Y(:,1:n)+Y(:,n+1:2*n)+Y(:,2*n+1:3*n)+Y(:,3*n+1:end);
    Yall=abs(Yall);%****cheat ;)
    %
    %Unlogged plots:
    hold on
    %plot(tout,Y1,'--','linewidth',lw,'color',[.165,.31,.431]);%[.165,.31,.431][.447,.553,.647]
    %plot(tout,Y2,'-','linewidth',lw,'color',[.165,.31,.431]);
    plot(tout,Ysum,'k','linewidth',lw);
    plot(tout,Yall);
    %}
    %{
    %Logged plots:
    semilogy(tout,Ysum,'k','linewidth',lw);
    hold on
    semilogy(tout,Yall);
    %}
    xlabel('Time (days)','FontSize',fs);
    ylabel('Prevalence','FontSize',fs);
    set(gca,'FontSize',fs);
    maxY=max(Ysum);
    %axis([0,tend,0,maxY])
    axis tight%([0,tend,0,maxY])
    %legend('Min','Max','location','NE')
    grid on
    grid minor
    box on
    hold off
end
%{
Vsum=sum(Vec,1);%/sum(NN);
maxV=max(Vsum);
fs=12; lw=2;
Vtot=Vsum(factor:factor:end);
Vall=Vec(1:n,factor:factor:end)+Vec(n+1:2*n,factor:factor:end)+Vec(2*n+1:3*n,factor:factor:end)+Vec(3*n+1:end,factor:factor:end);
%
%Unlogged plots:
figure
plot(1:tend/factor,Vtot,'-','linewidth',2,'color',[.447,.553,.647]);
plot(1:tend/factor,Vall);
axis([0,tend/factor,0,maxV]);
xlabel('Time (days)'); ylabel('Prevalence'); set(gca,'fontsize',fs)
grid on
grid minor
%}
%
%{
%Logged plots:
figure
semilogy(1:tend/factor,Vtot,'k','linewidth',lw)
hold on
semilogy(1:tend/factor,Vall);
axis([0,tend/factor,0,maxV]);
xlabel('Time (days)'); ylabel('Prevalence'); set(gca,'fontsize',fs)
grid on
grid minor
%}
f=R;
if sum(isnan(R))>0
    fuck=1;
    print('Somenting is NaN')
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

figure
fs=15;
Y=yout(:,nbar+1:2*nbar);
Y=sum(Y,2);
plot(tout,Y,'-','linewidth',2,'color',[.447,.553,.647]);
xlabel('Time (days)','FontSize',fs);
ylabel('Infectives','FontSize',fs);
set(gca,'FontSize',fs);
axis([0,tend,0,max(Y)])
%}