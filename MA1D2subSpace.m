function [f,g]=MA1D2subSpace(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,isdual)
%Parameters:
R0=1.8;%XXXX Feed in
eps=.5;
reduction=.2;%.61915634;
if isdual==0%heterogeneity
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
tauend=20;
time=(1:tauend);
lt=length(time);
t0=0; tend=360;
mu=0;
%
NN0=NNrep; NN0(NNrep==0)=1;%h
Nages=NNbar./NN0;%h
%
phi1=1; phi2=0;
numseed=10^(-8); seed=numseed;%*NNprob;
%
%Theta:
Prow1=[0,0,0,0,eps,1-eps];
%Prow1=[0,.1*ones(1,10)];
lp1=length(Prow1);
Prow1=Prow1/sum(Prow1);
%eps=.1;
Prow2=Prow1;%[0,0,0,0,0,1-eps,eps];
%Prow2=[0,.1*ones(1,10)];
lp2=length(Prow2);
Prow2=Prow2/sum(Prow2);
v=[5,14,45,16]; vover=1./v; V=kron(vover',ones(n,1)); V=repmat(V,1,lp1);%h
w=[0,5,14,45]; woverv=w.*vover; W=kron(woverv',ones(n,1)); W=repmat(W,1,lp1);%h
%%
%h{
Ni=repmat(NNrep,1,nbar); Nj=Ni';
Niover=1./Ni; Niover(Ni==0)=1; Njover=Niover';
Mj=(Kbar')*NNbar;
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
%h}
%%
A4=[1,1,1,1;1,0,1,0;0,1,0,1]; B4=[1,0,1,0;0,1,0,1;1,0,1,0;0,1,0,1];%
%h{
repD=repmat(D,4,4);
Nnew=4*n*na;%Age/space/immune class
A=kron(A4,D);%.*repD; is repD necessary here? Also, got dimension mismatch error - odd
B=kron(B4,D);%.*repD;
%h}
%Z0=zeros(3,1);%Initial condition of immunity
S0=[1;0;0];
S0=kron(S0,NNbar)/sum(S0);%Number sus to each distinct set of subtypes %h: bar
Shat=[S0(1:nbar);S0];%This is Shat(0)
zn=zeros(length(Shat),1);

A1=zeros(4*nbar,lt);%h
A2=zeros(2*nbar,lt); A2(:,1)=zeros(2*nbar,1);
%b1=zeros(1,lp1); b2=zeros(1,lp2);
b1=S0(2)*[Prow1(2:end),0]; b2=S0(3)*[Prow2(2:end),0];%h

amod=1/amax;%(amax-1)/amax;
options=odeset('refine',1);
xoverFcn=@(t,y)evZero(t,y,nbar); 
options=odeset(options,'Events',xoverFcn);
thresh=1-1/R0;
cross=1-reduction;
NNrep3=repmat(NNrep,3,1); NNrep4=[NNrep3;NNrep];
for tau=1:lt
%Only seed if below herd immunity threshold - otherwise integrators are
%temperamental:%h - taken out
%{
s1=sum(S0(1:nbar)); s2=sum(S0(nbar+1:2*nbar)); s3=sum(S0(2*nbar+1:3*nbar)); Ntot=sum(NN);
thresh1=1-(s1+cross*s2)/Ntot-thresh; thresh1(thresh1>0)=0; thresh1(thresh1<0)=1;
thresh2=1-(s1+cross*s3)/Ntot-thresh; thresh2(thresh2>0)=0; thresh2(thresh2<0)=1;
seeddot=[max(thresh1+thresh2);thresh1;thresh2];
seeddot=kron(seeddot,ones(nbar,1));
%}
seeddot=1;
y0=[S0;zn;zn];%XXXX
[tout,yout]=ode113(@(t,y)integr8all(t,y,beta,gamma,nbar,A,B,NN,NNrep3,NNrep4,seed,phi1,phi2,tau,seeddot,cross),[t0,tend],y0,options);
    Rt=yout(end,7*nbar+1:end);%XXXX
    
    %Rt(Rt<0)=0;
    
    Rt=Rt';
    Rt1=Rt(1:nbar); Rt2=Rt(nbar+1:2*nbar); Rt3=Rt(2*nbar+1:3*nbar); Rt4=Rt(3*nbar+1:end);
    A1tx=S0-[Rt1+Rt2;Rt3;Rt4]+[zeros(nbar,1);Rt2;Rt1]; 
    toSum=reshape(A1tx,nbar,3); sumA1tx=sum(toSum,2);
    
    A1t=[A1tx;NNbar-sumA1tx];
    A2t=[Rt1+Rt3;Rt2+Rt4];
    A2(:,tau)=A2t;
    
    b1=b1+(Rt(1)+Rt(3))*Prow1;
    b2=b2+(Rt(2)+Rt(4))*Prow2;
    w1=b1(1); b1=[b1(2:end),0];
    w2=b2(1); b2=[b2(2:end),0];
    %%
    %digits(100)
    %
    %A1t1=A1t(1:nbar); 
    A1t2=Rt(nbar+1:2*nbar); A1t3=Rt(2*nbar+1:3*nbar); A1t4=Rt(3*nbar+1:end);
    A1t34div=A1t3+A1t4; A1t34div(A1t34div==0)=1;%XXXX Added
    A1t24div=A1t2+A1t4; A1t24div(A1t24div==0)=1;%XXXX Added
    %Waning - new method - get negative S_00:
    from01=A1t3./A1t34div*w1; %from00s1=w1-from01;
    from00s1=A1t4./A1t34div*w1;
    from10=A1t2./A1t24div*w2; %from00s2=w2-from10;
    from00s2=A1t4./A1t24div*w2;
    %from00both=from00s1*from00s2;
    A1t4div=A1t4; A1t4div(A1t4==0)=1;%XXXX Added
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
    %Age the populations:
    if demog==1
        rem=A1t(nbar+1:end)*amod;
        toSum=reshape(rem,nbar,3); sumrem=sum(toSum,2);
        A1t=A1t+[sumrem;-rem];
        b1=b1*(1-amod); b2=b2*(1-amod);
    end
    %Update S0 and Shat:
    S0=A1t(1:3*nbar);%*NN;%Assume duration of immunity unrelated to cross infection
    Shat=[S0(1:nbar);S0];
    %Optional different output - susceptibility at start of season:
    %Need to comment out "A1(:,t)=..." above
    toSum=reshape(S0,nbar,3); sumS0=sum(toSum,2);
    A1(:,tau)=[S0;NNbar-sumS0];
end

f=A1;
g=A2;
end

%%

function f=integr8all(t,y,beta,gamma,nbar,A,B,NN,NNrep3,NNrep4,seed,phi1,phi2,tau,seeddot,cross)
mu=1/1800;%5*360=1800
phi=1;%phi1-phi2*cos(pi*t/180);
%S=y(1:3).*[1;cross;cross];%XXXX possibility of cross12 and cross21
S=[y(1:nbar);cross*y(nbar+1:2*nbar);cross*y(2*nbar+1:3*nbar)];
Shat=[S(1:nbar);S];
I=y(3*nbar+1:7*nbar);%XXXX
seedf=seed*S./NNrep3*exp(-t).*seeddot;%.* if seeddot is a vector (non-trivial case)
seedfhat=[seedf(1:nbar);seedf];
Iscaled=I./NNrep4;
Sdot=-beta*S.*(A*Iscaled)*phi-seedf;%+mu*R;
Idot=beta*Shat.*(B*Iscaled)*phi+seedfhat-gamma*I./NNrep4;
Rdot=gamma*I./NNrep4;%-mu*R;
f=[Sdot;Idot;Rdot];
end

%%

function [value,isterminal,direction]=evZero(t,y,nbar)
value=y(3*nbar+1:7*nbar);
isterminal=ones(4*nbar,1);
direction=-ones(4*nbar,1);
end

%% Junk

%{
    %Modify exidting Bs to account for new infections by other strain:
    b1move=Rt(4)/NN*B1(1,:);
    B1=B1+[-b1move;b1move];
    b2move=Rt(3)/NN*B2(1,:);
    B2=B2+[-b2move;+b2move];
    %
    %%Assign immunity:
    %Normalise to deal with proportions
    %A1t=A1t/NN;
    %A2t=A2t/NN;
    %Add new infections to B:
    rec=[Rt(1)*Prow1;Rt(2)*Prow2;Rt(3)*Prow1;Rt(3)*Prow2];
    B1=B1+[rec(4);-rec(4)+rec(1)];
    B2=B2+[rec(3);-rec(3)+rec(2)];
    b1=B1(:,1); b2=B2(:,1);
    B1=[B1(:,2:end),zeros(2,1)];
    B2=[B2(:,2:end),zeros(2,1)];
    %Wain from 1:
    A1t(1)=A1t(1)+b1(1);%S_11
    A1t(3)=A1t(3)+b1(2)
    %}

%Optional plot:%h - out ********16/04 to here
%{
if tau<=10
figure
%fs=15; lw=2
fs=12; lw=2;
Y=yout(:,4:7); %Z=yout(:,2*nbar+1:3*nbar);%XXXX
%Y=sum(Y,2);
hold on
plot(tout,Y,'-','linewidth',lw);%[.165,.31,.431][.447,.553,.647]
xlabel('Time (days)','FontSize',fs);
ylabel('Infectives','FontSize',fs);
set(gca,'FontSize',fs);
maxY=max(max(max(Y)))+.01;
axis([0,tend,-maxY,maxY])
legend('I_{11}^1','I_{11}^2','I_{10}^1','I_{01}^2','location','NE')
hold off
end
%}