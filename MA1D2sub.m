<<<<<<< HEAD
function [f,g]=MA1D2sub(eps)
%Parameters:
%solvetype: ODE only
R0=1.8; gamma=1/2.6;
NN=1; n=length(NN);
demog=1;%Ageing: 1=on
amax=80;
tauend=200;
time=(1:tauend);
lt=length(time);
t0=0; tend=360;
mu=0;
phi1=1; phi2=0;
%eps=.5;
%NNprob=NNbar/sum(NN); %NNprob=ones(nbar,1)/sum(NNbar);
numseed=10^(-8); seed=numseed;%*NNprob;
%
%Theta:
Prow1=[0,0,0,0,eps,1-eps];
%Prow1=[0,.1*ones(1,10)];
lp1=length(Prow1);
Prow1=Prow1/sum(Prow1);
%eps=.1;
Prow2=[0,0,0,0,0,1-eps,eps];
%Prow2=[0,.1*ones(1,10)];
lp2=length(Prow2);
Prow2=Prow2/sum(Prow2);
%%

beta=R0*gamma;
reduction=.3;%.61915634;
%%

A=[1,1,1,1;1,0,1,0;0,1,0,1]; B=[1,0,1,0;0,1,0,1;1,0,1,0;0,1,0,1];%XXXX
%Z0=zeros(3,1);%Initial condition of immunity
S0=[1;0;0];
%S0=[0.1168
%    0.3232
%    0.2614];
S0=S0/sum(S0)*NN;%Number sus to each distinct set of subtypes
Shat=[S0(1);S0];%This is Shat(0)
zn=zeros(length(Shat),1);

A1=zeros(4,lt); %A1(:,1)=[S0;NN-sum(S0)]; %If include keep intial condition
A2=zeros(2,lt); A2(:,1)=[0;0];
%b1=zeros(1,lp1); b2=zeros(1,lp2);
b1=S0(2)*[Prow1(2:end),0]; b2=S0(3)*[Prow2(2:end),0];

amod=1/amax;%(amax-1)/amax;
options=odeset('refine',1);
options=odeset(options,'Events',@evZero);
thresh=1-1/R0;
cross=1-reduction;
for tau=1:lt
%Only seed if below herd immunity threshold - otherwise ibntegrators are
%temperamental:
thresh1=1-(S0(1)+cross*S0(2))/NN-thresh; thresh1(thresh1>0)=0; thresh1(thresh1<0)=1;
thresh2=1-(S0(1)+cross*S0(3))/NN-thresh; thresh2(thresh2>0)=0; thresh2(thresh2<0)=1;
seeddot=[max(thresh1+thresh2);thresh1;thresh2];

y0=[S0;zn;zn];%XXXX
[tout,yout]=ode113(@(t,y)integr8all(t,y,beta,gamma,A,B,NN,seed,phi1,phi2,tau,seeddot,cross),[t0,tend],y0,options);
%Optional plot:
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
    Rt=yout(end,8:end);%XXXX
    Rt=Rt';
    A1tx=S0-[Rt(1)+Rt(2);Rt(3);Rt(4)]+[0;Rt(2);Rt(1)]; 
    
    A1t=[A1tx;NN-sum(A1tx)];
    A2t=[Rt(1)+Rt(3);Rt(2)+Rt(4)];
    %A1(:,tau)=A1t;
    A2(:,tau)=A2t;
    
    b1=b1+(Rt(1)+Rt(3))*Prow1;
    b2=b2+(Rt(2)+Rt(4))*Prow2;
    w1=b1(1); b1=[b1(2:end),0];
    w2=b2(1); b2=[b2(2:end),0];
    %%
    digits(100)
    %
    %Waning - new method - get negative S_00:
    from01=A1t(3)/(A1t(3)+A1t(4))*w1; %from00s1=w1-from01;
    from00s1=A1t(4)/(A1t(3)+A1t(4))*w1;
    from10=A1t(2)/(A1t(2)+A1t(4))*w2; %from00s2=w2-from10;
    from00s2=A1t(4)/(A1t(2)+A1t(4))*w2;
    %from00both=from00s1*from00s2;
    frac1=from00s1/A1t(4); frac2=from00s2/A1t(4);
    from00both=frac1*frac2*A1t(4);
    from00both(isnan(from00both)==1)=0;
    from00s1=from00s1-from00both;
    from00s2=from00s2-from00both;
    A1t=A1t+[from10+from01+from00both;
             from00s1-from10;
             from00s2-from01;
             -from00both-from00s1-from00s2];
    %}
    %{
    from01=A1t(3)/(A1t(3)+A1t(4))*w1; %from00s1=w1-from01;
    from00s1=A1t(4)/(A1t(3)+A1t(4))*w1;
    from10=A1t(2)/(A1t(2)+A1t(4))*w2; %from00s2=w2-from10;
    from00s2=A1t(4)/(A1t(2)+A1t(4))*w2;
    %}
    %Waning: first way round:
    %{
    %Wain from 1:
    from01=A1t(3)/(A1t(3)+A1t(4))*w1; from00=w1-from01;
    A1t=A1t+[from01;from00;-from01;-from00];
    %Wain from 2:
    from10=A1t(2)/(A1t(2)+A1t(4))*w2; from00=w2-from10;
    A1t=A1t+[from10;-from10;from00;-from00];
    %}
    %To try other way round:
    %{
    %Wain from 2:
    from10=A1t(2)/(A1t(2)+A1t(4))*w2; from00=w2-from10;
    A1t=A1t+[from10;-from10;from00;-from00];
    %Wain from 1:
    from01=A1t(3)/(A1t(3)+A1t(4))*w1; from00=w1-from01;
    A1t=A1t+[from01;from00;-from01;-from00];
    %}
    digits(32)
    %%
    %Age the populations:
    if demog==1
        rem=A1t(2:4)*amod;
        A1t=A1t+[sum(rem);-rem];
        b1=b1*(1-amod); b2=b2*(1-amod);
    end
    %Update S0 and Shat:
    S0=A1t(1:3);%*NN;%Assume duration of immunity unrelated to cross infection
    Shat=[S0(1);S0];
    %Optional different output - susceptibility at start of season:
    %Need to comment out "A1(:,t)=..." above
    A1(:,tau)=[S0;NN-sum(S0)];
end

f=A1;
g=A2;

%%

fs=15; lw=3; ms=3;
%figure

end

%%

function f=integr8all(t,y,beta,gamma,A,B,NN,seed,phi1,phi2,tau,seeddot,cross)
mu=1/1800;%5*360=1800
phi=1;%phi1-phi2*cos(pi*t/180);
S=y(1:3).*[1;cross;cross];%XXXX
Shat=[S(1);S];
I=y(4:7);%XXXX
seedf=seed*S./NN*exp(-t).*seeddot;
seedfhat=[seedf(1);seedf];
Sdot=-beta*S.*(A*(I./NN))*phi-seedf;%+mu*R;
Idot=beta*Shat.*(B*(I./NN))*phi+seedfhat-gamma*I./NN;
Rdot=gamma*I./NN;%-mu*R;
f=[Sdot;Idot;Rdot];
end

%%

function [value,isterminal,direction]=evZero(t,y)
value=y(4:7);
isterminal=zeros(4,1);
direction=-ones(4,1);
end

%%

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
=======
function [f,g]=MA1D2sub(eps)
%Parameters:
%solvetype: ODE only
R0=1.8; gamma=1/2.6;
NN=1; n=length(NN);
demog=1;%Ageing: 1=on
amax=80;
tauend=200;
time=(1:tauend);
lt=length(time);
t0=0; tend=360;
mu=0;
phi1=1; phi2=1;
%eps=.5;
%NNprob=NNbar/sum(NN); %NNprob=ones(nbar,1)/sum(NNbar);
numseed=10^(-8); seed=numseed;%*NNprob;
%
%Theta:
Prow1=[0,0,0,0,eps,1-eps];
%Prow1=[0,.1*ones(1,10)];
lp1=length(Prow1);
Prow1=Prow1/sum(Prow1);
%eps=.1;
Prow2=[0,0,0,0,0,1-eps,eps];
%Prow2=[0,.1*ones(1,10)];
lp2=length(Prow2);
Prow2=Prow2/sum(Prow2);
%%

beta=R0*gamma;
reduction=.3;%.61915634;
%%

A=[1,1,1,1;1,0,1,0;0,1,0,1]; B=[1,0,1,0;0,1,0,1;1,0,1,0;0,1,0,1];%XXXX
%Z0=zeros(3,1);%Initial condition of immunity
S0=[1;0;0];
%S0=[0.1168
%    0.3232
%    0.2614];
S0=S0/sum(S0)*NN;%Number sus to each distinct set of subtypes
Shat=[S0(1);S0];%This is Shat(0)
zn=zeros(length(Shat),1);

A1=zeros(4,lt); %A1(:,1)=[S0;NN-sum(S0)]; %If include keep intial condition
A2=zeros(2,lt); A2(:,1)=[0;0];
%b1=zeros(1,lp1); b2=zeros(1,lp2);
b1=S0(2)*[Prow1(2:end),0]; b2=S0(3)*[Prow2(2:end),0];

amod=1/amax;%(amax-1)/amax;
options=odeset('refine',1);
options=odeset(options,'Events',@evZero);
thresh=1-1/R0;
cross=1-reduction;
for tau=1:lt
%Only seed if below herd immunity threshold - otherwise ibntegrators are
%temperamental:
thresh1=1-(S0(1)+cross*S0(2))/NN-thresh; thresh1(thresh1>0)=0; thresh1(thresh1<0)=1;
thresh2=1-(S0(1)+cross*S0(3))/NN-thresh; thresh2(thresh2>0)=0; thresh2(thresh2<0)=1;
seeddot=[max(thresh1+thresh2);thresh1;thresh2];

y0=[S0;zn;zn];%XXXX
[tout,yout]=ode113(@(t,y)integr8all(t,y,beta,gamma,A,B,NN,seed,phi1,phi2,tau,seeddot,cross),[t0,tend],y0,options);
%Optional plot:
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
    Rt=yout(end,8:end);%XXXX
    Rt=Rt';
    A1tx=S0-[Rt(1)+Rt(2);Rt(3);Rt(4)]+[0;Rt(2);Rt(1)]; 
    
    A1t=[A1tx;NN-sum(A1tx)];
    A2t=[Rt(1)+Rt(3);Rt(2)+Rt(4)];
    %A1(:,tau)=A1t;
    A2(:,tau)=A2t;
    
    b1=b1+(Rt(1)+Rt(3))*Prow1;
    b2=b2+(Rt(2)+Rt(4))*Prow2;
    w1=b1(1); b1=[b1(2:end),0];
    w2=b2(1); b2=[b2(2:end),0];
    %%
    digits(100)
    %
    %Waning - new method - get negative S_00:
    from01=A1t(3)/(A1t(3)+A1t(4))*w1; %from00s1=w1-from01;
    from00s1=A1t(4)/(A1t(3)+A1t(4))*w1;
    from10=A1t(2)/(A1t(2)+A1t(4))*w2; %from00s2=w2-from10;
    from00s2=A1t(4)/(A1t(2)+A1t(4))*w2;
    %from00both=from00s1*from00s2;
    frac1=from00s1/A1t(4); frac2=from00s2/A1t(4);
    from00both=frac1*frac2*A1t(4);
    from00both(isnan(from00both)==1)=0;
    from00s1=from00s1-from00both;
    from00s2=from00s2-from00both;
    A1t=A1t+[from10+from01+from00both;
             from00s1-from10;
             from00s2-from01;
             -from00both-from00s1-from00s2];
    %}
    %{
    from01=A1t(3)/(A1t(3)+A1t(4))*w1; %from00s1=w1-from01;
    from00s1=A1t(4)/(A1t(3)+A1t(4))*w1;
    from10=A1t(2)/(A1t(2)+A1t(4))*w2; %from00s2=w2-from10;
    from00s2=A1t(4)/(A1t(2)+A1t(4))*w2;
    %}
    %Waning: first way round:
    %{
    %Wain from 1:
    from01=A1t(3)/(A1t(3)+A1t(4))*w1; from00=w1-from01;
    A1t=A1t+[from01;from00;-from01;-from00];
    %Wain from 2:
    from10=A1t(2)/(A1t(2)+A1t(4))*w2; from00=w2-from10;
    A1t=A1t+[from10;-from10;from00;-from00];
    %}
    %To try other way round:
    %{
    %Wain from 2:
    from10=A1t(2)/(A1t(2)+A1t(4))*w2; from00=w2-from10;
    A1t=A1t+[from10;-from10;from00;-from00];
    %Wain from 1:
    from01=A1t(3)/(A1t(3)+A1t(4))*w1; from00=w1-from01;
    A1t=A1t+[from01;from00;-from01;-from00];
    %}
    digits(32)
    %%
    %Age the populations:
    if demog==1
        rem=A1t(2:4)*amod;
        A1t=A1t+[sum(rem);-rem];
        b1=b1*(1-amod); b2=b2*(1-amod);
    end
    %Update S0 and Shat:
    S0=A1t(1:3);%*NN;%Assume duration of immunity unrelated to cross infection
    Shat=[S0(1);S0];
    %Optional different output - susceptibility at start of season:
    %Need to comment out "A1(:,t)=..." above
    A1(:,tau)=[S0;NN-sum(S0)];
end

f=A1;
g=A2;

%%

fs=15; lw=3; ms=3;
%figure

end

%%

function f=integr8all(t,y,beta,gamma,A,B,NN,seed,phi1,phi2,tau,seeddot,cross)
mu=1/1800;%5*360=1800
phi=1;%phi1-phi2*cos(pi*t/180);
S=y(1:3).*[1;cross;cross];%XXXX
Shat=[S(1);S];
I=y(4:7);%XXXX
seedf=seed*S./NN*exp(-t).*seeddot;
seedfhat=[seedf(1);seedf];
Sdot=-beta*S.*(A*(I./NN))*phi-seedf;%+mu*R;
Idot=beta*Shat.*(B*(I./NN))*phi+seedfhat-gamma*I./NN;
Rdot=gamma*I./NN;%-mu*R;
f=[Sdot;Idot;Rdot];
end

%%

function [value,isterminal,direction]=evZero(t,y)
value=y(4:7);
isterminal=zeros(4,1);
direction=-ones(4,1);
end

%%

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
>>>>>>> acfad3e4da850e525286e08138beef9f06a8c2cc
    %}