function [f,g]=CCalpha(A,B)%(A,B) %(D,R0)
alpha=(0:.05:1);
a=(.1:.1:2);
p=(.2:.2:4);
other=p;%********

lother=length(other);
%{
[a,b]=size(D);
%n=a*b;
la=length(alpha);
A=zeros(la,lother);
B=A;
for i=1:la
    aplhai=alpha(i);
    parfor j=1:lother
    otherj=other(j);
    [gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,ages0]=prepFlu(D,1.8,0,aplhai,otherj);
    [fs,x]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,0,2,10^(-8));
    [fi,x]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaI,2,2,10^(-8));
    fs1=trimbyk(fs,a,b);
    fi1=trimbyk(fi,a,b);
    nvec=trimbyk(NN,a,b); fs1(nvec==0)=[]; fi1(nvec==0)=[]; nvec(nvec==0)=[];
    nvec=log10(nvec);
    ccs=corrcoef(nvec,fs1); ccs=ccs(1,2);
    cci=corrcoef(nvec,fi1); cci=cci(1,2);
    A(i,j)=ccs;
    B(i,j)=cci;
    end
end
f=A;
g=B;
%}

%{
fs=18; ms=5;
lwx=1; lw=2;%1.5;
cmap=redblue(100);

figure
colormap(cmap)
imagesc(alpha,other,A')
hcb=colorbar;
caxis([0-1,1])
xlabel('Population power \alpha','FontSize',fs);
ylabel('Offset a')%Offset a Distance power p
set(gca,'FontSize',fs);
axis tight%([0,1,-1,1])
box on
figure
colormap(cmap)
imagesc(alpha,other,B')
hcb=colorbar;
caxis([0-1,1])
xlabel('Population power \alpha','FontSize',fs);
ylabel('Offset a')
set(gca,'FontSize',fs);
axis tight%([0,1,-1,1])
box on
%}

%{
[a,b]=size(D);
%n=a*b;
la=length(alpha);
A=zeros(la,1);
B=A;
for i=1:la
    ai=alpha(i);
    [gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,ages0]=prepFlu(D,1.8,0,ai);
    [fs,x]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,0,2,10^(-8));
    [fi,x]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaI,2,2,10^(-8));
    fs1=trimbyk(fs,a,b);
    fi1=trimbyk(fi,a,b);
    nvec=trimbyk(NN,a,b); fs1(nvec==0)=[]; fi1(nvec==0)=[]; nvec(nvec==0)=[];
    nvec=log10(nvec);
    ccs=corrcoef(nvec,fs1); ccs=ccs(1,2);
    cci=corrcoef(nvec,fi1); cci=cci(1,2);
    A(i)=ccs;
    B(i)=cci;
end
f=A;
g=B;
%}

%
figure
fs=18; ms=5;
lwx=1; lw=2;%1.5;
col1=[.165,.31,.431];
col2=[.447,.553,.647];

hold on
h1=plot(alpha,A,'--','color',col1,'LineWidth',lw);
h2=plot(alpha,B,'--','color',col2,'linewidth',lw);
plot([0,1],[0,0],'k-','linewidth',lwx);
xlabel('Population power \alpha','FontSize',fs);
ylabel('Corr. coef.')
set(gca,'FontSize',fs);
axis ([0,1,-1,1])
legend([h1,h2],'S-mob','I-mob','location','W')
grid on
grid minor
box on
hold off
%}