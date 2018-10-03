function f=plotDEbar(D,A)

del=(0:.05:1);
%eps=1-del;

ld=length(del);
[a,b]=size(D);
lz=(a-8)*(b-8);

%{
A=zeros(lz,ld);
parfor i=1:ld
    deli=del(i); 
    epsi=0;%del(i);%-deli;
    [gamma,NN,n,nbar,na,NNbar,NNrep,Kout,Kin,K1,Cbar,beta]=prepFluDE(D,R0,deli,epsi);
    zi=finalSize1YDE(gamma,n,nbar,na,NN,NNbar,NNrep,Kout,Kin,K1,Cbar,beta);
    A(:,i)=trimbyk(zi,a,b);
end
f=A;
%NN=trimbyk(NN,a,b);
%}

[a,b]=size(D);
NN=reshape(D,a*b,1);
NN=trimbyk(NN,a,b);
%NN=reshape(D(5:end-4,5:end-4),lz,1);
repN=repmat(NN,1,ld);
meanA=sum(A.*repN,1)./sum(NN);

nvec=log10(NN);
rem=find(NN==0);
nvec(rem)=[];
A(rem,:)=[];

cc=zeros(ld,1);
for i=1:ld
    ccMat=corrcoef(nvec,A(:,i));
    cc(i)=ccMat(1,2);
end

figure
fs=18; ms=5;
lwx=1; lw=1.5;
col1=[.165,.31,.431];%[1,0,0];%[0.0512,0.4600,0.8633];%[.165,.31,.431];
col2=[.447,.553,.647];
labels=cellstr(num2str(del'));
hold on

%yyaxis left
h1=plot(del,meanA','--k','linewidth',lw);%,'color',[.071,.208,.322]
bh=boxplot(A,'labels',labels,'position',del,'colors',[0,0,0],'symbol','.');%[.165,.31,.431][.447,.553,.647]
set(gca,'xtickmode','auto','xticklabelmode','auto')
axis([-.02,1.02,0,1]);%maxA
set(bh,'linewidth',1.5);

%ylabel('Attack rate','FontSize',fs); set(gca,'FontSize',fs)

yyaxis right
h2=plot(del,cc,'-','linewidth',lw);
axis([-.02,1.02,-1,1])
%ylabel('Corr. coef.','FontSize',fs);

legend([h1,h2],'Attack rate','Corr. coef.','location','SE')
grid on
grid minor
xlabel('\delta (\epsilon=1-\delta)','FontSize',fs); 
%\epsilon (\delta=1-\epsilon) %\epsilon (\delta=1-\epsilon)
set(gca,'fontsize',fs)