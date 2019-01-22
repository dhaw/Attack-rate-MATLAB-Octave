function f=plotARages(F,NN,NNbar,a,b)%,tit)
%Trim:
%{
NN=trimbyk(NN,a,b); n=length(NN); logN=log10(NN);
num=size(F,2); F2=zeros(n,num);
%{
for i=1:num
    F2(:,i)=trimbyk(F(:,i),a,b);
end
F=F2;
%}
meanA=nansum(F2(:,1).*NN)./sum(NN);
%}
n=length(NN);

meanA=sum(F.*repmat(NN,4,1))/sum(NN);
f1=F(1:n,:).*NN./NNbar(1:n);
f2=F(n+1:2*n,:).*NN./NNbar(n+1:2*n);
f3=F(2*n+1:3*n,:).*NN./NNbar(2*n+1:3*n);
f4=F(3*n+1:end,:).*NN./NNbar(3*n+1:end);
maxN=max(NN);

figure
%suptitle(tit)
%fs=40; ms=30; lwx=3; lw=3; %For presentations
fs=12; ms=5;%12; 20 %fs=18 for paper figures
lwx=1; lw=1.5; %lwx=1.5
col1=[.165,.31,.431];%[1,0,0];%[0.0512,0.4600,0.8633];%[.165,.31,.431];
col2=[.447,.553,.647];
col3=[0,0,0.5];
%col4=[zeros(4,2),[0;.25;.5;.75]];
col4=lines(4);
h1=semilogx(NN,f1,'o','color',col4(1,:),'markersize',ms,'LineWidth',lwx);%,'markerfacecolor',col3);
hold on
h2=semilogx(NN,f2,'o','color',col4(2,:),'markersize',ms,'LineWidth',lwx);%,'markerfacecolor',col3);
h3=semilogx(NN,f3,'o','color',col4(3,:),'markersize',ms,'LineWidth',lwx);%,'markerfacecolor',col3);
h4=semilogx(NN,f4,'o','color',col4(4,:),'markersize',ms,'LineWidth',lwx);%,'markerfacecolor',col3);
%semilogx([1,maxN],[meanA,meanA],'--','linewidth',lw,'color',[.5,0,0]);
legend([h1,h2,h3,h4],'0-4','5-19','20-64','65+','location','SE')

Aall=[f1;f2;f3;f4];
maxY=max(Aall); minY=min(min(Aall));

xlabel('Pop. density','FontSize',fs);
ylabel('Attack rate')
%xlabel('Proportion aged 65+','FontSize',fs);
%ylabel(strcat('Attack rate (',tit,')'),'FontSize',fs);
set(gca,'FontSize',fs,'xtick',[1,10,100,1000,10000,100000]);
minY=max(minY-.1,0); maxY=min(maxY+.1,1);
axis([0,maxN,0,maxY])%minY,maxY])
%set(gca,'yticklabels','')
grid on
grid minor
box on
hold off
end

%%
%{
figure; hold on; plot([meanfc,meanfc],[1000,1000],'--','color',[.5,0,0],'linewidth',2); h1=hist(fc); xlabel('Attack rate'), ylabel('Frequency'); set(gca,'fontsize',12); grid on; grid minor; box on
hist(fc,20)
plot([meanfc,meanfc],[0,1000],'--','color',[.5,0,0],'linewidth',2);
axis([.37,.45,0,800]))
%}