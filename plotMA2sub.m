%function f=plotMAs2sub(Asus,Arec,NN,NNbar,tmin,tmax)
function f=plotMA2sub(Asus,Arec,NN,NNbar,tmin,tmax,fp)

Asus=Asus(:,tmin:tmax);
Arec=Arec(:,tmin:tmax);
%ar=1;

tauend=size(Arec,2);
nbar=length(NNbar); n=length(NN);
R1=Arec(1:nbar,:);
R2=Arec(nbar+1:end,:);
%repNNbar=repmat(NNbar,1,tauend);
repNN=repmat(NN,1,tauend);
%R1=R1.*repNNbar; R2=R2.*repNNbar;
R1=R1(1:n,:)+R1(n+1:2*n,:)+R1(2*n+1:3*n,:)+R1(3*n+1:end,:);
R2=R2(1:n,:)+R2(n+1:2*n,:)+R2(2*n+1:3*n,:)+R2(3*n+1:end,:);
attack=R1+R2; attack=sum(attack,1)/sum(NN);
R1=R1./repNN; R1(repNN==0)=0;
R2=R2./repNN; R2(repNN==0)=0;
f=R1+R2;
maxAttack=max(max(f));
%%
%Comment out if don't want ccpand:
%{
thresh=0.005;
yy=1:tauend;
X=zeros(1,tauend); Y=X; Z=X;
fx=R1;
fsum=max(fx,[],1);
fx(:,fsum<thresh)=[];
yyi1=yy; yyi1(fsum<thresh)=[];
lf=size(fx,2);
for k=1:lf
    cck=corrcoef(fx(:,k),fp);
    X(yyi1(k))=cck(2);
end
%
gx=R2;
gsum=max(gx,[],1);
gx(:,gsum<thresh)=[];
yyi1=yy; yyi1(gsum<thresh)=[];
lg=size(gx,2);
for k=1:lg
    cck=corrcoef(gx(:,k),fp);
    Y(yyi1(k))=cck(2);
end
%
hx=f;
hsum=max(hx,[],1);
hx(:,hsum<thresh)=[];
yyi1=yy; yyi1(hsum<thresh)=[];
lh=size(hx,2);
for k=1:lh
    cck=corrcoef(hx(:,k),fp);
    Z(yyi1(k))=cck(2);
end
ymin=min([X,Y,Z,0]);
%}
%%
fs=12; lw=1;%30;%Font size
%
figure
%{
    T=1:tauend;
    plot(T,[R1;R2],'-','linewidth',lw)%,'color',[.5,0,0])%'o-'
    xlabel('Time (years)','FontSize',fs)
    ylabel('R (N1N1)','FontSize',fs)
    %axis ([0,tauend,0,1]);
    set(gca,'FontSize',fs);
    grid on
    grid minor
%}
%
    cmod=.3;
    T=tmin:tmax;%1:tauend;
    subplot(2,1,1)
    hold on
    plot(T,R1','-','linewidth',lw)%,'color',[.5,0,0])%'o-'
    hold off
    %xlabel('Time (years)','FontSize',fs)
    ylabel('R (H1N1)','FontSize',fs)
    axis ([tmin-1,tmax,0,maxAttack]);
    %axis tight
    set(gca,'FontSize',fs);
    grid on
    grid minor
    box on
    subplot(2,1,2)
    hold on
    plot(T,R2','-','linewidth',lw)%,'color',[0,0,.5])%'o-'
    hold off
    xlabel('Time (years)','FontSize',fs)
    ylabel('R (H3N2)','FontSize',fs)
    axis ([tmin-1,tmax,0,maxAttack]);
    %axis tight
    set(gca,'FontSize',fs);
    grid on
    grid minor
    box on
%Comment out if don't want ccpand (or total):
%{
    figure
    subplot(2,1,1)
    hold on
    s=plot(T,f','-','linewidth',lw);%,'color',[0,0,.5])%'o-'
    plot(T,attack,'-','color',cmod*[1,1,1],'linewidth',2)
    hold off
    %alpha(s,.5)
    %xlabel('Time (years)','FontSize',fs)
    ylabel('R (total)','FontSize',fs)
    axis ([tmin-1,tmax,0,maxAttack]);
    %axis tight
    set(gca,'FontSize',fs);
    grid on
    grid minor
    box on
    subplot(2,1,2)
    hold on
    plot([0,tauend],[0,0],'k-','linewidth',1)
    h1=plot(T,X,'-','linewidth',2);
    h2=plot(T,Y,'-','linewidth',2);
    h3=plot(T,Z,'-','color',cmod*[1,1,1],'linewidth',2);
    hold off
    %alpha(s,.5)
    xlabel('Time (years)','FontSize',fs)
    ylabel('cc_{pand}','FontSize',fs)
    axis ([tmin-1,tmax,-1,1]);
    %axis tight
    set(gca,'FontSize',fs);
    legend([h1,h2,h3],'H1','H3','Total','location','SW')
    grid on
    grid minor
    box on
%}
%{
    T=tmin:tmax;%1:tauend;
    subplot(3,1,1)
    hold on
    plot(T,R1','-','linewidth',lw)%,'color',[.5,0,0])%'o-'
    hold off
    xlabel('Time (years)','FontSize',fs)
    ylabel('R (H1N1)','FontSize',fs)
    axis ([tmin-1,tmax,0,maxAttack]);
    %axis tight
    set(gca,'FontSize',fs);
    grid on
    grid minor
    box on
    subplot(3,1,2)
    hold on
    plot(T,R2','-','linewidth',lw)%,'color',[0,0,.5])%'o-'
    hold off
    xlabel('Time (years)','FontSize',fs)
    ylabel('R (H3N2)','FontSize',fs)
    axis ([tmin-1,tmax,0,maxAttack]);
    %axis tight
    set(gca,'FontSize',fs);
    grid on
    grid minor
    box on
    subplot(3,1,3)
    hold on
    s=plot(T,f','-','linewidth',lw);%,'color',[0,0,.5])%'o-'
    plot(T,attack,'-','color',.2*[1,1,1],'linewidth',2)
    hold off
    %alpha(s,.5)
    xlabel('Time (years)','FontSize',fs)
    ylabel('R (total)','FontSize',fs)
    axis ([tmin-1,tmax,y0,maxAttack]);
    axis tight
    set(gca,'FontSize',fs);
    grid on
    grid minor
    box on
    %}