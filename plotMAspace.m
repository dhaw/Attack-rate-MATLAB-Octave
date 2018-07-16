function f=plotMAspace(A,Acum,NN,NNbar)
%ar=1;

tauend=size(Acum,2);
nbar=length(NNbar); n=length(NN);
R1=Acum(1:nbar,:);
R2=Acum(nbar+1:end,:);
%repNNbar=repmat(NNbar,1,tauend);
repNN=repmat(NN,1,tauend);
%R1=R1.*repNNbar; R2=R2.*repNNbar;
R1=R1(1:n,:)+R1(n+1:2*n,:)+R1(2*n+1:3*n,:)+R1(3*n+1:end,:);
R1=R1./repNN; R1(repNN==0)=0;
R2=R2(1:n,:)+R2(n+1:2*n,:)+R2(2*n+1:3*n,:)+R2(3*n+1:end,:);
R2=R2./repNN; R2(repNN==0)=0;

fs=15; lw=2;%30;%Font size
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
    T=1:tauend;
    subplot(2,1,1)
    plot(T,R1','-','linewidth',lw)%,'color',[.5,0,0])%'o-'
    xlabel('Time (years)','FontSize',fs)
    ylabel('R (N1N1)','FontSize',fs)
    %axis ([0,tauend,0,1]);
    axis tight
    set(gca,'FontSize',fs);
    grid on
    grid minor
    subplot(2,1,2)
    plot(T,R2','-','linewidth',lw)%,'color',[0,0,.5])%'o-'
    xlabel('Time (years)','FontSize',fs)
    ylabel('R (H3N2)','FontSize',fs)
    %axis ([0,tauend,0,1]);
    axis tight
    set(gca,'FontSize',fs);
    grid on
    grid minor
%}