function f=stochDom(A,f,NN,isdual)%Where does stochasticity dominate?
[a,b]=size(A);
sdf=var(f); sdf=sqrt(sdf);
B=A-repmat(f,1,b);
sdB=var(B,0,2); sdB=sqrt(sdB);
C=sdB/sdf; C=[NN,C]; C=sortrows(C,1);%In case
maxN=max(NN); maxY=max(C(:,2));
xmin=1;

if isdual==0%SM
    x=[.6,.8]; y=[.5,.3];
    col1=[.071,.208,.322];
else%IM
    x=[.75,.6]; y=[.55,.4];
    col1=[.447,.553,.647];
end
fs=24; lw=3; ms=10; col2=[.7,0,0];
figure
semilogx((C(:,1)),C(:,2),'.','color',col1,'linewidth',lw,'markersize',ms)
hold on
semilogx([xmin,maxN],[1,1],'--','color',col2,'linewidth',lw)
axis([xmin,maxN,0,maxY])
xlabel('Pop. density','FontSize',fs); ylabel('Ratio R')
grid on
grid minor
annotation('textarrow',x,y,'String','R=1','linewidth',lw,'fontsize',fs,'color',col2)%Threshold 
set(gca,'fontsize',fs)
hold off