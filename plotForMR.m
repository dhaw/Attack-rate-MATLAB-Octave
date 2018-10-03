function f=plotForMR
x=(-10:.1:10);
y=[2.^x;3.^x;1.2.^x;2.^(-x)];
%y=[2*x;3-x];

%x=25*rand(1000,1)+1;
%y=2.^(x)+1000*x.*(rand(1000,1)-.5);
%y(y<0)=NaN;

maxY=max(max(y));
%maxY=max([y1,y2]);
%minY=min([0,min([y1,y2])]);

figure
fs=20;%18; ms=5;%12; 20
lwx=1; lw=2; %lwx=1.5
%hold on
%semilogy(x,y,'ok','LineWidth',lw);%,'markersize',ms
%
hold on
plot([0,0],[-maxY,maxY],'-k','LineWidth',lwx)
plot([-10,10],[0,0],'-k','LineWidth',lwx)
h=plot(x,y,'-','LineWidth',lw);%,'markersize',ms
hold off
%}
hold on
xlabel('x','FontSize',fs);
ylabel('y','FontSize',fs,'rot',0);
set(gca,'FontSize',fs)%,'ytick',[10,100,1000,10000,100000,1000000,10000000,100000000]);
axis([-5,5,0,60])%([x(1),x(end),0,maxY])%minY,maxY])
%axis([-10,10,-20,20])
legend(h,'y=2^x','y=3^x','y=1.2^x','y=0.5^x','location','N')
grid on
grid minor