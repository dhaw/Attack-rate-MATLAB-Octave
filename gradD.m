function f=gradD(D,z)
[a,b]=size(D);
z=reshape(z,a,b); z=z(2:end-1,2:end-1); z=reshape(z,(a-2)*(b-2),1);

D=log10(D);

A=D(2:end-1,2:end-1);
a1=D(1:end-2,1:end-2);
a2=D(1:end-2,2:end-1);
a3=D(1:end-2,3:end);
a4=D(2:end-1,3:end);
a5=D(3:end,3:end);
a6=D(3:end,2:end-1);
a7=D(3:end,1:end-2);
a8=D(2:end-1,1:end-2);
N=a1+a2+a3+a4+a5+a6+a7+a8; N=N/8;
%G=A./N; G(N==0)=NaN; G(N==-inf)=NaN; G(A==0)=NaN;
G=A-N; G(N==0)=NaN; G(N==-inf)=NaN; G(A==0)=NaN;

%G(G==0)=NaN;

maxG=max(max(G));
%{
figure
imagesc(G)
xlabel('longitude')
ylabel('latitude')
hcb=colorbar;
caxis([0 maxG])
set(hcb,'YTick',[0,maxG])
%}

fs=18; ms=5;%12; 20
lwx=1; lw=1.5; %lwx=1.5
col1=[.165,.31,.431];%[1,0,0];%[0.0512,0.4600,0.8633];%[.165,.31,.431];
col2=[.447,.553,.647];
col3=[0,0,0];

figure
maxG2=maxG;%10000;
minG2=-maxG2;
GG=reshape(G,(a-2)*(b-2),1);
Gb=GG; Gb(Gb>maxG2)=NaN; Gb(Gb<minG2)=NaN;
%G(G==0)=NaN;
%
Gc=Gb; Gc(isnan(Gb)==1)=[]; zc=z; zc(isnan(Gb)==1)=[];
Gc(zc==0)=[]; zc(zc==0)=[];
cc=corrcoef(Gc,zc); cc=cc(2);
%
X=[ones(size(Gc)),Gc];
r=regress(zc,X);
maxGb=max(Gb); minGb=min(Gb); Gvec=[minGb,maxGb]; zvec=[r(1)+r(2)*minGb,r(1)+r(2)*maxGb];
%
grad=r(2);%(zvec(end)-zvec(1))/(Gvec(end)-Gvec(1));
f=[cc,grad];%[GG,z];%cc;%length(G(G==0));

A10=10.^A;
%A(isnan(A)==1)=0;
[aa,bb]=size(A10); pop=reshape(A10,aa*bb,1);
meanz=sum(pop.*z)/sum(pop);
%minGG=min(GG); maxGG=max(GG);

maxY=max(z); minY=min(z);
hold on
plot(GG,z,'o','color',col2,'markersize',ms,'LineWidth',lwx);
plot([minGb,maxGb],[meanz,meanz],'k--','linewidth',lw);
%plot(Gvec,zvec,'k-','LineWidth',lw);
xlabel('grad (log N)')
ylabel('Attack rate')
set(gca,'FontSize',fs);
minY=max(minY-.1,0); maxY=min(maxY+.1,1);
axis ([minGb,maxGb,0,maxY])%maxG
grid on
grid minor
box on
hold off