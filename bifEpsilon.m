<<<<<<< HEAD
function [F,G]=bifEpsilon%(F,G)
tauend=200;%As MA1d
solvetype=2;
burn=0;
%eps=[(0:.005:.05),(.051:.001:.3),(.305:.005:.6),(.601:.001:.9),(.905:.005:1)];2
eps=(0:.005:.5);
leps=length(eps);
%Comment out if input Fs and Gs
%
F=zeros(leps,tauend-burn); %F2=F1; G1=F1; G2=F1;
G=F;
for i=1:leps
    epsi=eps(i);
    %[f1,g1,other]=MA1D(solvetype,0,epsi);
    %F(i,:)=f1(burn+1:end); G(i,:)=g1(burn+1:end);
    [f1,g1]=MA1D2sub(epsi);
    F(i,:)=g1(1,burn+1:end); G(i,:)=g1(2,burn+1:end);
end
%}
fs=15; lw=3; ms=1; col1=[0,0,0]; col2=[.5,.5,.5];
figure
hold on
for i=1:leps
    ui1=unique(F(i,:));
    epsi1=eps(i)*ones(length(ui1),1);
    plot(epsi1,ui1,'o','color',col1,'markersize',ms);
    %ui2=unique(G(i,:));
    %epsi2=eps(i)*ones(length(ui2),1);
    %plot(epsi2,ui2,'o','color',col2,'markersize',ms);
end
maxF=max(max([F;G]));
axis([eps(1),eps(end),0,maxF])
xlabel('\epsilon','fontsize',fs)
ylabel('z','fontsize',fs,'rot',0)
set(gca,'fontsize',fs)
hold off

figure
hold on
for i=1:leps
    ui2=unique(G(i,:));
    epsi2=eps(i)*ones(length(ui2),1);
    plot(epsi2,ui2,'o','color',col2,'markersize',ms);
end
maxF=max(max(F));
axis([eps(1),eps(end),0,maxF])
xlabel('\epsilon','fontsize',fs)
ylabel('z','fontsize',fs,'rot',0)
set(gca,'fontsize',fs)
hold off
=======
function [F,G]=bifEpsilon%(F,G)
tauend=1000;%As MA1d
solvetype=2;
burn=600;
%eps=[(0:.005:.05),(.051:.001:.3),(.305:.005:.6),(.601:.001:.9),(.905:.005:1)];2
eps=(0:.005:.5);
leps=length(eps);
%Comment out if input Fs and Gs
%
F=zeros(leps,tauend-burn); %F2=F1; G1=F1; G2=F1;
G=F;
for i=1:leps
    epsi=eps(i);
    %[f1,g1,other]=MA1D(solvetype,0,epsi);
    %F(i,:)=f1(burn+1:end); G(i,:)=g1(burn+1:end);
    [f1,g1]=MA1D2sub(epsi);
    F(i,:)=g1(1,burn+1:end); G(i,:)=g1(2,burn+1:end);
end
%}
fs=15; lw=3; ms=1; col1=[0,0,0]; col2=[.5,.5,.5];
figure
hold on
for i=1:leps
    ui1=unique(F(i,:));
    epsi1=eps(i)*ones(length(ui1),1);
    plot(epsi1,ui1,'o','color',col1,'markersize',ms);
    %ui2=unique(G(i,:));
    %epsi2=eps(i)*ones(length(ui2),1);
    %plot(epsi2,ui2,'o','color',col2,'markersize',ms);
end
maxF=max(max([F;G]));
axis([eps(1),eps(end),0,maxF])
xlabel('\epsilon','fontsize',fs)
ylabel('z','fontsize',fs,'rot',0)
set(gca,'fontsize',fs)
hold off

figure
hold on
for i=1:leps
    ui2=unique(G(i,:));
    epsi2=eps(i)*ones(length(ui2),1);
    plot(epsi2,ui2,'o','color',col2,'markersize',ms);
end
maxF=max(max(F));
axis([eps(1),eps(end),0,maxF])
xlabel('\epsilon','fontsize',fs)
ylabel('z','fontsize',fs,'rot',0)
set(gca,'fontsize',fs)
hold off
>>>>>>> acfad3e4da850e525286e08138beef9f06a8c2cc
