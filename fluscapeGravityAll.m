function [Z1,Z2,Z3,Z4]=fluscapeGravityAll(Df,Ds,isdual,sample)
R0=1.8;
%isdual=1;%0=SM, 1=DM, 2=IM
solvetype=3;
runs=2;
[lscan,r]=fluscapeNNr(Df,Ds);
[K1,K2,K3,K4,NN]=makeK(lscan,r);
Z1=zeros(6*360,runs);
Z2=Z1; Z3=Z1; Z4=Z1;
%
for i=1:runs
[gamma,n,betaS,betaI,betaD,K]=prepFluGravityNoAge(R0,K1,NN);
%sample=randsample(n,5);
[f1,g1,z1]=finalSize1yearNoAge(gamma,n,NN,K,betaS,betaI,betaD,isdual,solvetype,sample,1);
%figure; imagesc(log(K));
Z1(1:length(z1),i)=z1;
%
[gamma,n,betaS,betaI,betaD,K]=prepFluGravityNoAge(R0,K2,NN);
[f2,g2,z2]=finalSize1yearNoAge(gamma,n,NN,K,betaS,betaI,betaD,isdual,solvetype,sample,1);
%figure; imagesc(log(K));
Z2(1:length(z2),i)=z2;
%
[gamma,n,betaS,betaI,betaD,K]=prepFluGravityNoAge(R0,K3,NN);
[f3,g3,z3]=finalSize1yearNoAge(gamma,n,NN,K,betaS,betaI,betaD,isdual,solvetype,sample,1);
%figure; imagesc(log(K));
Z2(1:length(z3),i)=z3;
%
[gamma,n,betaS,betaI,betaD,K]=prepFluGravityNoAge(R0,K4,NN);
[f4,g4,z4]=finalSize1yearNoAge(gamma,n,NN,K,betaS,betaI,betaD,isdual,solvetype,sample,1);
%figure; imagesc(log(K));
Z2(1:length(z4),i)=z4;
%
end
%f=[z1,z2,z3,z4];
save('stochDM.mat','f1','g1','z1','f2','g2','z2','f3','g3','z3','f4','g4','z4')
%
figure
fs=12; lw=2;
hold on
plot(f1,g1,'linewidth',lw);
plot(f2,g2,'linewidth',lw);
plot(f3,g3,'linewidth',lw);
plot(f4,g4,'linewidth',lw);
xlabel('Time (days)','FontSize',fs);
ylabel('Incidence','FontSize',fs);
set(gca,'FontSize',fs);
maxY=max([g1;g2;g3;g4]);
tstart=0;
tend=120;%min([f1(end),f2(end),f3(end),f4(end)]);
axis([tstart,tend,0,maxY])
grid on
grid minor
box on
legend('G','OG','R','OR','location','NE')
hold off