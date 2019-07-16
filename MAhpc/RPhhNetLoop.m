function f=RPhhNetLoop
alpha=0:6;
la=length(alpha);
%for i=1:la
    %load.data
    
H0=readtable('hhNetalpha0B.csv');
H1=readtable('hhNetalpha1B.csv');
H2=readtable('hhNetalpha2B.csv');
H3=readtable('hhNetalpha3B.csv');
H4=readtable('hhNetalpha4B.csv');
H5=readtable('hhNetalpha5B.csv');
H6=readtable('hhNetalpha6B.csv');
    
H0=RPcombine(table2array(H0(:,2:3)));
H1=RPcombine(table2array(H1(:,2:3)));
H2=RPcombine(table2array(H2(:,2:3)));
H3=RPcombine(table2array(H3(:,2:3)));
H4=RPcombine(table2array(H4(:,2:3)));
H5=RPcombine(table2array(H5(:,2:3)));
H6=RPcombine(table2array(H6(:,2:3)));

psave('H');