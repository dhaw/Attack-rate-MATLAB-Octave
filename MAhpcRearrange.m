function f=MAhpcRearrange(out)
lout=length(out);
[a,b,c]=size(out{1});%Assumes all cells same size
Y=cell(a,b);
for k=1:lout
    y=out{k};
    for i=1:a
        for j=1:b
            yijk=y(i,j,:);
            Y{i,j}(k,:)=yijk;
        end
    end
end
f=Y;
            