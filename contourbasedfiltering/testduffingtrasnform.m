% close all
figure
plot3(X(:,1),X(:,2),probs,'ro')

Z= zeros(size(X));
for i=1:size(X,1)
    DY=diag([1/6,1/20]);
    y=DY*X(i,:)';
   Z(i,1) = norm(y);
   Z(i,2) = atan(y(2)/y(1));
end

figure
plot3(Z(:,1),Z(:,2),probs,'b+')
xlabel('r')
ylabel('th')

figure
 Y=sort(X,1);
plot(1:size(Y,1),Y(:,1),'r',1:size(Y,1),Y(:,2),'b')


figure
Z = linkage(X,'ward');
c = cluster(Z,'Maxclust',2);
scatter(X(:,1),X(:,2),13,c)