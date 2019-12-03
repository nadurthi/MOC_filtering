function h=duff_meas_model_th(x)
[r,c] = size(x);
if r==1 || c==1
    x=x(:)';
end
a=[-20,-20];
hn=1;

h=zeros(size(x,1),hn);
for i=1:size(x,1)
%     h(i,:) = [x(i,1),x(i,2)];
%     h(i,1) = sqrt((x(i,1)-a(1))^2+(x(i,2)-a(2))^2);
%     h(i,2) = atan2(x(i,2)-a(2),x(i,1)-a(1));

        h(i,1) = atan2(x(i,2)-a(2),x(i,1)-a(1));
end


[r,c] = size(h);
if r==1 || c==1
    h=h(:);
end