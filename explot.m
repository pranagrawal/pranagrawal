% cl=length(msh.TRIANGLES);
% for i=1:cl
% n1=msh.TRIANGLES(i,1);
% n2=msh.TRIANGLES(i,2);
% n3=msh.TRIANGLES(i,3);
% n1_pos=msh.POS(n1,:);
% n2_pos=msh.POS(n2,:);
% n3_pos=msh.POS(n3,:);
points;
fopen("points.m",'r');
x=msh.POS(:,1);
y=msh.POS(:,2);
% P(1,:)=n1_pos;
% P(2,:)=n2_pos;
% P(3,:)=n3_pos;
% P(:,3)=[];
T = msh.TRIANGLES(:,1:3);
triplot(T,x,y)
% end