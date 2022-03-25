clear all
format rat
clc

%% INPUT PARAMETER
A=[1 2; 1 1; 1 -2];
B=[10; 6; 1];

%% PLOTTING LINES
y1=0:1:max(B);
x12=(B(1)-A(1,1)*y1)/A(1,2);
x22=(B(2)-A(2,1)*y1)/A(2,2);
x32=(B(3)-A(3,1)*y1)/A(3,2);
x12=max(0,x12);
x22=max(0,x22);
x32=max(0,x32);
plot(y1,x12,'r',y1,x22,'g',y1,x32,'b');
grid on;
title('x1 vs x2');
xlabel('value of x1');
ylabel('value of x2');

%% FINDING CORNER POINTS
cx1=find(y1==0);
c1=find(x12==0);
Line1=[y1(:,[c1 cx1]) ; x12(:,[c1 cx1])]';     
c2=find(x22==0);
Line2=[y1(:,[c2 cx1]) ; x22(:,[c2 cx1])]';
c3=find(x32==0);
Line3=[y1(:,[c3 cx1]) ; x32(:,[c3 cx1])]';
cornerpts=unique([Line1;Line2;Line3],'rows');

%% INTERSECTION POINTS
pt=[];
for i=1:size(A,1)      
    A1=A(i,:);
    B1=B(i,:);
    for j=i+1:size(A,1)
        A2=A(j,:);
        B2=B(j,:);
        P=[A1;A2];
        Q=[B1;B2];
        X=P\Q;
        pt=[pt X];
    end
end
ptt=pt';

%% ALL CORNER AND INTERSECTION POINTS
allpts=[ptt;cornerpts];
points=unique(allpts,'rows');

%% CONSTRAINTS OR FEASIBLE REGION
x1=points(:,1);
x2=points(:,2);
cons1=x1+2*x2-10;
h1=find(cons1>0);
points(h1,:)=[];

x1=points(:,1);
x2=points(:,2);
cons2=x1+x2-6;
h2=find(cons2>0);
points(h2,:)=[];

x1=points(:,1);
x2=points(:,2);
cons3=x1-2*x2-1;
h3=find(cons3>0);
points(h3,:)=[];

%% OBJECTIVE FUNCTION
C=[2;1];
for i=1:size(points,1)
    fX(i,:)=points(i,:)*C;
end
obj=[points fX];
maximum=max(fX);
indexOfMax=find(fX==maximum);
optimalSol=points(indexOfMax,:)
zOptimal=maximum