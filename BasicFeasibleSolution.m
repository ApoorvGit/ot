clear all
clc
format rat
%% PHASE 1: INPUT PARAMETER
A=[1 1 1 0;2 1 0 1];          
B=[450;600];
C=[3 4 0 0];

%% PHASE 2: SET OF ALL BASIC SOLUTIONS
m=size(A,1);   
n=size(A,2);
if(n>m)
    nCm=nchoosek(n,m)
    pair=nchoosek(1:n,m)
    sol=[];
    for i=1:nCm
        y=zeros(n,1)
        x=A(:,pair(i,:))\B
        if(x>=0 & x~=inf & x~=-inf)
            y(pair(i,:))=x
            sol=[sol,y]
        end
    end
else
    error('nCm does not exist')
end

%% BASIC FEASIBLE SOLUTION
z=C*sol
[Zmax,Zindex]=max(z)
bfs1=sol(:,Zindex)
bfs3=bfs1'
optimal_value=[bfs3 Zmax]
optimal_bfs=array2table(optimal_value);
optimal_bfs.Properties.VariableNames(1:size(optimal_bfs,2))={'x1','x2','s1','s2','Z'}