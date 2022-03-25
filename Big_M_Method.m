clc
clear
format rat

% MIN Z=2x1+x2
% 3x1+x2=3;
% 4x1+3x2>=6
% x1+2x2<=3

%% PHASE 1: INPUT PARAMETER

A=[3 1; 4 3; 1 2];
B=[3; 6; 3];
Z=[-2 -1];

%% PHASE 2: STANDARD FORM

s=eye(size(A,1));
m=size(A,1);
n=size(A,2);
col=size(A,2);
M=1000;
I=[2 1 0];

greater_than=find(I==1);
for i=1:size(greater_than,2)
    mat=zeros(1,m);
    mat(greater_than(i))=-1;
    mat=mat';
    A=[A mat];
end
artificial_var=find(I==2 | I==1);
artificial_var_in_table(1:size(artificial_var,2))=(artificial_var+size(A,2));
A=[A s B];
Cj=zeros(1,size(A,2));
Cj(1:n)=Z;

for i=1:size(artificial_var_in_table,2)
    Cj(artificial_var_in_table(i))=-M;
end

%% PHASE 3: FIRST TABLE

bv=(n+size(greater_than,2)+1):size(A,2)-1;
zjcj=Cj(bv)*A-Cj;
fprintf("INITIAL TABLE:\n");
ZjC=[zjcj; A];
simpTable=array2table(ZjC);
simpTable.Properties.VariableNames(1:size(ZjC,2))={'x1','x2','s1','a1','a2','s2','Sol'};
disp(simpTable);

%% PHASE 4: BIG-M METHOD- SIMPLEX TABLES

table=1;
RUN=true;
zc=zjcj(1:size(A,2)-1);

while RUN
    if any(zc<0)
        [minvalzjcj, minindzjcj]=min(zc);
        pivot_col_ind=minindzjcj;
        pivot_col=A(:,pivot_col_ind);
        if all(pivot_col<=0)
            print('LPP is unbounded');
        else
            for i=1:size(pivot_col,1)
                if pivot_col(i)>0
                    ratio(i)=B(i)./pivot_col(i);
                else
                    ratio(i)=inf;
                end
            end
            [min_ratio, ratio_ind]=min(ratio);
            pivot_row_ind=ratio_ind;
            bv(ratio_ind)=pivot_col_ind;
            pivot_value=A(pivot_row_ind,pivot_col_ind);
            A(pivot_row_ind,:)=A(pivot_row_ind,:)./pivot_value;
            for i=1:size(A,1)
                if i~=pivot_row_ind
                    A(i,:)=A(i,:)-(pivot_col(i)*A(pivot_row_ind,:));
                end
            end
            zjcj=zjcj-(zjcj(pivot_col_ind)*A(pivot_row_ind,:));
        end
        ZjC=[zjcj; A];
        zc=zjcj(1:size(A,2)-1);
        B(1:m)=(A(:,size(A,2)))';
        simpTable=array2table(ZjC);
        simpTable.Properties.VariableNames(1:size(ZjC,2))={'x1','x2','s1','a1','a2','s2','Sol'};
        fprintf("TABLE %d:\n",table);
        table=table+1;
        disp(simpTable);
    else
        RUN=false;
    end
end

%% PHASE 5: OPTIMAL SOLUTION

fprintf("FINAL TABLE:\n");
disp(simpTable);

optimal=true;

for i=1:size(bv,2)
    index=find(artificial_var_in_table==bv(i));
    if index~=0
        optimal=false;
    end
end

if optimal==false
    disp('INFEASIBLE SOLUTION');
else
    fprintf("OPTIMAL SOLUTION:\n");
    for i=1:col
        index=find(bv==i);
        if(index>0)
            fprintf("x%d = %.3f\n",i,A(index,size(A,2)));
        else
            fprintf("x%d = 0\n",i);
        end
    end
fprintf("Min Z = %f",ZjC(1,size(ZjC,2)));
end
