clc
clear
format rat


% Max Z=3x1+2x2
% 2x1+x2<=18
% 2x1+3x2<=42
% 3x1+x2<=24

%% PHASE 1: INPUT PARAMETER

Z=[3 2];
A=[2 1; 2 3; 3 1];
B=[18; 42; 24];

%% PHASE 2: COMPLETE MATRIX AND COST MATRIX

s=eye(size(A,1));
m=size(A,1);
n=size(A,2);
col=size(A,2);
A=[A s B];
Cj=zeros(1,size(A,2));
% Cost Matrix
Cj(1:n)=Z;

%% PHASE 3: FIRST TABLE

bv=n+1:size(A,2)-1;
zjcj=Cj(bv)*A-Cj;
fprintf("INITIAL TABLE:\n");
ZjC=[zjcj; A];
simpTable=array2table(ZjC);
simpTable.Properties.VariableNames(1:size(ZjC,2))={'x1','x2','s1','s2','s3','Sol'};
disp(simpTable);

%% PHASE 4: SIMPLEX TABLES

table=1;
optimal=true;
RUN=true;
zc=zjcj(1:size(A,2)-1);

while RUN
    if any(zc<0)
        zc=zjcj(1:size(A,2)-1);
        [minvalzjcj, minindzjcj]=min(zc);
        pivot_col_ind=minindzjcj;
        pivot_col=A(:,pivot_col_ind);
        if all(pivot_col<=0)
            print('LPP is unbounded');
            optimal=false;
            break;
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
        B(1:m)=(A(:,size(A,2)))';
        simpTable=array2table(ZjC);
        simpTable.Properties.VariableNames(1:size(ZjC,2))={'x1','x2','s1','s2','s3','Sol'};
        fprintf("TABLE %d:\n",table);
        table=table+1;
        disp(simpTable);
    else
        RUN=false;
    end
end

%% PHASE 5: OPTIMAL SOLUTION

if optimal==true
    fprintf("FINAL TABLE:\n");
    disp(simpTable);
    fprintf("OPTIMAL SOLUTION:\n");
    for i=1:col
        index=find(bv==i);
        if(index==0)
            fprintf("x%d = 0\n",i);
        else
            fprintf("x%d = %.3f\n",i,A(index,size(A,2)));
        end
    end
    fprintf("Max Z = %f",ZjC(1,size(ZjC,2)));
end



