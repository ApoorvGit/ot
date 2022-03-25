clc
clear
format rat

% MIN Z=7.5x1-3x2
% 3x1-x2-x3>=3
% x1-x2+x3>=2

%% PHASE 1: INPUT PARAMETER

A=[3 -1 -1; 1 -1 1];
B=[3; 2];
Z=[-7.5 3 0];

%% PHASE 2: STANDARD FORM

s=eye(size(A,1));
m=size(A,1);
n=size(A,2);
col=size(A,2);
I=[1 1];
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

for i=1:size(artificial_var_in_table,2)
    Cj(artificial_var_in_table(i))=-1;
end

%% PHASE I OF TWO-PHASE METHOD

% FIRST TABLE

bv=(n+size(greater_than,2)+1):size(A,2)-1;
zjcj=Cj(bv)*A-Cj;
fprintf("----------PHASE I----------\n\n");
fprintf("INITIAL TABLE:\n");
ZjC=[zjcj; A];
simpTable=array2table(ZjC);
simpTable.Properties.VariableNames(1:size(ZjC,2))={'x1','x2','x3','s1','s2','a1','a2','Sol'};
disp(simpTable);

% PHASE I- SIMPLEX TABLES

table=1;
RUN=true;
zc=zjcj(1:size(A,2)-1);
optimal=true;

while RUN
    if any(zc<0)
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
        zc=zjcj(1:size(A,2)-1);
        B(1:m)=(A(:,size(A,2)))';
        simpTable=array2table(ZjC);
        simpTable.Properties.VariableNames(1:size(ZjC,2))={'x1','x2','x3','s1','s2','a1','a2','Sol'};
        fprintf("TABLE %d:\n",table);
        table=table+1;
        disp(simpTable);
    else
        RUN=false;
    end
end


%% PHASE II OF TWO-PHASE METHOD

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
    fprintf("\n\n----------PHASE II----------\n");
    artificial_variables_number=size(artificial_var_in_table,2);
    A=A(1:end,1:end-artificial_variables_number-1);
    A=[A B];
    Cj=zeros(1,size(A,2));
    Cj(1:n)=Z;
    
    % FIRST TABLE
    zjcj=Cj(bv)*A-Cj;
    fprintf("INITIAL TABLE:\n");
    ZjC=[zjcj; A];
    simpTable=array2table(ZjC);
    simpTable.Properties.VariableNames(1:size(ZjC,2))={'x1','x2','x3','s1','s2','Sol'};
    disp(simpTable);
    
    % PHASE II SIMPLEX TABLES
    table=1;
    RUN=true;
    zc=zjcj(1:size(A,2)-1);
    optimal=true;
    
    while RUN
        if any(zc<0)
            [minvalzjcj, minindzjcj]=min(zc);
            pivot_col_ind=minindzjcj;
            pivot_col=A(:,pivot_col_ind);
            if all(pivot_col<=0)
                disp('LPP is unbounded');
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
            zc=zjcj(1:size(A,2)-1);
            B(1:m)=(A(:,size(A,2)))';
            simpTable=array2table(ZjC);
            simpTable.Properties.VariableNames(1:size(ZjC,2))={'x1','x2','x3','s1','s2','Sol'};
            fprintf("TABLE %d:\n",table);
            table=table+1;
            disp(simpTable);
        else
            RUN=false;
        end
    end
    
    % OPTIMAL SOLUTION
    
    if optimal==true
        fprintf("FINAL TABLE:\n");
        disp(simpTable);
        fprintf("OPTIMAL SOLUTION:\n");
        for i=1:col
            index=find(bv==i);
            if(index>0)
                fprintf("x%d = %.3f\n",i,A(index,size(A,2)));
            else
                fprintf("x%d = %d\n",i,0);
            end
        end
        fprintf("Min Z = %f",ZjC(1,size(ZjC,2)));
    end
end
