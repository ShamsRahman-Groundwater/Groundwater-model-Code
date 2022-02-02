% Author: Mostaquimur Rahman, University of Bristol, UK
% (ar15645@bristol.ac.uk)
% This is the groundwater model code

function [Hnew,runoff] = groundwater(model_topography,topo_resolution, ...
                         R,S,T,Hold,dt)
    
    dx         = topo_resolution; 
    dy         = topo_resolution; 
    topo_size  = size(model_topography);
    x_dim = topo_size(2)*dx; % total length in x direction (m)
    y_dim = topo_size(1)*dy; % total length in y direction (m)
    x = 0:dx:(x_dim-1); % x grid
    y = 0:dy:(y_dim-1); % y grid

    % indexing for sparse matrix
    LHS_sparse_row = zeros(length(y)*length(x)*5,1);
    LHS_sparse_col = zeros(length(y)*length(x)*5,1);
    LHS_sparse_val = zeros(length(y)*length(x)*5,1);
    RHS_sparse_row = zeros(length(y)*length(x),1);
    RHS_sparse_col = zeros(length(y)*length(x),1);
    RHS_sparse_val = zeros(length(y)*length(x),1);

    c = 1;
    cc = 1;
    for i = 1:length(x)
        for j = 1:length(y)
            % deal with the boundary nodes
            if (i==1 & j==1)
                row_idx = j; col_idx = i;
                left_row_idx = j; left_col_idx = i+1;
                right_row_idx = j; right_col_idx = i+1;
                up_row_idx = j+1; up_col_idx = i;
                down_row_idx = j+1; down_col_idx = i;
            elseif (i==1 & j==length(y))
                row_idx = j; col_idx = i;
                left_row_idx = j; left_col_idx = i+1;
                right_row_idx = j; right_col_idx = i+1;
                up_row_idx = j-1; up_col_idx = i;
                down_row_idx = j-1; down_col_idx = i;
            elseif (i==length(x) & j==1)
                row_idx = j; col_idx = i;
                left_row_idx = j; left_col_idx = i-1;
                right_row_idx = j; right_col_idx = i-1;
                up_row_idx = j+1; up_col_idx = i;
                down_row_idx = j+1; down_col_idx = i;
            elseif (i==length(x) & j==length(y))
                row_idx = j; col_idx = i;
                left_row_idx = j; left_col_idx = i-1;
                right_row_idx = j; right_col_idx = i-1;
                up_row_idx = j-1; up_col_idx = i;
                down_row_idx = j-1; down_col_idx = i;
            elseif (i==1 & 1<j<length(y))
                row_idx = j; col_idx = i;
                left_row_idx = j; left_col_idx = i+1;
                right_row_idx = j; right_col_idx = i+1;
                up_row_idx = j-1; up_col_idx = i;
                down_row_idx = j+1; down_col_idx = i;
            elseif (i==length(x) & 1<j<length(y))
                row_idx = j; col_idx = i;
                left_row_idx = j; left_col_idx = i-1;
                right_row_idx = j; right_col_idx = i-1;
                up_row_idx = j-1; up_col_idx = i;
                down_row_idx = j+1; down_col_idx = i;
            elseif (1<i<length(x) & j==1)
                row_idx = j; col_idx = i;
                left_row_idx = j; left_col_idx = i-1;
                right_row_idx = j; right_col_idx = i+1;
                up_row_idx = j+1; up_col_idx = i;
                down_row_idx = j+1; down_col_idx = i;
            elseif (1<i<length(x) & j==length(y))
                row_idx = j; col_idx = i;
                left_row_idx = j; left_col_idx = i-1;
                right_row_idx = j; right_col_idx = i+1;
                up_row_idx = j-1; up_col_idx = i;
                down_row_idx = j-1; down_col_idx = i;
            else
                row_idx = j; col_idx = i;
                left_row_idx = j; left_col_idx = i-1;
                right_row_idx = j; right_col_idx = i+1;
                up_row_idx = j-1; up_col_idx = i;
                down_row_idx = j+1; down_col_idx = i;
            end
            
            % create sparse matrices with the coefficient
            % LHS
            % coefficient of j,i
            LHS_sparse_row(cc) = c;
            LHS_sparse_col(cc) = (col_idx-1)*length(y)+row_idx;
            LHS_sparse_val(cc) = harmmean([T(left_row_idx,left_col_idx) T(row_idx,col_idx)]) ...
                                 + harmmean([T(right_row_idx,right_col_idx) T(row_idx,col_idx)]) ...
                                 + harmmean([T(up_row_idx,up_col_idx) T(row_idx,col_idx)]) ...
                                 + harmmean([T(down_row_idx,down_col_idx) T(row_idx,col_idx)]) ...
                                 + (dx^2*S(row_idx,col_idx)/dt);
            cc = cc + 1;
            % coefficient of j-1,i (up)
            LHS_sparse_row(cc) = c;
            LHS_sparse_col(cc) = (col_idx-1)*length(y)+up_row_idx;
            LHS_sparse_val(cc) = -1.*harmmean([T(up_row_idx,up_col_idx) T(row_idx,col_idx)]);
            cc = cc + 1;
            % coefficient of j+1,i (down)
            LHS_sparse_row(cc) = c;
            LHS_sparse_col(cc) = (col_idx-1)*length(y)+down_row_idx;
            LHS_sparse_val(cc) = -1.*harmmean([T(down_row_idx,down_col_idx) T(row_idx,col_idx)]);
            cc = cc + 1;
            % coefficient of j,i-1 (left)
            LHS_sparse_row(cc) = c;
            LHS_sparse_col(cc) = (left_col_idx-1)*length(y)+row_idx;
            LHS_sparse_val(cc) = -1.*harmmean([T(left_row_idx,left_col_idx) T(row_idx,col_idx)]);
            cc = cc + 1;
            % coefficient of j,i+1 (right)
            LHS_sparse_row(cc) = c;
            LHS_sparse_col(cc) = (right_col_idx-1)*length(y)+row_idx;
            LHS_sparse_val(cc) = -1.*harmmean([T(right_row_idx,right_col_idx) T(row_idx,col_idx)]);
            cc = cc + 1;        

            % RHS
            RHS_sparse_row(c) = c;
            RHS_sparse_col(c) = 1;
            RHS_sparse_val(c) = ((dx^2*S(row_idx,col_idx))/dt)*Hold(row_idx,col_idx) ...
                                + R(row_idx,col_idx)*dx^2;    
            c = c + 1;
        end
    end
    LHS = sparse(LHS_sparse_row,LHS_sparse_col,LHS_sparse_val);
    RHS = sparse(RHS_sparse_row,RHS_sparse_col,RHS_sparse_val);
%     sol = LHS\RHS;
    sol = pcg(LHS,RHS,10^-9,100);
    
    Hnew = zeros(topo_size);
    ck = 1;
    for i = 1:length(x)
        for j = 1:length(y)
            Hnew(j,i) = sol(ck);
            ck = ck + 1;
        end
    end
   
    runoff = max(0,(Hnew - model_topography).*S); 
    % adjust head
    Hnew = Hnew - (runoff./S);
end
