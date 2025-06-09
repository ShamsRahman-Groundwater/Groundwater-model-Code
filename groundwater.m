function [Hnew,runoff] = groundwater(model_topography,topo_resolution, ...
                         R,S,T,Hold,r_inac,c_inac,dt,sea_level)
    
    topo_size  = size(model_topography);
    x_dim = topo_size(2);
    y_dim = topo_size(1);
    
    dx = topo_resolution;
    dy = topo_resolution;

    % indexing for sparse matrix
    LHS_sparse_row = zeros(y_dim*x_dim*5,1);
    LHS_sparse_col = zeros(y_dim*x_dim*5,1);
    LHS_sparse_val = zeros(y_dim*x_dim*5,1);
    RHS_sparse_row = zeros(y_dim*x_dim,1);
    RHS_sparse_col = zeros(y_dim*x_dim,1);
    RHS_sparse_val = zeros(y_dim*x_dim,1);

    % LHS and RHS matrices
    c = 1;
    cc = 1;
    for i = 1:x_dim
        for j = 1:y_dim
            
            % set the cell indecies
            row_idx = j; col_idx = i;
            left_col_idx = i-1; right_col_idx = i+1;
            up_row_idx = j-1; down_row_idx = j+1;
            
            % handle boundary cell indecies
            if (i==1)
                left_col_idx = i;
            end
            if (i==x_dim)
                right_col_idx = x_dim;
            end
            if (j==1)
                up_row_idx = 1;
            end
            if (j==y_dim)
                down_row_idx = y_dim;
            end
             
            % calculate inter-cell transmissivity (harmonic mean)
            T_left = harmmean([T(row_idx,col_idx) T(row_idx,left_col_idx)]);
            T_right = harmmean([T(row_idx,col_idx) T(row_idx,right_col_idx)]);
            T_up = harmmean([T(row_idx,col_idx) T(up_row_idx,col_idx)]);
            T_down = harmmean([T(row_idx,col_idx) T(down_row_idx,col_idx)]);

            % set the coefficient matrix for constructing Ax = B
            % set the left hand side (A)
            % coefficient of j,i
            LHS_sparse_row(cc) = c;
            LHS_sparse_col(cc) = (col_idx-1)*y_dim+row_idx;
            pp = (T_left+T_right)/dx^2;
            qq = (T_up+T_down)/dy^2;
            rr = S(row_idx,col_idx)/dt;
            LHS_sparse_val(cc) = pp+qq+rr;
            cc = cc + 1;
            % coefficient of j-1,i (up)
            LHS_sparse_row(cc) = c;
            LHS_sparse_col(cc) = (col_idx-1)*y_dim+up_row_idx;
            LHS_sparse_val(cc) = -T_up/dy^2;
            cc = cc + 1;
            % coefficient of j+1,i (down)
            LHS_sparse_row(cc) = c;
            LHS_sparse_col(cc) = (col_idx-1)*y_dim+down_row_idx;
            LHS_sparse_val(cc) = -T_down/dy^2;
            cc = cc + 1;
            % coefficient of j,i-1 (left)
            LHS_sparse_row(cc) = c;
            LHS_sparse_col(cc) = (left_col_idx-1)*y_dim+row_idx;
            LHS_sparse_val(cc) = -T_left/dx^2;
            cc = cc + 1;
            % coefficient of j,i+1 (right)
            LHS_sparse_row(cc) = c;
            LHS_sparse_col(cc) = (right_col_idx-1)*y_dim+row_idx;
            LHS_sparse_val(cc) = -T_right/dx^2;
            cc = cc + 1;        

            % set the right hand side (B)
            RHS_sparse_row(c) = c;
            RHS_sparse_col(c) = 1;
            RHS_sparse_val(c) = R(row_idx,col_idx) + ...
                                (S(row_idx,col_idx)/dt)*Hold(row_idx,col_idx);
            c = c + 1;
        end
    end
    
    % pack everything in a sparse matrix to save memory
    LHS = sparse(LHS_sparse_row,LHS_sparse_col,LHS_sparse_val);
    RHS = sparse(RHS_sparse_row,RHS_sparse_col,RHS_sparse_val);
    
    % solve the linear system of equations
    sol = pcg(LHS,RHS,10^-6,300);
    
    % post-processing
    % reshape-to (row,col) format
    Hnew = zeros(topo_size);
    ck = 1;
    for i = 1:x_dim
        for j = 1:y_dim
            Hnew(j,i) = sol(ck);
            ck = ck + 1;
        end
    end
    
    % calculate groundwater runoff (anything above topography)
    runoff = max(0,(Hnew - model_topography).*S); 

    % adjust head after taking runoff out
    Hnew = Hnew - (runoff./S);
    
    % adjust for inactive (sea) cells
    for i = 1:length(r_inac)
        Hnew(r_inac(i),c_inac(i)) = sea_level;
    end
end
