function [A1, A2, A3, D1, D2, D3] = TDMA_v4(Tr, BC, farm)
% Tridiagonal matrices generated once per transport function
% 
% Input: (Bx,By,Bz,farm)
%            Bx, By, Bz are diffusion constants
%            farm = calls on dx, dy, dz_tr
%
% Output: 3 structures of tridiagonal matrices
%         1. x-implicit {y,z}
%         2. y-implicit {x,z}
%         3. z-implicit {x,y}
%
%

% Temporary variables for matrix sizes
dim_x = size(farm.gridx_tr,1);
dim_y = size(farm.gridx_tr,2);
dim_z = size(farm.gridx_tr,3);


%% X-IMPLICIT
% Tridiagonal matrix

    % Diffusion fully implicit (advection fully explicit)
    % Left:     -Bx
    % Center:   1 + 2Bx
    % Right:    -Bx

    % preallocate space
    A1 = cell(dim_y,1); 

    for j = 1:dim_y

        % TRIDIAGONAL MATRIX

        % Analyze as sheets of x,z with a for loop for y
        Bx_Plus  = reshape(Tr.BxPlus(:,j,:),dim_x*dim_z,1);
        Bx_Minus = reshape(Tr.BxMinus(:,j,:),dim_x*dim_z,1);

        A1{j} = spdiags([1 + Bx_Plus + Bx_Minus, -Bx_Minus, -Bx_Plus], ... % center, left, right
            [0 -1 1], ... % column diagonal index
            dim_x*dim_z, dim_x*dim_z+1); % size of matrix -> necessary to have the +1 on second dimension so that diagonals are organized horizontally across
        A1{j} = A1{j}(1:dim_x*dim_z,1:dim_x*dim_z); % get rid of extra column...
        
            % Repeating rows
            % x = 1; end
               next = dim_x:dim_x:dim_x*dim_z-dim_x;
               A1{j}(next+1,next) = 0;
               A1{j}(next,next+1) = 0;

        clear next Bx_Plus Bx_Minus

    end
    clear j


    % Boundary Contition: D matrix
    % D1
    
    D1 = cell(dim_y,1); % preallocate space

        % Identify which cells correspond to boundaries
        bcl = 1:dim_x:dim_z*dim_x;
        bcr = dim_x:dim_x:dim_z*dim_x;
        BC_t = reshape(repmat(BC,1,dim_x)',dim_z*dim_x,1);

    for j = 1:dim_y

        D1{j} = zeros(dim_x*dim_z,1);

        Bx_t = reshape(Tr.Bx(:,j,:),dim_z*dim_x,1);

        D1{j}(bcl) = BC_t(bcl) .* Bx_t(bcl);
        D1{j}(bcr) = BC_t(bcr) .* Bx_t(bcr);

    end
                
                
%% Y-IMPLICIT
% Tridiagonal matrix

    % preallocate space  
    A2 = cell(dim_z,1); 

    for k = 1:dim_z

        % TRIDIAGONAL MATRIX

        % Analyze as sheets of x,y with a for loop for z
        By_Plus = reshape(Tr.ByPlus(:,:,k),dim_x*dim_y,1);
        By_Minus = reshape(Tr.ByMinus(:,:,k),dim_x*dim_y,1);
        
        A2{k} = spdiags([1 + By_Minus + By_Plus, -By_Minus, -By_Plus], ... % center, left, right
            [0 -1 1], ... % column diagonal index
            dim_x*dim_y, dim_x*dim_y+1); % size of matrix -> necessary to have the +1 on second dimension so that diagonals are organized horizontally across
        A2{k} = A2{k}(1:dim_x*dim_y,1:dim_x*dim_y); % get rid of extra column...
       
            % Repeating rows
            % x = 1; end

                   next = dim_y:dim_y:dim_y*dim_x-dim_y;
                   A2{k}(next+1,next) = 0;
                   A2{k}(next,next+1) = 0;

        clear next By_Plus By_Minus

    end
    clear k


    % Boundary Contition: D matrix
    % D2
    
    D2 = cell(dim_z,1);

            bcl = 1:dim_y:dim_y*dim_x;
            bcr = dim_y:dim_y:dim_y*dim_x;

    for k = 1:dim_z

        % Boundary Condition
        D2{k} = zeros(dim_x*dim_y,1);

        BC_t = reshape(repmat(BC(k),dim_y,dim_x),dim_y*dim_x,1);
        By_t = reshape(Tr.By(:,:,k)',dim_y*dim_x,1);

        D2{k}(bcl) = BC_t(bcl) .* By_t(bcl);
        D2{k}(bcr) = BC_t(bcr) .* By_t(bcr);

    end
                
                
%% Z-IMPLICIT
% Tridiagonal matrix

        % preallocate space
        A3 = cell(dim_y,1); 

        for j = 1:dim_y

            % TRIDIAGONAL MATRIX

            % Analyze as sheets of z,x with a for loop for y
            Bz_Plus  = reshape(squeeze(Tr.BzPlus(:,j,:))',dim_x*dim_z,1);
            Bz_Minus = reshape(squeeze(Tr.BzMinus(:,j,:))',dim_x*dim_z,1);

            A3{j} = spdiags([1 + Bz_Minus + Bz_Plus, -Bz_Minus, -Bz_Plus], ... % center, left, right
                [0 -1 1], ... % column diagonal index
                dim_x*dim_z, dim_x*dim_z+1); % size of matrix -> necessary to have the +1 on second dimension so that diagonals are organized horizontally across
            A3{j} = A3{j}(1:dim_x*dim_z,1:dim_x*dim_z); % get rid of extra column...
       
                % Repeating rows
                % x = 1; end

                       next = dim_z:dim_z:dim_z*dim_x-dim_z;
                       A3{j}(next+1,next) = 0;
                       A3{j}(next,next+1) = 0;

            clear next Bz_z

        end
        clear j

                
    % Boundary Contition: D matrix
    % D3
    
        D3 = cell(dim_y,1);

                bcl = 1:dim_z:dim_z*dim_x; % surface
                bcr = dim_z:dim_z:dim_z*dim_x; % bottom
                BC_t = reshape(repmat(BC,1,dim_x),dim_z*dim_x,1);

        for j = 1:dim_y

            % Boundary Condition
            D3{j} = zeros(dim_x*dim_z,1);

            Bz_tp = reshape(squeeze(Tr.BzPlus(:,j,:))',dim_z*dim_x,1);
            Bz_tm = reshape(squeeze(Tr.BzMinus(:,j,:))',dim_z*dim_x,1);

            D3{j}(bcl) = BC_t(bcl+1) .* Bz_tm(bcl); % Neumann boundary condition; -z = +z; replace BC_t with C at z=2
            D3{j}(bcr) = ((BC_t(bcr)-BC_t(bcr-1))+BC_t(bcr)) .* Bz_tp(bcr); % interpolate to +1 beyond via interpolation
           
        end

 
end