function [Nutrient_out] = ADI_3D_v4(Nutrient_in,Tr)
% This function solves three-dimensional advection-diffusion equation in a
% non-homogenous flow environment using Alternating Direction Implicit
% (ADI) method.
%
% Inputs:
%       Nutrient_in -> 3-D nutrient field
%       Ax, Ay, Az, Bx, By, Bz -> Advection and diffusion field constants
%       A1, A2, A3 -> sets of tridiagonal matrices to solve LHS
%       FARM_DIM -> dimensions of farm for reshaping arrays
%
% Outputs:
%       Nutrient_out -> 3-D nutrient field for one transport step
%

%% Nutrient Field

           % Nt is the ghost matrix already generated in transport_setup
           % with nutrient-specific boundary conditions
           % Insert nutrient field into ghost matrix
           
                Nt = Tr.Nt;
                Nt(2:end-1,2:end-1,1:end-1) = Nutrient_in;
            
                % Discretization schemes are centered difference, so create
                % explicit matrices for at x+1, x-1, y+1, y-1, z+1, z-1
                
                    % For subsurface layer
                    Ntin = Nt(2:end-1,2:end-1,2:end-1);
                    NtxL = Nt(1:end-2,2:end-1,2:end-1);
                    NtxR = Nt(3:end,2:end-1,2:end-1);
                    NtyL = Nt(2:end-1,1:end-2,2:end-1);
                    NtyR = Nt(2:end-1,3:end,2:end-1);
                    NtzD = Nt(2:end-1,2:end-1,1:end-2); % Up is equivalent to z-1
                    NtzU = Nt(2:end-1,2:end-1,3:end); % Down is equivalent to z+1 
            
                    % For surface layer
                    Ntin_s = Nt(2:end-1,2:end-1,1);
                    NtxL_s = Nt(1:end-2,2:end-1,1);
                    NtxR_s = Nt(3:end,2:end-1,1);
                    NtyL_s = Nt(2:end-1,1:end-2,1);
                    NtyR_s = Nt(2:end-1,3:end,1);
                    NtzD_s = Nt(2:end-1,2:end-1,2); % Neumann Boundary Condition
                    NtzU_s = Nt(2:end-1,2:end-1,2);  
                    
                    
%% Explicit Advection

            % Advection is explicit (RHS of equation only) and so the
            % following is used at each implicit step. More efficient to
            % solve once and insert solution into each implicit step.
            
                Adv_explicit_sub = ... % subsurface matrix
                        + Tr.absAx_pos(:,:,2:end) .* NtxL ...
                        + Tr.absAx_neg(:,:,2:end) .* NtxR ...
                        + Tr.absAy_pos(:,:,2:end) .* NtyL ...
                        + Tr.absAy_neg(:,:,2:end) .* NtyR ...
                        + Tr.absAz_pos(:,:,2:end) .* NtzU ...
                        + Tr.absAz_neg(:,:,2:end) .* NtzD;

                Adv_explicit_sur = ... % surface matrix
                        + Tr.absAx_pos(:,:,1) .* NtxL_s ...
                        + Tr.absAx_neg(:,:,1) .* NtxR_s ...
                        + Tr.absAy_pos(:,:,1) .* NtyL_s ...
                        + Tr.absAy_neg(:,:,1) .* NtyR_s ...
                        + Tr.absAz_pos(:,:,1) .* NtzU_s ...
                        + Tr.absAz_neg(:,:,1) .* NtzD_s;
                            
 
% Transport dimensions
dimx = size(Nutrient_in,1);
dimy = size(Nutrient_in,2);
dimz = size(Nutrient_in,3);

%% ADI Method: First in X-direction
                
                % RHS solution
                
                N_RHS1 = NaN(size(Nutrient_in)); % preallocate space
                
                    % Subsurface
                    N_RHS1(:,:,2:end) = ...
                        Adv_explicit_sub ...
                        + Tr.ByMinus(:,:,2:end) .* NtyL ...
                        + Tr.ByPlus(:,:,2:end) .* NtyR ...
                        + Tr.BzMinus(:,:,2:end) .* NtzU ...
                        + Tr.BzPlus(:,:,2:end) .* NtzD ...
                        + (1 - Tr.absAx(:,:,2:end) - Tr.absAy(:,:,2:end) - Tr.absAz(:,:,2:end) - Tr.ByPlus(:,:,2:end) - Tr.ByMinus(:,:,2:end) - Tr.BzPlus(:,:,2:end) - Tr.BzMinus(:,:,2:end)) .* Ntin;


                    % Surface; z = 1; Neumann (-z = +z)
                    N_RHS1(:,:,1) = ...
                        Adv_explicit_sur ...
                        + Tr.ByMinus(:,:,1) .* NtyL_s ...
                        + Tr.ByPlus(:,:,1) .* NtyR_s ...
                        + Tr.BzMinus(:,:,1) .* NtzU_s ...
                        + Tr.BzPlus(:,:,1) .* NtzD_s ...
                        + (1 - Tr.absAx(:,:,1) - Tr.absAy(:,:,1) - Tr.absAz(:,:,1) - Tr.ByPlus(:,:,1) - Tr.ByMinus(:,:,1) - Tr.BzPlus(:,:,1) - Tr.BzMinus(:,:,1)) .* Ntin_s;
     
                  
           % LHS Solution
           
                % N_LHS1 = A1\(N_RHS1+D1)
            
                N_LHS1 = NaN(size(N_RHS1)); % create temporary matrix
               
                for j = 1:dimy
                   
                    % TDMA solver
                    temp = reshape(N_RHS1(:,j,:),dimx*dimz,1);
                    temp = Tr.A1{j}\(temp+Tr.D1{j});
                    N_LHS1(:,j,:) = reshape(temp,dimx,dimz);
                    
                end
                
                                
% Y-IMPLICIT

           % Nt is the ghost matrix already generated in transport_setup
           % with nutrient-specific boundary conditions
           % Insert nutrient field into ghost matrix
           
                N_LHS1t = Tr.N_LHS1t;
                N_LHS1t(2:end-1,2:end-1,1:end-1) = N_LHS1;
                   
           % RHS Solution     

                N_RHS2 = NaN(size(Nutrient_in)); % preallocate space
                
                    % Subsurface
                    N_RHS2(:,:,2:end) = ...
                        Adv_explicit_sub ...
                        + Tr.BzMinus(:,:,2:end) .* NtzU ...
                        + Tr.BzPlus(:,:,2:end) .* NtzD ...
                        + (1 - Tr.absAx(:,:,2:end) - Tr.absAy(:,:,2:end) - Tr.absAz(:,:,2:end) - Tr.BzPlus(:,:,2:end) - Tr.BzMinus(:,:,2:end)) .* Ntin ...
                        + Tr.BxMinus(:,:,2:end) .* N_LHS1t(1:end-2,2:end-1,2:end-1) ...
                        + (-Tr.BxMinus(:,:,2:end) -Tr.BxPlus(:,:,2:end)) .* N_LHS1t(2:end-1,2:end-1,2:end-1) ...
                        + Tr.BxPlus(:,:,2:end) .* N_LHS1t(3:end,2:end-1,2:end-1);


                    % Surface, z = 1; Neumann (-z = +z)
                    N_RHS2(:,:,1) = ...
                        Adv_explicit_sur ...
                        + Tr.BzMinus(:,:,1) .* NtzU_s ...
                        + Tr.BzPlus(:,:,1) .* NtzD_s ...
                        + (1 - Tr.absAx(:,:,1) - Tr.absAy(:,:,1) - Tr.absAz(:,:,1) - Tr.BzPlus(:,:,1) - Tr.BzMinus(:,:,1)) .* Ntin_s ...
                        + Tr.BxMinus(:,:,1) .* N_LHS1t(1:end-2,2:end-1,1) ...
                        + (-Tr.BxMinus(:,:,1) -Tr.BxPlus(:,:,1)) .* N_LHS1t(2:end-1,2:end-1,1) ...
                        + Tr.BxPlus(:,:,1) .* N_LHS1t(3:end,2:end-1,1);
                
                     
            % LHS Solution
            
                N_LHS2 = NaN(size(N_RHS2));
                                        
                for k = 1:dimz
                   
                    % TDMA solver
                    temp = reshape(N_RHS2(:,:,k)',dimx*dimy,1);
                    temp = Tr.A2{k}\(temp+Tr.D2{k}); % add B.C.
                    N_LHS2(:,:,k) = reshape(temp,dimy,dimx)';
                    
                end
                
                
% Z-IMPLICIT

           % N_LHS2t is the ghost matrix already generated in transport_setup
           % with nutrient-specific boundary conditions
           % Insert nutrient field into ghost matrix
           
            N_LHS2t = Tr.N_LHS2t;
            N_LHS2t(2:end-1,2:end-1,1:end-1) = N_LHS2;
                   
                    
           % Solve RHS; everything but surface layer           

                N_RHS3 = NaN(size(Nutrient_in)); % temporary matrix
                N_RHS3(:,:,2:end) = ...
                    Adv_explicit_sub ...
                    + (1 - Tr.absAx(:,:,2:end) - Tr.absAy(:,:,2:end) - Tr.absAz(:,:,2:end)) .* Ntin ...
                    + Tr.BxMinus(:,:,2:end) .* N_LHS1t(1:end-2,2:end-1,2:end-1) ...
                    + (-Tr.BxMinus(:,:,2:end) -Tr.BxPlus(:,:,2:end)) .* N_LHS1t(2:end-1,2:end-1,2:end-1) ...
                    + Tr.BxPlus(:,:,2:end) .* N_LHS1t(3:end,2:end-1,2:end-1) ...
                    + Tr.ByMinus(:,:,2:end) .* N_LHS2t(2:end-1,1:end-2,2:end-1) ...
                    + (-Tr.ByMinus(:,:,2:end) -Tr.ByPlus(:,:,2:end)) .* N_LHS2t(2:end-1,2:end-1,2:end-1) ...
                    + Tr.ByPlus(:,:,2:end) .* N_LHS2t(2:end-1,3:end,2:end-1);

                
                % z = 1; Surface boundary condition; Neumann (-z = +z)
                N_RHS3(:,:,1) = ...
                    Adv_explicit_sur ...
                    + (1 - Tr.absAx(:,:,1) - Tr.absAy(:,:,1) - Tr.absAz(:,:,1)) .* Ntin_s ...
                    + Tr.BxMinus(:,:,1) .* N_LHS1t(1:end-2,2:end-1,1) ...
                    +(-Tr.BxMinus(:,:,1) -Tr.BxPlus(:,:,1)) .* N_LHS1t(2:end-1,2:end-1,1) ...
                    + Tr.BxPlus(:,:,1) .* N_LHS1t(3:end,2:end-1,1) ...
                    + Tr.ByMinus(:,:,1) .* N_LHS2t(2:end-1,1:end-2,1) ...
                    + (-Tr.ByMinus(:,:,1) -Tr.ByPlus(:,:,1)) .* N_LHS2t(2:end-1,2:end-1,1) ...
                    + Tr.ByPlus(:,:,1) .* N_LHS2t(2:end-1,3:end,1);

                              
            % Solve LHS
                  
                N_LHS3 = NaN(size(Nutrient_in));
                      
                for j = 1:dimy
                    
                    % TDMA solver
                    temp = reshape(squeeze(N_RHS3(:,j,:))',dimx*dimz,1);
                    temp = Tr.A3{j}\(temp+Tr.D3{j});
                    N_LHS3(:,j,:) = reshape(temp,dimz,dimx)';
                    
                end
                
Nutrient_out = N_LHS3;  


end     