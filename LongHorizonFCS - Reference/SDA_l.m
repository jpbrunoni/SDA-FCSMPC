function[U_op,Niter,Rs2]=SDA_l(U_bar_unc, H, u_min, u_max, p, U_ini,Niter_max) 

%% Input Parameter defination
%
% U_bar_unc is the Unconstrained solution in the transformed space
% (hypersphere center)
% H = Lattice generator matrix (Lower Triangular)
%
% x(k+1)=Ax(k)+Bu(k) 
% x\in \R^m; u\in \R^n_u
% 
% u_i \in \V = \{ u_{\min},...,u_{\max} \} \subset \Z
% u \in \V^n_u \subset \Z^n_u
%
% \U_k = \V^n_u Np \subset \Z^p; p=n_u Np
%
% min |Yc - H*U_k|^2
%
% U_ini is the initial vector from the lattice
% d = An one dimensional array to store the calculated distance; d=zeros(1,(p+1));
% U = Input vector; U=zeros(1,p);
% U_op = The optimal input vector; U_op=zeros(1,p);

Niter=0; %variable that shows the number of iterations
d=zeros((p+1),1); 
U_k=zeros(p,1); 
U_op=U_ini; 

%% Initial Radius of Sphere (Squared value)
Rs2 = ceil((norm(U_bar_unc - H*U_ini))^2); %rounds the squared value to the nearest integer towards infinity
Rs2_ini = Rs2;


%%%%%%%%%%%%%%%%%% Sphere Decoding  Algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic
i = 0; C=0;
% i is the index of input vector dimension (i = 1 -> phase A; i = 2 -> B; i = 3 -> C)
% C is the index of input vector value (C = 0 or 1, i.e., the switching states)

while (i >=0) &&(Niter < Niter_max) % Optimal or early termination
    i=i+1; %i goes to the first layer
    j = u_min + C; %as u_min is 0 and C is 0(at the beggining), j starts in 0 - That is the value that will be considered in the i-th layer
    
    if j <= u_max %if j <= 1
        U_k(i) = j; %The vector U_k in the layer i is j (1st: i = 1, j = 0, U_k(1) = 0)
        d(i+1) = ((H(i,1:i)*U_k(1:i) - U_bar_unc(i))^2 + d(i)); %compute the accumulated distance (among the layers) between the computed U_k and the unconstrained solution, both in the transformed space, and stores in i+1                                 
        Niter=Niter+1; %increase in 1 the number of iterations
        if d(i+1) > Rs2 % If the distance is greater than the radius of the sphere
            i=i-1; C=C+1; %i decreases by 1 and C increases by 1 (this way, the U(i) = 1 will be tested, once that, in the next iteration, i=i+1)
        else            % If the distance is lower than the radius of the sphere
            if i==p           % If all the entire sequence has been computed
               U_op = U_k;    % Update the optimal solution
               Rs2 = d(i+1);  % Save the current squared radius
                    
               i=i-1; C=C+1; %i decreases by 1 and C increases by 1 in order to test U(i) = 1
            else        % Keep exploring sequence
               C=0;         %resets C
            end % if i==p         
        end % if d(i+1) > Rs2      
    else  %if j > 1
        i=i-2; % i decreases by 2 (if j>1 it means that U(i)=0 and U(i)=1 was computed, so it has to go back to the previous layer. Note that -2 is here because i=i+1 is in the next iteration)
        if i>=0 %if i >= 0, there are options to compute
            C=U_k(i+1) - u_min + 1; 
        end %if i>=0
    end % j <= u_max  
end % While
%toc
%%%%%%%%%%%%%%%%%%%%%%%% End of Sphere Decoding  Algorithm %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end