%convergeElastic function 
function [u,u_exact,r] = Elasticconvergence(a,b, E, Nu, yieldStress, p_max,xe,final_t,del_t,GP)

params = [E, Nu, yieldStress];
node = xe + 1 ;

r = zeros(xe, 2);   %Defining r vector
delta_r = (b-a)/xe;
temp_r = a;
for i = 1 : xe
    r(i,1) = temp_r;
    r(i,2) = r(i,1) + delta_r;
    temp_r = r(i,2);
end 

% allocate memory for matrices
% U_r displacement value allocates memory which depends on number of global iteration
u = zeros(node,1) ;
u_exact = zeros(node,1) ;

svarsGP     = zeros(3, xe); %svarsGP plastic strains as internal variables = (eps_prr, eps_pphiphi, eps_pzz)

%Allocate for every global matrix
k_global   = zeros(node,node) ;  %Tangent stiffness global matrix 
f_ext_global = zeros(node,1);   %F external global matrix
f_int_global = zeros(node,1);   %F internal global matrix
del_u     = zeros(node,1);

%allocate history 
u_history   = zeros(node, final_t/del_t +1);

%%____________________________________ITERATION____________________________ %%
t = 0 ; % time initiation  
iter=1;

p=0;
p_init = yieldStress*(b^2-a^2)/a^2 *(1/(sqrt(1-4*Nu + 4*Nu^2 + 3*(b/a)^4))); %pressure of plasticity begin
while p < p_init
    t = t + del_t; %time increment 
    p = p_max * t ; % pressure of function of time
    f_ext_global(1,1) = p*a; 
    con_iter=-1;
    %checking convergence by take the given convergence criterion
    while ~((norm(f_int_global-f_ext_global,inf)<.005*norm(f_int_global,inf)) && (norm(del_u,inf)<.005*norm(u,inf)))
        f_int_global = 0*f_int_global;
        k_global = 0*k_global;
        %assembly part from elements to global matrix
        for  i = 1 : xe % loop first element to El element
            [Kt_e, Fint_e,svarsGPNew] = Elementroutine(r(i,:)', svarsGP(:,i)', params, u(i:i+1),GP);
            svarsGP(:,i) = svarsGPNew;
            
            %Build Area Connectivity matrix AA for global tangent matrix
            Ae = connectivity(elem,ne);

            % assemblying element tangent stiffness into global tangent
            k_global = k_global + Ae'* Kt_e * Ae ; 
            % assemblying internal force into F internal global
            f_int_global = f_int_global + Ae'*Fint_e ;
        end
       
        
        G = f_int_global - f_ext_global; %calculating G matrix
        del_u = linsolve(k_global,-G); % calculating delta_u
        u = u + del_u; %increment of displacement by delta_u
        con_iter = con_iter+ 1;
        %exact solution calculation
        for i=1:node
            R=a+(i-1)*delta_r;
            [u_exact(i,1)] = Exact(p, R, a, b,params);
        end 
        
       
        
    end 
    %keeping displacement to history
    con_iter
    iter=iter + 1;
    u_history(:,iter) = u; 
   
    
end



