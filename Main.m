%Main function
function [u,u_exact,u_history,r, stress_rr_history, stress_phiphi_history] = Main(a,b, E, Nu, yieldStress, p_max,xe,final_t,del_t,GP)

params = [E, Nu, yieldStress];
ne = xe + 1 ;   % no of nodes


r = zeros(xe, 2);       %Defining r vector 
delta_r = (b-a)/xe;
temp_r = a;

for elem = 1 : xe
    r(elem,1) = temp_r;
    r(elem,2) = r(elem,1) + delta_r;
    temp_r = r(elem,2);
end 

% allocate memory for matrices
% U_r displacement value allocates memory which depends on number of global iteration
u = zeros(ne,1) ;
u_exact = zeros(ne,1) ;
svarsGP     = zeros(3, xe); %svarsGP plastic strains as internal variables = (eps_prr, eps_pphiphi, eps_pzz)

%Allocate for every global matrix
k_global   = zeros(ne,ne);  %Tangent stiffness global matrix 
f_ext_global = zeros(ne,1);   %F external global matrix
f_int_global = zeros(ne,1);   %F internal global matrix
del_u     = zeros(ne,1);

%allocate history 
u_history   =zeros(ne, final_t/del_t+1 );
stress_rr_history = zeros(xe, 1);
stress_phiphi_history = zeros(xe, 1);

%%____________________________________ITERATION____________________________ %%
t = 0 ; % time initiation  
iter=1;
p_init = pl_init(yieldStress,a,b,Nu); %pressure of plasticity begin
while t < final_t

    t = t + del_t
    p = p_max * t ;
    
    if p>p_init
        fprintf('The current pressure is larger than the plasticity-initiation pressure. The model is currently in plastic regime\n')
    end
    
    f_ext_global(1,1) = p*a;
    converge_iter=-1;
    
    while ~((norm(f_int_global-f_ext_global,inf)<.005*norm(f_int_global,inf)) && (norm(del_u,inf)<.005*norm(u,inf)))
        
        f_int_global = 0*f_int_global;
        k_global = 0*k_global;
        
        %assembly part from elements to global matrix
        for  elem = 1 : xe % loop first element to El element
            [k_el, f_int_el,svarsGPNew,stress2D] = Elementroutine(r(elem,:)',svarsGP(:,elem)',params,u(elem:elem+1),GP);
            svarsGP(:,elem) = svarsGPNew;
            stress_rr_history(elem,1) = stress2D(1,1);
            stress_phiphi_history(elem,1) = stress2D(2,1);
            
            %Build Area Connectivity matrix AA for global tangent matrix
            Ae = connectivity(elem,ne);            
            % assemblying element tangent stiffness into global tangent
            k_global = k_global + (Ae'* k_el * Ae) ; 
            % assemblying internal force into F internal global
            f_int_global = f_int_global + (Ae'*f_int_el) ;
        end       
        
        G = f_int_global - f_ext_global;
        del_u = linsolve(k_global,-G);
        u = u + del_u;
        converge_iter= converge_iter+ 1;
        
        %exact solution calculation
        for elem=1:ne
            R=a+(elem-1)*delta_r;
            [u_exact(elem,1)] = Exact(p, R, a, b,params);
        end                
        
    end 
    %converge_iter
    %keeping displacement to history
    iter=iter + 1;
    u_history(:,iter) = u; 
       
end



