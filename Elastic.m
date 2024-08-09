%convergeElastic function 
function [u,u_exact,r] = ElasticConvergence(r_in, r_out, E, nu, yield_stress, p_max, elems, t_final, dt, GP)

params = [E, nu, yield_stress];
nodes = elems + 1 ;

r = zeros(elems, 2);   %Defining r vector
dr = (r_out - r_in) / elems;
r_temp = r_in;
for i = 1 : elems
    r(i,1) = r_temp;
    r(i,2) = r(i,1) + dr;
    r_temp = r(i,2);
end 

% allocate memory for matrices
% U_r displacement value allocates memory which depends on number of global iteration
u = zeros(nodes, 1);
u_exact = zeros(nodes, 1);

svarsGP = zeros(3, elems); %svarsGP plastic strains as internal variables = (eps_prr, eps_pphiphi, eps_pzz)

%Allocate for every global matrix
k_global = zeros(nodes, nodes) ;  %Tangent stiffness global matrix 
f_ext_global = zeros(nodes, 1);   %F external global matrix
f_int_global = zeros(nodes, 1);   %F internal global matrix
del_u = zeros(nodes, 1);

%allocate history 
u_history   = zeros(nodes, t_final / dt + 1);

%%____________________________________ITERATION____________________________ %%
t = 0 ; % time initiation  
iter = 1;

p = 0;
p_init = yield_stress * (r_out^2 - r_in^2) / r_in^2 * (1 / (sqrt(1 - 4 * nu + 4 * nu^2 + 3 * (r_out / r_in)^4))); %pressure of plasticity begin
while p < p_init
    t = t + dt; %time increment 
    p = p_max * t ; % pressure of function of time
    f_ext_global(1, 1) = p * r_in; 
    con_iter = -1;
    %checking convergence by take the given convergence criterion
    while ~((norm(f_int_global - f_ext_global, inf) < .005 * norm(f_int_global, inf)) && (norm(del_u, inf) < .005 * norm(u,inf)))
        f_int_global = 0 * f_int_global;
        k_global = 0 * k_global;
        %assembly part from elements to global matrix
        for  elem = 1 : elems % loop first element to El element
            [Kt_e, Fint_e, svarsGPNew] = ElementRoutine(r(elem, :)', svarsGP(:, elem)', params, u(elem : elem + 1), GP);
            svarsGP(:, elem) = svarsGPNew;
            
            %Build Area Connectivity matrix AA for global tangent matrix
            Ae = connectivity(elem, nodes);

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
        for node=1:nodes
            R = r_in + (node - 1) * dr;
            [u_exact(node, 1)] = Exact(p, R, r_in, r_out, params);
        end 
        
       
        
    end 
    %keeping displacement to history
    con_iter
    iter = iter + 1;
    u_history(:, iter) = u; 
   
    
end



