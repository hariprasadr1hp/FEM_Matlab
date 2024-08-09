%Main function
function [u, u_exact,u_history, r, stress_rr_history, stress_phiphi_history] = Main(r_in, r_out, E, nu, yield_stress, p_max, elems, t_final, dt, GP)

params = [E, nu, yield_stress];
nodes = elems + 1;   % no of nodes


r = zeros(elems, 2);       %Defining r vector 
dr = (r_out - r_in) / elems;
r_temp = r_in;

for elem = 1 : elems
    r(elem, 1) = r_temp;
    r(elem, 2) = r(elem, 1) + dr;
    r_temp = r(elem, 2);
end 

% allocate memory for matrices
% U_r displacement value allocates memory which depends on number of global iteration
u = zeros(nodes, 1);
u_exact = zeros(nodes, 1);
svarsGP = zeros(3, elems); %svarsGP plastic strains as internal variables = (eps_prr, eps_pphiphi, eps_pzz)

%Allocate for every global matrix
k_global = zeros(nodes, nodes);  %Tangent stiffness global matrix 
f_ext_global = zeros(nodes, 1);   %F external global matrix
f_int_global = zeros(nodes, 1);   %F internal global matrix
du = zeros(nodes, 1);

%allocate history 
u_history = zeros(nodes, t_final / dt + 1 );
stress_rr_history = zeros(elems, 1);
stress_phiphi_history = zeros(elems, 1);

%%____________________________________ITERATION____________________________ %%
t_current = 0 ; % time initiation  
iter = 1;
p_init = pl_init(yield_stress, r_in, r_out, nu); %pressure of plasticity begin
while t_current < t_final

    t_current = t_current + dt
    p_current = p_max * t_current;
    
    if p_current > p_init
        fprintf('The current pressure is larger than the plasticity-initiation pressure. The model is currently in plastic regime\n')
    end
    
    f_ext_global(1, 1) = p_current * r_in;
    converge_iter = -1;
    
    while ~((norm(f_int_global - f_ext_global, inf) < .005 * norm(f_int_global, inf)) && (norm(du, inf) < .005 * norm(u, inf)))
        
        f_int_global = 0 * f_int_global;
        k_global = 0 * k_global;
        
        %assembly part from elements to global matrix
        for  elem = 1 : elems % loop first element to El element
            [k_el, f_int_el, svarsGPNew, stress2D] = ElementRoutine(r(elem, :)', svarsGP(:, elem)', params, u(elem : elem+1), GP);
            svarsGP(:, elem) = svarsGPNew;
            stress_rr_history(elem, 1) = stress2D(1, 1);
            stress_phiphi_history(elem, 1) = stress2D(2, 1);
            
            %Build Area Connectivity matrix AA for global tangent matrix
            Ae = connectivity(elem, nodes);            
            
            % assemblying element tangent stiffness into global tangent
            k_global = k_global + (Ae' * k_el * Ae); 
            
            % assemblying internal force into F internal global
            f_int_global = f_int_global + (Ae' * f_int_el);
        end       
        
        G = f_int_global - f_ext_global;
        du = linsolve(k_global,-G);
        u = u + du;
        converge_iter = converge_iter+ 1;
        
        %exact solution calculation
        for elem=1:nodes
            R = r_in + (elem - 1) * dr;
            [u_exact(elem, 1)] = Exact(p_current, R, r_in, r_out, params);
        end                
        
    end 
    %converge_iter
    %keeping displacement to history
    iter = iter + 1;
    u_history(:, iter) = u; 
       
end


