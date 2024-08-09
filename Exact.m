%___________________________EXACT CALCULATION_________________
function[u_exact] = exact(p, R, a, b, params)
    E=params(1); %Young's modulus
    Nu=params(2); %Poisson ratio
    u_exact = (1+Nu)*(p/E)*(a^2/(b^2-a^2))*((1-2*Nu)*R+b^2/R); %exact solution for displacement
end
