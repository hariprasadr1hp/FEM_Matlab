% Element routine for the problem 
% Input : xe - number of elements
 

function [k_el_GP,f_int_el_GP,svarsGPNew,stress2D] = Elementroutine(r,svarsGP,params,u,GP)


[xi,weight] = lgwt(GP,-1,1);
le = r(2) - r(1);               %distance between two nodes in an element
J = le/2 ;                      % Jacobian value
f_int_el_GP = zeros(2,1);
k_el_GP = zeros(2,2);


for i = 1:GP
    N = [0.5*(1-xi(i,1)) 0.5*(1+xi(i,1))]'; %linear shape function
    B = [-1/2*J^-1, 1/2*J^-1 ; N(1)/(N(1)*r(1)+N(2)*r(2)), N(2)/(N(1)*r(1)+N(2)*r(2))] ;% B matrix = operatorB * N
    eps2D = B * u ;% 2D strain
    [stress2D,matstiff,svarsGPNew] = Materialroutine(eps2D,0,0,svarsGP,params);
    k_el =  weight(i,1) * B' * matstiff * B * (N'*r )' *  J ;% Stiffness tangent matrix  
    f_int_el = weight(i,1) * B' * stress2D * N'* r * J ;% element F internal
    f_int_el_GP = f_int_el_GP + f_int_el;
    k_el_GP = k_el_GP + k_el;
end

end
