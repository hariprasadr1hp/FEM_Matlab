% %Assignment for Selected Topics of FEM in summer term 2017
% % lecturer: Dr. Geralf Huetter, TU Bergakademie Freiberg
% %
% %Material routine for elastic ideal-plastic material
% %  - plain strain, axisymmetric (i.e. without shear component)
% %  - elastic behavior: isotropic linear
% %  - plastic behavior: ideal Mises plasticity
% %  - Euler backward time integration
% %  - algorithmically consitent material tangent stiffness
% %
% % reference: 
% %  Belytschko, Liu, Moran: "Nonlinear Fnite Elements fo Continua and
% %  Structures", Sect. 5.9.5
% %
% % INPUT: 
% % - eps2D  : strain in directions (rr, phiphi)
% % - deps2D : corresponding strain increment
% % - dt     : time increment 
% % - svarsGP: plastic strains as internal variables from last time increment = (eps_prr,eps_pphiphi,eps_pzz)
% % - params : material parameters = (E, nu, yield_stress)

% % OUTPUT: 
% % - stress2D  : updated stress in plane (sigma_rr,sigma_phiphi)
% % - matstiff  : consistent material tangent stiffness
% % - svarsGPnew: updated values of state variables svarsGP


function [stress2D,matstiff,svarsGPnew] = Materialroutine(eps2D, deps2D, dt, svarsGP, params)
%function [stress2D,matstiff,svarsGPnew] = Materialroutine(eps2D,svarsGP,params)
%extract parameter values
E = params(1); %Young's modulus
nu = params(2); %Poisson ratio
yield_stress = params(3); %yield stress

%extract old values of state variables
eps_p_old = svarsGP(1:3)'; %last values of plastic strain in (rr, phiphi, zz)
eps = [eps2D; 0]; % strain in all 3 directions (rr, phiphi, zz)
eps_el = eps - eps_p_old; %predictor of elastic strain as difference between total and plastic strain
epsTrace = sum(eps_el); %trace of strain = volumetric strain
% stress predictor according to Hooke's law for isotropic material (in terms of Lame constants)
%   ("predictor"=assuming that yield condition is absent and thus plastic
%   deformations do not change
Lamee1 = E / (1 + nu);
Lamee2 = Lamee1 * nu / (1 - 2 * nu);
%stress predictor ("trial stress")
stress_pred = Lamee1 * eps_el + Lamee2 * epsTrace; 
stress_dev = stress_pred - sum(stress_pred) / 3; %stress deviator
stress_eq = sqrt(1.5 * stress_dev' * stress_dev); %v. Mises equivalent stress

%check plasticity:
if (stress_eq < yield_stress) %elastic case:
    stress2D = stress_pred(1:2); %stress identical to predictor
    svarsGPnew(1:3) = svarsGP(1:3);        % no update of plastic strains
    %tangent stiffness identical to elastic stiffness
    matstiff(1,1) = Lamee1 + Lamee2;
    matstiff(2,2) = matstiff(1,1);
    matstiff(1,2) = Lamee2;
    matstiff(2,1) = Lamee2;
else %plastic case: return-mapping of stress deviator to yield surface
    scaleparam = yield_stress / stress_eq; %scaling parameter for stress deviator
    stress2D = stress_pred(1:2) - (1 - scaleparam) * stress_dev(1:2); %update stress
    % update state variables (plastic strain)
    svarsGPnew(1:3) = svarsGP(1:3) + (1 - scaleparam) * stress_dev' / Lamee1; %plastic strains    
    % algorithmically consistent tangent stiffness
    matstiff = ones(2) * (Lamee2 + Lamee1 / 3) + scaleparam * ([2 -1; -1 2] * Lamee1 / 3 - (1.5 * Lamee1 / (stress_eq^2)) * stress_dev(1:2) * stress_dev(1:2)');
end
