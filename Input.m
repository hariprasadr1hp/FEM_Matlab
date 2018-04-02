%Assigment Code for FEM

%by Hari Prasad Radhakrishnan - 61492
tic()

%_________________________________INPUT________________________________________

%selecting mode 
%mode 1 = calculating elastic convergence
%mode 2 = calculating elastic + plastic strains
mode = 2;    %Enter 'mode' value either 1 or 2 accordingly

%Defining Parameters for Variant #1
p_max = 50; % maximum pressure
E = 70000; %young modulus
Nu = 0.25 ; % Poisson's ratio
yieldStress = 70; % yield stress


%Defining  final time, and delta time step (deltaT) 
final_t = 1;
del_t = 0.05;
GP = 2;

%Defining inner and outer radius
a = 40 ; % inner radius in meter
b = 80;  % outer radius in meter

%Defining how many elements we want to have
xe = 8 ;% takes arbitrary value


%params = [E, Nu, yieldStress];



%_____________________________END OF INPUT_________________________________

if mode == 1
    disp("You have chosen mode 1. The behavior of the model in the elastic regime is studied");
    [u,u_exact,r] = Elasticconvergence(a,b, E, Nu, yieldStress, p_max,xe,final_t, del_t,GP);
    f1=figure;
    plot([r(:,1) ;r(end,2)], u(1:(end))', [r(:,1) ;r(end,2)], u_exact(1:(end))','o')
    legend('u_{numerical}','u_{exact}')
    title('Converge Study in Elastic Regime')
    ylabel('u(r) in mm')
    xlabel('r in mm')
    

elseif mode == 2    
    fprintf("You have chosen mode 2. The behavior of the model in both elastic and plastic regime are studied");
    [u,u_exact,u_history,r,stress_rr_history, stress_phiphi_history] = Main(a,b, E, Nu, yieldStress, p_max,xe,final_t, del_t,GP);
    f1=figure;
    plot((0:del_t:final_t),u_history(end,:));
    title('Widening of the Outer Radius (b) over time');
    ylabel('u(r=b) mm');
    xlabel('time t');
    f2=figure;
    plot(r(:,2),stress_rr_history(:))
    title('\sigma_{rr} distribution');
    ylabel('\sigma_{rr} MPa');
    xlabel('r  mm');
    f3=figure;
    plot(r(:,2),stress_phiphi_history(:));
    title('\sigma_{\phi \phi} distribution');
    ylabel('\sigma_{\phi \phi} MPa');
    xlabel('r in mm');
    f4=figure;
    plot([r(:,1) ;r(end,2)], u_history(:,end));
    title('Displacement Distribution of u(r) at t=1');
    ylabel('u in mm');
    xlabel('r in mm');
else
    fprintf("WARNING: A wrong mode is chosen!. Please enter the mode value either 1 or 2");
end
    
        
toc() 
        
