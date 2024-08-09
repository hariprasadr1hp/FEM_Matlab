%Assigment Code for FEM
tic()

%_________________________________INPUT___________________________________

%selecting mode 
%mode 1 :: calculating elastic convergence
%mode 2 :: calculating elastic + plastic strains
mode = 2; %Enter 'mode' value either 1 or 2 accordingly

%Defining Parameters for Variant #1
p_max = 50; % maximum pressure
E = 70000; %young modulus
nu = 0.25 ; % Poisson's ratio
yield_stress = 70; % yield stress


%Defining  final time, and delta time step (deltaT) 
t_final = 1;
dt = 0.05;
GP = 2;

%Defining inner and outer radius
r_in = 40 ; % inner radius in meter
r_out = 80;  % outer radius in meter

%Defining how many elements we want to have
elems = 8 ;% takes arbitrary value

%params = [E, Nu, yieldStress];

%_____________________________END OF INPUT________________________________

if mode == 1
    disp("You have chosen mode 1. The behavior of the model in the elastic regime is studied");
    [u, u_exact, r] = Elastic(r_in, r_out, E, nu, yield_stress, p_max, elems, t_final, dt, GP);
    f = figure('MenuBar', 'none', 'Name', 'Elastic Convergence Based on Number of Elements', 'NumberTitle', 'off', 'Position', [200, 200, 1000, 600]); %,'Toolbar','figure');
    defaultBackground = get(0,'defaultUicontrolBackgroundColor');
    set(f, 'Color', defaultBackground)
    
    plot([r(:, 1); r(end, 2)], u(1: (end))', [r(:, 1); r(end, 2)], u_exact(1:(end))', 'o')
    legend('u_{numerical}', 'u_{exact}')
    title('Converge Study in Elastic Regime')
    ylabel('u(r) in mm')
    xlabel('r in mm')
    
    xeText = uicontrol('Style', 'Text', 'String', 'No of Elements:', 'Position', [10,500,100,20], 'HorizontalAlignment', 'left');
    xeEdit = uicontrol('Style', 'Edit', 'String', '10', 'Position', [120, 500, 50, 20], 'HorizontalAlignment', 'center');
    xeSlider = uicontrol('Style', 'Slider', 'Position', [180, 500, 100, 20], 'CallBack', @xeSliderCallBack, 'Value', 10, 'Min', 0, 'Max', 500, 'SliderStep', [0.01 0.1]);

elseif mode == 2
    fprintf("You have chosen mode 2. The behavior of the model in both elastic and plastic regime are studied");
    [u, u_exact, u_history, r, stress_rr_history, stress_phiphi_history] = Elastoplastic(r_in, r_out, E, nu, yield_stress, p_max, elems, t_final, dt, GP);
    
    f1 = figure;
    plot((0:dt:t_final), u_history(end, :));
    title('Widening of the Outer Radius (b) over time');
    ylabel('u(r=b) mm');
    xlabel('time t');
    
    f2 = figure;
    plot(r(:,2), stress_rr_history(:))
    title('\sigma_{rr} distribution');
    ylabel('\sigma_{rr} MPa');
    xlabel('r  mm');
    
    f3 = figure;
    plot(r(:, 2), stress_phiphi_history(:));
    title('\sigma_{\phi \phi} distribution');
    ylabel('\sigma_{\phi \phi} MPa');
    xlabel('r in mm');
    
    f4 = figure;
    plot([r(:, 1); r(end, 2)], u_history(:, end));
    title('Displacement Distribution of u(r) at t=1');
    ylabel('u in mm');
    xlabel('r in mm');

else
    fprintf("WARNING: A wrong mode is chosen!. Please enter the mode value either 1 or 2");
end


toc() 


function xeSliderCallBack(varargin)
        num = get(xeSlider, 'Value');
        set(xeEdit, 'String', num2str(num));
end
