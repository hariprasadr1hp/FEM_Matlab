function p = pl_init(yieldStress,inner_radius,outer_radius,Nu)

p = yieldStress*(((outer_radius^2)-(inner_radius^2))/(inner_radius^2))*(1/sqrt(1-(4*Nu)+(4*(Nu^2))+3*((outer_radius^4)/(inner_radius^4)))); %pressure of plasticity begin

end