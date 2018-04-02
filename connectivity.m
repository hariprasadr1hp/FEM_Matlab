% connectivity.m
%
% This script is for computing the area connectivity matrix for
% each element in the system

% Define the element number and no.of nodes in the system

function Ae = connectivity(element,nodes)
Ae = zeros(2,nodes);
Ae(1,element) = 1;
Ae(2,element+1) = 1;
end