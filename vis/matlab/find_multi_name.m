function [ output] = find_multi_name(species_name, name )
% This function is a part of Species Number Density Analysis Code.
% Used to find the index of all selected species name ( defined as 
% species_name) in full chemistry network (name).
% Written by Rui Xu. Oct. 2014



% total number of species in chemistry network
num = length(name);
%If the species is in network, error_label =1, otherwise error_label=0;
error_label = zeros(num);
% instore the species index
name_label  = zeros(num);
k=1;
for j=1:num
    for i=1:length(species_name)
        if strcmp(species_name(i),name(j)) ==1
            error_label(k) =1;
            name_label(k) =i;
            k = k+1;
        break;
        end
    end   
end

for i =1:num
% If the species not in the network, name_label set to be negative.
    if error_label(i,1) < 1
        name_label =-1;
    end
end
output = name_label;
end

