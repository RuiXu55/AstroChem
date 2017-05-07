function [ output] = find_name(species_name, find_name)
% Find species index from full chemistry network
% Written by Rui Xu. Oct. 2014.

error_label =-1;
for i=1:length(species_name)
    if strcmp(species_name(i),find_name) ==1
        error_label =1;
        str_num =i;
        break;
    end
end

if error_label <0
    output = -1;
    return; 
end
output = str_num;
end

