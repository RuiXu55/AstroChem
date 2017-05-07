function [ abundance ] = Find_element_abundance(species_name,sigma,k, find_element)
%Calculate the abundance of a specific species (find_element) at k-th iteration.
% Written by Rui Xu. Oct. 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_abundance = 0;       
num_label       = 0;         % indicate the position of the species
num             = 0;         % the num of the element in a specific species.
for i=1:length(species_name)
    %  e.g. H can be found in both H20 and H3O+.
    position = strfind(species_name{i},find_element);
    if length(position) < 1 
        continue;       % no such kind of element being found
    end
    
    current_species = species_name{i};   % get the name of current species
    %the total times of an element exist in species, e.g. CH3CH2, both C
    %and H show up twice.
    for j=1:length(position)           
        % The element only has one caracter, like H
        if length(find_element)==1
            if position(j)<length(current_species)
                if  current_species(position(j)+1)>'a' && current_species(position(j)+1)<'z'
                    continue;   % same first character with other species. like H and He
                                % otherwise it should followed by number or
                                % capital letter.
                end
                if current_species(position(j)+1)>'1' && current_species(position(j)+1)<'9'
                    num = str2double(current_species(position(j)+1));
                    % How many times the element has here. e.g. CH3, triple
                    % H here.
                else
                    num = 1; % The element is not followed by a number.
                end
            else  % The element is the last letter of the species.
                num = 1;
            end
        else  % The element contains two  caracter, e.g. Si
            if position(j)<length(current_species)-(length(current_species)-1)
                if current_species(position(j)+length(current_species))>'1' && current_species(position(j)+length(current_species))<'9'
                    num = str2double(current_species(position(j)+length(current_species)));
                else
                    num =1;
                end
            else
                num = 1;
            
            end
        end
        total_abundance = total_abundance+num*sigma(k,i);
        num_label = num_label+1;
    end
end
abundance=total_abundance;
