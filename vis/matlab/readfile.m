function [ data,name ] = readfile(Net_name,Num,file_tag)
% This function read the data from Data file 
% Net_name is a string identify the different network.
% Num for the total iteration in a network
% file_tag shows the radii of the disk
% Written by Rui Xu. Oct. 2014

init_label=1.0; % INIT_LABEL IS GANRENTEE JUST READ FROM THE FIRST FILE.
	            % NO USE IN THIS FUNCTION
                
%%%%%%% Read data from the file %%%%%%%%%%%%%%
Radius = num2str(file_tag); % str is the radius unit in AU
filename  = strcat('./Data/',Net_name,'-R',Radius,'.nsp-0.dat')
% change cell to string
filename  = filename{1}
fid   = fopen(filename,'r');% open certain file

% three lines of comments
tline =	fgetl(fid);	% zeta or rho 
tline =	fgetl(fid);	% comment
tline =	fgetl(fid);	% species in total


str3         = regexp(tline,'\s+','split');
species_num  = str2double(str3(2));	%total number of species

% init the value matrix sigma.
sigma 	= ones(Num,species_num);

tline = fgetl(fid);	%get the line of species
str3  = regexp(tline,'\s+','split');
% n_s_line is the total line of species
n_s_line        = str2double(str3(5));

% store all species_name in species_name
name_ind	= 1.0;
for j=1:n_s_line
    tline	= fgetl(fid);
    str4  	= regexp(tline,'\s+','split');	%str4 instore the name of species
    % analyse the name of the species
    for k=2:length(str4)-1
        species_name(name_ind)=str4(k);
        name_ind=name_ind+1;
    end
end

% store the data in matrix sigma
for i=1:Num  % Num is the total iteration times
    iter  = 1;
    tline = fgetl(fid);           %customised value, not included in our species.
    while iter<species_num+1      % till the end of the file    
            tline	=fgetl(fid);
            str3    =regexp(tline,'\s+','split');
        for k=1:length(str3)-1
                sigma(i,iter)=str2double(str3(k));
                iter=iter+1;
        end
    end
end
sta	= fclose(fid);
data = sigma;
name = species_name;
end

