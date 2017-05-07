% Plot multiple species number densities as a function of Z.
% Written by Rui Xu. Oct. 2014

clear all;

% Network name, only dont need to include **-R1-nsp.dat part of the name.
Net_name = {'cplx'};
plot_species ='H2O CO NH3 H2';
plot_species = 'S S+ SO SO+';
%plot_species = 'e-  gr1(2-) gr1(-) gr1 gr1(+)';
% length is total points, height is vertical height from midplane
length = 100;
height = 4;
r = 1;

fig = fopen(strcat('cplx',int2str(r)),'w');
z = linspace(0,height,length);
[sigma4,species_name4] = readfile(Net_name,length,r);

% Get the total abundance of Hydrogen
for k=1:length
   abundance_H4(k) = Find_element_abundance(species_name4,sigma4,k,'H');
end

% extract the species name.
select_name = regexp(plot_species,'\s+','split');
name_label2 = find_multi_name(species_name4, select_name);
len2 = size(name_label2)
full_name = regexp(species_name4,'\s+','split');


%figure(1);
%plot electron Abn
%semilogy(z,sigma4(:,name_label2(1))./abundance_H4(:),'r--','LineWidth',2);
%hold all;

% calculate total ion abundance (grain max charge is set to be 2)
num1=size(full_name);
num = num1(2);

Num_ion(1:length) = 0.0;
for i=1:num
    index1 = strfind(species_name4{i},'+');
    index2 = strfind(species_name4{i},'gr1(+)');
    index3 = strfind(species_name4{i},'gr1(2+)');
    index4 = strfind(species_name4{i},'gr1(-)');
    index5 = strfind(species_name4{i},'gr1(2-)');

    if index2 >0  | index3 >0 | index4 >0 | index5>0
        % The species is grains, left blank for other use
    elseif index1 >0
        for j =1:length
            Num_ion(j) = Num_ion(j) + sigma4(j,i);
        end
    end
end
%plot ions
% semilogy(z,Num_ion(:)./abundance_H4(:), ...
%            'b--','LineWidth',2);
% %plot grs
% semilogy(z,sigma4(:,name_label2(2))./abundance_H4(:),'--','color',[0,1.0,0],'LineWidth',1.5);  
% semilogy(z,sigma4(:,name_label2(3))./abundance_H4(:),'--','color',[0,0.5,0],'LineWidth',1.5);
% hold on;
% xlabel('Z/H','fontsize',18);
% ylabel('n/n_H','fontsize',18);
% 
% set(gca,'fontsize',18);
% box on;

if len2(1)==1
     fprintf(fig,'%7.20f %7.20f \n',[abundance_H4(:),sigma4(:,name_label2(1))]');
else
    size1 = size(abundance_H4);
    for x =1:size1(2)
        fprintf(fig,'%7.20f ',abundance_H4(x));
        for i =1:len2(1)
            fprintf(fig,'%7.20f  ',sigma4(x,name_label2(i)));
        end
        fprintf(fig,'\n');
    end
end
fclose(fig);


