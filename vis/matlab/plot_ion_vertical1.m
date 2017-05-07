% Plot multiple species number densities as a function of Z.
% Written by Rui Xu. Oct. 2014

clear all;

% Network name, only dont need to include **-R1-nsp.dat part of the name.
Net_name = {'cplx'};
plot_species ='e- gr1(-) gr1(2-) gr1 gr1(+) gr1(2+)';
% length is total points, height is vertical height from midplane
length = 16;
height = 4;
z = linspace(0,height,length);
[sigma4,species_name4] = readfile(Net_name,length,1);


% Get the total abundance of Hydrogen
for k=1:length
   abundance_H4(k) = Find_element_abundance(species_name4,sigma4,k,'H');
end

% extract the species name.
select_name = regexp(plot_species,'\s+','split');
name_label2 = find_multi_name(species_name4, select_name);

full_name = regexp(species_name4,'\s+','split');


figure(1);
%plot electron Abn
semilogy(z,sigma4(:,name_label2(1))./abundance_H4(:),'r-','LineWidth',2);
sigma4(:,name_label2(1))
hold all;

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
semilogy(z,Num_ion(:)./abundance_H4(:), ...
           'b-','LineWidth',2);


%plot grs
semilogy(z,sigma4(:,name_label2(2))./abundance_H4(:),'-','LineWidth',2);  
semilogy(z,sigma4(:,name_label2(3))./abundance_H4(:),'-','LineWidth',2);
semilogy(z,sigma4(:,name_label2(4))./abundance_H4(:),'-','LineWidth',2);  
semilogy(z,sigma4(:,name_label2(5))./abundance_H4(:),'-','LineWidth',2);
semilogy(z,sigma4(:,name_label2(6))./abundance_H4(:),'-','Linewidth',2);
hold on;
xlabel('Z/H','fontsize',18);
ylabel('n/n_H','fontsize',18);
legend('e-', 'cation','gr(-)', 'gr(2-)', 'gr', 'gr(+)', 'gr(2+)')
%ylim([1e-22,1e-8])
set(gca,'fontsize',18);
box on;





