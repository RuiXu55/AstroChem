% Plot multiple species number densities as a function of time
% comparison between full and reduced network
% Written by Rui Xu. Oct. 2014
clear all;
linestyles = cellstr(char('-','--','-.',':'));
colors = cellstr(char('k','r','b','g','c','y','m'));
% Network name, only dont need to include **-R1-nsp.dat part of the name.
% full network first, followed by reduced network
Net_name = {'temps', 'tempc'};
plot_species ='e- SO S H3+ H3O+';
radius = [1,10,100];


[sigma1,species_name1] = readfile1(Net_name{1},1,radius(1)); % full network
% extract the species name.
select_name = regexp(plot_species,'\s+','split');
ind = find_multi_name(species_name1, select_name);
len = size(select_name);
data1 = zeros(3,len(2));
data2 = zeros(3,len(2));

for x=1:3
    % start time/end time/# of data points
    [sigma1,species_name1] = readfile1(Net_name{1},1,radius(x)); % full network
    [sigma2,species_name2] = readfile1(Net_name{2},1,radius(x)); % reduced network
    
    % Get the total abundance of Hydrogen
    H1 = Find_element_abundance(species_name1,sigma1,1,'H');
    H2 = Find_element_abundance(species_name2,sigma2,1,'H');

    % extract the species name.
    select_name = regexp(plot_species,'\s+','split');
    ind1 = find_multi_name(species_name1, select_name);
    ind2 = find_multi_name(species_name2, select_name);
    len = size(select_name);
    
    for k=1:len(2)
          data1(x,k) = sigma1(:,ind1(k))/H1;
          data2(x,k) = sigma2(:,ind2(k))/H1;
    end
end


for k=1:len(2)
    loglog(radius,data1(:,k),[linestyles{1},colors{k}],'LineWidth',2.0,'DisplayName',select_name{k}); 
    hold on; 
end

xlabel('Radius');
ylabel('n_e/n_H');
set(gca,'fontsize',18);
legend('show')
for k=1:len(2)
    loglog(radius,data2(:,k),[linestyles{2},colors{k}],'LineWidth',2.0); 
end

% set(gca,'fontsize',28);
% set(gcf,'PaperPosition',[0,0,15,10]);
% set(gcf,'PaperSize',[15,10]);% %box on;
% box on;
% saveas(fig1,'evo.pdf','pdf');



%% relative error

% figure(2)
% len = size(select_name);
% for k=1:len(2)          
%     %fig=semilogx(time,100.*(sigma1(:,ind(k))-sigma2(:,ind(k)))./sigma1(:,ind(k)),[linestyles{1},colors{k}],'LineWidth',2.0);  
%     fig=semilogx(time,sigma1(:,ind(k)),[linestyles{1},colors{k}],'LineWidth',2.0);  
%     hold on;
% end
% 
% xlabel('Time/Year');
% ylabel('Relative Error %');
% axis([ts,te,0,40])
% set(gca,'fontsize',28);
% set(gcf,'PaperPosition',[0,0,15,10]);
% set(gcf,'PaperSize',[15,10]);% %box on;
% %saveas(fig,'plot.pdf','pdf');




