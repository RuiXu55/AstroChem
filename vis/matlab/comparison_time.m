% Plot multiple species number densities as a function of time
% comparison between full and reduced network
% Written by Rui Xu. Oct. 2014
clear all;
linestyles = cellstr(char('-','--','-.',':'));
colors = cellstr(char('k','r','b','g'));
% Network name, only dont need to include **-R1-nsp.dat part of the name.
% full network first, followed by reduced network
Net_name = {'temps', 'tempc'};
plot_species ='e-';
radius = [1,10,100];
ts = 1e0;  te = 1e6;  nt = 1;
for i=1:nt
   time(i) = ts*exp((i-1)/(nt-1)*log(te/ts));
end

for x=1:3
    % start time/end time/# of data points
    [sigma1,species_name1] = readfile1(Net_name{1},nt,radius(x)); % full network
    [sigma2,species_name2] = readfile1(Net_name{2},nt,radius(x)); % reduced network
    % Get the total abundance of Hydrogen
    for k=1:nt
      H1(k) = Find_element_abundance(species_name1,sigma1,k,'H');
      H2(k) = Find_element_abundance(species_name2,sigma2,k,'H');
    end
    % extract the species name.
    c
        if x==1
          h1 =loglog(time,sigma1(:,ind(k))/H1(k),[linestyles{1},colors{x}],'LineWidth',2.0); 
        elseif x==2
          h2 =loglog(time,sigma1(:,ind(k))/H1(k),[linestyles{1},colors{x}],'LineWidth',2.0);
        else
          h3 =loglog(time,sigma1(:,ind(k))/H1(k),[linestyles{1},colors{x}],'LineWidth',2.0);
        end
        hold on;
  
        fig1 =loglog(time,sigma2(:,ind(k))/H2(k),[linestyles{2},colors{x}],'LineWidth',2.0);
    end
end

xlabel('Time/Year');
ylabel('n_e/n_H');
set(gca,'fontsize',18);
legend([h1,h2,h3],'1AU','10AU','100AU');


set(gca,'fontsize',28);
set(gcf,'PaperPosition',[0,0,15,10]);
set(gcf,'PaperSize',[15,10]);% %box on;
box on;
saveas(fig1,'evo.pdf','pdf');



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




