function [] = SelectedWeek_SOC(linew,font,Time,start,finish,P_PV_opt,P_e_opt,P_imp_opt,P_exp_opt,P_e_nom,SOC_opt,P_b_ch_opt,P_b_disch_opt)

year=365*2023+126;
start_date=start/24+year;
finish_date=finish/24+year;
tt=datetime([start_date:(1/24):finish_date], 'ConvertFrom', 'datenum');

colors = [1 44 86 129 172 214];
cmap = crameri('batlow');

figure
plot(tt,P_PV_opt(start:finish),'LineWidth',linew,'Color',cmap(colors(2), :))
hold on
plot(tt,-P_e_opt(start:finish),'LineWidth',linew,'Color',cmap(colors(3), :))
hold on
plot(tt,-P_imp_opt(start:finish),'LineWidth',linew,'Color',cmap(colors(4), :))
hold on 
plot(tt,P_exp_opt(start:finish),'LineWidth',linew,'Color',cmap(colors(5), :))
hold on

xlim([min(tt) max(tt)])
% ylim([-1.2*max([max(P_imp_opt) max(P_e_opt)]) 1.2*max(P_PV_opt) ])
% title('Power from PV and PEMWE');
ylabel ('Power [kW]','fontweight','bold');
xlabel('Time','fontweight','bold');


% yyaxis right
% % plot(Time(start:finish),SOC_opt(start:finish),'-.','LineWidth',linew);
% plot(tt,SOC_opt(start:finish),'-.','LineWidth',linew);
% hold on 
% plot([min(tt) max(tt)],[0 0],'-','color', [.5 .5 .5],'LineWidth',linew);
% ylabel ('SOC [-]','fontweight','bold');
% rate=1.2*(max([max(P_imp_opt) max(P_e_opt)])/(max(P_PV_opt)));
% ylim([-rate 1.2])
% yticks([0:0.25:1])
% set(gca, 'YGrid', 'off', 'XGrid', 'on')

legend('P_{PV}','P_{PEM}','P_{imp}','P_{exp}','location','eastoutside')

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, 7, 4.5]); % Width=3.5in, Height=2.5in
set(gca, 'FontSize', font); % Set axis font size and font name
set(findall(gcf, 'Type', 'line'), 'LineWidth', 1); % Set line width

end