function [] = SelectedWeek_WHR_ext(linew,font,Time,start,finish,P_PV_opt,P_e_opt,P_imp_opt,P_exp_opt,P_th_HT,P_th_LT, ElyOn)

year=365*2023+126;
start_date=start/24+year;
finish_date=finish/24+year;
tt=datetime([start_date:(1/24):finish_date], 'ConvertFrom', 'datenum');

figure
plot(tt,P_PV_opt(start:finish),'LineWidth',linew)
hold on
plot(tt,P_imp_opt(start:finish),'LineWidth',linew)
hold on
area(tt,P_e_opt(start:finish),'FaceAlpha',.2)
% hold on
% area(tt,P_th_LT(start:finish),'LineWidth',linew)
hold on 
area(tt,P_th_HT(start:finish),'FaceColor',[0 100/255 0],'FaceAlpha',.3)
hold on 

xlim([min(tt) max(tt)])
% ylim([0 1500])
% title('Power from PV and PEMWE');
ylabel ('Power [kW]','fontweight','bold');
xlabel('Time','fontweight','bold');


% hold on
% yyaxis right
% plot(tt,T_out(start:finish),'-.','LineWidth',linew);
% hold on 
% plot([min(tt) max(tt)],[0 0],'-','color', [.7 .7 .7],'LineWidth',linew);

% legend('P_{PV}','P_{imp}','P_{PEM}','Q_{HEX}','Q_{HP}','location','eastoutside')
legend('P_{PV}','P_{imp}','P_{PEM}','Q_{HP}','location','eastoutside')

set(gcf, 'Units', 'inches');
% set(gcf, 'Position', [0, 0, 3.5, 2.5]); % Width=3.5in, Height=2.5in
set(gca, 'FontSize', font); % Set axis font size and font name
set(findall(gcf, 'Type', 'line'), 'LineWidth', 1); % Set line width


end