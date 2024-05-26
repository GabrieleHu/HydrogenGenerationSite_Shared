function [] = SelectedWeek_WHR(linew,font,Time,start,finish,P_PV_opt,P_e_opt,P_imp_opt,P_exp_opt,P_th_HT,P_th_LT, ElyOn, P_th,P_th_av)

year=365*2023+126;
start_date=start/24+year;
finish_date=finish/24+year;
tt=datetime([start_date:(1/24):finish_date], 'ConvertFrom', 'datenum');

P_th_loss = P_th - P_th_av;
P_th_amb  = P_th_av - P_th_LT - P_th_HT;

P_TH = [P_th_loss(start:finish) P_th_amb(start:finish) P_th_LT(start:finish) P_th_HT(start:finish)];
ncat = size(P_TH,2);

figure
ar = area(tt, P_TH,'LineWidth',linew);
hold on 

% get specific colour scheme
colors = [1 44 86 129 172 214];
cmap = crameri('batlow');

for i = 1:ncat
    ar(i).FaceColor = cmap(colors(i+1), :);
end

xlim([min(tt) max(tt)])
% ylim([-0.2*max([max(P_th_HT) max(P_th_LT)]) 1.2*max([max(P_th_HT) max(P_th_LT)]) ])
% title('Power from PV and PEMWE');
ylabel ('Power [kW]','fontweight','bold');
xlabel('Time','fontweight','bold');
legend('Q_{loss}','Q_{amb}','Q_{HEX}','Q_{HP}','location','eastoutside')


set(gcf, 'Units', 'inches');
% set(gcf, 'Position', [0, 0, 3.5, 2.5]); % Width=3.5in, Height=2.5in
set(gca, 'FontSize', font); % Set axis font size and font name
set(findall(gcf, 'Type', 'line'), 'LineWidth', 1); % Set line width



end