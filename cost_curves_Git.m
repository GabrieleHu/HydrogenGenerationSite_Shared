%% flags

plot_flag=0;
error_flag=0;


%% Heat exchanger cost curve

A = linspace(0.5,30);
c = log(A) - 0.6395*A.^2 + 947.2*A + 227.9;

N_bp2 = 3;
idx = round(linspace(1,100,N_bp2));
x_val  = A(idx)';
y_val  = c(idx)';
a_etae   = zeros(1,length(x_val)-1);
b_etae   = zeros(1,length(x_val)-1);
for i = 1 : length(x_val) - 1
    a_etae(i) = (y_val(i+1) - y_val(i)) / (x_val(i+1) - x_val(i));
    b_etae(i) = y_val(i+1) - a_etae(i) * x_val(i+1);
end
aa_HEX = a_etae;
bb_HEX = b_etae;

%% Heat exchanger cost curve - plot from Gabriele
 
if plot_flag==1
    
    A = linspace(0.5,30);
    c = log(A) - 0.6395*A.^2 + 947.2*A + 227.9;
    
    N_bp2 = 3;
    idx = round(linspace(1,length(A),N_bp2)); % Ensure linspace generates indices for A
    x_val  = A(idx)';
    y_val  = c(idx)';
    a_etae   = zeros(1,length(x_val)-1);
    b_etae   = zeros(1,length(x_val)-1);
    for i = 1 : length(x_val) - 1
        a_etae(i) = (y_val(i+1) - y_val(i)) / (x_val(i+1) - x_val(i));
        b_etae(i) = y_val(i+1) - a_etae(i) * x_val(i+1);
    end
    aa_HEX = a_etae;
    bb_HEX = b_etae;
    
    % Plotting
    figure;
    hold on;
    grid on;
    
    % Plot the non-linear cost curve
    plot(A, c, 'LineWidth', 2, 'DisplayName', 'Non-linear Cost Curve');
    
    % Plot the PWA approximation
    for i = 1:length(x_val)-1
        % Define the range for each linear segment
        x_range = linspace(x_val(i), x_val(i+1), 100);
        % Calculate the y values using the slope and intercept
        y_range = aa_HEX(i)*x_range + bb_HEX(i);
        % Plot the segment
        plot(x_range, y_range, '--', 'LineWidth', 2, 'DisplayName', ['PWA Segment ' num2str(i)]);
    end
    
    xlabel('Area (m^2)');
    ylabel('Cost');
    title('Heat Exchanger Cost Curve and PWA Approximation');
    legend('show');

end

%% error Estimation

if error_flag==1

    % Error calculation
    % Choose a dense set of points for error calculation to cover the entire range
    dense_A = linspace(min(A), max(A), 1000);
    dense_c = log(dense_A) - 0.6395*dense_A.^2 + 947.2*dense_A + 227.9; % Non-linear curve at dense points
    
    % Initialize the PWA approximation values at dense points
    pwa_approx = zeros(size(dense_A));
    
    % Calculate PWA values
    for i = 1:length(x_val)-1
        % Find indices within the current segment range
        idx_range = dense_A >= x_val(i) & dense_A <= x_val(i+1);
        % Calculate PWA approximation for this segment
        pwa_approx(idx_range) = aa_HEX(i)*dense_A(idx_range) + bb_HEX(i);
    end
    
    % Calculate errors
    errors = abs(dense_c - pwa_approx); % Absolute error
    MAE = mean(errors); % Mean Absolute Error
    RMSE = sqrt(mean(errors.^2)); % Root Mean Squared Error
    
    % Display error metrics
    disp(['Mean Absolute Error (MAE): ', num2str(MAE)]);
    disp(['Root Mean Squared Error (RMSE): ', num2str(RMSE)]);
    
    % Optionally, plot the error over the range
    figure;
    plot(dense_A, errors, 'LineWidth', 2);
    xlabel('Area (m^2)');
    ylabel('Absolute Error');
    title('Error between Non-linear Cost Curve and PWA Approximation');
    grid on;


    % Assuming dense_c contains the non-linear cost curve values
    average_cost = mean(dense_c);
    
    % Calculate percentage errors
    MAE_percent = (MAE / average_cost) * 100;
    RMSE_percent = (RMSE / average_cost) * 100;
    
    % Display percentage errors
    disp(['MAE as a percentage of average cost: ', num2str(MAE_percent), '%']);
    disp(['RMSE as a percentage of average cost: ', num2str(RMSE_percent), '%']);

end

%% Compressor cost curve

xx_c1    = 5;
xx_c2    = 12.5;
xx_c3    = 20;
% % Sensitivity analysis mass hydrogen:
% xx_c1    = 3;
% xx_c2    = 18;
% xx_c3    = 33;


yy_c1    = 43872 * (xx_c1)^0.5861;
yy_c2    = 43872 * (xx_c2)^0.5861;
yy_c3    = 43872 * (xx_c3)^0.5861;


aa_c  = [(yy_c2-yy_c1) / (xx_c2-xx_c1); (yy_c3-yy_c2) / (xx_c3-xx_c2)];
bb_c  = [yy_c1 - aa_c(1)*xx_c1; yy_c2 - aa_c(2)*xx_c2];

aa_c1 = (yy_c3-yy_c1) / (xx_c3-xx_c1);
bb_c1 = yy_c1 - aa_c1*xx_c1;

%% compressor cost curve - plot

S_c   = linspace(5,20);           % Compressor size                        [kW]
C_c   = 43872 * (S_c).^0.5861;    % Compressor cost curve                  [EUR/kW]

if plot_flag==1

    figure
    plot(S_c, C_c)
    hold on
    plot(S_c(1:67), aa_c(1)*S_c(1:67) + bb_c(1))
    hold on
    plot(S_c(67:end), aa_c(2)*S_c(67:end) + bb_c(2))
    % hold on
    % plot(S_c, aa_c1*S_c + bb_c1)
    legend('Actual curve', 'PWA pt1', 'PWA pt2')
    xlabel('Compressor size [kW]')
    ylabel('Compressor cost [EUR]')

end

%% Electrolyser cost curve: Reksten et al.

alpha = 0.622;
beta  = -158.9;
k_0   = 585.85;
k     = 9458.2;
V_0   = 2020;                          % Reference installation year
V     = 2023;                          % Installation year
SE    = 510;

xx_e1    = 250;
xx_e2    = 400;
xx_e3    = 550;
yy_e1    = (k_0 + k./xx_e1.*(xx_e1).^(alpha)) * (V/V_0)^(beta) .* xx_e1;
yy_e2    = (k_0 + k./xx_e2.*(xx_e2).^(alpha)) * (V/V_0)^(beta) .* xx_e2;
yy_e3    = (k_0 + k./xx_e3.*(xx_e3).^(alpha)) * (V/V_0)^(beta) .* xx_e3;

aa_e  = [(yy_e2-yy_e1) / (xx_e2-xx_e1); (yy_e3-yy_e2) / (xx_e3-xx_e2)];
bb_e  = [yy_e1 - aa_e(1)*xx_e1; yy_e2 - aa_e(2)*xx_e2];

% Sensitivity analysis mass hydrogen:
% xx_e1    = 100;
% xx_e2    = 250;
% xx_e3    = 400;
% xx_e4    = 550;
% xx_e5    = 850;
% yy_e1    = (k_0 + k./xx_e1.*(xx_e1).^(alpha)) * (V/V_0)^(beta) .* xx_e1;
% yy_e2    = (k_0 + k./xx_e2.*(xx_e2).^(alpha)) * (V/V_0)^(beta) .* xx_e2;
% yy_e3    = (k_0 + k./xx_e3.*(xx_e3).^(alpha)) * (V/V_0)^(beta) .* xx_e3;
% yy_e4    = (k_0 + k./xx_e4.*(xx_e4).^(alpha)) * (V/V_0)^(beta) .* xx_e4;
% yy_e5    = (k_0 + k./xx_e5.*(xx_e5).^(alpha)) * (V/V_0)^(beta) .* xx_e5;
% 
% yy_e1    = e_surplus * yy_e1;
% yy_e2    = e_surplus * yy_e2;
% yy_e3    = e_surplus * yy_e3;
% yy_e4    = e_surplus * yy_e4;
% yy_e5    = e_surplus * yy_e5;
% 
% aa_e  = [(yy_e2-yy_e1) / (xx_e2-xx_e1); (yy_e3-yy_e2) / (xx_e3-xx_e2); (yy_e4-yy_e3) / (xx_e4-xx_e3); (yy_e5-yy_e4) / (xx_e5-xx_e4)];
% bb_e  = [yy_e1 - aa_e(1)*xx_e1; yy_e2 - aa_e(2)*xx_e2; yy_e3 - aa_e(3)*xx_e3; yy_e4 - aa_e(4)*xx_e4]; 


% Plotting the cost curve/approximation and relative error of the
% approximation

% S_el   = linspace(0,2000,2000);         % Electrolyser size                        [kW]
% C_el   = (k_0 + k./S_el.*(S_el).^(alpha)) * (V/V_0)^(beta) .* S_el;  % Electrolyser cost

% figure
% plot(S_el, C_el/1000)
% hold on
% plot(S_el(281:401), (aa_e(1)*S_el(281:401) + bb_e(1))/1000, 'color',[0.8500 0.3250 0.0980])
% hold on
% plot(S_el(401:511), (aa_e(2)*S_el(401:511) + bb_e(2))/1000, 'color',[0.8500 0.3250 0.0980])
% xlim([260 530])
% legend('Cost curve', 'PWA')
% xlabel('Electrolyser size [kW]')
% ylabel('Electrolyser cost [kEUR/kW]')
% 
% 
% figure
% plot( S_el(281:401), (C_el(281:401)-(aa_e(1)*S_el(281:401) + bb_e(1)))./C_el(281:401)*100, 'color',[0.8500 0.3250 0.0980] )
% hold on
% plot( S_el(401:511), (C_el(401:511)-(aa_e(2)*S_el(401:511) + bb_e(2)))./C_el(401:511)*100, 'color',[0.8500 0.3250 0.0980] )
% xlim([260 530])
% % title('Relative error in electrolyser cost curve')
% xlabel('Electrolyser size [kW]')
% ylabel('Relative error [%]')
