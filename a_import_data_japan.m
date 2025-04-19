
%% import and clean data from Japan

% import quarterly data relevant to estimate our model
[output_table, ~, T] = call_dbnomics( ...
    'OECD/QNA/JPN.B1_GE.CQRSA.Q', ...        % Real GDP, seasonally adjusted
    'OECD/QNA/JPN.P31S14_S15.CQRSA.Q', ...   % Household final consumption, seasonally adjusted
    'OECD/QNA/JPN.P51.CQRSA.Q', ...          % Gross fixed capital formation, seasonally adjusted
    'OECD/QNA/JPN.B1_GE.DNBSA.Q', ...        % GDP deflator, not seasonally adjusted
    'OECD/MEI/JPN.LRUNTTTT.STSA.Q', ...      % Harmonized unemployment rate, seasonally adjusted
    'OECD/KEI/IR3TIB01.JPN.ST.Q');           % 3-month interest rate


% variables : Time, Output, Consumption, Investment, Deflator, Unemployment, Nominal Rate

% select non-NaN observations
idx = find(~isnan(sum(output_table(:,2:end),2)));
output_table = output_table(idx,:);
T = T(idx);

% normalize prices to one in 2015
id2015 = find(T == 2015);
def = output_table(:,5) / output_table(id2015,5);

%% taking in real growth rates per capita
gy_obs  = diff(log(output_table(:,2)./(def)));
gc_obs  = diff(log(output_table(:,3)./(def)));
gi_obs  = diff(log(output_table(:,4)./(def)));

%%

%unemployment rate 
u_obs	= output_table(2:end,6);

% inflation rate
pi_obs  = diff(log(def));
% quarterly interest rate
r_obs	= output_table(2:end,7)/400;

T = T(2:end);

% save into myobs.mat
save myobs gy_obs gc_obs gi_obs u_obs T pi_obs r_obs;

%% plots
% First set of plots (3 subplots for growth and unemployment)
figure;
% First column of plots
subplot(3,2,1)
plot(T, gy_obs)  % GDP growth rate (log difference)
xlim([min(T) max(T)]);
title('GDP Growth (Log Difference)')

subplot(3,2,2)
plot(T, gc_obs)  % Consumption growth rate (log difference)
xlim([min(T) max(T)]);
title('Consumption Growth (Log Difference)')

% Second column of plots
subplot(3,2,3)
plot(T, gi_obs)  % Investment growth rate (log difference)
xlim([min(T) max(T)]);
title('Investment Growth (Log Difference)')

subplot(3,2,4)
plot(T, u_obs)  % Unemployment rate (no transformation)
xlim([min(T) max(T)]);
title('Unemployment Rate')

subplot(3,2,5)
plot(T, pi_obs)  % Inflation rate (log difference of deflator)
xlim([min(T) max(T)]);
title('Inflation Rate (Log Difference of Deflator)')

subplot(3,2,6)
plot(T, r_obs)  % Nominal interest rate (adjusted quarterly rate)
xlim([min(T) max(T)]);  % X-axis limits remain the same
ylim([0, 0.01]);  % Set Y-axis limits to be between 0% and 5% (i.e., 0 to 0.05)
title('Nominal Interest Rate')

% Save the figure as a PNG file
saveas(gcf, 'japan_economic_time_series.png');
