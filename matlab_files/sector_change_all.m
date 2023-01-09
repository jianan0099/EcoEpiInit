function sector_change_all(raw_epidemic_data_info, I_thre, Re_thre, phi, k, CHI_thre, rho, varphi, NPI_policy_scenario)
% sector change all
font_size = 11;
% ------- basic situation -------------------------
main_file_name = strcat('main_rho_', rho, '_varphi_', varphi);
base_scenario_key = strcat(rho, '_', varphi, '_', NPI_policy_scenario);
eco_epi_hyper_paras_info = strcat( 'I_thre_', I_thre,'_Re_thre_',Re_thre,'_phi_',phi,'_k_',k,'_CHI_thre_',CHI_thre);
% path
results_path = strcat('results/', raw_epidemic_data_info, '/', eco_epi_hyper_paras_info,  '/', main_file_name,'/full_market_change.xlsx');
heatmap_x_labels = table2cell(readtable(results_path, 'Sheet', 'heat_x_labels', 'ReadVariableNames', false));
heatmap_y_labels = table2cell(readtable(results_path, 'Sheet', 'heat_y_labels', 'ReadVariableNames', false));
figure('Position', [239,156,1086,745])

time_stamp_settings = {'_re', ''};
keep_year_settings = {'1.0', '1.5', '2.0'};
reopen_scenarios = {'Early full-reopening', 'Moderate full-reopening', 'Late full-reopening'};
time_stamps_info = {'At the end of the restriction period','At the end of the projection horizon'};

for time_stamp_index=1:2
    time_stamp = time_stamp_settings(time_stamp_index);
for keep_year_index=1:3
    keep_year = keep_year_settings(keep_year_index);
for part_num=1:2
    
figure_save_path = strcat('results/', raw_epidemic_data_info, '/', eco_epi_hyper_paras_info,'/', main_file_name, '/main_plot/full_market_change_results',...
    keep_year, time_stamp,'_',string(part_num), '.eps');

ax = subplot(1,1,1);
country_num_each_part = 50;
if part_num<=2
country_index_range = country_num_each_part*(part_num-1)+1:country_num_each_part*part_num;
else
    country_index_range = country_num_each_part*(part_num-1)+1:141;
end

heatmap_sheet_name = string(strcat(base_scenario_key, keep_year, time_stamp, '_heat'));
full_table = table2array(readtable(results_path, 'Sheet', heatmap_sheet_name)) * 100;
final_full_table = full_table(country_index_range, :);
disp(min(min(final_full_table)))
disp(max(max(final_full_table)))
imagesc(final_full_table, [-100, 100])
colormap(ax, cbrewer2('RdBu'))
xticks = linspace(1, 65, 65);
yticks = linspace(1, country_num_each_part, country_num_each_part);
set(gca, 'XTick', xticks, 'XTickLabel', heatmap_x_labels(1:65), 'TickLabelInterpreter','none')
xtickangle(90)
set(gca, 'YTick', yticks, 'YTickLabel', heatmap_y_labels(country_index_range), 'TickLabelInterpreter','none')
colorbar
clb = colorbar;
ylabel(clb, 'Relative change in market share (%)','FontSize',font_size, 'Rotation', 270,...
     'Units', 'normalized','Position', [4.35377779006958,0.500000476837158,0])
 
title({strcat(string(reopen_scenarios(keep_year_index)), ';', ' Part-', string(part_num)); string(time_stamps_info(time_stamp_index))});
saveas(gcf,figure_save_path,'epsc')
end
end
end
end