function sector_change_modify(raw_epidemic_data_info, I_thre, Re_thre, phi, k, CHI_thre, rho, ...
    varphi, NPI_policy_scenario, if_base, supple_name)

% ------- basic situation -------------------------
if strcmp(NPI_policy_scenario, 'keep_curr_')
    main_file_name = strcat('main_rho_', rho, '_varphi_', varphi);
end
if strcmp(NPI_policy_scenario, 'linear_decrease_')
    main_file_name = strcat('linear_', 'main_rho_', rho, '_varphi_', varphi);
end
base_scenario_key = strcat(rho, '_', varphi, '_', NPI_policy_scenario);
eco_epi_hyper_paras_info = strcat( 'I_thre_', I_thre,'_Re_thre_',Re_thre,'_phi_',phi,'_k_',k,'_CHI_thre_',CHI_thre);
% path
results_path = strcat('results/', raw_epidemic_data_info, '/', eco_epi_hyper_paras_info,  '/', main_file_name,'/market_change.xlsx');
if if_base
figure_save_path = strcat('results/', raw_epidemic_data_info, '/', eco_epi_hyper_paras_info,'/', main_file_name, '/main_plot/market_change_main.eps');
else
    figure_save_path = strcat('results/', raw_epidemic_data_info, '/', eco_epi_hyper_paras_info,'/', main_file_name, '/main_plot/', supple_name,'market_change_main.eps');
end
% --------- plot settings --------------------------
row_info = {'', '_re', ''};
col_info = {'1.0', '1.5', '2.0'};
titles = {'Early full-reopening', 'Moderate full-reopening', 'Late full-reopening'};
texts = char(97:120);   
font_size = 11;
figure('Position', [91,219,789,741])

% ---- set basic information -------
heatmap_min_max = table2array(readtable(results_path, 'Sheet', 'heat_min_max')) * 100;
heatmap_min = -60;
heatmap_max = 60;

heatmap_x_labels = table2cell(readtable(results_path, 'Sheet', 'heat_x_labels', 'ReadVariableNames', false));
heatmap_y_labels = table2cell(readtable(results_path, 'Sheet', 'heat_y_labels', 'ReadVariableNames', false));


% ---- general -------
text(0.356791962174941,0.967698519515478, 'market share',  'Units', 'normalized', 'FontSize',font_size,'FontWeight','bold');
figs = 1;

for row=1:1
    for col=1:3
        ax = subplot(3,3,figs,'Units','normalized','Position', [0.0712+(col-1)*0.285, 0.785, 0.24, 0.163]); 
        result_info = string(strcat(base_scenario_key, col_info(col), '_mar'));
        sector_structure_change = table2array(readtable(results_path, 'Sheet', result_info))*100;
        bar(sector_structure_change)
        if col==1
            ylabel('$\Delta \mathcal{M}\,(\times10^2)$','interpreter','latex', 'Units', 'normalized','Position', [-0.1737,0.5,-1])
        end
        ylim([0, 1.2])
        xticks = linspace(1, 10, 10);
        set(gca, 'XTick', xticks, 'XTickLabel', heatmap_x_labels,'FontSize',font_size, 'XTickLabelRotation',90)
        title(titles(col),'FontSize',font_size,'Units', 'normalized', 'Position', [0.5,1.1,0])
        set(gca,'FontSize',font_size)
        text(-0.1, 1.17, texts(figs),  'Units', 'normalized', 'FontSize',font_size,'FontWeight','bold');
        figs = figs +1;
    end
end

for row=2:3
    for col=1:3
        ax = subplot(3,3,figs );
        set(gca,'FontSize',font_size,'Position', [0.0712+(col-1)*0.285, 0.433-(row-2)*0.358, 0.24, 0.22])
        heatmap_sheet_name = string(strcat(base_scenario_key, col_info(col), row_info(row), '_heat'));
        C = table2array(readtable(results_path, 'Sheet', heatmap_sheet_name)) * 100;
        disp(min(min(C)))
        disp(max(max(C)))
        imagesc(C, [heatmap_min heatmap_max])

        xticks = linspace(1, 10, 10);
        yticks = linspace(1, 10, 10);
        set(gca, 'XTick', xticks, 'XTickLabel', heatmap_x_labels,'FontSize',font_size)
        xtickangle(90)
        set(gca, 'YTick', yticks, 'YTickLabel', heatmap_y_labels,'FontSize',font_size)

        colormap(ax, cbrewer2('RdBu'))

        if col>2
            colorbar
            clb = colorbar;
            ylabel(clb, {'Relative change in', 'market share (%)'},'FontSize',font_size, 'Rotation', 270, ...
                'Units', 'normalized','Position', [5.8,0.49,0])
            clb.Position = [0.8937, 0.433-(row-2)*0.358,0.017494089934515,0.22];
        end

        if col>1
        set(gca,'YTickLabel',[]);
        end

        %title(titles(col),'FontSize',font_size)

        text(-0.1, 1.1, texts(figs),  'Units', 'normalized', 'FontSize',font_size,'FontWeight','bold');
        figs = figs +1;
        
        if col==1
            if row==2
            text(1.03, 1.22,'At the end of the restriction period',  'Units', 'normalized', 'FontSize',font_size,'FontWeight','bold');
            else
            text(1.03, 1.22,'At the end of the projection horizon',  'Units', 'normalized', 'FontSize',font_size,'FontWeight','bold');
            end
        end
    end
end


saveas(gcf,figure_save_path,'epsc')   
end
