function market_structure_change_fit(raw_epidemic_data_info, I_thre, Re_thre, phi, k, CHI_thre, rho, varphi, NPI_policy_scenario)

main_file_name = strcat('main_rho_', rho, '_varphi_', varphi);
base_scenario_key = strcat(rho, '_', varphi, '_', NPI_policy_scenario);
eco_epi_hyper_paras_info = strcat( 'I_thre_', I_thre,'_Re_thre_',Re_thre,'_phi_',phi,'_k_',k,'_CHI_thre_',CHI_thre);
sector_structure_change_results_path = strcat('results/', raw_epidemic_data_info, '/', eco_epi_hyper_paras_info,  '/', main_file_name,'/sector_regression_results.xlsx');
sector_structure_change = readtable(sector_structure_change_results_path, 'Sheet', 'sector_structure_change', 'PreserveVariableNames',true);
% log_shapiro = readtable(sector_structure_change_results_path, 'Sheet', 'sector_stru_dis_shapiro', 'PreserveVariableNames',true);
log_likelihood = readtable(sector_structure_change_results_path, 'Sheet', 'sector_stru_log_norm', 'PreserveVariableNames',true);
figure_save_path = strcat('results/', raw_epidemic_data_info, '/', eco_epi_hyper_paras_info,'/', main_file_name, '/main_plot/market_change_fit.eps');

col_info = {'1.0', '1.5', '2.0'};
font_size = 11;
titles = {'Early full-reopening', 'Moderate full-reopening', 'Late full-reopening'};
figs = 1;
figure('Position', [596,310,789,253])
texts = char(97:120);   
for col=1:3
    ax = subplot(1,3,figs ,'Position', [0.1+(col-1)*0.3, 0.2, 0.25, 0.65],'Units','normalized'); 
    result_info = string(strcat(base_scenario_key, col_info(col), '_original'));
    result = sector_structure_change.(result_info); 
    histogram(result,20,'Normalization','probability')

   ylim([0 0.9])


    if col==1
    xlabel('$\Delta \mathcal{M}$','interpreter','latex')
    end
    if col==2
    xlabel('$\Delta \mathcal{M}$','interpreter','latex')
    end
    if col==3
    xlabel('$\Delta \mathcal{M}$','interpreter','latex')
    end


    if col>1
    set(gca,'YTickLabel',[]);
    end

    if col==1
        ylabel('Distribution', 'Units', 'normalized','Position', [-0.22,0.5,-1])
    end


    %R_info = string(strcat(base_scenario_key, col_info(col), '_R'));
    ks_p_info = string(strcat(base_scenario_key, col_info(col), '_p_log'));
    %log_R = log_likelihood.(R_info);
    ks_p = string(log_likelihood.(ks_p_info));
    %text(0.55, 0.8, strcat('$\mathcal{R} = ', log_R, '$'), 'interpreter','latex', 'Units', 'normalized', 'FontSize',font_size);
    text(0.55, 0.65, strcat('$P =\,', ks_p, '$'), 'interpreter','latex','Units', 'normalized', 'FontSize',font_size);

    title(titles(col),'FontSize',font_size)

%     if col==1
%         text(0.9, 1.36,'Distribution of market structure changes',  'Units', 'normalized', 'FontSize',font_size,'FontWeight','bold');
%      end
    set(gca,'FontSize',font_size)
    text(-0.1, 1.17, texts(figs),  'Units', 'normalized', 'FontSize',font_size,'FontWeight','bold');
    figs = figs +1;
end
saveas(gcf,figure_save_path,'epsc')  
end