% regression variables norm data distribution
function reg_var_dis(raw_epidemic_data_info, I_thre, Re_thre, phi, k, CHI_thre, rho, varphi)
% ------- basic situation -------------------------
main_file_name = strcat('main_rho_', rho, '_varphi_', varphi);
% base_scenario_key = strcat(rho, '_', varphi, '_', NPI_policy_scenario);
eco_epi_hyper_paras_info = strcat( 'I_thre_', I_thre,'_Re_thre_',Re_thre,'_phi_',phi,'_k_',k,'_CHI_thre_',CHI_thre);
% path
results_path = strcat('results/', raw_epidemic_data_info, '/', eco_epi_hyper_paras_info,  '/', main_file_name,'/sector_regression_results.xlsx');
figure_save_path = strcat('results/', raw_epidemic_data_info, '/', eco_epi_hyper_paras_info,'/', main_file_name, '/main_plot/regression_var_dis.eps');
results_all = readtable(results_path, 'Sheet', 'norm_data', 'ReadVariableNames', true);

% --------- plot settings --------------------------
metrics = {'Impact_to_labourMultiplier__vartheta_P__','RiskDiversity',  ...
   'Impact_to_demandMultiplier__vartheta_D__',  'FractionOfSelf_suppliedDemand', ...
   'BetweennessCentrality', 'Degree', 'EigenvectorCentrality'};

texts = char(97:120);   
font_size = 11;
figure('Position', [91,132,846,428])
figs = 1;
for row=1:2
    if row==1
        col_max = 4;
    else
        col_max = 3;
    end
     for col=1:col_max
        subplot(2,4,4 * (row-1) + col ,'Position', [0.08+(col-1)*0.23, 0.62-(row-1)*0.5, 0.18, 0.28],'Units','normalized') %
        h_all = [];
        metric = string(metrics(figs));
        histogram(results_all.(metric),20,'Normalization','probability')
        
        if col==1
        ylabel({'Fraction of sectors'},  'Units', 'normalized','Position', [-0.26,0.53,0])
        end
        set(gca,'FontSize',font_size)
        text(-0.1, 1.17, texts(figs),  'Units', 'normalized', 'FontSize',font_size,'FontWeight','bold');
        

        
        if figs==1
            xlabel({'Impact-to-labour', 'multiplie'},'FontSize',font_size)
            ylim([0 0.7])
        end
                
        if figs==2
            xlabel({'Impact-to-demand', 'multiplie'},'FontSize',font_size)
        end
        
       if figs==3
            xlabel({'Risk diversity'},'FontSize',font_size)
        end
        
       if figs==4
            xlabel({'Fraction of self-supplied', 'demand'},'FontSize',font_size)
            ylim([0 0.12])
       end
        
       if figs==5
           xlabel({'Betweenness centrality'},'FontSize',font_size)
            ylim([0 0.5])
       end
        
       if figs==6
            xlabel({'Degree'},'FontSize',font_size)
       end
       
       if figs==7
            xlabel({'Eigenvector centrality'},'FontSize',font_size)
            ylim([0 0.5])
       end

        figs = figs +1;
    end
end
saveas(gcf,figure_save_path,'epsc')
end
