function risk_diversity_sample(raw_epidemic_data_info, I_thre, Re_thre, phi, k, CHI_thre)
% risk diversity sample
% ------- basic situation -------------------------

eco_epi_hyper_paras_info = strcat( 'I_thre_', I_thre,'_Re_thre_',Re_thre,'_phi_',phi,'_k_',k,'_CHI_thre_',CHI_thre);
results_path = strcat('results/', raw_epidemic_data_info, '/', eco_epi_hyper_paras_info, '/risk_diversity_sample.xlsx');
figure_save_path = strcat('results/', raw_epidemic_data_info, '/', eco_epi_hyper_paras_info,'/risk_diversity_sample.eps');

% ---------- results ------------------
risk_diversity_data = table2array(readtable(results_path,'Sheet', 'risk','PreserveVariableNames',true));
production_data = table2array(readtable(results_path,'Sheet','production','PreserveVariableNames',true));
CHI_data = table2array(readtable(results_path,'Sheet', 'CHI','PreserveVariableNames',true));

% --------- draw -----------------------
font_size = 11;
figure('Position', [367,342,701,566])

% ---- production distribution --
for col=2:4
    production = production_data(col-1, :);
    ax = subplot(4,4, col, 'Position', [0.338545541164731+ 0.21*(col-2),0.8,0.19,0.153714522865626]);
    bar(production, 'BarWidth', 0.4)
    ylim([0 1]);
    if col==2
    ylabel('Market share');
        else
        set(gca,'ytick',[])
    end
    set(gca, 'XTick', linspace(1, 2, 2), 'XTickLabel', {'Firm 1', 'Firm 2'},'FontSize',font_size,'XTickLabelRotation',0)
end
% ---- CHI distribution --
for row=2:4
    CHI = CHI_data(row-1, :);
    ax = subplot(4,4, 4 * (row-1) + 1, 'Position', [0.065381444281858,0.54 - 0.22*(row-2),0.2,0.153714522865626]);
    bar(CHI, 'BarWidth', 0.4)
    ylim([0 12]);
    ylabel('CHI');
    if row==4
    set(gca, 'XTick', linspace(1, 2, 2), 'XTickLabel', {'Firm 1', 'Firm 2'},'FontSize',font_size,'XTickLabelRotation',0)
    else
        set(gca,'xtick',[])
    end
    set(gca,'FontSize',font_size)
end
% ------ risk diversity --
ax = subplot(4,4, 4 + 2, 'Units', 'normalized',  'Position', [0.336117021276596,0.095406360424028,0.626793107111421,0.610507618070595]);
%risk_diversity_data = flip(risk_diversity_data);
heatmap(risk_diversity_data)
colorbar off
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
set(gca,'FontSize',font_size)
saveas(gcf,figure_save_path,'epsc')    
end



