function com_explain(raw_epidemic_data_info, I_thre, Re_thre, phi, k, CHI_thre, rho, ...
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
results_path = strcat('results/', raw_epidemic_data_info, '/', eco_epi_hyper_paras_info, '/', main_file_name, '/competition_explain.xlsx');
if if_base
    figure_save_path = strcat('results/', raw_epidemic_data_info, '/', eco_epi_hyper_paras_info,'/', main_file_name, '/main_plot/competition_explain.eps');
else
    figure_save_path = strcat('results/', raw_epidemic_data_info, '/', eco_epi_hyper_paras_info,'/', main_file_name, '/main_plot/', supple_name, 'competition_explain.eps');
end
results = readtable(results_path,'Sheet', 'results','PreserveVariableNames',true);
init_x_info = readtable(results_path,'Sheet', 'init_x_info','PreserveVariableNames',true);
CHI_scores = readtable(results_path,'Sheet', 'CHI_scores','PreserveVariableNames',true);
CHI_levels = readtable(results_path,'Sheet', 'CHI_levels','PreserveVariableNames',true);
vartheta = readtable(results_path,'Sheet', 'vartheta','PreserveVariableNames',true);
largest_producer_info = readtable(results_path,'Sheet', 'large_reg','PreserveVariableNames',true);

% --- figs info ------
sectors = {'eeq','pdr', 'edu'};
metrics = {'re','var_pro','ka'}; 
ylabels = {'Fraction',  '\Theta', '\Delta\xi','\kappa'};
baseline_values = [1, 0, 1];

year_info = {'1.0', '1.5', '2.0'};
keep_weeks_info = [52*1, 52*1.5, 52*2];
specific_year = 2;


point_length = 2; % plot every [point_length] points
texts = char(97:120);
marker_size = 2;
line_width = 1.5;
font_size = 11;
line_color_all = [%8/255 29/255 88/255;
                37/255 52/255 148/255;
                34/255,94/255,168/255;
                65/255,182/255,196/255;
                127/255,205/255,187/255;
                199/255,233/255,180/255;];
face_color_all = [251/255 200/255 196/255;
            255/255,223/255,191/255;
            204/255,242/255,217/255;];
y_lim_min = [0, -1, 0.5];
y_lim_max = [1.65, 1.7, 1.5];
          
% ---- main ---------------------------
figure('Position', [91,132,846,763])
figs = 1;

for row=1:1
     for col=1:3
         
         sector = string(sectors(col));
         vartheta_sec = vartheta.(sector);
         CHI_score = CHI_scores.(sector);
         CHI_level = CHI_levels.(sector);        
         largest_producer =largest_producer_info.(sector);
         
         subplot(4,3,figs,'Units','normalized','Position', [0.1+(col-1)*0.32, 0.76-(row-1)*0.23, 0.25, 0.16],'Units','normalized')
         for k = 1:length(CHI_score)   
            b=bar(k,CHI_score(k));
            if strcmp(string(CHI_level(k)),'stringent_npi')
                set(b,'FaceColor',face_color_all(1,:));
            end
            if strcmp(string(CHI_level(k)),'moderate_npi')
                set(b,'FaceColor',face_color_all(2,:));
            end
            if strcmp(string(CHI_level(k)),'mild_npi')
                set(b,'FaceColor',face_color_all(3,:));
            end
            hold on
         end
         ylim([0 80])
         %b = bar(init_x);
         xticks = linspace(1, 5, 5);
         text(-0.2, 1.2, texts(figs), 'Units', 'Normalized','FontSize',font_size,'FontWeight','bold');
         
         if col==1
           ylabel({'CHI score'}, 'Units', 'normalized', 'Position', [-0.18, 0.5, 0]);
         end
         
        if col==1
            % strcat(' (\vartheta^P=', string(vartheta_sec(1)), ', \vartheta^D=',string(vartheta_sec(2)),')')
             title({'Manufacture of', 'electrical equipment'},'Units', 'normalized', 'Position', [0.5, 1.02, 0]); % 
         end
         
         if col==2
             % strcat(' (\vartheta^P=', string(vartheta_sec(1)), ', \vartheta^D=',string(vartheta_sec(2)),')')
             title({'Rice'},'Units', 'normalized' , 'Position', [0.5, 1.02, 0]); %
         end
         
         if col==3
             % strcat(' (\vartheta^P=',string(vartheta_sec(1)), ', \vartheta^D=',string(vartheta_sec(2)),')')
             title({'Education'},'Units', 'normalized', 'Position', [0.5, 1.02, 0]);
         end
         
         set(gca, 'XTick', xticks, 'XTickLabel',largest_producer,'FontSize',font_size, 'XTickLabelRotation',90)
         figs = figs+1;
    end
end

for row=2:2
     for col=1:3
         
         sector = string(sectors(col));
         
         init_x = init_x_info.(sector);         
         largest_producer =largest_producer_info.(sector);
         
         subplot(4,3,figs,'Units','normalized','Position', [0.1+(col-1)*0.32, 0.76-(row-1)*0.23, 0.25, 0.16],'Units','normalized')
         for k = 1:length(init_x)   
            b=bar(k,init_x(k));
            set(b,'FaceColor',line_color_all(k,:));
            hold on
         end

         xticks = linspace(1, 5, 5);
         text(-0.2, 1.2, texts(figs), 'Units', 'Normalized','FontSize',font_size,'FontWeight','bold');
         
         if col==1
           ylabel({'Pre-pandemic', 'market share (%)'}, 'Units', 'normalized', 'Position', [-0.18, 0.5, 0]);
         end
        

         set(gca, 'XTick', xticks, 'XTickLabel',largest_producer,'FontSize',font_size, 'XTickLabelRotation',90)
         figs = figs+1;
    end
end


        
for row=3:4
    metric = metrics(row-1);
    
     for col=1:3
         sector = string(sectors(col));
         largest_producer =largest_producer_info.(sector);
         reopen_date = 1+ keep_weeks_info(specific_year)/point_length; 
         
         subplot(4,3,figs,'Units','normalized','Position', [0.1+(col-1)*0.32, 0.76-(row-1)*0.23, 0.25, 0.16])
         
         patch('XData',[1 reopen_date reopen_date 1],'YData',[y_lim_min(row-1) y_lim_min(row-1)  y_lim_max(row-1) y_lim_max(row-1)],...
          'FaceColor', [240/255,240/255,240/255],'EdgeColor','none'); % [240/255,240/255,240/255]
         hold on
         xline(reopen_date,'LineWidth', 1, 'color', [82/255,82/255,82/255])
         hold on

         h_all = [];

         for i=0:4
            line_info = string(strcat(base_scenario_key, year_info(specific_year),'_', sector, '_', string(i), metric)); % specific_year
            line_ = results.(line_info);
            h = plot(line_(1:point_length:length(line_)), '-', 'LineWidth', line_width, 'color', line_color_all(i+1,:));
            h_all = [h_all, h];
            hold on
         end
        
         yline(baseline_values(row-1),'--', 'LineWidth', line_width, 'color', 'black');

         
         if col==1
             if row==2
                  ylabel(string(strcat('Reliablity (', ylabels(row), ')')), 'Units', 'normalized', 'Position', [-0.18, 0.5, 0]);
             end
            if row==3
                  ylabel({'Variation in', string(strcat('profitability (', ylabels(row), ')'))}, 'Units', 'normalized', 'Position', [-0.18, 0.5, 0]);
            end
            if row==4
                  ylabel({'Production', string(strcat('adjustment ratio (', ylabels(row), ')'))}, 'Units', 'normalized', 'Position', [-0.18, 0.5, 0]);
            end


         end

         if row==3
            legend(h_all,largest_producer, 'Units', 'normalized','FontSize',font_size-3, 'NumColumns', 2, 'Position', [0.213+(col-1) * 0.32,0.4, 0.1, 0.06])
            legend boxoff 
         end
        text(-0.2, 1.2, texts(figs), 'Units', 'Normalized','FontSize',font_size,'FontWeight','bold');
         
        xlim([0 260/point_length])
        ylim([y_lim_min(row-1) y_lim_max(row-1)])
        if row==4
            xlabel('t [years]','FontSize',font_size);
        end
        set(gca, 'XTick', linspace(1, (length(line_)-1)/point_length+1, 6), 'XTickLabel',{'0','1','2','3','4','5'},'FontSize',font_size, 'Layer', 'top')
        box on
        
        figs = figs+1;
    end
end

saveas(gcf,figure_save_path,'epsc')
end
