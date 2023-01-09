function eco_draw_define_position(raw_epidemic_data_info, I_thre, Re_thre, phi, k, CHI_thre, rho, ...
    varphi, NPI_policy_scenario,y_lim_min, y_lim_max, if_base, supple_name, ...
    legend_position1, legend_position2, legend_position3, legend_position4, ...
    text1, text2, text3,text4)

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
results_path = strcat('results/', raw_epidemic_data_info, '/', eco_epi_hyper_paras_info, '/', main_file_name, '/main_results.csv');
if if_base
    figure_save_path = strcat('results/', raw_epidemic_data_info, '/', eco_epi_hyper_paras_info,'/', main_file_name, '/main_plot/eco_results.eps');
else
    figure_save_path = strcat('results/', raw_epidemic_data_info, '/supple_figs_summary/', supple_name,'_eco_results.eps');
end

%% ---------- read results -----------------------
results = readtable(results_path,'PreserveVariableNames',true);

%% --------- plot settings --------------------------
row_data = {'World', 'High-CHI countries/regions', 'Mid-CHI countries/regions', 'Low-CHI countries/regions'};
row_info = {'world', 'stringent_npi', 'moderate_npi',  'mild_npi'};
line_info_world = {'_di_GVA_', '_ind_propGVA_', '_ind_evoGVA_'};
line_info = {'_ave_di_GVA_', '_ave_ind_propGVA_', '_ave_ind_evoGVA_'};
col_info = {'1.0', '1.5', '2.0'};
keep_weeks_info = [52*1, 52*1.5, 52*2];
point_length = 4; % plot every [point_length] points
texts = char(97:120);
point_interval = 2;
marker_size = 2;
line_width = 1.5;
font_size = 11;

reopen_text_x = [0.23, 0.313, 0.4259];
face_color_all = [203/255, 230/255, 243/255;
            251/255 200/255 196/255;
            255/255,223/255,191/255;
            204/255,242/255,217/255;];
        
 line_color_all = [[43/255 140/255 190/255;
                    3/255,78/255,123/255;
                    3/255,78/255,123/255;];
             
             [239/255 59/255 44/255;
             153/255 0/255 13/255;
             153/255 0/255 13/255;];
             
             [255/255,127/255,0/255;
             204/255,76/255,2/255;
             204/255,76/255,2/255;];
             
             [35/255,139/255,69/255;
              0/255,88/255,36/255;
              0/255,88/255,36/255;]];
            
line_style = {'-', ':', '-.'};   
figure('Position', [91,132,846,828])

%% ------- plot ----------------------------------

figs = 1;
for row=1:4
    % ---- define legend and text positions --------------
    if row==1
        legend_position = legend_position1;
        text_position = text1;
    end
    
    if row==2
        legend_position = legend_position2;
        text_position = text2;
    end
    
    if row==3
        legend_position = legend_position3;
        text_position = text3;
    end
    
    if row==4
        legend_position = legend_position4;
        text_position = text4;
    end

% ----------------------------------------------------
     for col=1:3
         
        subplot(4,3,3*(row-1)+col,'Position', [0.08+(col-1)*0.32, 0.8-(row-1)*0.243, 0.24, 0.16],'Units','normalized')
        h_all = [];
        
        reopen_date = 1 + keep_weeks_info(col)/point_length;

        patch('XData',[0 reopen_date reopen_date 0],'YData',[y_lim_min(row) y_lim_min(row)  y_lim_max(row) y_lim_max(row)],...
              'FaceColor', [240/255,240/255,240/255],'EdgeColor','none');
        hold on
        
        %if row == 1
            line_info_specific = line_info_world;
            total_loss_name_specific = '_tot_GVA_';
%         else
%             line_info_specific = line_info;
%             total_loss_name_specific = '_ave_tot_GVA_';
%         end
        
        total_loss_line_ = string(strcat(base_scenario_key, col_info(col), total_loss_name_specific, row_info(row)));
        line_result = results.(total_loss_line_) * 100;
        a = area(line_result(1:point_length:length(line_result))', 'FaceColor', face_color_all(row, :), 'EdgeColor','none', 'ShowBaseLine', 'off');
        h_all = [h_all, a];
        hold on
        
        xline(reopen_date,'LineWidth', 1, 'color', [82/255,82/255,82/255])
        rotation_text_reopen = 90;
        y_text_reopen = 0.02;
        if col>1
            rotation_text_reopen = 0;
            y_text_reopen = 0.1;
        end
        
        if (col==1 && row==1)||(col==1 && row==4)
            rotation_text_reopen = 0;
            y_text_reopen = 0.9072;
        end

        text(text_position(1)+(col-1)*0.1, text_position(2),0, 'Full reopening','FontSize',font_size-2, 'color',[82/255,82/255,82/255],'Rotation',text_position(3), ...
             'Units','normalized'); % , 'Rotation',90
        
        yline(0,'--', 'LineWidth', line_width, 'color',[150/255,150/255,150/255]);
        hold on
        
   
        for i=1:3
           	line_ = string(strcat(base_scenario_key, col_info(col),line_info_specific(i), row_info(row)));
            line_result = results.(line_) * 100;
            h = plot(line_result(1:point_length:length(line_result)),string(line_style(i)), 'LineWidth', line_width, 'color', line_color_all((row-1)*3+i,:));
            h_all = [h_all, h];
            hold on
        end
        
        set(gca,'FontSize',font_size, 'Layer', 'top')
        box off
        

        xlim([1 (length(line_result)-1)/point_length+1])
        ylim([y_lim_min(row) y_lim_max(row)])
        xticks(1:52/point_length:(length(line_result)-1)/point_length+1)
        xticklabels({'0','1','2','3','4','5'})

        text(-0.2, 1.2, texts(figs), 'Units', 'Normalized','FontSize',font_size,'FontWeight','bold');
        figs = figs +1;
        
        title(row_data(row), 'FontSize',font_size)
        
        if row==4
            xlabel('t [years]','FontSize',font_size);
        end
        
        if col==1
            ylabel('GVA change (%)','FontSize',font_size, 'Units', 'normalized', 'Position', [-0.16, 0.5, 0]);
        end
        
        if col==1
            legend(h_all,'Total','Pandemic','Supply-chain', 'Competition', 'Units', 'normalized',...
                'Position', legend_position,'FontSize',font_size-3)
            legend boxoff 
        end 
    end
end
saveas(gcf,figure_save_path,'epsc')
end