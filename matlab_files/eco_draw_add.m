function eco_draw_add(raw_epidemic_data_info, I_thre, Re_thre, phi, k, CHI_thre, rho, ...
    varphi, NPI_policy_scenario,y_lim_min, y_lim_max, if_base, supple_name, legend_position)

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
timing_path = strcat('results/', raw_epidemic_data_info, '/', eco_epi_hyper_paras_info, '/', main_file_name, '/timing_of_reopen.xlsx');
if if_base
    figure_save_path = strcat('results/', raw_epidemic_data_info, '/', eco_epi_hyper_paras_info,'/', main_file_name, '/main_plot/eco_results.eps');
else
    figure_save_path = strcat('results/', raw_epidemic_data_info, '/supple_figs_summary/', supple_name,'_eco_results.eps');
end

%% ---------- read results -----------------------
results = readtable(results_path,'PreserveVariableNames',true);
timing = readtable(timing_path,'PreserveVariableNames',true);
%% --------- plot settings --------------------------
row_data = {'World', 'High-CHI countries', 'Mid-CHI countries', 'Low-CHI countries'};
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
figure_width_range = 0.24;
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
figure('Position', [91,132,969,828])

%% ------- plot ----------------------------------


for row=1:4
     for col=1:3
        subplot(4,4,4*(row-1)+col,'Position', [0.06+(col-1)*figure_width_range, 0.81-(row-1)*0.243, 0.19, 0.15],'Units','normalized')
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
        text(reopen_text_x(col),y_text_reopen,0, 'Full reopening','FontSize',font_size-2, 'color',[82/255,82/255,82/255],'Rotation',rotation_text_reopen, ...
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
        
        figs = 4*(row-1)+col;
        text(-0.1, 1.2, texts(figs), 'Units', 'Normalized','FontSize',font_size,'FontWeight','bold');
        
        
        title(row_data(row), 'FontSize',font_size)
        
        if row==4
            xlabel('t [years]','FontSize',font_size);
        end
        
        if col==1
            ylabel('GVA change (%)','FontSize',font_size, 'Units', 'normalized', 'Position', [-0.2, 0.5, 0]);
        end
        
        if col==1
            legend(h_all,'Total','Pandemic','Supply disruption', 'Competition', 'Units', 'normalized',...
                'Position', [legend_position(1)+(col-1) * figure_width_range, legend_position(2)-(row-1)* 0.243, legend_position(3), legend_position(4)],'FontSize',font_size-1.5)
            legend boxoff 
        end 
    end
end

col=4;
for row=1:4
     subplot(4,4,4*(row-1)+col,'Position', [0.06+(col-1)*figure_width_range, 0.81-(row-1)*0.243, 0.19, 0.15],'Units','normalized')
     
     timing_line_ = string(row_info(row));
     timing_result = timing.(timing_line_)*100;
     h = plot(timing_result, 'LineWidth', line_width,'color', line_color_all((row-1)*3+1,:),'Marker', 'o','MarkerSize', 5, 'MarkerFaceColor', line_color_all((row-1)*3+1,:));
     
     set(gca,'FontSize',font_size, 'Layer', 'top')
     box off
     
     title(row_data(row), 'FontSize',font_size)
     if row==4
            xlabel({'Timing of full-reopening', '[months]'},'FontSize',font_size);
     end
             figs = 4*(row-1)+col;
        text(-0.1, 1.2, texts(figs), 'Units', 'Normalized','FontSize',font_size,'FontWeight','bold');
    xlim([1 13])
    xticks(1:3:13)
    xticklabels({'12','15','18','21','24'})
    
    if row==1
        legend_position = [0.781517029493831,0.816884057971014,0.16202270055826,0.05];
    end
   
        if row==2
        legend_position = [0.781517029493831,0.5705,0.16202270055826,0.05];
        end
       
        if row==3
        legend_position = [0.824860682744601,0.324975846058504,0.16202270055826,0.048309177448208];
        end
       
        if row==4
        legend_position = [0.830836319498749,0.064009661835749,0.155830750278398,0.087898550724636];
    end
     legend(h,strcat('At the end of the', string(newline), 'projection horizon'), 'Units', 'normalized',...
                'Position', legend_position, 'FontSize',font_size-1.5)
     legend boxoff 
end
saveas(gcf,figure_save_path,'epsc')
end