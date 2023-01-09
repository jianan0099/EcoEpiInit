% ------- baseline --------------------
raw_epidemic_data_info = '2022-03-01_IP_7_2022-11-11';
I_thre = '0.001';
Re_thre = '1';
phi = '0.01';
k = '100';
CHI_thre = '25_50';
rho = '0.01';
varphi = '0.7';
NPI_policy_scenario = 'keep_curr_';
% -------------------------------------

% main text
epi_draw(raw_epidemic_data_info, I_thre, Re_thre, phi, k, CHI_thre, rho, ...
    varphi, NPI_policy_scenario,  [0 4], [0.75 1.05], [0.75 0.85 0.95 1.05], true, '')
eco_draw_add(raw_epidemic_data_info, I_thre, Re_thre, phi, k, CHI_thre, rho, ...
    varphi, NPI_policy_scenario, [-50, -80, -60, -40], [20, 10, 10, 25], true, '', [0.15, 0.83, 0.08, 0.05])
com_explain(raw_epidemic_data_info, I_thre, Re_thre, phi, k, CHI_thre, rho, ...
    varphi, NPI_policy_scenario, true, '')
sector_change_modify(raw_epidemic_data_info, I_thre, Re_thre, phi, k, CHI_thre, rho, ...
    varphi, NPI_policy_scenario, true, '')


% supple only for base
market_structure_change_fit(raw_epidemic_data_info, I_thre, Re_thre, phi, k, CHI_thre, rho, varphi, NPI_policy_scenario)
reg_var_dis(raw_epidemic_data_info, I_thre, Re_thre, phi, k, CHI_thre, rho, varphi)
risk_diversity_sample(raw_epidemic_data_info, I_thre, Re_thre, phi, k, CHI_thre)
sector_change_all(raw_epidemic_data_info, I_thre, Re_thre, phi, k, CHI_thre, rho, varphi, NPI_policy_scenario)

% supple sensitivity analysis
% change rho
eco_draw(raw_epidemic_data_info, I_thre, Re_thre, phi, k, CHI_thre, '0.002', ...
     varphi, NPI_policy_scenario, [-50, -80, -60, -40], [20, 10, 10, 25], false, 'rho_sensi_0.002', [0.19, 0.82, 0.08, 0.05])
eco_draw(raw_epidemic_data_info, I_thre, Re_thre, phi, k, CHI_thre, '0.02', ...
     varphi, NPI_policy_scenario, [-50, -80, -60, -40], [30, 10, 10, 45], false, 'rho_sensi_0.02', [0.19, 0.81, 0.08, 0.05])
eco_draw_define_position(raw_epidemic_data_info, I_thre, Re_thre, phi, k, CHI_thre, '0.05', ...
    varphi, NPI_policy_scenario, [-50, -80, -60, -40], [60, 10, 40, 80], false, 'rho_sensi_0.05', ...
    [0.184491726860369,0.881533817432929,0.133569737768624,0.074275360303225],...
    [0.183309693763442,0.591432368157567,0.133569737768624,0.074275360303225], ...
    [0.163215131115688,0.321862319848388,0.133569737768624,0.074275360303225], ...
    [0.180945627569589,0.15640338265032,0.133569737768624,0.074275360303225], ...
    [0.223,0.088, 0], [0.223, 0.088, 0], [0.223, 0.9072, 0],[0.223, 0.088, 0])

% change varphi
eco_draw(raw_epidemic_data_info, I_thre, Re_thre, phi, k, CHI_thre, rho, ...
     '0.3', NPI_policy_scenario, [-50, -80, -60, -40], [20, 10, 10, 25], false, 'varphi_sensi_0.3', [0.19, 0.82, 0.08, 0.05])
eco_draw(raw_epidemic_data_info, I_thre, Re_thre, phi, k, CHI_thre, rho, ...
    '0.9', NPI_policy_scenario,[-50, -80, -60, -40], [20, 10, 10, 25], false, 'varphi_sensi_0.9', [0.19, 0.82, 0.08, 0.05])

% change phi
epi_draw(raw_epidemic_data_info, I_thre, Re_thre, '0.005', k, CHI_thre, rho, ...
       varphi, NPI_policy_scenario,  [0 4], [1 1.5], [1 1.25 1.5], false, 'phi_sensi_0.005')
eco_draw(raw_epidemic_data_info, I_thre, Re_thre, '0.005', k, CHI_thre, rho, ...
     varphi, NPI_policy_scenario, [-50, -80, -60, -40], [20, 10, 10, 25], false, 'phi_sensi_0.005', [0.19, 0.82, 0.08, 0.05])
% change I_thre
eco_draw(raw_epidemic_data_info, '0.1', Re_thre, phi, k, CHI_thre, rho, ...
     varphi, NPI_policy_scenario, [-50, -80, -60, -40], [20, 10, 10, 25], false, 'I_thre_sensi_0.1', [0.19, 0.82, 0.08, 0.05])