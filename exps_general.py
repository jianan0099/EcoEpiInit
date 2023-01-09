import scipy
import statsmodels.api as sm
from scipy.stats import shapiro
import igraph as ig
from ARIO_modify import *
import pandas as pd
import numpy as np
import networkx as nx
import powerlaw


class FinalExps:
    def __init__(self, I_fraction_threshold=1e-3, Re_thre=1, phi=0.01, k=100, num_time_steps_in_inv_init=3,
                 NPI_policy_scenario='keep_curr'):
        self.eco_epi_model = ARIO(I_fraction_threshold=I_fraction_threshold,
                                  Re_thre=Re_thre, phi=phi, k=k,
                                  num_time_steps_in_inv_init=num_time_steps_in_inv_init)
        self.economic_region_npi_policy, self.economic_region_chi_value = \
            self.eco_epi_model.epidemic_model.define_country_npi_policy_class_in_economic()
        self.NPI_policy_scenario = NPI_policy_scenario

    def get_stringency_info_csv(self):
        stringency_info = [['ISO_A3', 'CHI_class']]
        for i in range(self.eco_epi_model.epidemic_model.c_num):
            stringency_info.append([self.eco_epi_model.epidemic_model.common[i].upper(),
                                    self.eco_epi_model.epidemic_model.stringency_info_current[i]])
        df = pd.DataFrame(stringency_info)
        df.to_csv(self.eco_epi_model.hyper_matlab_results_path + '/stringency_info.csv', index=False, header=None)

    def specify_timing(self, rho=0.01, varphi=0.7, timing_range=range(52 * 7 * 1, 52 * 7 * 2+1, 4 * 7)):
        main_results_file = self.eco_epi_model.hyper_matlab_results_path + '/' \
                            + 'main_rho_' + str(rho) + '_varphi_' + str(varphi)
        if not os.path.exists(main_results_file):
            os.mkdir(main_results_file)

        results = [['scenario', 'world', 'mild_npi', 'moderate_npi', 'stringent_npi']]
        for keep_time_step in timing_range:
            scenario_key, eco_data = \
                self.eco_epi_model.IO_trans_simplify(rho, varphi, self.NPI_policy_scenario, keep_time_step)
            value_added_actual_by_country = \
                self.eco_epi_model.cal_sum_by_regions_for_time_series_arr(eco_data['value_added_actual'])
            value_added_original_by_country = self.eco_epi_model.cal_sum_by_regions_for_time_series_arr(
                np.array(eco_data['value_added_original'])[None, :]).flatten()

            current_result = [scenario_key]
            for country_class in ['world', 'mild_npi', 'moderate_npi', 'stringent_npi']:
                corr_regions_index_list = self.economic_region_npi_policy[country_class]
                current_result.append(np.sum(value_added_actual_by_country[-1, corr_regions_index_list]) / \
                                      sum(value_added_original_by_country[corr_regions_index_list]) - 1)
            results.append(current_result)
        df = pd.DataFrame(results)
        df.to_excel(main_results_file + '/timing_of_reopen.xlsx', index=False, header=None)

    def main_results(self, rho=0.01, varphi=0.7):
        main_results_file = self.eco_epi_model.hyper_matlab_results_path + '/' \
                            + 'main_rho_' + str(rho) + '_varphi_' + str(varphi)
        if not os.path.exists(main_results_file):
            os.mkdir(main_results_file)

        results = {}
        for time_index, keep_time_step in enumerate([52 * 7 * 1, 52 * 7 * 1.5, 52 * 7 * 2]):
            print(time_index)
            scenario_key, eco_epi_data = self.eco_epi_model.load_economic_data(rho, varphi, self.NPI_policy_scenario,
                                                                               keep_time_step,
                                                                               if_check_economic_data=True)
            print('get eco data')

            prevalence = eco_epi_data['epidemic_data']['I_frac']
            cum_mortality_rate = eco_epi_data['epidemic_data']['unavailable_l_frac'] - eco_epi_data['epidemic_data'][
                'I_frac']

            value_added_actual_by_country = \
                self.eco_epi_model.cal_sum_by_regions_for_time_series_arr(eco_epi_data['value_added_actual'])
            value_added_without_evo_by_country = \
                self.eco_epi_model.cal_sum_by_regions_for_time_series_arr(
                    np.array(eco_epi_data['value_added_without_evo']))
            value_added_max_by_country = \
                self.eco_epi_model.cal_sum_by_regions_for_time_series_arr(np.array(eco_epi_data['value_added_max']))
            value_added_change_due_to_evo_by_country = value_added_max_by_country - value_added_without_evo_by_country
            value_added_original_by_country = self.eco_epi_model.cal_sum_by_regions_for_time_series_arr(
                np.array(eco_epi_data['value_added_original'])[None, :]).flatten()

            for country_class in self.economic_region_npi_policy:
                corr_regions_index_list = self.economic_region_npi_policy[country_class]
                results[scenario_key + '_prev_' + country_class] = \
                    np.mean(prevalence[:, corr_regions_index_list], axis=1)
                results[scenario_key + '_cum_mor_' + country_class] = \
                    np.mean(cum_mortality_rate[:, corr_regions_index_list], axis=1)

                # ----- total ----------
                results[scenario_key + '_ave_tot_GVA_' + country_class] = \
                    np.mean(value_added_actual_by_country[:, corr_regions_index_list] /
                            value_added_original_by_country[corr_regions_index_list][None, :], axis=1) - 1
                results[scenario_key + '_tot_GVA_' + country_class] = \
                    np.sum(value_added_actual_by_country[:, corr_regions_index_list], axis=1) / \
                    sum(value_added_original_by_country[corr_regions_index_list]) - 1

                # ----- direct ----------
                results[scenario_key + '_ave_di_GVA_' + country_class] = \
                    np.mean(value_added_without_evo_by_country[:, corr_regions_index_list] /
                            value_added_original_by_country[corr_regions_index_list][None, :], axis=1) - 1
                results[scenario_key + '_di_GVA_' + country_class] = \
                    np.sum(value_added_without_evo_by_country[:, corr_regions_index_list], axis=1) / \
                    sum(value_added_original_by_country[corr_regions_index_list]) - 1

                # ----- indirect by evo ----------
                results[scenario_key + '_ave_ind_evoGVA_' + country_class] = \
                    np.mean(value_added_change_due_to_evo_by_country[:, corr_regions_index_list] /
                            value_added_original_by_country[corr_regions_index_list][None, :], axis=1)
                results[scenario_key + '_ind_evoGVA_' + country_class] = \
                    np.sum(value_added_change_due_to_evo_by_country[:, corr_regions_index_list], axis=1) / \
                    sum(value_added_original_by_country[corr_regions_index_list])

                # ----- indirect by propagation ----------
                results[scenario_key + '_ave_ind_propGVA_' + country_class] = \
                    results[scenario_key + '_ave_tot_GVA_' + country_class] - \
                    results[scenario_key + '_ave_di_GVA_' + country_class] - \
                    results[scenario_key + '_ave_ind_evoGVA_' + country_class]
                results[scenario_key + '_ind_propGVA_' + country_class] = \
                    results[scenario_key + '_tot_GVA_' + country_class] - \
                    results[scenario_key + '_di_GVA_' + country_class] - \
                    results[scenario_key + '_ind_evoGVA_' + country_class]
        results_df = pd.DataFrame(results)
        results_df.to_csv(main_results_file + '/main_results.csv', index=False)

    def market_structure_change(self, rho=0.01, varphi=0.7):
        main_results_file = self.eco_epi_model.hyper_matlab_results_path + '/' + 'main_rho_' \
                            + str(rho) + '_varphi_' + str(varphi)

        # ------ for market share --------------------
        full_market_share_change_results_save_path = main_results_file + '/full_market_change.xlsx'
        all_comms_df = pd.DataFrame([comm.upper() for comm in self.eco_epi_model.activities])
        all_comms_df.to_excel(full_market_share_change_results_save_path,
                              sheet_name='heat_x_labels', index=False, header=None)

        all_regions_df = pd.DataFrame([reg.upper() for reg in self.eco_epi_model.regions])
        with pd.ExcelWriter(full_market_share_change_results_save_path, engine="openpyxl", mode='a') as writer:
            all_regions_df.to_excel(writer, sheet_name='heat_y_labels', index=False, header=None)

        market_share_change_results_save_path = main_results_file + '/market_change.xlsx'
        market_share_change_results = {}

        # output sector
        target_sector = [
            'pdr',  # Paddy rice production=0.5, demand=1
            'rmk',  # Raw milk production=0.5, demand=1

            'afs',  # Accommodation, Food and service activities production=1, demand=0.25
            'tex',  # Textiles  production=1, demand=0.5
            'eeq',  # Electrical equipment  production=1, demand=0.9
            'coa',  # Coal production=1, demand=1
            'trd',  # Trade production=1, demand=1.1

            'edu',  # Education production=0.1, demand=0.85
            'gdt',  # Gas manufacture production=0.1, demand=1
            'osg',  # Public Administration and defense production=0.1, demand=1.25
        ]
        reshape_init_x = np.reshape(self.eco_epi_model.x, (self.eco_epi_model.num_regions, -1))
        target_region_index = \
            np.argsort(np.sum(reshape_init_x, axis=1))[::-1][:10]  # top 10
        target_region = [self.eco_epi_model.regions[index] for index in target_region_index]
        min_value_in_heatmap, max_value_in_heatmap = np.inf, -np.inf
        example_sector = ['pdr',  # Paddy rice production=0.5, demand=1
                          'eeq',  # Electrical equipment  production=1, demand=0.9
                          'tex',  # Textiles  production=1, demand=0.5
                          'edu',  # Education production=0.1, demand=0.85
                          ]
        regions_CHI_class = {}
        items_with_largest_three_regions_in_example_sector = []
        for sector in example_sector:
            example_sector_index = self.eco_epi_model.activities.index(sector)
            for region_index in np.argsort(reshape_init_x[:, example_sector_index])[::-1][:3]:
                regions_CHI_class[self.eco_epi_model.regions[region_index]] = \
                    self.eco_epi_model.get_eco_region_CHI_class(self.economic_region_npi_policy, region_index)
                items_with_largest_three_regions_in_example_sector.append(
                    self.eco_epi_model.regions[region_index] + '_' + sector)
        # --------------------------------------------

        for time_index, keep_time_step in enumerate([52 * 7 * 1, 52 * 7 * 1.5, 52 * 7 * 2.0]):
            target_sectors_market_structure_change = []
            scenario_key, eco_epi_data = self.eco_epi_model.load_economic_data(rho, varphi, self.NPI_policy_scenario,
                                                                          keep_time_step, if_check_economic_data=True)

            # ----- market share ------------
            # heat_map
            init_production_distribution = eco_epi_data['production_distribution'][0]
            final_production_distribution = eco_epi_data['production_distribution'][-1]
            sector_distribution_change = np.sum(
                np.reshape(np.power(final_production_distribution - init_production_distribution, 2),
                           (self.eco_epi_model.num_regions, -1)), axis=0)
            product_distribution_change = \
                np.reshape(self.eco_epi_model.division_define(
                    final_production_distribution - init_production_distribution, init_production_distribution, 0),
                    (self.eco_epi_model.num_regions, -1))
            needed_production_distribution_change = \
                product_distribution_change[target_region_index][:, [self.eco_epi_model.activities.index(sector)
                                                                     for sector in target_sector]]
            df = pd.DataFrame(needed_production_distribution_change)
            if time_index > 0:
                with pd.ExcelWriter(market_share_change_results_save_path, engine="openpyxl", mode='a') as writer:
                    df.to_excel(writer, sheet_name=scenario_key + '_heat', index=False, header=None)
            else:
                df.to_excel(market_share_change_results_save_path, sheet_name=scenario_key + '_heat',
                            index=False, header=None)

            last_day_in_epidemic_control = \
                self.eco_epi_model.get_last_day_in_epidemic_control_in_eco_time_step(
                    self.NPI_policy_scenario, keep_time_step)
            reopen_time_production_distribution = eco_epi_data['production_distribution'][last_day_in_epidemic_control]
            reopen_time_product_distribution_change = \
                np.reshape(self.eco_epi_model.division_define(
                    reopen_time_production_distribution - init_production_distribution, init_production_distribution,
                    0),
                    (self.eco_epi_model.num_regions, -1))
            reopen_time_needed_production_distribution_change = \
                reopen_time_product_distribution_change[target_region_index][:,
                [self.eco_epi_model.activities.index(sector) for sector in target_sector]]

            for sector in target_sector:
                target_sectors_market_structure_change.append(
                    sector_distribution_change[self.eco_epi_model.activities.index(sector)])

            df = pd.DataFrame(target_sectors_market_structure_change)
            with pd.ExcelWriter(market_share_change_results_save_path, engine="openpyxl", mode='a') as writer:
                df.to_excel(writer, sheet_name=scenario_key + '_mar', index=False, header=None)

            df = pd.DataFrame(reopen_time_needed_production_distribution_change)
            with pd.ExcelWriter(market_share_change_results_save_path, engine="openpyxl", mode='a') as writer:
                df.to_excel(writer, sheet_name=scenario_key + '_re_heat', index=False, header=None)

            min_value_in_heatmap, max_value_in_heatmap = min(min_value_in_heatmap,
                                                             np.min(needed_production_distribution_change)), \
                                                         max(max_value_in_heatmap,
                                                             np.max(needed_production_distribution_change))

            for example_item in items_with_largest_three_regions_in_example_sector:
                market_share_change_results[scenario_key + '_' + example_item] = \
                    np.array(eco_epi_data['production_distribution'])[:,
                    list(self.eco_epi_model.item_info).index(example_item)]

            # ------ full market share -------------------
            temp_market_share_change_results = np.reshape(product_distribution_change,
                                                          (self.eco_epi_model.num_regions, -1))
            df = pd.DataFrame(temp_market_share_change_results)
            with pd.ExcelWriter(full_market_share_change_results_save_path, engine="openpyxl", mode='a') as writer:
                df.to_excel(writer, sheet_name=scenario_key + '_heat', index=False, header=None)

            reopen_time_temp_market_share_change_results = np.reshape(reopen_time_product_distribution_change,
                                                                      (self.eco_epi_model.num_regions, -1))
            df = pd.DataFrame(reopen_time_temp_market_share_change_results)
            with pd.ExcelWriter(full_market_share_change_results_save_path, engine="openpyxl", mode='a') as writer:
                df.to_excel(writer, sheet_name=scenario_key + '_re_heat', index=False, header=None)

        # ---- heatmap info ------
        df = pd.DataFrame([min_value_in_heatmap, max_value_in_heatmap])
        with pd.ExcelWriter(market_share_change_results_save_path, engine="openpyxl", mode='a') as writer:
            df.to_excel(writer, sheet_name='heat_min_max', index=False, header=None)

        df = pd.DataFrame([sector.upper() for sector in target_sector])
        with pd.ExcelWriter(market_share_change_results_save_path, engine="openpyxl", mode='a') as writer:
            df.to_excel(writer, sheet_name='heat_x_labels', index=False, header=None)

        df = pd.DataFrame([region.upper() for region in target_region])
        with pd.ExcelWriter(market_share_change_results_save_path, engine="openpyxl", mode='a') as writer:
            df.to_excel(writer, sheet_name='heat_y_labels', index=False, header=None)

        # ---- market share change info ---------
        df = pd.DataFrame(market_share_change_results)
        with pd.ExcelWriter(market_share_change_results_save_path, engine="openpyxl", mode='a') as writer:
            df.to_excel(writer, sheet_name='market_share_example', index=False)
        df = pd.DataFrame(regions_CHI_class, index=[0]).T
        with pd.ExcelWriter(market_share_change_results_save_path, engine="openpyxl", mode='a') as writer:
            df.to_excel(writer, sheet_name='regions_CHI_class')

    @staticmethod
    def define_significance_level(p_value):
        if p_value < 0.1:
            if p_value < 0.05:
                if p_value < 0.01:
                    return '***'
                else:
                    return '**'
            else:
                return '*'
        else:
            return ''

    def regression(self, rho=0.01, varphi=0.7):
        main_results_file = self.eco_epi_model.hyper_matlab_results_path + '/' \
                            + 'main_rho_' + str(rho) + '_varphi_' + str(varphi)

        sector_structure_change = {}
        sector_structure_distribution = {}
        sector_structure_distribution_power_law_compare = {}
        sector_regression_results = [{}, {}, {}]
        sector_regression_results_save_path = main_results_file + '/sector_regression_results.xlsx'

        for keep_time_step in [52 * 7 * 1, 52 * 7 * 1.5, 52 * 7 * 2.0]:  # 52 * 7 * 1, 52 * 7 * 2,
            # ---- get results -------
            scenario_key, eco_epi_data = self.eco_epi_model.load_economic_data(rho, varphi, self.NPI_policy_scenario,
                                                                          keep_time_step)
            before_distribution = np.reshape(eco_epi_data['production_distribution'][0],
                                             (self.eco_epi_model.num_regions, -1))
            final_distribution = np.reshape(eco_epi_data['production_distribution'][-1],
                                            (self.eco_epi_model.num_regions, -1))

            trade_by_sector = self.eco_epi_model.group_Z_by_sector(self.eco_epi_model.Z)
            self_feed = trade_by_sector.diagonal().copy() / np.sum(trade_by_sector, axis=0)
            np.fill_diagonal(trade_by_sector, 0)
            trade_by_sector_G = nx.from_numpy_array(trade_by_sector, create_using=nx.DiGraph())
            trade_by_sector_ig = ig.Graph.from_networkx(trade_by_sector_G)

            CHI_region = np.array([self.economic_region_chi_value[region.upper()]
                                   for region in self.eco_epi_model.regions])
            sector_distribution_change = []
            sector_diversity = []

            for i in range(self.eco_epi_model.num_commodities):
                weighted_x = before_distribution[:, i] * CHI_region
                sector_diversity.append(-1 * scipy.stats.entropy(weighted_x / sum(weighted_x)))
                sector_distribution_change.append(
                    np.sum(np.power(before_distribution[:, i] - final_distribution[:, i], 2)))
            log_sector_distribution_change = np.log(sector_distribution_change)

            sector_structure_change[scenario_key + '_original'] = sector_distribution_change
            sector_structure_change[scenario_key + '_log'] = log_sector_distribution_change


            ks_log, p_log = scipy.stats.kstest(sector_distribution_change, "lognorm",
                                               scipy.stats.lognorm.fit(sector_distribution_change))
            ks_pare, p_pare = scipy.stats.kstest(sector_distribution_change, "pareto",
                                                 scipy.stats.pareto.fit(sector_distribution_change))
            R, p = powerlaw.Fit(sector_distribution_change).distribution_compare('power_law', 'lognormal')
            sector_structure_distribution_power_law_compare[scenario_key + '_R'] = [str(round(R, 3))]
            sector_structure_distribution_power_law_compare[scenario_key + '_p'] = ['p < 0.1' if p < 0.1 else '']
            sector_structure_distribution_power_law_compare[scenario_key + '_ks_log'] = [str(round(ks_log, 3))]
            sector_structure_distribution_power_law_compare[scenario_key + '_p_log'] = [str(round(p_log, 3))]
            sector_structure_distribution_power_law_compare[scenario_key + '_ks_pare'] = [str(round(ks_pare, 3))]
            sector_structure_distribution_power_law_compare[scenario_key + '_p_pare'] = [str(round(p_pare, 3))]


            sector_regression_data = {
                # # containment policy influence
                # 'CHI': CHI,
                'Impact-to-labour multiplier ($vartheta^P$)': self.eco_epi_model.NPI_sector_multiplier,
                'Impact-to-demand multiplier ($vartheta^D$)': self.eco_epi_model.demand_multiplier_sector,
                'Risk diversity': sector_diversity,
                'Fraction of self-supplied demand': self_feed,

                # choose 1 4/5/6
                'Betweenness centrality': trade_by_sector_ig.betweenness(directed=False),
                'Degree': [trade_by_sector_G.degree(node) for node in trade_by_sector_G.nodes()],
                'Eigenvector centrality': trade_by_sector_ig.eigenvector_centrality(directed=False),

                'sector_distribution_change': log_sector_distribution_change}

            sector_regression_df = pd.DataFrame(sector_regression_data)
            normalized_sector_regression_df = (sector_regression_df - sector_regression_df.mean()) \
                                              / sector_regression_df.std()

            coff_sector_model = normalized_sector_regression_df.iloc[:-1, :-1].corr().round(3)

            y = sector_regression_df[
                'sector_distribution_change']  # normalized_sector_regression_df['sector_distribution_change']
            sector_structure_distribution[scenario_key] = [str(round(shapiro(y)[1], 3))]

            keep_info = str(keep_time_step / (52 * 7))
            for coff_index in [4, 5, 6]:
                selected_coff_name = list(sector_regression_data.keys())[coff_index]
                final_coff_index = list(range(4)) + [coff_index]
                x = normalized_sector_regression_df.iloc[:, final_coff_index]
                x = sm.add_constant(x)  # adding a constant
                model = sm.OLS(y, x).fit()
                # print(model.summary())
                # coff_sector_model = np.corrcoef(normalized_sector_regression_df.iloc[:, final_coff_index].T)

                if keep_info == '1.0':
                    sector_regression_results[coff_index - 4]['Variables'] = []
                    for item in list(sector_regression_data.keys())[:4]:
                        sector_regression_results[coff_index - 4]['Variables'].append(item)
                    sector_regression_results[coff_index - 4]['Variables'].append(selected_coff_name)
                    sector_regression_results[coff_index - 4]['Variables'].append('No. Observations')
                    sector_regression_results[coff_index - 4]['Variables'].append('$R^2$')
                    sector_regression_results[coff_index - 4]['Variables'].append('Adjusted $R^2$')
                    sector_regression_results[coff_index - 4]['Variables'].append('F-statistic')
                    sector_regression_results[coff_index - 4]['Variables'].append('Residual Shapiro-Wilks p')
                regression_coff = model.params.values[1:]
                regression_std = model.bse.values[1:]
                regression_p = model.pvalues.values[1:]
                shapiro_p = shapiro(model.resid)[1]

                final_regression_coff = []
                for i in range(len(regression_coff)):
                    coff, std = round(regression_coff[i], 3), round(regression_std[i], 3)
                    final_regression_coff.append('$' + str(coff) + '\,(' + str(std) + ')'
                                                 + '^{' + self.define_significance_level(regression_p[i]) + '}$')
                final_regression_coff.append(model.nobs)
                final_regression_coff.append(round(model.rsquared, 3))
                final_regression_coff.append(round(model.rsquared_adj, 3))
                final_regression_coff.append('$' + str(round(model.fvalue, 3))
                                             + '^{' + self.define_significance_level(model.f_pvalue) + '}$')
                final_regression_coff.append('$' + str(round(shapiro_p, 3)) + '^{'
                                             + self.define_significance_level(shapiro_p) + '}$')

                sector_regression_results[coff_index - 4][keep_info] = final_regression_coff

        for i in range(len(sector_regression_results)):
            selected_coff_name = list(sector_regression_data.keys())[i + 4]
            df = pd.DataFrame(sector_regression_results[i])
            if i == 0:
                df.to_excel(sector_regression_results_save_path, sheet_name=selected_coff_name, index=False)
            else:
                with pd.ExcelWriter(sector_regression_results_save_path, engine="openpyxl", mode='a') as writer:
                    df.to_excel(writer, sheet_name=selected_coff_name, index=False)
        df = pd.DataFrame(sector_structure_change)
        with pd.ExcelWriter(sector_regression_results_save_path, engine="openpyxl", mode='a') as writer:
            df.to_excel(writer, sheet_name='sector_structure_change', index=False)
        df = pd.DataFrame(sector_structure_distribution)
        with pd.ExcelWriter(sector_regression_results_save_path, engine="openpyxl", mode='a') as writer:
            df.to_excel(writer, sheet_name='sector_stru_dis_shapiro', index=False)
        df = pd.DataFrame(sector_structure_distribution_power_law_compare)
        with pd.ExcelWriter(sector_regression_results_save_path, engine="openpyxl", mode='a') as writer:
            df.to_excel(writer, sheet_name='sector_stru_log_norm', index=False)
        with pd.ExcelWriter(sector_regression_results_save_path, engine="openpyxl", mode='a') as writer:
            coff_sector_model.to_excel(writer, sheet_name='corr')
        with pd.ExcelWriter(sector_regression_results_save_path, engine="openpyxl", mode='a') as writer:
            normalized_sector_regression_df.to_excel(writer, sheet_name='norm_data', index=False)

    def competition_explain(self, rho=0.01, varphi=0.7):
        main_results_file = self.eco_epi_model.hyper_matlab_results_path + '/' + 'main_rho_' + str(rho) + '_varphi_' \
                            + str(varphi)
        results_save_path = main_results_file + '/competition_explain.xlsx'

        reshape_init_x = np.reshape(self.eco_epi_model.x, (self.eco_epi_model.num_regions, -1))
        sector_list = ['pdr', 'eeq', 'edu']
        competition_results = {}
        init_x_info = {}
        vartheta_info = {}
        corr_max_producer_name = {}
        CHI_scores = {}
        CHI_levels = {}
        for keep_time_step in [52 * 7 * 1, 52 * 7 * 1.5, 52 * 7 * 2]:
            scenario_key, eco_epi_data = self.eco_epi_model.load_economic_data(rho, varphi, self.NPI_policy_scenario,
                                                                          keep_time_step, if_check_economic_data=True)
            for sec in sector_list:
                sec_index = self.eco_epi_model.activities.index(sec)

                sector_ave_profit = \
                    np.sum(np.array(eco_epi_data['profit'])[:, self.eco_epi_model.product_index_list_in_items(sec_index)] *
                           np.array(eco_epi_data['production_distribution'])[:,
                           self.eco_epi_model.product_index_list_in_items(sec_index)], axis=1)

                target_region_list_index = np.argsort(reshape_init_x[:, sec_index])[::-1][:5]  # top 5
                if keep_time_step == 52 * 7 * 1:
                    init_x_info[sec] = 100 * reshape_init_x[target_region_list_index, sec_index] / sum(
                        reshape_init_x[:, sec_index])
                    vartheta_info[sec] = [self.eco_epi_model.NPI_sector_multiplier[sec_index],
                                          self.eco_epi_model.demand_multiplier_sector[sec_index]]
                    corr_max_producer_name[sec] = [self.eco_epi_model.regions[reg_index].upper() for reg_index in
                                                   target_region_list_index]
                    CHI_scores[sec] = [self.economic_region_chi_value[self.eco_epi_model.regions[reg_index].upper()] for reg_index
                                       in target_region_list_index]
                    CHI_levels[sec] = [self.eco_epi_model.get_eco_region_CHI_class(self.economic_region_npi_policy, reg_index) for
                                       reg_index in target_region_list_index]
                for rank, reg_index in enumerate(target_region_list_index):
                    item_index = reg_index * 65 + sec_index
                    reliability = np.array(eco_epi_data['reliability'])[:, item_index]
                    kappa = np.array(eco_epi_data['kappa'])[:, item_index]
                    profit = np.array(eco_epi_data['profit'])[:, item_index]
                    var_profit = profit - sector_ave_profit
                    total_order = np.array(eco_epi_data['total_order'])[:, item_index]

                    competition_results[scenario_key + '_' + sec + '_' + str(rank) + 're'] = reliability
                    competition_results[scenario_key + '_' + sec + '_' + str(rank) + 'ka'] = kappa
                    competition_results[scenario_key + '_' + sec + '_' + str(rank) + 'pro'] = profit
                    competition_results[scenario_key + '_' + sec + '_' + str(rank) + 'sec_pro'] = sector_ave_profit
                    competition_results[scenario_key + '_' + sec + '_' + str(rank) + 'var_pro'] = var_profit
                    competition_results[scenario_key + '_' + sec + '_' + str(rank) + 'or'] = total_order

        df = pd.DataFrame(competition_results)
        df.to_excel(results_save_path, sheet_name='results', index=False)
        df = pd.DataFrame(init_x_info)
        with pd.ExcelWriter(results_save_path, engine="openpyxl", mode='a') as writer:
            df.to_excel(writer, sheet_name='init_x_info', index=False)
        df = pd.DataFrame(CHI_scores)
        with pd.ExcelWriter(results_save_path, engine="openpyxl", mode='a') as writer:
            df.to_excel(writer, sheet_name='CHI_scores', index=False)
        df = pd.DataFrame(CHI_levels)
        with pd.ExcelWriter(results_save_path, engine="openpyxl", mode='a') as writer:
            df.to_excel(writer, sheet_name='CHI_levels', index=False)
        df = pd.DataFrame(corr_max_producer_name)
        with pd.ExcelWriter(results_save_path, engine="openpyxl", mode='a') as writer:
            df.to_excel(writer, sheet_name='large_reg', index=False)
        df = pd.DataFrame(vartheta_info)
        with pd.ExcelWriter(results_save_path, engine="openpyxl", mode='a') as writer:
            df.to_excel(writer, sheet_name='vartheta', index=False)

    def output_change(self, rho=0.01, varphi=0.7):
        main_results_file = self.eco_epi_model.hyper_matlab_results_path + '/' + 'main_rho_' + str(rho) + '_varphi_' \
                            + str(varphi)
        results_save_path = main_results_file + '/output_info.xlsx'
        output_info = {}
        for keep_time_step in [52 * 7 * 1, 52 * 7 * 1.5, 52 * 7 * 2]:
            scenario_key, eco_epi_data = self.eco_epi_model.load_economic_data(rho, varphi, self.NPI_policy_scenario,
                                                                          keep_time_step, if_check_economic_data=True)
            actual_output = np.array(eco_epi_data['value_added_actual']) / self.eco_epi_model.q_L[None, :]
            max_output = np.array(eco_epi_data['value_added_max']) / self.eco_epi_model.q_L[None, :]
            total_order = np.array(eco_epi_data['total_order'])

            actual_output_by_country = self.eco_epi_model.cal_sum_by_regions_for_time_series_arr(actual_output)
            max_output_by_country = self.eco_epi_model.cal_sum_by_regions_for_time_series_arr(max_output)
            total_order_by_country = self.eco_epi_model.cal_sum_by_regions_for_time_series_arr(total_order)

            for country_class in self.economic_region_npi_policy:
                corr_regions_index_list = self.economic_region_npi_policy[country_class]
                actual_output_by_country_class = np.mean(actual_output_by_country[:, corr_regions_index_list], axis=1)
                max_output_by_country_class = np.mean(max_output_by_country[:, corr_regions_index_list], axis=1)
                total_order_by_country_class = np.mean(total_order_by_country[:, corr_regions_index_list], axis=1)

                output_info[scenario_key + '_actual_' + country_class] = actual_output_by_country_class
                output_info[scenario_key + '_max_' + country_class] = max_output_by_country_class
                output_info[scenario_key + '_total_' + country_class] = total_order_by_country_class

        df = pd.DataFrame(output_info)
        df.to_excel(results_save_path, index=False)

    @staticmethod
    def risk_diversity_cal(before_distribution_, CHI_region_):
        # before_distribution = np.array([0.1, 0.9])
        # CHI_region = np.array([1, 10])
        weighted_x = before_distribution_ * CHI_region_
        return -scipy.stats.entropy(weighted_x / sum(weighted_x))

    def risk_diversity_sample(self):
        results_save_path = self.eco_epi_model.hyper_matlab_results_path + '/risk_diversity_sample.xlsx'
        risk_diversity_info = []
        production_distribution_info = [np.array([0.5, 0.5]), np.array([0.3, 0.7]), np.array([0.1, 0.9])]
        CHI_info = [np.array([10, 10]), np.array([1, 10]), np.array([10, 1])]
        for CHI_region in CHI_info:
            risk_diversity_info_temp = []
            for before_distribution in production_distribution_info:
                risk_diversity_info_temp.append(self.risk_diversity_cal(before_distribution, CHI_region))
            risk_diversity_info.append(risk_diversity_info_temp)
        df = pd.DataFrame(risk_diversity_info)
        df.to_excel(results_save_path, sheet_name='risk', index=False, header=None)
        df = pd.DataFrame([info.tolist() for info in production_distribution_info])
        with pd.ExcelWriter(results_save_path, engine="openpyxl", mode='a') as writer:
            df.to_excel(writer, sheet_name='production', index=False, header=None)
        df = pd.DataFrame([info.tolist() for info in CHI_info])
        with pd.ExcelWriter(results_save_path, engine="openpyxl", mode='a') as writer:
            df.to_excel(writer, sheet_name='CHI', index=False, header=None)

    def all_experiments(self, rho, varphi):
        self.main_results(rho=rho, varphi=varphi)
        self.market_structure_change(rho=rho, varphi=varphi)
        self.competition_explain(rho=rho, varphi=varphi)
        self.output_change(rho=rho, varphi=varphi)
        self.regression(rho=rho, varphi=varphi)


# main results
exp_main = FinalExps()
exp_main.all_experiments(rho=0.01, varphi=0.7)
exp_main.risk_diversity_sample()
exp_main.specify_timing(rho=0.01, varphi=0.7, timing_range=range(52 * 7 * 1, 52*2 * 7 * +1, 4 * 7))
exp_main.get_stringency_info_csv()

# supplementary results
# change rho
# rho=0.002
exp_supple_change_phi = FinalExps()
exp_supple_change_phi.main_results(rho=0.002)
# rho=0.02
exp_supple_change_phi = FinalExps()
exp_supple_change_phi.main_results(rho=0.02)
# rho=0.05
exp_supple_change_phi = FinalExps()
exp_supple_change_phi.main_results(rho=0.05)

# change varphi
# varphi=0.3
exp_supple_change_phi = FinalExps()
exp_supple_change_phi.main_results(varphi=0.3)
# varphi=0.9
exp_supple_change_phi = FinalExps()
exp_supple_change_phi.main_results(varphi=0.9)

# change phi
# phi=0.005
exp_supple_change_phi = FinalExps(phi=0.005)
exp_supple_change_phi.main_results()

# change I_threshold
exp_supple_change_I_fraction_threshold = FinalExps(I_fraction_threshold=1e-1)
exp_supple_change_I_fraction_threshold.main_results()

