import numpy as np
import pickle
import json
import pandas as pd
import eco_epi_hyper
import os
from collections import *


class DiscreteSVEIRD:
    """
    for disease transmission
    """

    def __init__(self, Re_thre, phi, k, hyper_file_path):

        hyper_paras = eco_epi_hyper.load_hyper()

        # -------------- self setting need sensitivity analysis -----
        self.Re_thre = Re_thre
        self.phi = phi
        self.k = k
        self.hyper_file_path = hyper_file_path
        # -----------------------------------------------------------

        # ---------------- not change ----------------------
        with open('processed_gtap/regions_list.pkl', 'rb') as f:
            self.gtap_regions_list = pickle.load(f)
        self.gtap_regions_list = [reg.upper() for reg in self.gtap_regions_list]
        with open("processed_key_paras/final_regs_gtap_map.json", "r") as f:
            self.final_gtap_details_region_dict = json.load(f)
        # --------------------------------------------------

        # ----------------- pre define ---------------------
        self.economic_time_step_length = hyper_paras['days_in_each_eco_time_step']
        self.projection_start_date = hyper_paras['projection_start_date_owid']
        self.c_stringent_npi_ref = hyper_paras['c_stringent_npi_ref']
        self.c_mild_npi_ref = hyper_paras['c_mild_npi_ref']
        self.gamma_ref = hyper_paras['gamma_ref']
        self.tau = hyper_paras['tau']  # simulation period
        self.upsilon = hyper_paras['upsilon']  # natural death rate
        self.beta = hyper_paras['R0'] * hyper_paras['alpha']  # transmissibility multiplied by number of contacts
        self.eta = hyper_paras["eta"]  # vaccine efficacy against infection
        self.epsilon = hyper_paras["epsilon"]  # vaccine efficacy against death
        self.nu = hyper_paras['nu']  # infection fatality rate
        self.sigma = hyper_paras["sigma"]  # 1/sigma incubation period
        self.alpha = hyper_paras['alpha']  # 1/alpha infectious period
        self.varepsilon = hyper_paras['varepsilon']  # rate losing vaccinal immunity
        self.psi = hyper_paras['psi']  # rate at which lossing natural immunity
        self.CHI_thre_mild = hyper_paras['CHI_thre_mild']
        self.CHI_thre_moderate = hyper_paras['CHI_thre_moderate']
        self.t0_date_owid, self.active_case_cal_length = \
            hyper_paras['t0_date_owid'], hyper_paras['active_case_cal_length']
        # --------------------------------------------------

        # ------------ read modified country info ----------
        country_info = \
            eco_epi_hyper.read_country_info(self.hyper_file_path, self.t0_date_owid, self.active_case_cal_length,
                                            self.projection_start_date,
                                            self.CHI_thre_mild, self.CHI_thre_moderate)
        self.G_before_pandemic = country_info['G']
        self.N0 = country_info['N']
        self.V0 = country_info['V']
        self.R0 = country_info['R']
        self.D0 = country_info['D']
        self.I0 = country_info['IS']
        self.S0 = self.N0 - self.V0 - self.I0 - self.R0 - self.D0
        self.V_raw0 = country_info['V_raw']
        self.common = country_info['common']
        self.CHI_current_known = country_info['stringency_index_current']
        self.CHI_current = self.CHI_current_known[:, 0]
        self.CHI_CHN_ref = self.CHI_current_known[self.common.index('CHN'), hyper_paras["CHI_ref_CHN_date_index"]]
        self.sum_popu = np.sum(self.N0)
        self.c_num = len(self.N0)
        self.stringency_info_current = country_info['stringency_info_current']
        self.init_stringency_class_details = defaultdict(list)
        for i in range(self.c_num):
            self.init_stringency_class_details[self.stringency_info_current[i]].append(i)
        # --------------------------------------------------

    @staticmethod
    def set_array_bounds(array_v, min_v, max_v):
        array_v[array_v > max_v] = max_v
        array_v[array_v < min_v] = min_v
        return array_v

    def transfer_CHI_to_c_bounds_and_gamma(self, CHI_value_t):
        c_max = self.set_array_bounds(CHI_value_t / self.CHI_CHN_ref * self.c_stringent_npi_ref, 0, 0.9)
        c_min = self.set_array_bounds(CHI_value_t / self.CHI_CHN_ref * self.c_mild_npi_ref, 0, 0.9)
        gamma = self.set_array_bounds(CHI_value_t / self.CHI_CHN_ref * self.gamma_ref, 0, 0.9)
        return c_max, c_min, gamma

    def cal_Re(self, c_level, S, V):
        return (1 - c_level) * self.beta / (self.alpha * self.N0) * ((1 - self.eta) * V + S)

    def save_stringency_class_to_csv(self, csv_path):
        stringency_dict = {'ISO_A3': 'CHI'}
        stringency_dict.update({self.common[i]: self.stringency_info_current[i] for i in range(self.c_num)})
        df = pd.DataFrame(stringency_dict, index=[0]).T
        df.to_csv(csv_path, header=False)

    def set_c_and_gamma_matrix(self, S, V, c_current, CHI_value_t):
        c_max, c_min, gamma = self.transfer_CHI_to_c_bounds_and_gamma(CHI_value_t)
        Re = self.cal_Re(c_current, S, V)
        model_gamma_each_country = np.tile(1 - gamma, (self.c_num, 1))
        new_c = c_min + (c_max - c_min) / (1 + np.exp((1 - Re) * self.k))
        new_gamma_matrix = np.minimum(model_gamma_each_country, model_gamma_each_country.T)
        return new_c, new_gamma_matrix

    @staticmethod
    def international_mobility_cal(E, all_mobility, could_move_nodes):
        E_frac = (E / np.tile(could_move_nodes, (len(E), 1))).T
        out_ = np.diag(np.sum(all_mobility, axis=1)).dot(E_frac)
        in_ = (all_mobility / np.tile(np.sum(all_mobility, axis=1), (all_mobility.shape[0], 1))).dot(out_)
        return out_, in_

    def define_country_npi_policy_class_in_economic(self):
        economic_region_npi_policy = {'world': [],
                                      'mild_npi': [],
                                      'moderate_npi': [],
                                      'stringent_npi': []}
        economic_region_chi_value = {}

        for r, country_code in enumerate(self.gtap_regions_list):
            economic_region_npi_policy['world'].append(r)

            if len(self.final_gtap_details_region_dict[country_code]) == 1:
                index_in_epidemic_model = self.common.index(self.final_gtap_details_region_dict[country_code][0])
                stringency_level = self.stringency_info_current[index_in_epidemic_model]
                economic_region_chi_value[country_code] = self.CHI_current[index_in_epidemic_model]
                if stringency_level == 'mild':
                    economic_region_npi_policy['mild_npi'].append(r)
                elif stringency_level == 'moderate':
                    economic_region_npi_policy['moderate_npi'].append(r)
                elif stringency_level == 'stringent':
                    economic_region_npi_policy['stringent_npi'].append(r)

            elif len(self.final_gtap_details_region_dict[country_code]) == 0:
                economic_region_npi_policy['mild_npi'].append(r)
                economic_region_chi_value[country_code] = 10

            else:
                corr_countries_list_index_in_epidemic_model = [self.common.index(detailed_country_code)
                                                               for detailed_country_code in
                                                               self.final_gtap_details_region_dict[country_code]]

                economic_region_chi_value[country_code] = \
                    np.mean(self.CHI_current[corr_countries_list_index_in_epidemic_model])
                stringency_level = self.get_country_stringency_class_given_CHI(economic_region_chi_value[country_code])


                if stringency_level == 'mild':
                    economic_region_npi_policy['mild_npi'].append(r)
                elif stringency_level == 'moderate':
                    economic_region_npi_policy['moderate_npi'].append(r)
                elif stringency_level == 'stringent':
                    economic_region_npi_policy['stringent_npi'].append(r)

        with open(self.hyper_file_path + '/economic_region_npi_policy.json', 'w') as f:
            json.dump(economic_region_npi_policy, f)
        with open(self.hyper_file_path + '/economic_region_chi_value.json', 'w') as f:
            json.dump(economic_region_chi_value, f)
        # -------------
        return economic_region_npi_policy, economic_region_chi_value

    def trans_SVEIRD(self, NPI_policy_scenario, keep_time_step):
        self.active_cases_num = []
        self.cumulative_cases_num = []
        self.cumulative_deaths_num = []
        self.dynamic_c_info = []
        self.population_info = []
        self.vaccinated_num = []
        self.mortality_num = []
        self.unable_to_work_popu = []
        self.recovered_num = []
        self.vaccination_coverage_num = []

        S = self.S0.copy()
        V = self.V0.copy()
        ES = np.zeros(self.c_num)
        EV = np.zeros(self.c_num)
        IS = self.I0.copy()
        IV = np.zeros(self.c_num)
        R = self.R0.copy()
        D = self.D0.copy()
        N = self.N0.copy()
        vaccine_coverage = self.V_raw0.copy()
        c_current = np.ones(self.c_num) * self.c_mild_npi_ref



        cumulative_cases_temp = IS + IV + R + D
        self.active_cases_num.append((IS + IV).tolist())
        self.cumulative_cases_num.append(cumulative_cases_temp.tolist())
        self.cumulative_deaths_num.append(D.tolist())
        self.dynamic_c_info.append(c_current.tolist())
        self.population_info.append(N.tolist())
        self.vaccinated_num.append(V.tolist())
        self.mortality_num.append([0] * len(D))
        self.unable_to_work_popu.append((IS + IV + D).tolist())
        self.recovered_num.append(R.tolist())
        self.vaccination_coverage_num.append(vaccine_coverage.tolist())

        for t in range(1, self.tau+1):
            print('pure epidemic data cal----t=', t, '--------')
            if NPI_policy_scenario == 'keep_curr':
                if t > keep_time_step:
                    CHI_value_t = np.zeros(self.c_num)
                else:
                    if t < len(self.CHI_current_known[0]):
                        CHI_value_t = self.CHI_current_known[:, t]
                    else:
                        CHI_value_t = self.CHI_current_known[:, -1]
            c_current, gamma_matrix_current \
                = self.set_c_and_gamma_matrix(S, V, c_current, CHI_value_t)

            all_mobility = gamma_matrix_current * self.G_before_pandemic
            international_out_in = \
                self.international_mobility_cal(np.array([S, V, ES, EV, IS, IV, R]), all_mobility, N - D)

            I_sum = IS + IV
            Lambda = self.upsilon * (N - D)  # new birth
            S_ES_trans = I_sum * S / N * (1 - c_current) * self.beta
            S_V_trans = np.minimum(S, self.phi * N)
            V_S_trans = self.varepsilon * V
            V_EV_trans = I_sum * V / N * (1 - c_current) * self.beta * (1 - self.eta)
            ES_IS_trans = self.sigma * ES
            EV_IV_trans = self.sigma * EV
            IS_R_trans = self.alpha * IS * (1 - self.nu)
            IS_D_trans = self.alpha * IS * self.nu
            IV_R_trans = self.alpha * IV * (1 - self.nu * (1 - self.epsilon))
            IV_D_trans = self.alpha * IV * self.nu * (1 - self.epsilon)
            R_S_trans = R * self.psi


            dS = Lambda - S_V_trans - S_ES_trans + V_S_trans + R_S_trans - self.upsilon * S \
                 - international_out_in[0][:, 0] + international_out_in[1][:, 0]

            dV = S_V_trans - V_EV_trans - V_S_trans - self.upsilon * V \
                 - international_out_in[0][:, 1] + international_out_in[1][:, 1]

            dES = S_ES_trans - ES_IS_trans - self.upsilon * ES \
                  - international_out_in[0][:, 2] + international_out_in[1][:, 2]

            dEV = V_EV_trans - EV_IV_trans - self.upsilon * EV \
                  - international_out_in[0][:, 3] + international_out_in[1][:, 3]

            dIS = ES_IS_trans - IS_R_trans - IS_D_trans - self.upsilon * IS \
                  - international_out_in[0][:, 4] + international_out_in[1][:, 4]

            dIV = EV_IV_trans - IV_R_trans - IV_D_trans - self.upsilon * IV \
                  - international_out_in[0][:, 5] + international_out_in[1][:, 5]

            dR = IS_R_trans + IV_R_trans - R_S_trans - self.upsilon * R \
                 - international_out_in[0][:, 6] + international_out_in[1][:, 6]

            dD = IS_D_trans + IV_D_trans

            S = S + dS
            V = V + dV
            ES = ES + dES
            EV = EV + dEV
            IS = IS + dIS
            IV = IV + dIV
            R = R + dR
            D = D + dD
            N = S + V + ES + EV + IS + IV + R + D
            vaccine_coverage += S_V_trans

            IS_temp = IS.copy()
            IS_temp[IS_temp < 1] = 0

            IV_temp = IV.copy()
            IV_temp[IV_temp < 1] = 0

            ES_temp = ES.copy()
            ES_temp[ES_temp < 1] = 0

            EV_temp = EV.copy()
            EV_temp[EV_temp < 1] = 0


            cumulative_cases_temp = cumulative_cases_temp + ES_IS_trans + EV_IV_trans
            self.active_cases_num.append((IS + IV).tolist())
            self.cumulative_cases_num.append(cumulative_cases_temp.tolist())
            self.cumulative_deaths_num.append(D.tolist())
            self.dynamic_c_info.append(c_current.tolist())
            self.population_info.append(N.tolist())
            self.vaccinated_num.append(V.tolist())
            self.mortality_num.append(dD.tolist())
            self.unable_to_work_popu.append((IS + IV + D).tolist())
            self.recovered_num.append(R.tolist())
            self.vaccination_coverage_num.append(vaccine_coverage.tolist())


    def get_country_stringency_class_given_CHI(self, CHI_value):
        if CHI_value < self.CHI_thre_mild:
            return 'mild'
        elif CHI_value < self.CHI_thre_moderate:
            return 'moderate'
        else:
            return 'stringent'

    def get_country_index_list_given_country_class(self, country_class, **kwargs):
        if country_class == 'world':
            return list(range(self.c_num))
        elif country_class in ['mild', 'moderate', 'stringent']:
            return self.init_stringency_class_details[country_class]
        elif country_class == 'self_define':
            return [self.common.index(c_iso) for c_iso in kwargs['country_iso_list']]

    def get_result_given_gtap_country_mapping_and_time_step(self,
                                                            result_array,
                                                            if_divided_by_popu=True):
        mean_results = []

        needed_time_step_list = [0] + [temp * self.economic_time_step_length - 1
                                       for temp in range(1, int(self.tau / self.economic_time_step_length) + 1)]
        for country_code in self.gtap_regions_list:
            if len(self.final_gtap_details_region_dict[country_code]) == 1:
                corr_result = result_array[needed_time_step_list,
                                           self.common.index(self.final_gtap_details_region_dict[country_code][0])]
                if if_divided_by_popu:
                    mean_results.append(corr_result / self.N0[self.common.index(
                        self.final_gtap_details_region_dict[country_code][0])])
                else:
                    mean_results.append(corr_result)
            elif len(self.final_gtap_details_region_dict[country_code]) == 0:

                countries_adopting_mild_npi = self.get_country_index_list_given_country_class('mild')
                corr_result = result_array[needed_time_step_list, :][:, countries_adopting_mild_npi]
                if if_divided_by_popu:
                    mean_results.append(np.mean(corr_result / self.N0[countries_adopting_mild_npi][None, :], axis=1))
                else:
                    mean_results.append(np.mean(corr_result, axis=1))

            else:
                corr_countries_list = [self.common.index(detailed_country_code)
                                       for detailed_country_code in self.final_gtap_details_region_dict[country_code]]
                corr_result = result_array[needed_time_step_list, :][:, corr_countries_list]
                if if_divided_by_popu:
                    mean_results.append(np.mean(corr_result / self.N0[corr_countries_list][None, :], axis=1))
                else:
                    mean_results.append(np.mean(corr_result, axis=1))
        return np.transpose(np.array(mean_results))

    def get_mean_result_for_specific_countries(self, result_array, country_class,
                                               average_type='ave_country_frac_with_popu', **kwargs):

        country_index_list = self.get_country_index_list_given_country_class(country_class, **kwargs)
        specific_country_result = result_array[:, country_index_list]
        if average_type == 'ave_country_frac_with_popu':
            specific_country_popu = self.N0[country_index_list]

            specific_country_frac = specific_country_result / specific_country_popu[None, :]
            return np.mean(specific_country_frac, axis=1)
        elif average_type == 'ave_total_with_popu':
            specific_country_popu = self.N0[country_index_list]
            return np.sum(specific_country_result, axis=1) / np.sum(specific_country_popu[None, :], axis=1)
        elif average_type == 'ave_country':

            return np.mean(specific_country_result, axis=1)

    def trans_length(self, li):

        li = list(li)
        li.extend([0] * (int(self.tau / self.economic_time_step_length) - len(li)))
        return li

    # ----- trans length for CUM Values -------------
    def trans_length1(self, li):

        li = list(li)
        li.extend([li[-1]] * (int(self.tau / self.economic_time_step_length) - len(li)))
        return li

    # ------------------------------------------------

    @staticmethod
    def trans_length_matrix(mat, target_length):
        if len(mat) >= target_length:
            return mat
        mat = np.array(mat)
        return np.concatenate([mat, np.array([[0] * len(mat[0])] * (target_length - len(mat)))], axis=0)


    @staticmethod
    def trans_length_matrix1(mat, target_length):
        if len(mat) >= target_length:
            return mat
        # matrix size[0] -> target_length
        mat = np.array(mat)
        return np.concatenate([mat, np.array([list(mat[-1])] * (target_length - len(mat)))], axis=0)



    def load_epidemic_data(self, NPI_policy_scenario, keep_time_step):
        scenario_key = \
            NPI_policy_scenario + '_' + str(keep_time_step / (52 * 7))
        save_path = self.hyper_file_path + '/epidemic_data_' + scenario_key + '.json'
        if os.path.exists(save_path):
            print('------ found pure epidemic data -----------')
            with open(save_path, 'r') as f:
                epidemic_data = json.load(f)
        else:
            print('------ not found pure epidemic data -----------')
            self.trans_SVEIRD(NPI_policy_scenario, keep_time_step)
            epidemic_data = {'active_cases_num': self.active_cases_num,
                             'cumulative_cases_num': self.cumulative_cases_num,
                             'cumulative_deaths_num': self.cumulative_deaths_num,
                             'vaccinated_num': self.vaccinated_num,
                             'mortality_num': self.mortality_num,
                             'unable_to_work_popu': self.unable_to_work_popu,
                             'reduction_in_contacts': self.dynamic_c_info,
                             'recovered_num': self.recovered_num,
                             'vaccination_coverage_num': self.vaccination_coverage_num}
            with open(save_path, 'w') as f:
                json.dump(epidemic_data, f)
        return epidemic_data

    def load_epidemic_data_for_economic_cal(self, NPI_policy_scenario, keep_time_step):
        scenario_key = NPI_policy_scenario + '_' + str(keep_time_step / (52 * 7))
        save_path = self.hyper_file_path + '/epid_data_for_eco_cal_' + scenario_key + '.json'
        if os.path.exists(save_path):
            print('--------- found epidemic data for economic cal ---------------------')
            with open(save_path, 'r') as f:
                needed_epidemic_data_in_economic_model = json.load(f)
        else:
            print('--------- not found epidemic data for economic cal ---------------------')
            needed_economic_cal_info = ['unable_to_work_popu', 'active_cases_num', 'reduction_in_contacts']
            key_in_transfer_result = {'unable_to_work_popu': 'unavailable_l_frac',
                                      'active_cases_num': 'I_frac',
                                      'reduction_in_contacts': 'reduction_in_contact'}
            cal_rule_if_divided_by_popu = {'unable_to_work_popu': True,
                                           'active_cases_num': True,
                                           'reduction_in_contacts': False}

            pure_epidemic_data = self.load_epidemic_data(NPI_policy_scenario, keep_time_step)
            needed_epidemic_data_in_economic_model = {}
            for key in needed_economic_cal_info:
                needed_epidemic_data_in_economic_model[key_in_transfer_result[key]] = \
                    self.get_result_given_gtap_country_mapping_and_time_step(np.array(pure_epidemic_data[key]),
                                                                             cal_rule_if_divided_by_popu[key]).tolist()
            with open(save_path, 'w') as f:
                json.dump(needed_epidemic_data_in_economic_model, f)
        return needed_epidemic_data_in_economic_model

    def load_epidemic_data_for_economic_cal_simplify(self, NPI_policy_scenario, keep_time_step):
        needed_economic_cal_info = ['unable_to_work_popu', 'active_cases_num', 'reduction_in_contacts']
        key_in_transfer_result = {'unable_to_work_popu': 'unavailable_l_frac',
                                  'active_cases_num': 'I_frac',
                                  'reduction_in_contacts': 'reduction_in_contact'}
        cal_rule_if_divided_by_popu = {'unable_to_work_popu': True,
                                       'active_cases_num': True,
                                       'reduction_in_contacts': False}

        self.trans_SVEIRD(NPI_policy_scenario, keep_time_step)
        pure_epidemic_data = {'active_cases_num': self.active_cases_num,
                              'unable_to_work_popu': self.unable_to_work_popu,
                              'reduction_in_contacts': self.dynamic_c_info}

        needed_epidemic_data_in_economic_model = {}
        for key in needed_economic_cal_info:
            needed_epidemic_data_in_economic_model[key_in_transfer_result[key]] = \
                self.get_result_given_gtap_country_mapping_and_time_step(np.array(pure_epidemic_data[key]),
                                                                         cal_rule_if_divided_by_popu[key]).tolist()
        return needed_epidemic_data_in_economic_model
