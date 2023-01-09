import numpy as np
import pandas as pd
import math
import datetime
import json
import os
import networkx as nx
import pickle


def raw_data_from_owid(t0_date_owid):
    owid = pd.read_csv('https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv')
    owid = owid[owid['date'] == t0_date_owid]
    People_fully_vaccinated = {}
    Population = {}

    for _, data in owid.iterrows():
        country = data['iso_code']
        if len(country) == 3 and not np.isnan(data['population']):
            Population[country] = data['population']
            People_fully_vaccinated[country] = math.floor(data['total_vaccinations'] / 2)
    return People_fully_vaccinated, Population


def datetime_delete_owid_format(end_date, n_days_ago):
    date_n_days_ago = datetime.datetime.strptime(end_date, '%Y-%m-%d') - datetime.timedelta(days=n_days_ago)
    date_n_days_ago = date_n_days_ago.strftime('%Y-%m-%d')
    return date_n_days_ago


def get_npi_stringency_of_target_date(owid, target_date_):
    info_before_target_date = owid[owid['date'] <= target_date_]
    npi_stringency_before_target_date = info_before_target_date[['iso_code', 'date', 'stringency_index']]
    npi_stringency_before_target_date_dict = {}
    for _, row in npi_stringency_before_target_date.iterrows():
        if row['stringency_index'] == row['stringency_index']:
            if row['iso_code'] not in npi_stringency_before_target_date_dict:
                npi_stringency_before_target_date_dict[row['iso_code']] = (row['stringency_index'],
                                                                           row['date'])
            else:
                if row['date'] > npi_stringency_before_target_date_dict[row['iso_code']][1]:
                    npi_stringency_before_target_date_dict[row['iso_code']] = (row['stringency_index'],
                                                                               row['date'])
    return {k: npi_stringency_before_target_date_dict[k][0] for k in npi_stringency_before_target_date_dict}

def cal_days_between_dates(date_start_owid, date_end_owid):
    d_end = datetime.datetime.strptime(date_end_owid, '%Y-%m-%d')
    d_start = datetime.datetime.strptime(date_start_owid, '%Y-%m-%d')
    delta = d_end - d_start
    return delta.days

def get_npi_stringency_data_range(owid, start_date_, end_date_):
    info_date_range = owid[owid['date'] <= end_date_]
    info_date_range = info_date_range[info_date_range['date'] >= start_date_]
    date_length = cal_days_between_dates(start_date_, end_date_)+1
    npi_stringency_date_range_row = info_date_range[['iso_code', 'date', 'stringency_index']]

    npi_stringency_date_range_dict = {}
    for _, row in npi_stringency_date_range_row.iterrows():
        date_index = cal_days_between_dates(start_date_, row['date'])
        if row['iso_code'] not in npi_stringency_date_range_dict:
            npi_stringency_date_range_dict[row['iso_code']] = [0]*date_length

        npi_stringency_date_range_dict[row['iso_code']][date_index] = row['stringency_index']

    return npi_stringency_date_range_dict, date_length


def replace_nan(list_):
    if list_[0] != list_[0]:
        return list_
    else:
        replaces_list = [list_[0]]
        for i in range(1, len(list_)):
            if list_[i] == list_[i]:
                replaces_list.append(list_[i])
            else:
                replaces_list.append(replaces_list[i-1])
        return replaces_list


def get_epidemic_data_for_specific_date(save_dir, target_t0_date_owid_, infectious_period_in_days_,
                                        projection_start_date_owid_):
    active_case_cal_date_begin = datetime_delete_owid_format(target_t0_date_owid_, infectious_period_in_days_)
    owid = pd.read_csv('https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv')
    owid_target_date = owid[owid['date'] == target_t0_date_owid_]
    owid_data_on_active_case_cal_date_begin = owid[owid['date'] == active_case_cal_date_begin]
    info_before_target_date = owid[owid['date'] <= target_t0_date_owid_]
    vaccination_data_before_target_date = info_before_target_date[['iso_code', 'date',
                                                                   'total_vaccinations']]

    total_cases_on_active_case_cal_date_begin = {}
    for _, row in owid_data_on_active_case_cal_date_begin.iterrows():
        total_cases_on_active_case_cal_date_begin[row['iso_code']] = row['total_cases']

    vaccination_data_before_target_date_dict = {}
    for _, row in vaccination_data_before_target_date.iterrows():
        if row['total_vaccinations'] == row['total_vaccinations']:
            if row['iso_code'] not in vaccination_data_before_target_date_dict:
                vaccination_data_before_target_date_dict[row['iso_code']] = (row['total_vaccinations'],
                                                                             row['date'])
            else:
                if row['date'] > vaccination_data_before_target_date_dict[row['iso_code']][1]:
                    vaccination_data_before_target_date_dict[row['iso_code']] = (row['total_vaccinations'],
                                                                                 row['date'])

    npi_info_date_range, date_length = \
        get_npi_stringency_data_range(owid, target_t0_date_owid_, projection_start_date_owid_)

    People_fully_vaccinated = {}
    total_cases = {}
    total_deaths = {}
    active_cases = {}
    NPI_stringency_index_range = {}
    for _, data in owid_target_date.iterrows():
        country = data['iso_code']
        if country in vaccination_data_before_target_date_dict:
            People_fully_vaccinated[country] = \
                math.ceil((vaccination_data_before_target_date_dict[country][0]) / 3)
        if country in npi_info_date_range:
            NPI_stringency_index_range[country] = replace_nan(npi_info_date_range[country])
        else:
            NPI_stringency_index_range[country] = [np.nan]*date_length
        total_cases[country] = data['total_cases']
        total_deaths[country] = data['total_deaths']
        active_cases[country] = total_cases[country] - total_cases_on_active_case_cal_date_begin[country]

    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    with open(save_dir + "/fully_vaccinated.json", "w") as f_:
        json.dump(People_fully_vaccinated, f_)
    with open(save_dir + "/total_cases.json", "w") as f_:
        json.dump(total_cases, f_)
    with open(save_dir + "/total_deaths.json", "w") as f_:
        json.dump(total_deaths, f_)
    with open(save_dir + "/active_cases.json", "w") as f_:
        json.dump(active_cases, f_)
    with open(save_dir + "/NPI_stringency_index.json", "w") as f_:
        json.dump(NPI_stringency_index_range, f_)


def read_raw_epidemic_data(target_t0_date_owid_, infectious_period_in_days_, projection_start_date_owid_):
    save_dir = 'raw_epidemic_data/' + target_t0_date_owid_ + '_IP_' + str(infectious_period_in_days_) + \
               '_' + str(projection_start_date_owid_)
    if not os.path.exists(save_dir):
        get_epidemic_data_for_specific_date(save_dir, target_t0_date_owid_, infectious_period_in_days_,
                                            projection_start_date_owid_)
    with open(save_dir + "/fully_vaccinated.json", "r") as f_:
        fully_vaccinated_ = json.load(f_)
    with open(save_dir + "/total_cases.json", "r") as f_:
        total_cases_ = json.load(f_)
    with open(save_dir + "/total_deaths.json", "r") as f_:
        total_deaths_ = json.load(f_)
    with open(save_dir + "/active_cases.json", "r") as f_:
        active_cases_ = json.load(f_)
    with open(save_dir + "/NPI_stringency_index.json", "r") as f_:
        NPI_stringency_index_ = json.load(f_)
    return fully_vaccinated_, total_cases_, total_deaths_, active_cases_, NPI_stringency_index_


def raw_dict_to_array(dict_, common):
    array = np.ones(len(common)) * np.nan
    for c in dict_:
        if c in common:
            array[common.index(c)] = dict_[c]
    return array


def complementary_popu_info():
    # iso code country name mapping source: https://www.nationsonline.org/oneworld/country_code_list.htm
    complementary_popu = {
        'ERI': 3645280,  # Eritrea
        'GGY': 63590,  # Guernsey
        'JEY': 103267,  # Jersey
        'SPM': 5794,  # Saint Pierre and Miquelon
        'ESH': 597339,  # Western Sahara
        'FLK': 3480,  # Falkland Islands
        'SHN': 6077,  # Saint Helena
        'AIA': 15003,  # Anguilla
        'GUF': 298682,  # French Guiana
        'MYT': 272815,  # Mayotte
        'COK': 17564,  # Cook Islands
        'MSR': 4992,  # Montserrat
        'TWN': 23816775,  # Taiwan
        'NIU': 1626,  # Niue
        'WLF': 11239  # Wallis and Futuna Islands / Wallis & Futuna
    }
    return complementary_popu


def get_actual_daily_mobility_matrix(comm):
    G_air = nx.read_gexf('raw_data/G_air_complete.gexf')
    update_G_air = nx.DiGraph(G_air.subgraph(comm))
    update_G_air.remove_edges_from(nx.selfloop_edges(update_G_air))  # remove self loop
    G_ij = np.array(nx.adjacency_matrix(update_G_air, comm, weight='weight').todense())
    G_ij_half = (np.tril(G_ij, -1) + np.transpose(np.triu(G_ij, 1))) / 2
    G_ij = G_ij_half + G_ij_half.T
    return G_ij / 365


def generate_country_info(hyper_file_path, target_t0_date_owid, infectious_period_in_days,  projection_start_date_owid,
                          CHI_thre_mild, CHI_thre_moderate):
    CHI_date_length = cal_days_between_dates(target_t0_date_owid, projection_start_date_owid)+1
    # target_t0_date_owid, infectious_period_in_days
    with open("processed_key_paras/final_considered_regs.json", "r") as f:
        comm_regs = json.load(f)

    # ------ popu -------------------------------------
    # complementary_popu_info()
    with open("raw_data/raw_popu_data_2020.json", "r") as f:
        raw_popu_data_2020_dict = json.load(f)
    raw_popu_data_2020_dict.update(complementary_popu_info())
    # --------------------------------------------------------

    # ----- epidemic ----------------------------------
    fully_vaccinated, total_cases, total_deaths, active_cases, NPI_stringency_index = \
        read_raw_epidemic_data(target_t0_date_owid, infectious_period_in_days,  projection_start_date_owid)
    # --------------------------------------------------------

    # ----- NPI ---------------------------------------
    stringency_info = [np.nan] * len(comm_regs)
    stringency_index_range_epi = np.zeros((len(comm_regs), CHI_date_length))
    for c in comm_regs:
        c_index = comm_regs.index(c)
        if c in NPI_stringency_index and NPI_stringency_index[c][0] == NPI_stringency_index[c][0]:
            stringency_index_range_epi[comm_regs.index(c)] = np.array(NPI_stringency_index[c])
            if NPI_stringency_index[c][0] < CHI_thre_mild:
                stringency_info[c_index] = 'mild'
            elif NPI_stringency_index[c][0] < CHI_thre_moderate:
                stringency_info[c_index] = 'moderate'
            else:
                stringency_info[c_index] = 'stringent'
        else:
            stringency_index_range_epi[comm_regs.index(c)] = \
                np.array([10] * (cal_days_between_dates(target_t0_date_owid, projection_start_date_owid)+1))
            stringency_info[c_index] = 'mild'
    # --------------------------------------------------------

    country_info = {
        'common': comm_regs,
        'G': get_actual_daily_mobility_matrix(comm_regs),
        'N': raw_dict_to_array(raw_popu_data_2020_dict, comm_regs),
        'stringency_info_current': np.array(stringency_info),
        'stringency_index_current': stringency_index_range_epi,
    }

    IS_array = raw_dict_to_array(active_cases, comm_regs)
    IS_array[np.isnan(IS_array)] = 0
    country_info['IS'] = IS_array

    D_array = raw_dict_to_array(total_deaths, comm_regs)
    D_array[np.isnan(D_array)] = 0
    country_info['D'] = D_array

    R_array = raw_dict_to_array(total_cases, comm_regs) - D_array - IS_array
    R_array[np.isnan(R_array)] = 0
    R_array = np.maximum(R_array, np.zeros_like(R_array))
    country_info['R'] = R_array

    V_array = raw_dict_to_array(fully_vaccinated, comm_regs)
    V_array[np.isnan(V_array)] = 0
    country_info['V_raw'] = V_array

    V_array = np.minimum(np.array(country_info['N']) - (IS_array + R_array + D_array), V_array)
    country_info['V'] = V_array

    save_path = hyper_file_path + '/clean_epidemic_data.pkl'
    with open(save_path, 'wb') as f:
        pickle.dump(country_info, f)
    # --------------------------------------------------------


def read_country_info(hyper_file_path, target_t0_date_owid, infectious_period_in_days, projection_start_date_owid,
                      CHI_thre_mild, CHI_thre_moderate):
    save_path = hyper_file_path + '/clean_epidemic_data.pkl'
    if not os.path.exists(save_path):
        generate_country_info(hyper_file_path, target_t0_date_owid, infectious_period_in_days,
                              projection_start_date_owid, CHI_thre_mild, CHI_thre_moderate)
    with open(save_path, 'rb') as f:
        country_info = pickle.load(f)
    country_info['common'] = list(country_info['common'])
    return country_info


def get_lat_lon_eco():
    lat_lon = pd.read_csv('raw_data/countries_codes_and_lat_lon.csv')
    with open('processed_gtap/regions_list.pkl', 'rb') as f:
        regions = pickle.load(f)
    lat_lon_dict = {}

    for _, row in lat_lon.iterrows():
        if row['Alpha-3 code'][2:-1].lower() in regions:
            lat_lon_dict[row['Alpha-3 code'][2:-1].lower()] = {'lat': float(row['Latitude (average)'][2:-1]),
                                                               'lon': float(row['Longitude (average)'][2:-1])}
    return lat_lon_dict


def generate_init_parameters(MRIOT, days_in_each_eco_time_step, num_regions, num_commodities, num_all_items):
    save_dir = 'processed_gtap/time_step_length_' + str(days_in_each_eco_time_step)
    if os.path.exists(save_dir):
        Z = np.load(save_dir + '/Z.npy')
        y = np.load(save_dir + '/y.npy')
        v = np.load(save_dir + '/v.npy')
        x = np.load(save_dir + '/x.npy')
        for_cal_sum_by_product = np.load(save_dir + '/for_cal_sum_by_product.npy')
        Z_Dis_sum = np.load(save_dir + '/Z_Dis_sum.npy')
        F_Dis_sum = np.load(save_dir + '/F_Dis_sum.npy')
        sum_input_product = np.load(save_dir + '/sum_input_product.npy')
    else:
        os.mkdir(save_dir)
        MRIOT_array = MRIOT.to_numpy() / 365 * days_in_each_eco_time_step
        Z = MRIOT_array[:num_all_items, :num_all_items]
        y = MRIOT_array[:num_all_items, num_all_items:-1]
        v = MRIOT_array[-2, :num_all_items]
        x = MRIOT_array[:num_all_items, -1]
        for_cal_sum_by_product = np.tile(np.diag([1] * num_commodities), (num_regions, num_regions))
        Z_Dis_sum = np.dot(for_cal_sum_by_product, Z)
        F_Dis_sum = np.dot(for_cal_sum_by_product, y)
        sum_input_product = np.dot(for_cal_sum_by_product, Z)[:num_commodities, :]

        np.save(save_dir + '/Z.npy', Z)
        np.save(save_dir + '/y.npy', y)
        np.save(save_dir + '/v.npy', v)
        np.save(save_dir + '/x.npy', x)
        np.save(save_dir + '/for_cal_sum_by_product.npy', for_cal_sum_by_product)
        np.save(save_dir + '/Z_Dis_sum.npy', Z_Dis_sum)
        np.save(save_dir + '/F_Dis_sum.npy', F_Dis_sum)
        np.save(save_dir + '/sum_input_product.npy', sum_input_product)
    return {'Z': Z,
            'y': y,
            'v': v,
            'x': x,
            'for_cal_sum_by_product': for_cal_sum_by_product,
            'Z_Dis_sum': Z_Dis_sum,
            'F_Dis_sum': F_Dis_sum,
            'sum_input_product': sum_input_product}


def load_hyper():
    with open('processed_key_paras/sector_supply_multiplier.pkl', 'rb') as f:
        sector_supply_multiplier = pickle.load(f)
    with open('processed_key_paras/sector_demand_multiplier.pkl', 'rb') as f:
        sector_demand_multiplier = pickle.load(f)
    with open('processed_gtap/regions_list.pkl', 'rb') as f:
        regions = pickle.load(f)
    with open('processed_gtap/activities_list.pkl', 'rb') as f:
        activities = pickle.load(f)
    with open('processed_gtap/adjusted_MRIOT_array.pkl', 'rb') as f:
        MRIOT = pickle.load(f)

    # 'NPI_sector_multiplier':
    hyper_paras = \
        {'regions': regions,
         'activities': activities,
         'NPI_sector_multiplier': [sector_supply_multiplier[sector] for sector in activities],
         'demand_multiplier_sector': [sector_demand_multiplier[sector] for sector in activities],
         'eco_lat_lon': get_lat_lon_eco(),
         'days_in_each_eco_time_step': 7,
         'items_info': MRIOT.index.values[:-2],
         'T': 260,  # simulation length in weeks
         'projection_start_date_owid': '2022-11-11',
         }

    init_parameters = generate_init_parameters(MRIOT, hyper_paras['days_in_each_eco_time_step'],
                                               len(regions),
                                               len(activities),
                                               len(regions) * len(activities))
    hyper_paras.update(init_parameters)

    hyper_epi_paras = {
        "t0_date_owid": '2022-03-01',
        "active_case_cal_length": 7,  # active cases [days]
        # simulation length 5 years * 52 weeks/year * 7 days/week
        "tau": hyper_paras['T'] * hyper_paras['days_in_each_eco_time_step'],

        # https://data.worldbank.org/indicator/SP.DYN.CDRT.IN
        "upsilon": 7.7 / 1000 / 365,  # natural death rate 7.7/1000/365 from world bank average

        # https://www.nature.com/articles/s41591-022-01855-7
        "alpha": 1 / 5.6,  # 1/alpha is the infectious period
        "sigma": 1 / 1.2,  # 1/sigma is the incubation period

        # booster info
        # https://www.nature.com/articles/s41591-022-01753-y
        # https://www.nature.com/articles/s41467-022-30895-3
        # https://jamanetwork.com/journals/jama/fullarticle/2787609
        # https://www.bmj.com/content/376/bmj-2021-069761.long
        # https://www.medrxiv.org/content/10.1101/2022.03.22.22272769v2.full.pdf+html
        # https://www.nature.com/articles/s41591-022-01699-1
        "eta": 0.75,
        "epsilon": 0.98,

        # https://www.nature.com/articles/s41467-022-30895-3
        "varepsilon": 1 / (7 * 30),  # 1/\varepsilon is the duration of vaccinal immunity

        # https://www.sciencedirect.com/science/article/pii/S2213260021005592?via%3Dihub
        "R0": 10,

        # https://www.nature.com/articles/s41591-022-01855-7#Sec23
        "psi": 1 / 900,

        # https://www.nature.com/articles/s41591-022-01855-7#Sec15
        # https://www.covidvaccine.gov.hk/pdf/death_analysis.pdf
        "nu": 0.03,  # infection fatality rate

        "gamma_ref": 0.8,  # oag 2020 relative to 2019

        # https://www.nature.com/articles/s41591-022-01855-7
        "c_stringent_npi_ref": 0.8,
        "c_mild_npi_ref": 0.6,

        # CHI<CHI_thre_mild -> mild; CHI<CHI_thre_moderate -> moderate; else -> stringent
        "CHI_thre_mild": 25,
        "CHI_thre_moderate": 50,

    }
    hyper_paras.update(hyper_epi_paras)
    hyper_paras["CHI_ref_CHN_date_index"] = cal_days_between_dates(hyper_paras["t0_date_owid"], '2022-04-15')
    return hyper_paras
