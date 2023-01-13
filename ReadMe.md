

# Supply chain loss from easing COVID-19 restrictions: an evolutionary epidemiological-economic modeling study

1 Modeling trade flows based on ARIO_modify.py

    parameter settings from eco_epi_hyper.load_hyper()
    
    initializing the epidemic model based on parameter settings [SVEIRD_meta.py]
    
    main function: load_epidemic_data_for_economic_cal()
    
        setting global reopening scenarios and calculating epidemic dynamics and trade flows
        
        Modeling epidemic dynamics based on SVEIRD_meta.py
        
            parameter settings from eco_epi_hyper.load_hyper()
            
            initial epidemic situations from eco_epi_hyper.read_country_info() [our world in data]
            
                raw data is automatically stored in raw_epidemic_data
                
                processed data used in the economic model is automatically stored in results_for_economic_cal/clean_epidemic_data.pkl
                
                epidemic data for a specific reopening scenario is automatically stored in results_for_economic_cal
                
            calculating corresponding trade flows based on epidemic data
            
            output data is stored in results_for_economic_cal

2 all experiments 

    all from exps_general.py
    
    getting simulation results from results_for_economic_cal
    
    processing data for visualization
    
    processed data is automatically stored in matlab_files/results


3 visualization 

using Matlab

summarizing in exp_all.m

4 important files 

    processed_key_paras
    
        sector_supply_multiplier.pkl -> impact-to-labour multiplier (\vartheta^P)
        
        sector_demand_multiplier.pkl -> impact-to-demand multiplier (\vartheta^D)
        
        final_regs_gtap_map.json -> epidemic model and economic model region code mapping
        
        final_considered_regs.json -> country list in the final model

    processed_gtap
    
        regions_list.pkl -> regions in GTAP dataset
        
        activities_list.pkl -> commodities and services in GTAP dataset
        
        adjusted_MRIOT_array.pkl -> raw MRIOT (sample data can be generated from generate_sample_data.py generate_sample_MRIOT())
        
        time_step_length_7 -> dir for specific time steps (weekly)
        
            Z.npy -> Z in MRIOT
            
            y.npy -> y in MRIOT
            
            v.npy -> v in MRIOT
            
            x.npy -> x in MRIOT
            
            for_cal_sum_by_product.npy -> auxiliary parameter in the economic model
            
            Z_Dis_sum.npy -> auxiliary parameter in the economic model
            
            F_Dis_sum.npy -> auxiliary parameter in the economic model
            
            sum_input_product.npy -> auxiliary parameter in the economic model

    raw_data
    
        countries_codes_and_lat_lon.csv -> latitude and longitude of countries/regions
        
        raw_popu_data_2020.json -> population data
        
        G_air_complete.gexf -> traffic data from OAG (sample data)







