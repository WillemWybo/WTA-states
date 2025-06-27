#imports
import yaml
import numpy as np
from yaml_io import *

def nest_parameters_preparation(times, config, is_verbose, nest_pms):
    print("in nest_parameters_preparation: verbose mode is", is_verbose)
    
    # Extract config subsets of params 
    brain_state = config['brain_state']
    use_nestml = config['use_nestml']
    use_single_compartment_environment = config['use_single_compartment_environment']
    exc_pms_file_name = config['exc_neu_params_filename']
    inh_pms_file_name = config['inh_neu_params_filename']
    network = config['network']
    if use_single_compartment_environment:
        cf = 1
    else:
        cf = config['cf']
    weights = config['weights']
    poisson = config['poisson']
    dc_exc = config['dc_exc']
    dc_inh = config['dc_inh']
    
    #reading neural parameters
    exc_pms = read_neural_parameters(exc_pms_file_name, is_verbose)
    inh_pms = read_neural_parameters(inh_pms_file_name, is_verbose)
    
    nest_pms["exc_pms"] = exc_pms
    nest_pms["inh_pms"] = inh_pms
    
    if 't_ref' not in exc_pms['equation_params']: 
        exc_t_ref_ms=0.0 
    else: 
        exc_t_ref_ms=exc_pms['equation_params']['t_ref']
    
    #setting the sinaptic delay
    if exc_t_ref_ms <= 0.5: #ms
        network['min_syn_delay_ms'] = 0.6
    else:
        network['min_syn_delay_ms'] = exc_t_ref_ms + 0.1

    

    # Calculate number of inhibitory neurons based on excitatory populations and neurons per pop
    network['num_inh_neu'] = int(network['num_exc_pop'] * network['num_exc_neu_per_pop'] / 4)
    num_inh_neu = network['num_inh_neu']
    if 'conn_rule' not in network:
        network['conn_rule'] = 'pairwise_bernoulli'     
    if 'p_conn_exc' not in network:
        network['p_conn_exc'] = 1.0
    if 'p_conn_inh' not in network:
        network['p_conn_exc'] = 1.0
    if 'inter_pop_conn' not in network:
        network['inter_pop_conn'] = False
   
    use_recurrency = network["use_exc_recurrency"]
    if is_verbose:
        print("use_exc_recurrency:",use_recurrency)
        print("weights:",weights)

    if brain_state == "awake":
        weights=weights['awake']
    else:
        weights=weights['NREM'] 

    recurrent_weight = weights['recurrent_weight'] * cf if use_recurrency else 0.0
    inh_to_exc_weight = weights["inh_to_exc_weight"] * cf
    exc_to_inh_weight = weights["exc_to_inh_weight"]
    inh_to_inh_weight = weights["inh_to_inh_weight"]
    if network['inter_pop_conn'] == True:
        inter_pop_weight =  weights["inter_pop_weight"] * cf
    
    if is_verbose:
        print("recurrent_weight:", recurrent_weight)
        print("inh_to_exc_weight", inh_to_exc_weight)
        print("exc_to_inh_weight", exc_to_inh_weight)
        print("inh_to_inh_weight", inh_to_inh_weight)
        if network['inter_pop_conn'] == True:
            print("inter_pop_weight", inter_pop_weight)

    if 'default_plasticity' not in network:
        network['default_plasticity']={}
        network['default_plasticity']['synapse_model'] = 'stdp_synapse'
        network['default_plasticity']['tau_plus'] = 20   
        network['default_plasticity']['lambda'] = 0.1   
        network['default_plasticity']['alpha'] = 1.2
        network['default_plasticity']['mu_plus'] = 1.0   
        network['default_plasticity']['mu_minus'] = 1.0 
        network['default_plasticity']['weight'] = 0.01      
        network['default_plasticity']['Wmax'] = recurrent_weight  
        network['default_plasticity']['delay'] = 'from-exc_to_exc_delay_ms'  
    
    use_poisson_generators = network["use_poisson_generators"]
    if use_poisson_generators:
        num_poisson_generators=network["num_poisson_generators"]
        # Set stimulation parameters based on brain state
        if brain_state == "awake":
            poisson=poisson['awake']
        else:
            if brain_state == "NREM":
                poisson=poisson['NREM']
            else:
                assert(false==true)
        poisson_spreading_factor=poisson['poisson_noise_spreading_factor']
        poisson_rate = poisson['poisson_basic_rate']/ poisson_spreading_factor
        poisson_weight = poisson['poisson_weight'] * poisson_spreading_factor * cf
        poisson_delta = poisson["poisson_delta"]
    
    # Create rate arrays for Poisson generators
    #poisson_rate_arrays = [poisson_rate + i * poisson_delta for i in range(network["num_exc_neu_per_pop"])]
    #print("number of poisson generators",len(poisson_rate_arrays))
    
    #dc inhjectors to excitatory neurons
    use_dc_exc_injectors = network["use_dc_exc_injectors"]
    if use_dc_exc_injectors: 
        dc_exc_amplitudes = dc["amplitudes"]
        dc_exc_weight = dc["weight"]
        dc_exc_delay = dc["delay"]
        dc_exc_start = 0.0
        dc_exc_stop = times["sim_pms"]["stop_ms"]
    
    #dc inhjectors to inhibitory neurons
    use_dc_inh_injector = network["use_dc_inh_injector"]  # DC injection for inhibitory neurons
    if use_dc_inh_injector:
        dc_inh_amplitude = dc_inh["amplitude"] # DC amplitude for inhibitory neurons
        dc_inh_weight = dc_inh["weight"]
        dc_inh_delay = dc_inh["delay"]
        dc_inh_start = 0.0
        dc_inh_stop = times["sim_pms"]["stop_ms"]
    
    if is_verbose:
        print ("exc neu params before brain-state specific tuning:", exc_pms)
    
    if is_verbose:
        print ("inh neu params before brain-state specific tuning:", inh_pms)
    
    assert(exc_pms['neuron_kind']=="excitatory")
    if exc_pms['multi_compartment']==False:
        if is_verbose:
            print("setting brain state INDEPENDENT neural params in single-comp neu")
        # Set brain state independent parameters for single-compartment neuron
        exc_neu_params={}
        exc_neu_params['receptors']={} 
        exc_neu_params['equation_params'] = {
            "a": exc_pms['equation_params']['a'],
            "t_ref": exc_pms['equation_params']['t_ref'],
            "Delta_T": exc_pms['equation_params']['Delta_T'],
            "C_m": exc_pms['equation_params']['C_m'],
            "g_L": exc_pms['equation_params']['g_L'],
            "tau_w": exc_pms['equation_params']['tau_w'],
            "V_th": exc_pms['equation_params']['V_th'],
            "V_peak": exc_pms['equation_params']['V_peak'],
        }
        if is_verbose:
            print("setting brain state DEPENDENT neural params in single-comp neu")
        # Check if brain_state variable is set correctly
        assert(brain_state in ["awake", "NREM"])
        # Set brain state dependent parameters for single-compartment neuron
        if brain_state == "awake":
            brain_state_dependent_params = {
                "b": exc_pms['equation_params'].get('b_awake', None),
                "E_L": exc_pms['equation_params'].get('E_L_awake', None),
                "V_reset": exc_pms['equation_params'].get('V_reset_awake', None)
            }
            if is_verbose:
                print(f"Brain state 'awake' parameters: {brain_state_dependent_params}")
            exc_neu_params['equation_params'].update(brain_state_dependent_params)
        elif brain_state == "NREM":
            brain_state_dependent_params = {
                "b": exc_pms['equation_params'].get('b_NREM', None),
                "E_L": exc_pms['equation_params'].get('E_L_NREM', None),
                "V_reset": exc_pms['equation_params'].get('V_reset_NREM', None)
            }
            exc_neu_params['equation_params'].update(brain_state_dependent_params)
    else:
        assert(exc_pms['multi_compartment']==True)
        assert(use_single_compartment_environment==False)
        if is_verbose:
            print("setting brain state dependent neural params in multi-comp neu")
        exc_neu_params = {}
        exc_neu_params['receptors']=exc_pms['receptors']
        exc_neu_params['equation_params']={}
        exc_neu_params['equation_params'].update(exc_pms['equation_params'])
        # Set brain state dependent parameters for multi-compartment neuron
        if brain_state == "awake":
            print("for multi-compartment neurons the default set of parameters is the awake one")
        elif brain_state == "NREM":
            brain_state_dependent_params = {
                'b': exc_pms['NREM_changers'].get('b_NREM', None),
                'e_L_s': exc_pms['equation_params'].get('e_L_s', None) + exc_pms['NREM_changers'].get('delta_e_L_s_NREM', 0),
                'e_L_d': exc_pms['equation_params'].get('e_L_d', None) + exc_pms['NREM_changers'].get('delta_e_L_d_NREM', 0),
                'g_C_d': exc_pms['equation_params'].get('g_C_d', None) * exc_pms['NREM_changers'].get('g_C_d_NREM_multiplier', 1)
            }
            if is_verbose:
                print(f"Brain state 'sleep' parameters for multi-compartment neuron: {brain_state_dependent_params}")
            exc_neu_params['equation_params'].update(brain_state_dependent_params)
    
    if is_verbose:
        print("exc neu params AFTER brain-state specific tuning:", exc_neu_params)
    
    assert(inh_pms['neuron_kind']=="inhibitory")
    assert(inh_pms['multi_compartment']==False)
    inh_neu_params=inh_pms['equation_params']
    
    if is_verbose:
        print ("inh neu params AFTER brain-state specific tuning:", inh_neu_params)


    # Add contextual poisson signal configuration to nest_pms
    if 'contextual_poisson' in config:
        nest_pms['contextual_poisson'] = config['contextual_poisson']
    # Add contextual inhibiton signal configuration to nest_pms
    if 'contextual_inhibition' in config:
        nest_pms['contextual_inhibition'] = config['contextual_inhibition']
    nest_pms["brain_state"]= brain_state
    nest_pms["sim_pms"]=times["sim_pms"]
    nest_pms["use_nestml"]= use_nestml
    nest_pms["use_single_compartment_environment"]= use_single_compartment_environment
    nest_pms["exc_neu_params"]=exc_neu_params
    nest_pms["inh_neu_params"]=inh_neu_params
    nest_pms["exc_t_ref_ms"]=exc_t_ref_ms
    nest_pms["recording_pms"]=times["recording_pms"]
    nest_pms["events_pms"]=times["events_pms"]
    
    #how many neurons
    nest_pms["network"]={}
    nest_pms["network"]["num_exc_neu_per_pop"]=network['num_exc_neu_per_pop']
    nest_pms["network"]["num_exc_pop"]=network['num_exc_pop']
    nest_pms["network"]["num_inh_neu"]=num_inh_neu
    
    #connection rules
    nest_pms["network"]["conn_rule"] = network['conn_rule']
    nest_pms["network"]["p_conn_exc"] = network['p_conn_exc']
    nest_pms["network"]["p_conn_inh"] = network['p_conn_inh']
    nest_pms["network"]["inter_pop_conn"] = network['inter_pop_conn']
    
    #synaptic weights
    nest_pms["network"]["use_exc_recurrency"]=use_recurrency
    nest_pms["network"]["weights"]={}
    nest_pms["network"]["weights"]["exc_to_exc"] = recurrent_weight 
    nest_pms["network"]["weights"]["inh_to_exc_weight"] = inh_to_exc_weight
    nest_pms["network"]["weights"]["exc_to_inh_weight"] = exc_to_inh_weight
    nest_pms["network"]["weights"]["inh_to_inh_weight"] = inh_to_inh_weight
    if network['inter_pop_conn'] == True:
        nest_pms["network"]["weights"]["inter_pop_weight"] = inter_pop_weight
    
    #synaptic delay
    nest_pms["network"]["min_syn_delay_ms"] = network['min_syn_delay_ms']

    #plasticity default parameter (if switched on)
    nest_pms["network"]["default_plasticity"] = network['default_plasticity']
    
    #poisson generators
    nest_pms["network"]["use_poisson_generators"]=use_poisson_generators
    if(use_poisson_generators):
        nest_pms["network"]["num_poisson_generators"]=num_poisson_generators
        nest_pms["poisson"]={}
        nest_pms["poisson"]["poisson_noise_spreading_factor"]=poisson_spreading_factor
        nest_pms["poisson"]["poisson_rate"]=poisson_rate
        nest_pms["poisson"]["poisson_weight"]=poisson_weight
        nest_pms["poisson"]["poisson_delta"]=poisson_delta
    
    
    #DC exc injectors
    nest_pms["use_dc_exc_injectors"] = use_dc_exc_injectors
    if(use_dc_exc_injectors):
        nest_pms["dc_exc"]={}
        nest_pms["dc_exc"]["start_ms"] = dc_exc_start
        nest_pms["dc_exc"]["stop_ms"] = dc_exc_stop
        nest_pms["dc_exc"]["dc_exc_amplitudes"] = dc_exc_amplitudes
        nest_pms["dc_exc"]["dc_exc_weight"] = dc_exc_weight
        nest_pms["dc_exc"]["dc_exc_delay_ms"] = dc_exc_delay
    
    #DC inh injectors
    nest_pms["use_dc_inh_injector"] = use_dc_inh_injector
    if use_dc_inh_injector:
        nest_pms["dc_inh"]={}
        nest_pms["dc_inh"]["start_ms"]=dc_inh_start 
        nest_pms["dc_inh"]["stop_ms"]=dc_inh_stop 
        nest_pms["dc_inh"]["dc_inh_amplitude"] = dc_inh_amplitude
        nest_pms["dc_inh"]["weight"]=dc_inh_weight
        nest_pms["dc_inh"]["delay_ms"]=dc_inh_delay 

    return nest_pms
