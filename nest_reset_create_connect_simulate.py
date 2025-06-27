import nest
from cm_neuron import create_cm_neuron
from cm_neuron_nestml import create_nestml_neuron, create_receptor_mapping
from synaptic_io import *

# 2024-1002
# List of receptors used in Ca-AdEx neuron
# It should be moved in a different file (eg. Ca-AdEx_neural_params.yaml)
GABA_soma, AMPA_soma, GABA_dist, AMPA_dist, REFRACT_soma, ADAPT_soma, BETA_dist, NMDA_soma, AMPA_NMDA_soma, ALPHAexc_soma, ALPHAinh_soma, ALPHAexc_dist, ALPHAinh_dist, NMDA_dist, AMPA_NMDA_dist = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14


def nest_reset_create_connect_simulate(nest_pms, num_threads, verbose):

    sim_completed=False
    sim_pms=nest_pms["sim_pms"]
    use_nestml=nest_pms["use_nestml"]
    nest.ResetKernel()
    if use_nestml: 
        nest.Install("ca_adex_2expsyn_module")

    nest.SetKernelStatus({"resolution": sim_pms["resolution_ms"]})
    nest.SetKernelStatus({'local_num_threads':num_threads})
    
    num_exc_neu_per_pop = nest_pms["network"]["num_exc_neu_per_pop"]
    num_exc_pop = nest_pms["network"]["num_exc_pop"]
    conn_rule = nest_pms["network"]["conn_rule"]
    p_conn_exc = nest_pms["network"]["p_conn_exc"]
    p_conn_inh = nest_pms["network"]["p_conn_exc"]
    default_plasticity = nest_pms["network"]["default_plasticity"]

    use_single_compartment_environment=nest_pms["use_single_compartment_environment"]
    
    print("IN nest_reset_create_connect_simulate: use_single_compartment_environment =", use_single_compartment_environment)
    print("IN nest_reset_create_connect_simulate: use_nestml =", use_nestml)
    exc_neu_params=nest_pms['exc_neu_params'] 
    if not use_nestml:
        if use_single_compartment_environment:
            # Create excitatory neuron populations
            params=exc_neu_params['equation_params']
            params["tau_minus"]={}
            params["tau_minus"]=20.0
            #nest.Create("iaf_psc_alpha", params={"tau_minus": 20.0})
            neurons = [
                nest.Create(
                    "aeif_cond_alpha",
                    num_exc_neu_per_pop,
                    params=params
                )
                for _ in range(num_exc_pop)       
            ]
        else:
            # Create excitatory Ca-AdEx neuron populations
            nest.CopyModel('static_synapse', 'REFRACT_syn')
            nest.CopyModel('static_synapse', 'ADAPT_syn')
            nest.CopyModel('static_synapse', 'BAP_syn')
            
            neurons = [
                create_cm_neuron(num_exc_neu_per_pop, params=exc_neu_params['equation_params']) for _ in range(num_exc_pop)
            ]
    else:
        neurons = [
            create_nestml_neuron(
                num_exc_neu_per_pop, 
                params=exc_neu_params['equation_params'], 
                multi_comp=not(use_single_compartment_environment)
            ) for _ in range(num_exc_pop)
        ]
        ri = create_receptor_mapping(multi_comp=not(use_single_compartment_environment))
        
    # Create inhibitory neurons
    inh_neu_params = nest_pms["inh_neu_params"]
    num_inh_neu = nest_pms["network"]["num_inh_neu"]
    inh_neurons = nest.Create(
        "aeif_cond_alpha",
        num_inh_neu,
        params=inh_neu_params,
    )

    # Spike recorders for excitatory and inhibitory populations
    recording_pms=nest_pms["recording_pms"]
    spike_recorders = [nest.Create("spike_recorder",
        {"start": recording_pms["start_ms"], "stop": recording_pms["stop_ms"]}) 
                   for _ in range(num_exc_pop)]
    inh_spike_recorder = nest.Create("spike_recorder", 
        {"start": recording_pms["start_ms"], "stop": recording_pms["stop_ms"]})

    # Connect spike detectors
    for i in range(num_exc_pop):
        nest.Connect(neurons[i], spike_recorders[i])
    nest.Connect(inh_neurons, inh_spike_recorder)


    # 2024-1002
    # Added multimeter on Ca-AdEx only (temporary solution)
    # Target pop and neu and all the parameters to be specified in a yaml
    if not use_nestml:
        if use_single_compartment_environment:
            multimeter = nest.Create('multimeter', 1, {'record_from': ['V_m'], 'interval': .1, 
                                'start': recording_pms['start_ms'], 'stop': recording_pms['stop_ms']})
            pop = 0
            neu = 0
            # nest.Connect(multimeter, neurons[pop])
            nest.Connect(multimeter, neurons[pop][neu])
        else:
            # create multimeter to record compartment voltages and various state variables
            rec_list = ['v_comp0', 'v_comp1', 'w_5','m_Ca_1','h_Ca_1','i_AMPA_9']
            pop = 0
            neu = 0
            multimeter = nest.Create('multimeter', 1, {'record_from': rec_list, 'interval': .1, 
                                'start': recording_pms['start_ms'], 'stop': recording_pms['stop_ms']})
            # nest.Connect(multimeter, neurons[pop])
            nest.Connect(multimeter, neurons[pop][neu])
    else:
        # create multimeter to record compartment voltages and various state variables
        rec_list = ['v_comp0', f'syn_2exp{ri.ALPHAexc_soma}']
        if not use_single_compartment_environment:
            rec_list.append('v_comp1')
        pop = 0
        neu = 0
        multimeter = nest.Create('multimeter', 1, {'record_from': rec_list, 'interval': .1, 
                            'start': recording_pms['start_ms'], 'stop': recording_pms['stop_ms']})
        # nest.Connect(multimeter, neurons[pop])
        nest.Connect(multimeter, neurons[pop][neu])



    # Connection specifications
    exc_t_ref_ms=nest_pms["exc_t_ref_ms"]
    min_syn_delay_ms=nest_pms["network"]["min_syn_delay_ms"]
    use_recurrency=nest_pms["network"]["use_exc_recurrency"]
    recurrent_weight = nest_pms["network"]["weights"]["exc_to_exc"] 
    inh_to_exc_weight = nest_pms["network"]["weights"]["inh_to_exc_weight"]
    exc_to_inh_weight = nest_pms["network"]["weights"]["exc_to_inh_weight"] 
    inh_to_inh_weight = nest_pms["network"]["weights"]["inh_to_inh_weight"]

     
    exc_to_exc_delay_ms = min_syn_delay_ms + 1.
    inh_to_inh_delay_ms = min_syn_delay_ms + exc_t_ref_ms + 1.0
    assert(conn_rule == 'pairwise_bernoulli')
    conn_spec_dict_exc = {"rule": conn_rule, "p": p_conn_exc}
    conn_spec_dict_inh = {"rule": conn_rule, "p": p_conn_inh}
    conn_spec_dict_exc = {"rule": 'all_to_all', "allow_autapses":False}
    conn_spec_dict_inh = {"rule": 'all_to_all', "allow_autapses":False}

    #set to false the book-keeping of intra-pop and inter pop syn types
    present_exc_conn={}
    for i in range(num_exc_pop):
        present_exc_conn[i]={}
        for j in range(num_exc_pop):
            present_exc_conn[i][j] = {'synapse_model': False}
        
    #set initial intra-pop connections    
    for i in range(num_exc_pop):
        static_synapse_i_i='static_synapse_'+str(i)+'_'+str(i)
        nest.CopyModel('static_synapse',static_synapse_i_i)
        syn_spec={"synapse_model": static_synapse_i_i, "weight": recurrent_weight, "delay": exc_to_exc_delay_ms}
        if not use_nestml:
            if(not use_single_compartment_environment):
                syn_spec.update({'receptor_type': ALPHAexc_soma})
        else:
            syn_spec.update({'receptor_type': ri.ALPHAexc_soma})         
            
        nest.Connect(neurons[i], neurons[i], conn_spec=conn_spec_dict_exc, syn_spec=syn_spec)
        present_exc_conn[i][i] = {'synapse_model': static_synapse_i_i} 

    #inter assemblies connections: initial set-up
    if nest_pms["network"]["inter_pop_conn"] == True:

        inter_pop_weight = nest_pms["network"]["weights"]["inter_pop_weight"]
        
        for i in range(num_exc_pop):
            for j in range(num_exc_pop):
                if i != j:
                    static_synapse_i_j='static_synapse_'+str(i)+'_'+str(j)
                    nest.CopyModel('static_synapse',static_synapse_i_j)
                    syn_spec={"synapse_model":static_synapse_i_j,"weight": inter_pop_weight, "delay": exc_to_exc_delay_ms}
                    if not use_nestml:
                        if(not use_single_compartment_environment):
                            syn_spec.update({'receptor_type': ALPHAexc_soma})
                    else:
                        syn_spec_dict.update({'receptor_type': ri.ALPHAexc_soma})
                    
                    nest.Connect(neurons[i], neurons[j], conn_spec=conn_spec_dict_exc, syn_spec=syn_spec)
                    present_exc_conn[i][j]['synapse_model'] = static_synapse_i_j

    #inh to inh connections: initial setup
    nest.CopyModel('static_synapse','static_synapse_inh_inh')
    nest.Connect(inh_neurons, inh_neurons, conn_spec=conn_spec_dict_inh,\
                syn_spec={"synapse_model": 'static_synapse_inh_inh', "weight": inh_to_inh_weight, "delay": inh_to_inh_delay_ms})
    
    # Connect inhibitory neurons to all excitatory neurons and vice versa
    #before 2024-1006 (QUESTA NON ROMPE IL KERNEL)
    inh_to_exc_delay_ms = min_syn_delay_ms + exc_t_ref_ms + 0.55
    exc_to_inh_delay_ms = min_syn_delay_ms + 0.5

    #IN PROVA QUESTA DOPO FUNZIONA ANCORA
    #inh_to_exc_delay_ms = min_syn_delay_ms + 0.55
    #exc_to_inh_delay_ms = min_syn_delay_ms + 0.5 

    #IN PROVA 
    #inh_to_exc_delay_ms = min_syn_delay_ms + 0.1
    #exc_to_inh_delay_ms = min_syn_delay_ms + 0.1   
    
    #QUESTA CONFIGURAZIONE FA CRASHARE IL KERNEL
    #inh_to_exc_delay_ms = min_syn_delay_ms 
    #exc_to_inh_delay_ms = min_syn_delay_ms
    
    for i in range(num_exc_pop):
        
        # inh to exc(i)
        new_syn_model_inh_exc_n = 'static_synapse_inh_exc_'+str(i)
        nest.CopyModel('static_synapse', new_syn_model_inh_exc_n)
        syn_spec={"synapse_model": new_syn_model_inh_exc_n, "weight": inh_to_exc_weight, "delay": inh_to_exc_delay_ms}

        if not use_nestml:
            if(not use_single_compartment_environment):
                syn_spec.update({"weight": -inh_to_exc_weight, 'receptor_type': ALPHAinh_soma})
        else: 
            syn_spec.update({
                "weight": -inh_to_exc_weight, 
                "receptor_type": ri.ALPHAinh_soma,
            })

        nest.Connect(inh_neurons, neurons[i], conn_spec=conn_spec_dict_inh, syn_spec=syn_spec)
        
        # exc(i) to inh    
        new_syn_model_exc_n_inh = 'static_synapse_exc_' + str(i) + '_inh'
        nest.CopyModel('static_synapse', new_syn_model_exc_n_inh)
        syn_spec={"synapse_model": new_syn_model_exc_n_inh, "weight": exc_to_inh_weight, "delay": exc_to_inh_delay_ms}
        nest.Connect(neurons[i], inh_neurons, conn_spec=conn_spec_dict_exc, syn_spec=syn_spec)      

    # Create and connect Poisson generators if enabled
    use_poisson_generators=nest_pms["network"]["use_poisson_generators"]
    if use_poisson_generators:
        num_poisson_generators=nest_pms["network"]["num_poisson_generators"]
        poisson_rate=nest_pms["poisson"]["poisson_rate"]
        pgs = nest.Create("poisson_generator", num_poisson_generators, params={"rate": poisson_rate})
        poisson_weight=nest_pms["poisson"]["poisson_weight"]
        
        syn_spec={"weight": poisson_weight, "delay": 1.0}
        for i in range(num_exc_pop):
            if not use_nestml:         
                if not use_single_compartment_environment:
                    syn_spec.update({'receptor_type': ALPHAexc_soma})
            else:
                syn_spec.update({'receptor_type': ri.ALPHAexc_soma}) 
            nest.Connect(pgs, neurons[i], syn_spec=syn_spec)
            
    # Add contextual poisson signal if configured
    if 'contextual_poisson' in nest_pms and \
            nest_pms['brain_state'] in nest_pms['contextual_poisson']: 
        # Extract the upper-level configuration parameters
        contextual_poisson_config = nest_pms['contextual_poisson'][nest_pms['brain_state']]
        
        # Ensure all required keys are present in the configuration
        if all(key in contextual_poisson_config for key in ['events', 'spreading_factor', 'basic_rate', 'poisson_weight']):
            events = contextual_poisson_config['events']
            spreading_factor = contextual_poisson_config['spreading_factor']
            basic_rate = contextual_poisson_config['basic_rate']
            poisson_weight = contextual_poisson_config['poisson_weight']
    
            # Iterate through each event in the events list
            for event in events:
                target_pop = event['target_population']
                start_time = event['start_time_ms']
                stop_time = event['stop_time_ms']
                
                # Create a Poisson generator for the event
                contextual_poisson_gen = nest.Create('poisson_generator', 1, {'rate': basic_rate})
                
                # Set the start and stop times for the generator
                nest.SetStatus(contextual_poisson_gen, {'start': start_time, 'stop': stop_time})
    
                # Define the synapse specifications
                syn_spec = {'weight': poisson_weight * spreading_factor, "delay": 1.0}
    
                # Connect the generator to the target population
                if not use_nestml:
                    if not use_single_compartment_environment:             
                        syn_spec.update({'receptor_type': AMPA_NMDA_dist})
                else:
                    # if use_single_compartment_environment:    
                    #     syn_spec.update({'receptor_type': ri.ALPHAexc_soma})
                    # else:    
                    #     syn_spec.update({'receptor_type': ri.AMPA_NMDA_dist})
                    syn_spec.update({'receptor_type': ri.PSC_context})

                nest.Connect(contextual_poisson_gen, neurons[target_pop], syn_spec=syn_spec)


    # DC current injection for all neurons if enabled
    use_dc_exc_injectors=nest_pms["use_dc_exc_injectors"]
    if use_dc_exc_injectors:
        dc_exc_start_ms=nest_pms["dc_exc"]["start_ms"]
        dc_exc_stop_ms=nest_pms["dc_exc"]["stop_ms"]
        dc_exc_amplitudes=nest_pms["dc_exc"]["dc_exc_amplitudes"]
        dc_exc_weight=nest_pms["dc_exc"]["dc_exc_weight"]
        dc_exc_delay_ms=nest_pms["dc_exc"]["dc_exc_delay_ms"]
        dcgs = [nest.Create("dc_generator", 1,\
                {"start": dc_exc_start_ms, "stop": dc_exc_stop_ms, "amplitude": amp})\
                for amp in dc_exc_amplitudes]
        for i in range(num_exc_pop):
            if use_single_compartment_environment:
                nest.Connect(dcgs[i], neurons[i],\
                    syn_spec={'weight': dc_exc_weight, 'delay': dc_exc_delay_ms})
            else:
                #assert("MC neuron not yet implemented")
                nest.Connect(dcgs[i], neurons[i], syn_spec={'weight': dc_exc_weight, 'delay': dc_exc_delay, 'receptor_type': 0})

    # Specific DC current injection for inhibitory neurons
    use_dc_inh_injector=nest_pms["use_dc_inh_injector"]
    if use_dc_inh_injector:
        dc_inh_start_ms = nest_pms["dc_inh"]["start_ms"]
        dc_inh_stop_ms = nest_pms["dc_inh"]["stop_ms"]
        dc_inh_amplitude = nest_pms["dc_inh"]["dc_inh_amplitude"]
        dc_inh_weight = nest_pms["dc_inh"]["weight"]
        dc_inh_delay_ms = nest_pms["dc_inh"]["delay_ms"]
        
        dc_inh_generator = nest.Create("dc_generator", 1,\
                {"start": dc_inh_start_ms, "stop": dc_inh_stop_ms, 
                 "amplitude": dc_inh_amplitude})
        nest.Connect(dc_inh_generator, inh_neurons,\
                     syn_spec={'weight': dc_inh_weight, 'delay': dc_inh_delay_ms})

    def make_plastic_conn(source_exc_pop, target_exc_pop, Wmax, verbose):
        if verbose:
            print("in make_plastic_conn: source and target pops", source_exc_pop, target_exc_pop)
        assert(source_exc_pop >= 0 and source_exc_pop < num_exc_pop)
        assert(target_exc_pop >= 0 and target_exc_pop < num_exc_pop)
        new_syn_model = default_plasticity['synapse_model']+'_'+str(source_exc_pop)+'_'+str(target_exc_pop)
        if verbose:
            print("in make_plastic_conn: new_syn_model = ", new_syn_model)
            
        if new_syn_model not in nest.Models(mtype='synapses'):
            if verbose:
                print("in make_plastic_conn: new_syn_model not already among Models")
            if present_exc_conn[source_exc_pop][target_exc_pop]['synapse_model'] != False:
                if verbose:
                    print("in make_plastic_conn: present_exc_conn =",present_exc_conn[source_exc_pop][target_exc_pop]['synapse_model'])
                    print("in make_plastic_conn: before GetConnections (old set)")  
                existing_conns = nest.GetConnections(neurons[source_exc_pop],
                            neurons[target_exc_pop], 
                            synapse_model=present_exc_conn[source_exc_pop][target_exc_pop]['synapse_model'])
                if verbose:
                    print("in make_plastic_conn: before Disconnect (old set)")
                nest.Disconnect(existing_conns)
                if verbose:
                    print("in make_plastic_conn: after Disconnect (old set)")
                present_exc_conn[source_exc_pop][target_exc_pop] = {'synapse_model': False}

            conn_spec_dict = conn_spec_dict_exc
            nest.CopyModel(default_plasticity['synapse_model'],new_syn_model)
            syn_spec_dict = {"synapse_model": new_syn_model, 
                             "Wmax": Wmax,
                             "tau_plus": default_plasticity['tau_plus'], 
                             "lambda": default_plasticity['lambda'], 
                             "alpha": default_plasticity['alpha'], 
                             "mu_plus": default_plasticity['mu_plus'], 
                             "mu_minus": default_plasticity['mu_minus'],
                             "weight": default_plasticity['weight'], 
                             "delay": exc_to_exc_delay_ms}
            if not use_nestml:            
                if not use_single_compartment_environment:
                        syn_spec_dict.update({'receptor_type': ALPHAexc_soma})
            else:
                syn_spec_dict.update({'receptor_type': ri.ALPHAexc_soma})
                
            if verbose:
                print("in make_plastic_conn: before connecting ",new_syn_model) 
            nest.Connect(neurons[source_exc_pop], neurons[target_exc_pop],
                         conn_spec=conn_spec_dict, syn_spec=syn_spec_dict)
            if verbose:
                print("in make_plastic_conn: after connecting ",new_syn_model) 
            present_exc_conn[source_exc_pop][target_exc_pop] = {'synapse_model': new_syn_model}
            if verbose:
                print("Wmax =", Wmax, "initial_weight", default_plasticity['weight'])
        return 

    def disconnect_conn(source_exc_pop, target_exc_pop, verbose):
        if verbose:
            print("in disconnect_conn: source_exc_pop -> target_exc_pop",source_exc_pop, target_exc_pop)
        assert(source_exc_pop >= 0 and source_exc_pop < num_exc_pop)
        assert(target_exc_pop >= 0 and target_exc_pop < num_exc_pop)
        if present_exc_conn[source_exc_pop][target_exc_pop]['synapse_model'] != False:
            synapse_model = present_exc_conn[source_exc_pop][target_exc_pop]['synapse_model']
            existing_conns = nest.GetConnections(neurons[source_exc_pop],
                            neurons[target_exc_pop], synapse_model=synapse_model)
            nest.Disconnect(existing_conns)
            present_exc_conn[source_exc_pop][target_exc_pop] = {'synapse_model': False}
        # Next lines reconnect with zero values 
        conn_spec_dict = conn_spec_dict_exc
        new_disconnected_static_syn_model = \
            'static_synapse_weight_zero_'+str(source_exc_pop)+'_'+str(target_exc_pop)
        nest.CopyModel('static_synapse',new_disconnected_static_syn_model)
        syn_spec_dict = {"synapse_model":new_disconnected_static_syn_model,"weight": 0.0, "delay": exc_to_exc_delay_ms}
        if not use_nestml:
            if not use_single_compartment_environment:
                syn_spec_dict.update({'receptor_type': ALPHAexc_soma})
        else:
            syn_spec_dict.update({'receptor_type': ri.ALPHAexc_soma})
        
        nest.Connect(neurons[source_exc_pop], neurons[target_exc_pop], 
                             conn_spec=conn_spec_dict, syn_spec=syn_spec_dict)
        present_exc_conn[source_exc_pop][target_exc_pop] = {'synapse_model': new_disconnected_static_syn_model, 'weight': 0.0}
        return

    def perform_event_action(requested_action, tot_simulated_time_ms, verbose):
        if verbose:
            print("in perform_event_action: requested action:", requested_action, "at time (ms)", tot_simulated_time_ms) 
        if requested_action['kind']=='store_syn':
            #STORE SYN MATRIX
            store_syn(requested_action['syn_file_name'], tot_simulated_time_ms, verbose)
        elif requested_action['kind']=='disconnect_exc_pop_syn':
            #DISCONNECT (several options to select target intra and inter pops)

            #can be an integer, a list of integers, or a string    
            target_exc_pop = requested_action['target_exc_pop']     
            if isinstance(target_exc_pop, int):
                # Perform actions for integer value
                disconnect_conn(target_exc_pop, target_exc_pop, verbose)
            elif isinstance(target_exc_pop, list):
                # The value is a list
                # Check that all elements in the list are integers
                if all(isinstance(item, int) for item in target_exc_pop):
                    # Perform actions for list of integers
                    for item in target_exc_pop:
                        disconnect_conn(item, item, verbose)
                else:
                    raise TypeError("All items in 'target_exc_pop' list must be integers.")
            elif isinstance(target_exc_pop, str):
                # The value is a string
                # Perform actions for string value
                if target_exc_pop == 'all_intra_exc':
                    for index_target_pop in range(0,num_exc_pop):
                        disconnect_conn(index_target_pop, index_target_pop, verbose)
                elif target_exc_pop == 'all_inter_exc':
                    for index_source_pop in range(0,num_exc_pop):
                        for index_target_pop in range(0,num_exc_pop):
                            if index_source_pop!=index_target_pop:
                                disconnect_conn(index_source_pop, index_target_pop, verbose)
                elif target_exc_pop == 'all_intra_inter_exc':
                    for index_source_pop in range(0,num_exc_pop):
                        for index_target_pop in range(0,num_exc_pop):
                               disconnect_conn(index_source_pop, index_target_pop, verbose)
                else:
                    print("unrecognized command in set plastic synapses")
                    assert(False) 
            else:
                # Handle unexpected types
                raise TypeError(f"Invalid type for 'target_exc_pop': {type(target_exc_pop)}")

        elif requested_action['kind']=='exc_plasticity_ON':
            #make PLASTIC exc synapses (several options to select target intra and inter pops)
                
            Wmax = recurrent_weight * requested_action['W_max_factor']                                       

            #can be an integer, a list of integers, or a string    
            target_exc_pop = requested_action['target_exc_pop'] 
            
            if isinstance(target_exc_pop, int):
                # Perform actions for integer value
                make_plastic_conn(target_exc_pop, target_exc_pop, Wmax, verbose)
            elif isinstance(target_exc_pop, list):
                # The value is a list
                # Check that all elements in the list are integers
                if all(isinstance(item, int) for item in target_exc_pop):
                    # Perform actions for list of integers
                    for item in target_exc_pop:
                        make_plastic_conn(item, item, Wmax, verbose)
                else:
                    raise TypeError("All items in 'target_exc_pop' list must be integers.")
            elif isinstance(target_exc_pop, str):
                # The value is a string
                # Perform actions for string value
                if target_exc_pop == 'all_intra_exc':
                    for index_target_pop in range(0,num_exc_pop):
                        make_plastic_conn(index_target_pop, index_target_pop, Wmax, verbose)
                elif target_exc_pop == 'all_inter_exc':
                    for index_source_pop in range(0,num_exc_pop):
                        for index_target_pop in range(0,num_exc_pop):
                            if index_source_pop!=index_target_pop:
                                make_plastic_conn(index_source_pop, index_target_pop, Wmax, verbose)
                elif target_exc_pop == 'all_intra_inter_exc':
                    for index_source_pop in range(0,num_exc_pop):
                        for index_target_pop in range(0,num_exc_pop):
                                make_plastic_conn(index_source_pop, index_target_pop, Wmax, verbose)
                else:
                    print("unrecognized command in set plastic synapses")
                    assert(False) 
            else:
                # Handle unexpected types
                raise TypeError(f"Invalid type for 'target_exc_pop': {type(target_exc_pop)}")
        elif requested_action['kind']=='set_learning':
            # Set learning rate
            learning_rate = default_plasticity['lambda']
            alpha = default_plasticity['alpha']
            if 'lambda' in requested_action:
                learning_rate = requested_action['lambda']
            if 'alpha' in requested_action:
                alpha = requested_action['alpha']
            target_exc_pop = requested_action['target_exc_pop']
            assert(target_exc_pop >= 0 and target_exc_pop < num_exc_pop)
            # add check if target_exc_pop is in the list of plastic populations (issue #34)
            existing_conns_plastic = nest.GetConnections(neurons[target_exc_pop],
                                               neurons[target_exc_pop], synapse_model='stdp_synapse')
            nest.SetStatus(existing_conns_plastic, {'lambda': learning_rate, 'alpha': alpha})
            print('Set lambda =', learning_rate, 'alpha =', alpha, 'in population', target_exc_pop)     
        else:
            #UNRECOGNIZED
            print("in perform_event_action: requested action:",requested_action['kind'],"NOT SUPPORTED")
            assert False
        return
        
    list_of_syn_matrix_file_names=[]
    
    def store_syn(file_name, tot_simulated_time_ms, verbose):
        if(verbose):
            print("in store_syn:")
            print("will save to file:", file_name)
            print("at time (ms):", tot_simulated_time_ms)
        array_of_dicts = []
        for i in range(num_exc_pop):

            synapse_model = present_exc_conn[i][i]['synapse_model']
            array_of_dicts.append(\
                {"time (ms)": tot_simulated_time_ms,\
                "conn_index": i,\
                "synapse_model":present_exc_conn[i][i]['synapse_model'],\
                "connections": nest.GetConnections(source=neurons[i],target=neurons[i],
                           synapse_model=synapse_model).get(("source", "target", "weight"), output="pandas")})
        with open(file_name, "wb") as f:
            pickle.dump(array_of_dicts, f, pickle.HIGHEST_PROTOCOL)
        list_of_syn_matrix_file_names.append(file_name)
        return

    def simulate(verbose):
        # Run the simulation
        if verbose:
            print("in nest_..._simulate: LIST OF PLANNED ACTIONS")
            for item in nest_pms['events_pms']:
                if 'action' in item:
                        print("ACTION",item['action'])
                if 'sim_period_ms' in item:    
                    if verbose:       
                        print("expected NEXT SIM PERIOD ms",item['sim_period_ms'])
            print("in nest_..._simulate: JUST BEFORE ACTUAL EXECUTION")
    
        tot_simulated_time_ms=0.0
        if 'events_pms' in nest_pms:
            for item in nest_pms['events_pms']:
                if 'action' in item:
                    if verbose:
                        print("ACTION",item['action'])
                    #performing the action
                    perform_event_action(item['action'], tot_simulated_time_ms, verbose)   
                if 'sim_period_ms' in item:    
                    if verbose:       
                        print("expected NEXT SIM PERIOD ms",item['sim_period_ms']) 
                    next_sim_period_ms = item['sim_period_ms']
                else:
                    next_sim_period_ms = 0.0
                if next_sim_period_ms + tot_simulated_time_ms > sim_pms["stop_ms"]:
                    next_sim_period_ms = sim_pms["stop_ms"] - tot_simulated_time_ms
                    store_syn('syn_matrix-end_of_sim',sim_pms["stop_ms"],verbose) 
                    if verbose:
                        print("reduced last sim time ms to", next_sim_period_ms) 
                        print("launching next period of ", next_sim_period_ms, "ms of simulation")
                #launching simulation period
                if next_sim_period_ms > 0:
                    nest.Simulate(next_sim_period_ms)
                    tot_simulated_time_ms += next_sim_period_ms
                    if verbose:
                        print("performed until now:",tot_simulated_time_ms, "TOTAL ms of simulation")

        if tot_simulated_time_ms < sim_pms["stop_ms"]:
            next_sim_period_ms = sim_pms["stop_ms"] - tot_simulated_time_ms
            if verbose:
                print("ADDING ", next_sim_period_ms, "ms sim to complete expected" ,sim_pms["stop_ms"], "TOTAL ms of simulation")
            #possible last simulation period
            assert(next_sim_period_ms>0)
            nest.Simulate(next_sim_period_ms)
            tot_simulated_time_ms += next_sim_period_ms
            store_syn('syn_matrix-end_of_sim', tot_simulated_time_ms, verbose) 
            if verbose:
                print("performed until now:",tot_simulated_time_ms, "TOTAL ms of simulation")
            assert(tot_simulated_time_ms==sim_pms["stop_ms"])

        return 
        
    #HERE IS THE SIMULATION
    verbose=True
    simulate(verbose)
    sim_completed=True

    return sim_completed, spike_recorders, inh_spike_recorder, multimeter, list_of_syn_matrix_file_names
