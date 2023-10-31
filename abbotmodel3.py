#Code for Song-Abbot-Miller 2001 paper


import numpy as np
from scipy.integrate import odeint

import matplotlib.pyplot as plt


#ODE that g_ex follows
def odeg_ex(g_ex, t, tau_ex):
    dg_ex_dt = -g_ex/tau_ex
    return dg_ex_dt


#ODE that g_in follows
def odeg_in(g_in, t, tau_in):
    dg_in_dt = -g_in/tau_in
    return dg_in_dt


#ODE that M follows
def odeM(M, t, tau_minus):
    dM_dt = -M/tau_minus
    return dM_dt


#ODE that Pa follows
def odePa(Pa, t, tau_plus):
    dPa_dt = -Pa/tau_plus
    return dPa_dt


#ODE that V follows
def odeV(V, t, params):
    
    tau_m = 20
    V_rest = -70
    E_ex = 0
    E_in = -70
    
    g_in = params[0]
    g_ex = params[1]
    
    
    dV_dt = (1/tau_m)*(V_rest-V + g_ex*(E_ex-V) + g_in*(E_in-V))
    
    return dV_dt


#changes to make when post synaptic neuron fires Action Potential (AP)
def post_fireAP(M, A_minus, g_a_bar, g_bar_max, Pa_sols, tval):
    for a in range(len(Pa_sols)):
        g_a_bar[a] += Pa_sols[a][tval]*g_bar_max
        if g_a_bar[a] > g_bar_max: g_a_bar[a] = g_bar_max
    M[tval] -= A_minus
    return M, g_a_bar
    

#changes to make when pre synaptic inhibitory neuron fires Action Potential (AP)
def pre_inhibitory_receiveAP(g_in, g_in_bar, tval):
    g_in[tval] += g_in_bar
    return g_in
    
    
#changes to make when pre synaptic excitatory neuron fires Action Potential (AP)
def pre_excitatory_receiveAP(Pa_sols, A_plus, g_a_bar, M, g_bar_max, g_ex, tval, a):
    g_ex[tval] += g_a_bar[a]
    Pa_sols[a][tval] += A_plus
    g_a_bar[a] += M[tval]*g_bar_max
    if g_a_bar[a] < 0: g_a_bar[a] = 0
    return Pa_sols, g_a_bar, g_ex






def main():
    
    
    # Parameters for the first set of spike trains (N = 1000 synapses)
    N = 1000  # Number of synapses labeled 'N'

    # Parameters for the second set of spike trains (Q = 200 synapses)
    Q = 200  # Number of synapses labeled 'Q'

    total_time = 100*1000  # Total time interval in milliseconds
    time_resolution = 1  # Smaller time resolution (0.001s)

    res = 0

    # Define the rates for both sets of spike trains
    rate_N = 10/1000  # Rate for 'N' synapses (10Hz)
    rate_Q = 10/1000  # Rate for 'Q' synapses (10Hz)


    # Initialize lists to store spike trains
    spike_trains = []
    
    #spike_time_arr is a list of all spike_times
    spike_time_arr = []
    spike_times = []
    
    #create Poisson spike train for each excitatory synapse
    for synapse_number in range(N):
        #generate Poisson process
        num_spikes = int(rate_N * total_time)
        
        spike_times = np.cumsum(np.random.exponential(1 / rate_N, num_spikes))
        spike_times = spike_times[spike_times < total_time]
        
        spike_trains.extend([(round(spike_time, res), 'E', synapse_number) for spike_time in spike_times])
        for elem in spike_times:
            spike_time_arr.append(round(elem, res))
        
        
    # Generate spike trains for 'Q' synapses with the rate_Q
    
    #create Poisson spike train for each inhibitory synapse
    for synapse_number in range(Q):
        #generate Poisson process
        num_spikes = int(rate_Q * total_time)
        spike_times = np.cumsum(np.random.exponential(1 / rate_Q, num_spikes))
        spike_times = spike_times[spike_times < total_time]
        
        spike_trains.extend([(round(spike_time, res), 'I', synapse_number + N) for spike_time in spike_times])
        for elem in spike_times:
            spike_time_arr.append(round(elem, res))

    # Sort all spike times
    spike_time_arr.sort()
    spike_trains.sort()


    # Initialize a dictionary to store spike times and associated information
    spike_time_dict = {}
    
    # Loop through each spike and add it to the dictionary
    for spike_time, synapse_type, synapse_number in spike_trains:
    # Check if the time key already exists in the dictionary
        if spike_time not in spike_time_dict:
            spike_time_dict[spike_time] = []
    
        # Append the current spike information to the list associated with the time key
        spike_time_dict[spike_time].append((synapse_type, synapse_number))




    # Create an array of time values
    t = np.arange(0, total_time, time_resolution)
    
    
    #initialize constants
    tau_ex = 5
    tau_in = 5
    g_in_bar = 0.05
    g_bar_max = 0.015
    A_plus = 0.005
    A_ratio = 1.05
    A_minus = A_ratio * A_plus

    tau_minus = 20
    tau_plus = 20
    

    
    #free to choose initial conditions as paper said this does not modify results
    g_a_bar = [0.015 for i in range(N)]      #choose value between 0 and g_bar_max
    M_0 = 0
    V0 = -70
    g_ex_0 = 1
    g_in_0 = 1
    Pa_0 = 0
    
    
    #solve ODEs for g_ex and g_in
    g_ex = odeint(odeg_ex, g_ex_0, t, args = (tau_ex,)) 
    g_in = odeint(odeg_in, g_in_0, t, args = (tau_in,))
    
    
    #initialize all N functions Pa
    Pa_sols = []
    Pa = odeint(odePa,Pa_0, t, args = (tau_plus,))
    
    for i in range(N):
        Pa_sols.append(Pa.copy())
    #print(Pa_sols)

    #solve ODE for M
    M = odeint(odeM, M_0, t, args = (tau_minus,))

    #solve ODE for V
    params = (g_in[0], g_ex[0])
    V = [V0]
    
    
    # Integrate the ODE for V
    for i in range(1, len(t)):
        params = (g_in[i-1], g_ex[i-1])
        V_new = odeint(odeV, V[-1], [t[i-1], t[i]], args=(params,))
        V.append(V_new[1])  # Take the second value from the integration result
    # V now contains the solution to the ODE for V, and g_ex and g_in are also computed over time    


    #run simulation iterating over each time index
    for tval in range(len(t)-1):
        current_time = t[tval]
        print(current_time)
        #tarr is equal to [tval, tval+1]
        tarr = t[tval:tval+2]
        
        if current_time in spike_time_arr:
            print("SPIKE INCOMING")
            spike_info_list = spike_time_dict[current_time]
            
            for spike_info in spike_info_list:
                synapse_type, synapse_number = spike_info
                
                if synapse_type == 'E':
                    a = synapse_number
                    Pa_sols, g_a_bar, g_ex = pre_excitatory_receiveAP(Pa_sols, A_plus, g_a_bar, M, g_bar_max, g_ex, tval, a)
                
                
                elif synapse_type =='I':
                    g_in = pre_inhibitory_receiveAP(g_in, g_in_bar, tval)

                    
        print(V[tval])
        #check for voltage threshold and fire post synaptic AP if needed
        if V[tval] >= -54:
            print("HERE")
            print()
            print()
            
            M, g_a_bar = post_fireAP(M, A_minus, g_a_bar, g_bar_max, Pa_sols, tval)
            V[tval] = -60
        
        
        #integrate all ODEs forward one step
        for a in range(len(Pa_sols)):
            Pa = Pa_sols[a]
            Pa_arr = odeint(odePa, Pa[tval], tarr, args = (tau_plus,))
            Pa[tval+1] = Pa_arr[-1]
            Pa_sols[a] = Pa 
        
        
        g_ex_arr = odeint(odeg_ex, g_ex[tval], tarr, args = (tau_ex,))
        g_ex[tval+1] = g_ex_arr[-1]
        
        g_in_arr = odeint(odeg_in, g_in[tval], tarr, args = (tau_in,))
        g_in[tval+1] = g_in_arr[-1]
        
        M_arr = odeint(odeM, M[tval], tarr, args = (tau_minus,))
        M[tval+1] = M_arr[-1]
        
        params = (g_in[tval], g_ex[tval])
        V_arr = odeint(odeV, V[tval], tarr, args=(params,))
        V[tval+1] = V_arr[-1]
        
        
        
    # Define the bin edges (from 0 to 1 in intervals of 0.05)
    bin_edges = np.arange(0, 1.05, 0.05)

    # Calculate the ratios and create a list of values within the specified bin edges
    ratios = [value / g_bar_max for value in g_a_bar]
    hist, bins = np.histogram(ratios, bins=bin_edges)
    
    #print(ratios)
    
    # Calculate the percentage values for each bin
    percentage_values = (hist / len(g_a_bar))
    
    # Create a bar plot for the histogram
    plt.bar(bins[:-1], percentage_values, width=0.05, align='edge')
    
    # Show the plot
    plt.show()
        
        
if __name__ == '__main__':
    main()