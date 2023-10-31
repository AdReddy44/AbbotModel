% Parameters for the first set of spike trains (N = 1000 synapses)
N = 1000;

% Parameters for the second set of spike trains (M = 200 synapses)
Q = 200;

total_time = 100 * 1000;  % Total time interval in milliseconds
time_resolution = 1;  % time resolution (1 ms)

res = 0;

% Define the rates for both sets of spike trains
rate_N = 10 / 1000;  % Rate for 'N' synapses (10Hz)
rate_Q = 10 / 1000;  % Rate for 'Q' synapses (10Hz)


spike_time_arr = [];
spike_trains = [];

% Generate spike trains for 'N' synapses with the rate_N
for synapse_number = 1:N
    num_spikes = round(rate_N * total_time);
    
    spike_times = cumsum(exprnd(1 / rate_N, 1, num_spikes));
    spike_times = spike_times(spike_times < total_time);

    % Round spike times individually
    rounded_spike_times = round(spike_times, res);
    
    % Add (spike_time, synapse_number) for each spike_time
    spike_trains = [spike_trains; rounded_spike_times(:), repmat(synapse_number, numel(rounded_spike_times), 1)];
    
    % Append spike times to spike_time_arr
    spike_time_arr = [spike_time_arr, rounded_spike_times];
end


% Generate spike trains for 'Q' synapses with the rate_Q
for synapse_number = 1:Q
    num_spikes = round(rate_Q * total_time);
    
    spike_times = cumsum(exprnd(1 / rate_Q, 1, num_spikes));
    spike_times = spike_times(spike_times < total_time);
    
    % Round spike times individually
    rounded_spike_times = round(spike_times, res);
    
    % Add (spike_time, synapse_number) for each spike_time
    spike_trains = [spike_trains; rounded_spike_times(:), repmat(synapse_number + N, numel(rounded_spike_times), 1)];
    
    % Append spike times to spike_time_arr
    spike_time_arr = [spike_time_arr, rounded_spike_times];
end



% Initialize a cell array to store synapse numbers for each spike time
spike_time_dict = containers.Map('KeyType', 'double', 'ValueType', 'any');

% Loop through each spike and add it to the dictionary
for idx = 1:length(spike_trains)
    spike_time = spike_trains(idx,1);
    synapse_number = spike_trains(idx,2);
    
    % Check if the time key already exists in the dictionary
    if ~isKey(spike_time_dict, spike_time)
        spike_time_dict(spike_time) = [];
    end
    
    % Append the current synapse number to the list associated with the time key
    spike_time_dict(spike_time) = [spike_time_dict(spike_time), synapse_number];
end




% Create an array of time values
tspan = 0:time_resolution:total_time;


% Initialize constants
tau_ex = 5;
tau_in = 5;
g_in_bar = 0.05;
g_bar_max = 0.015;
A_plus = 0.005;
A_ratio = 1.05;
A_minus = A_ratio * A_plus;

tau_minus = 20;
tau_plus = 20;


% Free to choose initial conditions
g_a_bar = repmat(0.015, 1, N); % Choose value between 0 and g_bar_max
M_0 = 0;
V_0 = -70;
g_ex_0 = 1;
g_in_0 = 1;
Pa_0 = 0;


% Initialize all N functions Pa
Pa_sols = cell(1, N);
[t, Pa] = ode45(@(t, Pa) ode_Pa(t, Pa, tau_plus), tspan, Pa_0);

for i = 1:N
    Pa_sols{i} = Pa;
end

disp(length(tspan));


% Solve ODE for M
[t, M] = ode45(@(t, M) ode_M(t, M, tau_minus), tspan, M_0);

[t, g_ex] = ode45(@(t, g_ex) ode_g_ex(t, g_ex, tau_ex), tspan, g_ex_0);

[t, g_in] = ode45(@(t, g_in) ode_g_ex(t, g_in, tau_ex), tspan, g_in_0);

params = [g_in_0, g_ex_0];

[t, V] = ode45(@(t, V) ode_V(t, V, params), tspan, V_0);


%run simulation iterating over each time index
for tval = 1:length(tspan) - 1
    current_time = tspan(tval);
    disp(current_time);
    
    tarr = tspan(tval:tval+1);
    
    %check for spike at current time
    if isKey(spike_time_dict, current_time)
        % Retrieve spike information list for the current time
        spike_info_list = spike_time_dict(current_time);
           
        %iterate over all synapses that have a spike at this time
        for k = 1:length(spike_info_list)
        
            synapse_number = spike_info_list(k);
    
            %check for excitatory synapse
            if synapse_number <= 1000
                a = synapse_number;
                [Pa_sols, g_a_bar, g_ex] = pre_excitatory_receiveAP(Pa_sols, A_plus, g_a_bar, M, g_bar_max, g_ex, tval, a);
                
            %check for inhibitory synapse
            elseif synapse_number > 1000
                g_in = pre_inhibitory_receiveAP(g_in, g_in_bar, tval);
            end
    
        end
   
    end


    ﻿%check for voltage threshold and fire post synaptic AP if needed
    if V(tval) >= -54
        [M, g_a_bar] = post_fireAP(M, A_minus, g_a_bar, g_bar_max, Pa_sols, tval);
        V(tval) = -60;
    end
        

    
    %below update all ODE solutions by one step

    for a = 1:length(Pa_sols)
        Pa = Pa_sols{a};
        [t, Pa_arr] = ode45(@(t, Pa) ode_Pa(t, Pa, tau_plus), tarr, Pa(tval));
        Pa(tval+1) = Pa_arr(length(Pa_arr));
        Pa_sols{a} = Pa;
    end
    
    
    [t, g_ex_arr] = ode45(@(t, g_ex) ode_g_ex(t, g_ex, tau_ex), tarr, g_ex(tval));
    g_ex(tval+1) = g_ex_arr(length(g_ex_arr));
    
    [t, g_in_arr] = ode45(@(t, g_in) ode_g_in(t, g_in, tau_in), tarr, g_in(tval));
    g_in(tval+1) = g_in_arr(length(g_in_arr));
    
    [t, M_arr] = ode45(@(t, M) ode_M(t, M, tau_minus), tarr, M(tval));
    M(tval+1) = M_arr(length(M_arr));
    
    params = [g_in(tval), g_ex(tval)];
    [t, V_arr] = ode45(@(t, V) ode_V(t, V, params), tarr, V(tval));
    V(tval+1) = V_arr(length(V_arr));
end



% Define the bin edges (from 0 to 1 in intervals of 0.05)
bin_edges = 0:0.05:1;

% Calculate the ratios and create an array of values within the specified bin edges
ratios = g_a_bar / g_bar_max;
histogram_values = histcounts(ratios, bin_edges);

% Calculate the percentage values for each bin
percentage_values = histogram_values / numel(g_a_bar);

% Create a bar plot for the histogram
bar(bin_edges(1:end-1), percentage_values, 'BarWidth', 0.05, 'EdgeColor', 'k');

% Show the plot
xlabel('Ratio');
ylabel('Percentage');
title('Histogram of Ratios');
grid on;




