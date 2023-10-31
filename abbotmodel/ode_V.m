function dV_dt = ode_V(t, V, params)
    % Parameters
    tau_m = 20;
    V_rest = -70;
    E_ex = 0;
    E_in = -70;
    
    % Extract parameters
    g_in = params(1);
    g_ex = params(2);
    
    % ODE
    dV_dt = (1/tau_m) * (V_rest - V + g_ex * (E_ex - V) + g_in * (E_in - V));
end
