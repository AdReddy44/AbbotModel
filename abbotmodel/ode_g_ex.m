% ODE for g_ex
function dg_ex_dt = ode_g_ex(t, g_ex, tau_ex)
    dg_ex_dt = -g_ex / tau_ex;
end
