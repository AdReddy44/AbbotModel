function [Pa_sols, g_a_bar, g_ex] = pre_excitatory_receiveAP(Pa_sols, A_plus, g_a_bar, M, g_bar_max, g_ex, tval, a)
    g_ex(tval) = g_ex(tval) + g_a_bar(a);
    Pa_sols{a}(tval) = Pa_sols{a}(tval) + A_plus;
    g_a_bar(a) = g_a_bar(a) + M(tval) * g_bar_max;
    
    % Ensure g_a_bar is non-negative
    if g_a_bar(a) < 0
        g_a_bar(a) = 0;
    end
    
    % Output updated values
end
