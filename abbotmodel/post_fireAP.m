function [M, g_a_bar] = post_fireAP(M, A_minus, g_a_bar, g_bar_max, Pa_sols, tval)
    for a = 1:length(Pa_sols)
        g_a_bar(a) = g_a_bar(a) + Pa_sols{a}(tval) * g_bar_max;
        
        % Ensure g_a_bar does not exceed g_bar_max
        if g_a_bar(a) > g_bar_max
            g_a_bar(a) = g_bar_max;
        end
    end
    
    M(tval) = M(tval) - A_minus;
    
    % Output updated values
end
