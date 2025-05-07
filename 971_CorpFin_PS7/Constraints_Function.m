% Constraints function

function [c, ceq] = Constraints_Function(x, a, parameters)

    % Unpacking the parameters
    beta    = parameters(1);
    delta   = parameters(2);
    L       = parameters(3);
    R       = parameters(4);
    lambda  = parameters(5);
    Y_L     = parameters(6);
    Y_H     = parameters(7);
    pi_H    = parameters(8);
    na      = parameters(9);
    mu      = parameters(10);
    a_max   = parameters(11);
    a_min   = parameters(12);
    a_grid  = parameters(13:end);

    % Extract decision variables
    d_H = x(1);
    d_L = x(2);
    p_H = x(3);
    p_L = x(4);
    a_H = x(5);
    a_L = x(6);

    Payoff_Agent_H  = (d_H + p_H * R) + delta * (1 - p_H) * a_H;
    Payoff_Agent_L  = (d_L + p_L * R) + delta * (1 - p_L) * a_L;
    
    % Promise-Keeping constraint
    ceq(1)          = (pi_H * Payoff_Agent_H + (1 - pi_H) * Payoff_Agent_L) - a;

    % Incentive Compatibility constraint
    IC_HH           = d_H + (1 - p_H) * delta * a_H + p_H * R;
    IC_HL           = lambda * (Y_H - Y_L) + d_L + (1 - p_L) * delta * a_L + p_L * R;
    ceq(2)          = IC_HH - IC_HL;
    
    % Since there are no inequality constraints, then we set c to empty.
    c = [];
end