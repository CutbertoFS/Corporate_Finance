% Solve the model

function [b, d_H, d_L, p_H, p_L, a_H, a_L] = Solve_Model(parameters, results)

    % Unpacking the parameters and policy functions
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
    a_grid  = parameters(13);
    b       = results(1);
    d_H     = results(2);
    d_L     = results(3);
    p_H     = results(4);
    p_L     = results(5);
    a_H     = results(6);
    a_L     = results(7);

    % Set tolerance 
    tol = 1e-3;

    % Initialize variables
    next_b = zeros(na, 1);    % Pre-allocate next_b
    old_b  = b;               % First best guess
    error  = 10;              % Initialize error
    n      = 0;               % Iteration counter

    % Iteration
    while error > tol
        for i = 1:na
            % Call the solver and update next_b[i]
            next_b(i) = Solver(a_grid(i), old_b, parameters);
        end

        % Calculate the error
        error = max(abs(next_b - old_b));   % Obtain the error
        old_b = next_b;                     % Update Investor's Value Function
        n = n + 1;                          % Increase counter
        fprintf("Solve_Model is in %d iteration with error %.6f\n", n, error);
    end

    % After convergence, fill in the results for the policy functions
    for i = 1:na
        res = Solver(a_grid(i), old_b, param, results, other_param);
        b(i)    = res.obj;   
        d_H(i)  = res.d_H;   
        d_L(i)  = res.d_L;   
        p_H(i)  = res.p_H;   
        p_L(i)  = res.p_L;   
        a_H(i)  = res.a_H;   
        a_L(i)  = res.a_L;   
    end
end


