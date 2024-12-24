function x_out = coleman_hodgon_model(x_in, t_span, alpha, beta, c)
    % Defines the Coleman Hodgdon Hysteresis model
    %    x_out = c * x_in + x_h
    %    x_h_p = beta * x_in_p - alpha * abs(x_in_p) * x_h
    %  Parameters
    %  ----------
    %  x_in:    the displacement at the proximal end
    %  t_span:  timesteps for the function to be evaluated at
    %  dt:      timestep
    %  alpha:   model paramater
    %  beta:    model paramater
    %  c:       model paramater
    %
    %  Return
    %  ----------
    %  x_out:   the displacement at the distal end
    
    x_in_ = zeros(length(t_span) + 1,1);
    x_in_(2:end) = x_in;
    dx_in = (x_in_(2:end) - x_in_(1:end-1));

    x_out = zeros(length(t_span) + 1,1);

    for i = 1:length(t_span)
        if dx_in(i) >= 0
            % Intermediary variable C'
            foo = (x_out(i) - c * x_in_(i) - beta / alpha) * exp(alpha * x_in_(i));
            x_out(i + 1) = c * x_in_(i+1) + beta / alpha + foo * exp(-alpha * x_in_(i+1));
        else
            % Intermediary variables C''
            foo = (x_out(i) - c * x_in_(i) + beta / alpha) * exp(-alpha * x_in_(i));
            x_out(i + 1) = c * x_in_(i+1) - beta / alpha + foo * exp(alpha * x_in_(i+1));
        end
    end
    x_out = x_out(2:end);
end


