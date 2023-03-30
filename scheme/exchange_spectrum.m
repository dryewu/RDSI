function exchange_spectrum()

% Morozov, Sergey, et al.
% "Diffusion processes modeling in magnetic resonance imaging." 
% Insights into Imaging 11.1 (2020): 1-9.

    adc_restricted = [];  % restricted with impermeable membranes
    adc_semi_restricted = [];    % restricted with semi-permeable membrances
    adc_hindered = [];    % hindered
    adc_isotropic = (0e-3 : 0.1e-3 : 1e-3)';  % isotropic
    adc_free = (2.5e-3 : 0.1e-3 : 3e-3)';  % free water

    for adc_long_diff = 0.5e-3 : 0.1e-3 : 2e-3
        for adc_trans_diff = 0e-3 : 0.2e-3 : 1e-3
            rate = adc_long_diff/adc_trans_diff;
            if isinf(rate) || rate > pi^2
                adc_restricted = [adc_restricted; adc_long_diff adc_trans_diff];
            elseif rate < (pi^2) && rate > (pi/2)^2
                adc_semi_restricted = [adc_semi_restricted; adc_long_diff adc_trans_diff];
            elseif rate > (pi/3)^2 && rate < (pi/2)^2
                adc_hindered = [adc_hindered; adc_long_diff adc_trans_diff];
            end
        end
    end

    num_restricted  = size(adc_restricted,1);
    num_semi_restricted  = size(adc_semi_restricted,1);
    num_hindered    = size(adc_hindered,1);
    num_isotropic   = size(adc_isotropic,1);
    num_free        = size(adc_free,1);
    save('exchange_spectrum.mat','adc_restricted','num_restricted',...
                                'adc_hindered','num_hindered',...
                                'adc_isotropic','num_isotropic',...
                                'adc_semi_restricted','num_semi_restricted',...
                                'adc_free','num_free');
end

