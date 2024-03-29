function default_spectrum()
    adc_restricted = [];
    adc_hindered = [];
    adc_isotropic = (0.5e-3 : 0.1e-3 : 3e-3)';

    for adc_long_diff = 0.5e-3 : 0.1e-3 : 2e-3
        for adc_trans_diff = 0e-3 : 0.2e-3 : 1e-3
            rate = adc_long_diff/adc_trans_diff;
            if isinf(rate) || rate > pi^2/2
                adc_restricted = [adc_restricted; adc_long_diff adc_trans_diff];
            elseif rate > 1.1 && rate < (pi/2)^2
                adc_hindered = [adc_hindered; adc_long_diff adc_trans_diff];
            end
        end
    end

    num_restricted  = size(adc_restricted,1);
    num_hindered    = size(adc_hindered,1);
    num_isotropic   = size(adc_isotropic,1);
    save('default_spectrum.mat','adc_restricted','num_restricted',...
                                'adc_hindered','num_hindered',...
                                'adc_isotropic','num_isotropic');
end