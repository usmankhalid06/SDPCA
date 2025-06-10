function x_firm = firm_thresholding_nonadaptive(x, lambda1, lambda2)
    x_firm = zeros(size(x));
    for i = 1:length(x)
        if abs(x(i)) <= lambda1
            x_firm(i) = 0;
        elseif abs(x(i)) > lambda1 && abs(x(i)) <= lambda2
            x_firm(i) = sign(x(i)) * (lambda2 * (abs(x(i)) - lambda1) / (lambda2 - lambda1));
        else
            x_firm(i) = x(i);
        end
    end
end
