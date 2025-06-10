function tc_event = create_event_related_signal(N, event_times, duration, amplitude, TR)
    t = (0:N-1)';
    tc_event = zeros(N, 1);
    hrf_duration = 20; 
    hrf_t = (0:hrf_duration-1)';
    a1 = 6; a2 = 16; b1 = 1; b2 = 1; c = 1/6;
    hrf = ((hrf_t/a1).^(a1) .* exp(-(hrf_t-a1)/b1) - ...
           c * (hrf_t/a2).^(a2) .* exp(-(hrf_t-a2)/b2));
    hrf = hrf / max(hrf); 
    for event_start = event_times
        if event_start + hrf_duration <= N
            tc_event(event_start:event_start+hrf_duration-1) = ...
                tc_event(event_start:event_start+hrf_duration-1) + amplitude * hrf;
        end
    end
    cardiac_noise = 0.05 * sin(2*pi*1.0*t*TR);  
    resp_noise = 0.03 * sin(2*pi*0.3*t*TR);    
    tc_event = tc_event + cardiac_noise + resp_noise + 0.02*randn(N,1);
end