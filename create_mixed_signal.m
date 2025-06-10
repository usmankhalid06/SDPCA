function tc_mixed = create_mixed_signal(N, task_events, rest_periods, TR)
    tc_mixed = zeros(N, 1);
    for i = 1:length(task_events)
        event_resp = create_event_related_signal(N, task_events{i}, 5, 1.0, TR);
        tc_mixed = tc_mixed + event_resp;
    end
    rest_component = create_resting_state_signal(N, [0.01, 0.1], 0.3, TR);
    tc_mixed = tc_mixed + rest_component;
    global_signal = 0.1 * sin(2*pi*(0:N-1)'/(N/4));
    tc_mixed = tc_mixed + global_signal;
end