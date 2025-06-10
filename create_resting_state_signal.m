function tc_resting = create_resting_state_signal(N, freq_band, power, TR)
    t = (0:N-1)';
    low_freq = freq_band(1);
    high_freq = freq_band(2);
    white_noise = randn(N, 1);
    Y = fft(white_noise);
    freq = (0:N-1) / N / TR; 
    filter_mask = (freq >= low_freq & freq <= high_freq) | ...
                  (freq >= (1/TR - high_freq) & freq <= (1/TR - low_freq));
    
    Y_filtered = Y .* filter_mask';
    tc_resting = real(ifft(Y_filtered));
    tc_resting = tc_resting * sqrt(power) / std(tc_resting);
    drift = 0.05 * sin(2*pi*t/(N/2));
    tc_resting = tc_resting + drift;
end
