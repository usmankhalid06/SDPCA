function tc_block = create_block_design_signal(N, block_length, rest_length, amplitude, TR)
    t = (0:N-1)';
    tc_block = zeros(N, 1);
    cycle_length = block_length + rest_length;
    for start_idx = 1:cycle_length:N
        block_end = min(start_idx + block_length - 1, N);
        tc_block(start_idx:block_end) = amplitude;
    end
    tc_block = smoothdata(tc_block, 'gaussian', 3);
    drift = 0.1 * sin(2*pi*t/(N/3));
    physio_noise = 0.03 * (sin(2*pi*1.0*t*TR) + 0.5*sin(2*pi*0.3*t*TR));
    tc_block = tc_block + drift + physio_noise + 0.02*randn(N,1);
end