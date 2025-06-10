function TC = generate_TC(onsets, dur, N)
    stim = get_stim(onsets, dur, N);
    TC = zscore(stim);  % Normalizes the stimulus
end

function stim = get_stim(onsets, dur, nStim)
    stim = zeros(1, nStim);
    for i = 1:length(onsets)
        stim(onsets(i) + 1 : onsets(i) + dur) = 1;
    end
end
