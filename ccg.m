% chatgpt translation of ccg.py from https://github.com/jiaxx/modular_network/tree/master/code
function ccgjitter = get_ccgjitter(spikes, FR, jitterwindow)
    % spikes: neuron*ori*trial*time
    % currently, time need to be devidible by jitterwindow
    if size(spikes, 1) ~= length(FR)
        error('Dimensions of spikes and FR do not match.');
    end

    n_unit = size(spikes, 1);
    n_t = size(spikes, 4);
    % triangle function
    t = -(n_t-1) : (n_t-1);
    theta = n_t - abs(t);
    NFFT = 2 .^ nextpow2(2*n_t);
    target = arrayfun(@(x) round(x), NFFT/2 + (-n_t+2:n_t));

    ccgjitter = [];
    pair = 0;

    for i = 1 : n_unit-1 % V1 cell
        for m = i+1 : n_unit % V2 cell
            % temp format: [rep, (ori), time], *time has to be the last dim
            if FR(i) > 2 && FR(m) > 2
                % temp format: [rep, (ori), time], *time has to be the last dim
                % input shape is neuron*ori*rep*time
                temp1 = squeeze(spikes(i, :, :, :));
                temp2 = squeeze(spikes(m, :, :, :));
                FR1 = squeeze(mean(sum(temp1, 3), 2));
                FR2 = squeeze(mean(sum(temp2, 3), 2));
                tempccg = xcorrfft(temp1, temp2, NFFT);

                % input shape is time*neuron*ori*rep
                temp1 = permute(temp1, [3 2 1]);
                temp2 = permute(temp2, [3 2 1]);
                ttemp1 = jitter(temp1, jitterwindow);
                ttemp2 = jitter(temp2, jitterwindow);
                tempjitter = xcorrfft(permute(ttemp1, [3 2 1]), permute(ttemp2, [3 2 1]), NFFT);
                tempjitter = squeeze(nanmean(tempjitter(:, :, target), 1));

                ccgjitter(end + 1) = (tempccg - tempjitter)' / (sqrt(FR(i) * FR(m)) * theta');
            end
        end
    end
end

function CCG = xcorrfft(a, b, NFFT)
    CCG = fftshift(ifft(fft(a, NFFT) .* conj(fft(b, NFFT))));
end

function p = nextpow2(n)
    p = ceil(log2(n));
end

function output = jitter(data, l)
    % Jittering multidemntational logical data where 
    % 0 means no spikes in that time bin and 1 indicates a spike in that time bin.

    if numel(size(data)) > 3
        flag = 1;
        sd = size(data);
        data = reshape(data, [size(data,1), size(data,2), numel(data) / (size(data,1) * size(data,2))]);
    else
        flag = 0;
    end

    psth = mean(data, 2);
    len = size(data, 1);

    if mod(size(data,1), l) ~= 0
        data(len+1:len+mod(-size(data,1),l), :, :) = 0;
        psth(len+1:len+mod(-size(data,1),l), :) = 0;
    end

    if size(psth,2) > 1
        dataj = squeeze(sum(reshape(data, [l, size(data,1) / l, size(data,2), size(data,3)]), 1));
        psthj = squeeze(sum(reshape(psth, [l, size(psth,1) / l, size(psth,2)]), 1));
    else
        dataj = squeeze(sum(reshape(data, [l, size(data,1) / l, size(data,2)])));
        psthj = sum(reshape(psth, [l, size(psth,1) / l]));
    end

    if size(data,1) == l
        dataj = reshape(dataj, [1, size(dataj,2), size(dataj,3)]);
        psthj = reshape(psthj, [1, size(psthj,2)]);
    end

    psthj = reshape(psthj, [size(psthj,1), 1, size(psthj,2)]);
    psthj(psthj == 0) = 10e-10;

    corr = dataj ./ repmat(psthj, [1, size(dataj,2), 1]);
    corr = reshape(corr, [1, size(corr,1), size(corr,2), size(corr,3)]);
    corr = repmat(corr, [l, 1, 1, 1]);
    corr = reshape(corr, [size(corr,1) * size(corr,2), size(corr,3), size(corr,4)]);

    psth = reshape(psth, [size(psth,1), 1, size(psth,2)]);
    output = repmat(psth, [1, size(corr,2), 1]) .* corr;

    output = output(1:len, :, :);
end
