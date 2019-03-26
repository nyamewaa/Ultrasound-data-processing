function[b, x, z] = dynamic_beam(rf, acq_params, bf_params)
%This function takes in pulse echo data as input and converts it
%to a focused beamformed image with pre-defined element and transmit
%pitches, sampling frequency and focus

    rx_pos = acq_params.rx_pos(:,1);
    tx_pos = acq_params.tx_pos(:,1);
    rf = rf(acq_params.t0:end,:,:); 
    num_samp = acq_params.samples-acq_params.t0+1;
    d = size(rf); 
    num_el = d(2);
    num_tx = d(3); 
    max_delay = 0;
    apod_rf = zeros(size(rf));
    apod_mat = ones(num_samp, num_el);
    depth = num_samp/acq_params.fs*acq_params.c/2;
    focus = (depth/bf_params.num_foci)/2:(depth/bf_params.num_foci):depth-(depth/bf_params.num_foci)/2;
    
    n_images = zeros(num_samp, num_tx, length(focus));

    rf_out = zeros(size(rf));
    t = (0:size(rf,1)-1).*1/acq_params.fs;
    z = t.*acq_params.c/2;
    left_side_transducer = min(rx_pos);

    for n = 1:length(focus)
        dz = focus(n);
        for i = 1:num_tx
            el_min_idx = ceil(((tx_pos(i)-bf_params.app_size/2)-left_side_transducer)/bf_params.pitch);
            el_max_idx = floor(((tx_pos(i)+bf_params.app_size/2)-left_side_transducer)/bf_params.pitch);
            for j = el_min_idx:el_max_idx
                dx = abs(tx_pos(i)-rx_pos(j));
                totdistance = sqrt(dz.^2+dx.^2)-dz; % total dist. traveled by signal
                tdelay = totdistance/acq_params.c; % time delay
                num_samples_delayed = round(tdelay*acq_params.fs)+1;
                if num_samples_delayed > max_delay
                    max_delay = num_samples_delayed;
                end
                s_length = num_samp-num_samples_delayed+1;
                rf_out(1:s_length,j,i) = rf(num_samples_delayed:end,j,i);
            end
        
            % apodization: user defined profile
            if strcmp(bf_params.apodization,'hamming')
                h = hamming(num_el);
            elseif strcmp(bf_params.apodization,'hann');
                h = hann(num_el);
            elseif strcmp(bf_params.apodization,'flat');
                h = flattopwin(num_el);
            else
                h = ones(num_el,1);
            end

            for k = 1:num_el
                apod_mat(:,k) = h(k);
            end
            apod_rf(:,:,i) = rf_out(:,:,i).*apod_mat;
        end

        rfsum = squeeze(sum(apod_rf,2));
        n_images(:,:,n) = rfsum;
    end
    
    % plotting apodization profiles
    if strcmp(bf_params.apodization,'hamming')
        wvtool(h);
    elseif strcmp(bf_params.apodization,'flat')
        wvtool(h);
    elseif strcmp(bf_params.apodization,'hann')
        wvtool(h);
    end
    
    combined = zeros(num_samp-max_delay, num_tx);
   
    % concatenate images of different foci
    step = ceil(num_samp/length(focus));
    for i = 1:length(focus)
        if i ~= length(focus)
            combined((i-1)*step+1:i*step,:) = n_images((i-1)*step+1:i*step,:,i);
        else
            combined((i-1)*step+1:end,:) = n_images((i-1)*step+1:end-max_delay,:,i);
        end
    end
        
    % extract envelope by squaring and lowpass filtering
    rfsumSquared = combined.^2;
    [lpB,lpA] = butter(32,2*acq_params.f0/(acq_params.fs/2)); % cutoff freq is 2*acq_params.f0
    env = filter(lpB,lpA,rfsumSquared);

    % filter and compress image
    Wn = [(acq_params.f0-acq_params.f0*.55/2)/(2*acq_params.f0) 
        (acq_params.f0+acq_params.f0*.55/2)/(2*acq_params.f0)];
    [B,A] = butter(32, Wn);
    filtered_env = abs(filter(B,A,env));
    compressed_env = exp(bf_params.compression*log(filtered_env));
    compressed_env = compressed_env./max(compressed_env(:));
    b = db(compressed_env);
    x = tx_pos;
end

