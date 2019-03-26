function [] =focused_beam(Data)
%This function takes in pulse echo data as input and converts it
%to a focused beamformed image with pre-defined element and transmit
%pitches, sampling frequency and focus

[rf,ntx,nel,nsamp] = readBinData(Data);

el_pitch = 0.201e-3; % inter element spacing (m)
trans_space = 0.177e-3; % inter beam spacing (m)
receive_pos = el_pitch*(nel-1)/2.*linspace(-1,1,nel); % element position
trans_pos = trans_space*(ntx-1)/2.*linspace(-1,1,ntx); % transmit position 
c = 1540; % speed of sound
fs = 40e6; % sampling frequency in Hz
focus = 2.5*1e-2; % in meters
fc = 7e6; % center frequency in Hz
compression = 0.6;

rfout = zeros(size(rf));
t = (0:size(rf,1)-1).*1/fs;
z = t.*c/2;
    

    for i = 1:ntx
        for j = 1:nel
            dz = focus;
            dx = abs(trans_pos(i)-receive_pos(j));
            totdistance = sqrt(dz.^2+dx.^2)-dz; % total dist. traveled by signal
            tdelay = totdistance/c; % time delay
            num_samples_delayed = ceil(tdelay*fs)+1;
            s_length = size(rf,1)-num_samples_delayed+1;
            rfout(1:s_length,j,i) = rf(num_samples_delayed:end,j,i);
        end
    end

    % extract envelope by squaring and lowpass filtering
    rfsum = squeeze(sum(rfout,2));
    rfsumSquared = rfsum.^2;
    [lpB,lpA] = butter(32,2*fc/(fs/2)); % cutoff freq is 2*fc
    env = filter(lpB,lpA,rfsumSquared);
%     compressed_env = log10(1 + 0.6*env)./log10(1 + 0.6);

    % filter and compress image
    Wn = [(fc-fc*.55/2)/(2*fc) (fc+fc*.55/2)/(2*fc)];
    [B,A] = butter(32, Wn);
    filtered_env = filter(B,A,env);
    compressed_env = log10(1 + 0.6*filtered_env)./log10(1 + 0.6);
    figure();
    imagesc(trans_pos,z,abs(compressed_env)); colormap gray; axis image;


end