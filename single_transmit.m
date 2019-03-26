function [] =single_transmit(Data)
%This function takes in pulse echo data as input and converts it
%to a focused beamformed  image with pre-defined element and transmit
%pitches, sampling frequency and focus

[rf,ntx,nel,nsamp] = readBinData(Data);
el_pitch = 0.201e-3; %inter element spacing (m)
trans_space = 0.177e-3; %inter beam spacing (m)
receive_pos = el_pitch*(nel-1)/2.*linspace(-1,1,nel); %(element position)
trans_pos = trans_space*(ntx-1)/2.*linspace(-1,1,ntx); %(transmit positionP
c = 1540; %speed of sound
fs = 40e6; %sampling frequency
focus = 2.5*1e-2; %(m)

rfout = zeros(size(rf));
t = (0:size(rf,1)-1).*1/fs;
z = t.*c/2;
    

    for i = 1:ntx
        for j = 1:nel
            dz = focus;
            dx = abs(trans_pos(i)-receive_pos(j));
            totdistance = dz+sqrt(dz.^2+dx.^2)-2*dz;
            tdelay = totdistance/c;
            nsample_delay = ceil(tdelay*fs)+1;
            s_length = size(rf,1)-nsample_delay+1;
            rfout(1:s_length,j,i) = rf(nsample_delay:end,j,i);
        end
    end

    
    rfsum = squeeze(sum(rfout,2));
    env = abs(hilbert(rfsum));
    env_comp=log(max(env,0.6));
    figure;
    env_single=env_comp(:,21);
    imagesc(trans_pos,z,env_single); colormap gray; axis image