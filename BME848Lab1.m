% Initialize variables: 
bf_params.pitch = 0.201e-3; % inter element spacing (m)
bf_params.compression = 0.6; % compression factor
bf_params.app_size = 0.03/2; % receive apperature size in m
bf_params.apodization = 'none';
bf_params.num_foci = 1;
load('s2000_hypo_phantom.mat');
%% Part A: A fixed receive focus beamformed image

focus_cm = 2.5;

data = 'imageData_Focused.bin';
[rf,num_tx,num_el,num_samp] = readBinData(data);
[b, x, z] = focused_beam(rf, acq_params, bf_params, focus_cm);
figure();
imagesc(x,z,b,[-40,0]); colormap gray; axis image;

data = 'imageData_PlaneWave.bin';
[rf,num_tx,num_el,num_samp] = readBinData(data);
[b, x, z] = focused_beam(rf, acq_params, bf_params, focus_cm);
figure();
imagesc(x,z,b,[-40,0]); colormap gray; axis image;

%% Part B: A dynamically focused image (on receive) with a minimum of 5 focal position updates

bf_params.num_foci = 5;

data = 'imageData_Focused.bin';
[rf,num_tx,num_el,num_samp] = readBinData(data);
[b, x, z] = dynamic_beam(rf, acq_params, bf_params);
figure();
imagesc(x,z,b,[-40,0]); colormap gray; axis image;

data = 'imageData_PlaneWave.bin';
[rf,num_tx,num_el,num_samp] = readBinData(data);
[b, x, z] = dynamic_beam(rf, acq_params, bf_params);
figure();
imagesc(x,z,b,[-40,0]); colormap gray; axis image;

%% Part C: Apodization 

bf_params.apodization = 'hamming';

data = 'imageData_Focused.bin';
[rf,num_tx,num_el,num_samp] = readBinData(data);
[b, x, z] = dynamic_beam(rf, acq_params, bf_params);
figure();
imagesc(x,z,b,[-40,0]); colormap gray; axis image;

data = 'imageData_PlaneWave.bin';
[rf,num_tx,num_el,num_samp] = readBinData(data);
[b, x, z] = dynamic_beam(rf, acq_params, bf_params);
figure();
imagesc(x,z,b,[-40,0]); colormap gray; axis image;
%%
bf_params.apodization = 'flat';

data = 'imageData_Focused.bin';
[rf,num_tx,num_el,num_samp] = readBinData(data);
[b, x, z] = focused_beam(rf, acq_params, bf_params);
figure();
imagesc(x,z,b,[-40,0]); colormap gray; axis image;

data = 'imageData_PlaneWave.bin';
[rf,num_tx,num_el,num_samp] = readBinData(data);
[b, x, z] = focused_beam(rf, acq_params, bf_params);
figure();
imagesc(x,z,b,[-40,0]); colormap gray; axis image;


%% Part D

bf_params.apodization = 'flat'; 
bf_params.num_foci = 6;
bf_params.compression = 0.6;

data = 'imageData_Focused.bin';
[rf,num_tx,num_el,num_samp] = readBinData(data);
[b, x, z] = dynamic_beam(rf, acq_params, bf_params);
figure();
h= imagesc(x,z,b,[-40,0]); colormap gray; axis image;

%% Part E

acq_params.num_foci = 1;
focus_cm = 2.5;

data = 'imageData_Focused.bin';
[rf,num_tx,num_el,num_samp] = readBinData(data);
[b, x, z] = focused_beam_singleTransmit(rf, acq_params, bf_params, focus_cm);
figure();
imagesc(x,z,b,[-40,0]); colormap gray; axis image;
%%
%CNR Calculation with best images
best_I=rgb2gray(imread('best.tif'));
roi_sig=roipoly(best_I);
ba_gr=roipoly(best_I);
roi_sig2=double(best_I(roi_sig));
ba_gr2=double(best_I(ba_gr));

mean_roi=mean(roi_sig2);
mean_baground=mean(ba_gr2);
stdev_noise=std(ba_gr2)/0.7;
CNR=abs(mean_roi-mean_baground)/stdev_noise


