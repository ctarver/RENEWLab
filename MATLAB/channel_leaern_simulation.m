%% Settings.
clear;clc;%close all;
addpath('modules');
p.n_users = 1;
p.n_ants = 6;

sim_mode = 1;
rf_freq = 3.56e9;

p.users.distance = 10;
p.users.theta = 100;

p.mod.name = 'OFDM';
p.mod.required_domain = 'freq';
p.mod.required_fs = 7.68e6;
p.mod.n_users = p.n_users;
p.mod.n_scs = 180;
p.mod.fft_size = 512;
p.mod.n_symbols = 6;
p.mod.make_cyclic = true;
p.mod.use_windowing = true;
p.mod.window_length = 32;
p.mod.use_random = false;


p.bs_array.name = 'Sim_Array'; % or 'ARGOS';
p.bs_array.n_antennas = p.n_ants;
p.bs_array.required_fs = p.mod.required_fs;
p.bs_array.required_domain = 'time';

p.ue_array.name = 'Sim_Array';
p.ue_array.n_antennas = p.n_users;
p.ue_array.required_fs = p.mod.required_fs;
p.ue_array.required_domain = 'freq';

% Only used for simulation arrays.
p.sim_channel.name = 'Quadriga';
p.sim_channel.required_fs = p.mod.required_fs;
p.sim_channel.required_domain = 'freq';
p.sim_channel.n_ants = p.n_ants;
p.sim_channel.n_users = p.n_users;
p.sim_channel.fft_size =  p.mod.fft_size;

% Real channel that we will learn.
p.real_channel.name = 'RealChannel';
p.real_channel.required_fs = p.mod.required_fs;
p.real_channel.required_domain = 'freq';
p.real_channel.n_scs = p.mod.fft_size;
p.real_channel.n_ants = p.n_ants;
p.real_channel.one_shot = true;

% Precoder.
p.precoder.name = 'MRT';
p.precoder.required_fs = p.mod.required_fs;
p.precoder.required_domain = 'freq';

%% Test
dataflow = ACLR_Sim(p);
dataflow.run()

dataflow.plot();
