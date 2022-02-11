%% Settings.
clear;clc;%close all;
addpath('modules');
p.n_users = 1;
p.n_ants = 6;

sim_mode = 0;
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

if sim_mode
    p.bs_array.name = 'Sim_Array'; % or 'ARGOS';
    p.bs_array.n_antennas = p.n_ants;
    p.bs_array.required_fs = p.mod.required_fs;
    p.bs_array.required_domain = 'time';
else
    p.bs_array.name = 'IRIS';
    p.bs_array.n_antennas = p.n_ants;
    p.bs_array.required_fs = p.mod.required_fs;
    p.bs_array.required_domain = 'time';
    p.bs_array.use_hub = true;
    p.bs_array.wired_ue = false;
    p.bs_array.tx_freq = rf_freq;
    p.bs_array.rx_freq = rf_freq;
    p.bs_array.tx_gain = 60;
    p.bs_array.rx_gain = 60;
    p.bs_array.chain_ids = 5;
    p.bs_array.node_ids = 1:6;
    p.bs_array.sched = "PPPP";
    p.bs_array.use_tdd = true;
    p.bs_array.n_frame = 10000;
    p.bs_array.max_amp = 1;
end
if sim_mode
    p.ue_array.name = 'Sim_Array';
    p.ue_array.n_antennas = p.n_users;
    p.ue_array.required_fs = p.mod.required_fs;
    p.ue_array.required_domain = 'freq';
else
    p.ue_array.name = 'IRIS';
    p.ue_array.n_antennas = 1;
    p.ue_array.required_fs = p.mod.required_fs;
    p.ue_array.required_domain = 'time';
    p.ue_array.wired_ue = true;
    p.ue_array.use_hub = false;
    p.ue_array.tx_freq = rf_freq;
    p.ue_array.rx_freq = rf_freq;
    p.ue_array.tx_gain = 75;
    p.ue_array.rx_gain = 45;
    p.ue_array.chain_ids = 7;
    p.ue_array.node_ids = 4;
    p.ue_array.sched = [];%"GGRG";
    p.ue_array.is_bs = false;
    p.ue_array.use_tdd = false;
    p.ue_array.n_frame = 1;
    p.ue_array.n_samp = 2*8192;
end
% Only used for simulation arrays.
p.sim_channel.name = 'Quadriga';
p.sim_channel.required_fs = p.mod.required_fs;
p.sim_channel.required_domain = 'freq';
p.sim_channel.n_ant = p.n_ants;
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
dataflow = ACLR_Dataflow(p);
dataflow.run()

dataflow.plot();
