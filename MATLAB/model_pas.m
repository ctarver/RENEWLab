%% Settings.
clear;clc;%close all;
p.n_users = 1;
p.n_ants = 1;

sim_mode = 0;

p.users.distance = 10;
p.users.theta = 100;

p.pa.name = 'GMP';
p.pa.required_domain = 'time';
p.pa.required_fs = 7.68e6/2;
p.pa.P = 7;
p.pa.M = 4;
p.pa.L = 0;

p.mod.name = 'OFDM';
p.mod.required_domain = 'freq';
p.mod.required_fs = 7.68e6/2;
p.mod.n_users = p.n_users;
p.mod.n_scs = 100;
p.mod.fft_size = 512/2;
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
    p.bs_array.tx_freq = 3.6e9;
    p.bs_array.rx_freq = 3.6e9;
    p.bs_array.tx_gain = 80;
    p.bs_array.rx_gain = 60;
    p.bs_array.chain_ids = 6;
    p.bs_array.node_ids = 8;
    p.bs_array.sched = "BGPG";
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
    p.ue_array.wired_ue = false;
    p.ue_array.use_hub = false;
    p.ue_array.tx_freq = 3.6e9;
    p.ue_array.rx_freq = 3.6e9;
    p.ue_array.tx_gain = 75;
    p.ue_array.rx_gain = 65;
    p.ue_array.chain_ids = 5;
    p.ue_array.node_ids = 8;
    p.ue_array.sched = "GGRR";
    p.ue_array.is_bs = false;
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
p.real_channel.one_shot = false;

% Precoder.
p.precoder.name = 'MRT';
p.precoder.required_fs = p.mod.required_fs;
p.precoder.required_domain = 'freq';

%% Test
dataflow = PA_Dataflow(p);
dataflow.run()

dataflow.plot();
