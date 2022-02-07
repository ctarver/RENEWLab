%% Settings.
p.n_users = 1;
p.n_ants = 16;

p.users.distance = 10;
p.users.theta = 100;

p.mod.name = 'OFDM';
p.mod.required_domain = 'freq';
p.mod.required_fs = 7.68e6;
p.mod.n_users = p.n_users;
p.mod.n_scs = 200;
p.mod.fft_size = 512;
p.mod.n_symbols = 14;
p.mod.make_cyclic = true;
p.mod.use_windowing = true;
p.mod.window_length = 64;
p.mod.use_random = false;

p.bs_array.name = 'Sim_Array'; % or 'ARGOS';
p.bs_array.n_antennas = p.n_ants;
p.bs_array.required_fs = p.mod.required_fs;
p.bs_array.required_domain = 'time';

p.bs_array.name = 'IRIS';
p.bs_array.n_antennas = p.n_ants;
p.bs_array.required_fs = p.mod.required_fs;
p.bs_array.required_domain = 'time';
p.bs_array.use_hub = true;
p.bs_array.wired_ue = false;
p.bs_array.tx_freq = 3.6e9;  
p.bs_array.rx_freq = 3.6e9;
p.bs_array.tx_gain = 75;
p.bs_array.rx_gain = 60;
p.bs_array.chain_ids = 2:3;
p.bs_array.node_ids = 1:8;
p.bs_array.sched = "BGPG";

p.ue_array.name = 'Sim_Array';
p.ue_array.n_antennas = p.n_users;
p.ue_array.required_fs = p.mod.required_fs;
p.ue_array.required_domain = 'freq';

p.ue_array.name = 'IRIS';
p.ue_array.n_antennas = 1;
p.ue_array.required_fs = p.mod.required_fs;
p.ue_array.required_domain = 'time';
p.ue_array.wired_ue = false;
p.ue_array.use_hub = false;
p.ue_array.tx_freq = 3.6e9;  
p.ue_array.rx_freq = 3.6e9;
p.ue_array.tx_gain = 75;
p.ue_array.rx_gain = 60;
p.ue_array.chain_ids = 7;
p.ue_array.node_ids = 1;
p.ue_array.sched = "GGRG";

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

% Precoder.
p.precoder.name = 'MRT';
p.precoder.required_fs = p.mod.required_fs;
p.precoder.required_domain = 'freq';

%% Test
dataflow = ACLR_Dataflow(p);
dataflow.run()

dataflow.plot();
