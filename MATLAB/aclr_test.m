%% Settings.
p.n_users = 1;
p.n_ants = 8;

p.mod.name = 'OFDM';
p.mod.required_domain = 'freq';
p.mod.required_fs = 7.68e6;
p.mod.n_users = p.n_users;
p.mod.n_scs = 200;
p.mod.fft_size = 512;
p.mod.n_symbols = 14;

p.bs_array.name = 'Sim_Array'; % or 'ARGOS';
p.bs_array.n_antennas = p.n_ants;
p.bs_array.required_fs = p.mod.required_fs;
p.bs_array.required_domain = 'time';

p.ue_array.name = 'Sim_Array';
p.ue_array.n_antennas = p.n_users;
p.ue_array.required_fs = p.mod.required_fs;
p.ue_array.required_domain = 'time';

% Only used for simulation arrays.
p.sim_channel.name = 'LOS';
p.sim_channel.required_fs = p.mod.required_fs;
p.sim_channel.required_domain = 'freq';

% Real channel that we will learn. 
p.real_channel.name = 'RealChannel';
p.real_channel.required_fs = p.mod.required_fs;
p.real_channel.required_domain = 'freq';

% Precoder.
p.precoder.name = 'MRT';
p.precoder.required_fs = p.mod.required_fs;
p.precoder.required_domain = 'freq';

%% Test
dataflow = ACLR_Dataflow(p);
dataflow.run()
dataflow.report();
dataflow.plot();
