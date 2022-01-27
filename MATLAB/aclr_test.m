%% Settings.
p.n_users = 1;
p.n_ants = 8;

p.mod.name = 'OFDM';
p.mod.required_domain = 'freq';
p.mod.required_fs = 7.92e6;
p.mod.n_users = p.n_users;
p.mod.n_scs = 200;
p.mod.fft_size = 528;
p.mod.n_symbols = 14;

p.bs_array.name = 'Simulation'; % or 'ARGOS';

p.ue_array.name = 'Simulation';

p.sim_channel.name = 'LOS';

p.real_channel.name = 'RealChannel';

%% Test
dataflow = ACLR_Dataflow(p);
dataflow.run()
dataflow.report();
dataflow.plot();
