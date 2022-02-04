%% Settings.
p.n_users = 1;

p.mod.name = 'OFDM';
p.mod.required_domain = 'freq';
p.mod.required_fs = 2*7.68e6;
p.mod.n_users = p.n_users;
p.mod.n_scs = 200;
p.mod.fft_size = 2*512;
p.mod.n_symbols = 14;
p.mod.make_cyclic = false;
p.mod.use_windowing = true;
p.mod.window_length = 32;
p.mod.use_random = false;

v0_downlink_data = Signal.make_ofdm(p.n_users, p.mod);

v0_downlink_data.plot_iq;
v0_downlink_data.plot_psd;
v0_downlink_data.plot_constellation;
v0_downlink_data.measure_channels([-6 -2],[-2 2],[2 6]);

v0_downlink_data.match_this('time');
v0_downlink_data.plot_constellation;
v0_downlink_data.match_this('freq');
ber = v0_downlink_data.calculate_bit_errors;
evm = v0_downlink_data.calculate_evm;
v0_downlink_data.plot_spectrogram;
fprintf('BER  = %d. EVM = %d\n', ber, evm);