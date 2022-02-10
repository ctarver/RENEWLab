clear
load('capture1.mat')

% Time AND phase align. 
dataflow.v3_ue_rx.align_to(dataflow.v2_bs_out)
dataflow.v3_ue_rx.normalize_to_this_amp(1)

p.pa.name = 'GMP';
p.pa.required_domain = 'time';
p.pa.required_fs = 7.68e6;
p.pa.P = 9;
p.pa.M = 4;
p.pa.L = 0;
pa = Module.create('pa', p);


this_pa_in_sig = dataflow.v2_bs_out.data;
this_pa_out = dataflow.v3_ue_rx.data;
pa.learn(this_pa_in_sig, this_pa_out);

sim_pa_out = pa.use(this_pa_in_sig);

error = sim_pa_out - this_pa_out.';
NMSE = norm(error)/norm(this_pa_out.');
pa.coeffs
fprintf('NMSE  = %d\n', NMSE);

