classdef PA_Dataflow < handle
    %ACLR_DATAFLOW
    
    properties
        p
        n_users
        n_ants
        bs
        ues
        sim_mode
        simulated_channel  % This is used when bs and ues are in sim mode.
        real_channel       % Will hold the RealChannel measured OTA.
        precoder
        pas   % Will hold the PA models measured.
        
        v0_downlink_data
        v1_pre_out
        v2_bs_out
        v3_ue_rx
        
        v10_pilot_signals
        v12_bs_out
        v13_ue_rx
        
        v21_pre_out
        v22_bs_out
        v23_ue_rx
    end
    
    methods
        function obj = PA_Dataflow(p)
            %ACLR_DATAFLOW
            obj.p = p;
            obj.n_users = p.n_users;
            obj.n_ants = p.n_ants;
            
            obj.bs = Module.create('bs_array', p);
            obj.ues = Module.create('ue_array', p);
            %obj.simulated_channel = Module.create('sim_channel', p);
            obj.real_channel = Module.create('real_channel', p);
            obj.precoder = Module.create('precoder', p);
            
            obj.pas = Module.create('pa', p, p.n_ants);
            
            obj.bs.subscribe_rx(obj.ues)
        end
        
        function run(obj)
            %% Step 1. No beamforming. See what RX power is.
            obj.v0_downlink_data = Signal.make_ofdm(obj.n_users, obj.p.mod);
            bypass_precoder = PrecoderBypass('name', 'PrecoderBypass', ...
                'required_fs', obj.v0_downlink_data.fs, ...
                'required_domain', 'bypass', ...
                'n_ant', obj.n_ants);
            obj.v1_pre_out = bypass_precoder.use(obj.v0_downlink_data);
            obj.v2_bs_out = obj.bs.tx(obj.v1_pre_out);        
            %ue_rx = obj.simulated_channel.use(obj.v2_bs_out); % Only used if array are sim.
            obj.v3_ue_rx = obj.ues.rx(obj.v2_bs_out);  % Arg is ignored if real array.
            
            %% Step 2. Learn PAs.
            % Make pilots.
            % Run Downlink from each TX to the UEs
            ofdm_settings = obj.p.mod;
            ofdm_settings.n_users = obj.n_ants;

            for i_ant = 1:obj.n_ants
                % Make data for all antennas then zero out the ones we
                % aren't using.
                this_pa_in_sig = Signal.make_ofdm(obj.n_ants, obj.p.mod);
                complement = setxor(i_ant, 1:obj.n_ants);
                this_pa_in_sig.data(complement, :, :)  = ...
                    zeros(size(this_pa_in_sig.data(complement, :, :)));
                
                bs_out = obj.bs.tx(this_pa_in_sig);
%                ue_rx = obj.simulated_channel.use(bs_out); % Only used if array are sim.
                this_pa_out = obj.ues.rx(bs_out); % Hopefully wired...
                obj.pas(i_ant).learn(this_pa_in_sig, this_pa_out);
            end
            % We should probably check the ACLR of each PA using the above.
            
            %% Step 3. Main Downlink
            obj.precoder.update(obj.real_channel.H);
            obj.v21_pre_out = obj.precoder.use(obj.v0_downlink_data);
            obj.v22_bs_out = obj.bs.tx(obj.v21_pre_out);
            %ue_rx = obj.simulated_channel.use(obj.v22_bs_out); % Only used if array are sim.
            obj.v23_ue_rx = obj.ues.rx(ue_rx);
            
        end
        
        function report(obj)
            %% What was the beamforming gain?
            % Compare the ue_rx_sigs to the final beamformed data.
            
        end
        
        function plot(obj)
            obj.v3_ue_rx.plot_psd();
            obj.v23_ue_rx.plot_psd();
        end
    end
end

