classdef ACLR_Dataflow < handle
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
        
        v0_downlink_data
        v1_pre_out
        v2_bs_out
        v3_ue_rx
        v3_ue_rx_raw
        
        v10_pilot_signals
        v11_pre_out
        v12_bs_out
        v13_ue_rx
        
        v21_pre_out
        v22_bs_out
        v23_ue_rx
    end
    
    methods
        function obj = ACLR_Dataflow(p)
            if nargin == 0
               aclr_test 
               return;
            end
            
            %ACLR_DATAFLOW
            obj.p = p;
            obj.n_users = p.n_users;
            obj.n_ants = p.n_ants;
            
            obj.bs = Module.create('bs_array', p);
            obj.ues = Module.create('ue_array', p);
            %obj.simulated_channel = Module.create('sim_channel', p);
            obj.real_channel = Module.create('real_channel', p);
            obj.precoder = Module.create('precoder', p);
            
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
            obj.v1_pre_out.match_this('time');
            obj.v1_pre_out.normalize_to_this_rms(3);
            obj.v2_bs_out = obj.bs.tx(obj.v1_pre_out);
            %ue_rx = obj.simulated_channel.use(obj.v2_bs_out); % Only used if array are sim.
            obj.v3_ue_rx_raw = obj.ues.rx(obj.v2_bs_out);  % Arg is ignored if real array.
            obj.v3_ue_rx = obj.v3_ue_rx_raw.copy()
            % Is my signal here?
            obj.v3_ue_rx.align_to(obj.v0_downlink_data)
            
            
            %% Step 2. Learn Channel.
            % Make pilots.
            % Run Downlink from each TX to the UEs
            for i_ant = 1:obj.n_ants
                obj.v10_pilot_signals = obj.real_channel.create_pilots(obj.p.mod, i_ant);
                
                if obj.real_channel.one_shot
                    obj.v11_pre_out = obj.v10_pilot_signals;
                else
                    obj.v11_pre_out = bypass_precoder.use(obj.v10_pilot_signals);
                
                    % Delete data for other antennas..
                    complement = setxor(i_ant, 1:obj.n_ants);
                    obj.v11_pre_out.data(complement, :, :)  = zeros(size(obj.v11_pre_out.data(complement, :, :)));
                end
                obj.v11_pre_out.match_this('time');
                obj.v11_pre_out.normalize_to_this_amp(0.9);
                obj.v12_bs_out = obj.bs.tx(obj.v11_pre_out);
                %ue_rx = obj.simulated_channel.use(obj.v12_bs_out); % Only used if array are sim.
                ue_rx_raw = obj.ues.rx(obj.v12_bs_out);
                obj.v13_ue_rx = ue_rx_raw.copy();
                obj.v13_ue_rx.align_to(obj.v10_pilot_signals)
                obj.real_channel.learn(obj.v10_pilot_signals, obj.v13_ue_rx, i_ant);
                
                if obj.real_channel.one_shot
                    break
                end
            end
            % We should probably check the ACLR of each PA using the above.
            
            
            obj.real_channel.normalize();
            %obj.precoder.P = fftshift(obj.precoder.P, 3);
            %% Step 3. Main Downlink
            obj.precoder.update(obj.real_channel.H);
            obj.v21_pre_out = obj.precoder.use(obj.v0_downlink_data);
            obj.v21_pre_out.match_this('time');
            obj.v21_pre_out.normalize_to_this_rms(3);
            obj.v22_bs_out = obj.bs.tx(obj.v21_pre_out);
            %ue_rx = obj.simulated_channel.use(obj.v22_bs_out); % Only used if array are sim.
            obj.v23_ue_rx = obj.ues.rx(obj.v22_bs_out);
            
        end
        
        function report(obj)
            %% What was the beamforming gain?
            % Compare the ue_rx_sigs to the final beamformed data.
            
        end
        
        function plot(obj)
            obj.v3_ue_rx_raw.plot_psd();
            obj.v23_ue_rx.plot_psd();
        end
    end
end

