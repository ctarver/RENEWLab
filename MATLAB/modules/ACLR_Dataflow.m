classdef ACLR_Dataflow < handle
    %ACLR_DATAFLOW
    
    properties
        p
        n_users
        n_ants
        bs
        ues
        simulated_channel  % This is used when bs and ues are in sim mode.
        real_channel       % Will hold the RealChannel measured OTA.
        precoder
        v0_downlink_data
        v1_pre_out
        v2_ue_rx
    end
    
    methods
        function obj = ACLR_Dataflow(p)
            %ACLR_DATAFLOW
            obj.p = p;
            obj.n_users = p.n_users;
            obj.n_ants = p.n_ants;
            
            obj.bs = Module.create('bs', p);
            obj.ues = Module.create('ue', p);
            
            if strcmp(obj.bs, 'Simulation')
                obj.simulated_channel = Module.create('sim_channel', p);
            end
            
            obj.real_channel = Module.create('channel', p);
            obj.precoder = Module.create('precoder', p);
        end
        
        function run(obj)
            %% Make Data
            obj.v0_downlink_data = Signal.make_ofdm(obj.n_users, obj.p.mod);
            
            
            %% Run Downlink from each TX to the UEs
            % How steady are these channels.
            
            % Make a placeholder signal.
            ue_rx_sigs = Signal.make_ofdm(obj.n_ants, obj.p.mod);
            ue_rx_sigs.match_this('domain', 'time');
            ue_rx_sigs.data = zeros(size(ue_rx_sigs.data));
            
            for i_bs = 1:obj.n_ants
                tx_sig = obj.v0_downlink_data.zero_all_but(i_bs);
                obj.bs.tx(tx_sig)  % Need to check on the trigger.
                ue_rx_sigs.data(:) = obj.ues.rx();
            end
            
            %% We should probably check the ACLR of each PA using the above.
            
            
            %% Xcorr and Calculate Downlink Channel
            obj.real_channel.learn(obj.v0_downlink_data, ue_rx_sigs);
            
            %% Main Downlink
            obj.precoder.update(obj.real_channel.H);
            obj.v1_pre_out = obj.precoder.use(obj.v0_downlink_data);
            obj.bs.tx(obj.v1_pre_out);
            obj.v2_ue_rx = obj.ue.rx();
            
        end
        
        function report(obj)
            %% What was the beamforming gain?
            % Compare the ue_rx_sigs to the final beamformed data.
            
        end
        
        function plot()
            
        end
    end
end

