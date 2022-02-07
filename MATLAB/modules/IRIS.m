classdef IRIS < Array
    
    
    properties
        chain_ids  % Which IRIS chains to use
        node_ids   % Which nodes in those iris chains.
        n_antennas
        serials
        tx_freq
        rx_freq
        tx_gain
        rx_gain
        use_hub
        wired_ue
        hub_id
        node
        is_bs
    end
    
    methods
        function obj = IRIS(varargin)
            vars = inputParser;
            valid_constellations = {'BPSK', 'QPSK', '16QAM','64QAM'};
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            validVectorPosNum = @(x) all(isnumeric(x)) && all(x > 0);
            validBool = @(x) islogical(x);
            
            addParameter(vars, 'name', 'IRIS', @(x) any(validatestring(x,{'IRIS'})));
            addParameter(vars, 'required_domain', 'time', @(x) any(validatestring(x,{'time'})));
            addParameter(vars, 'required_fs', 20e6, validScalarPosNum);
            addParameter(vars, 'n_antennas', 8, validScalarPosNum);
            addParameter(vars, 'index', 1, validScalarPosNum);
            addParameter(vars, 'chain_ids', 2, validVectorPosNum);
            addParameter(vars, 'node_ids', 1:8, validVectorPosNum);
            addParameter(vars, 'tx_gain', 80, validScalarPosNum);
            addParameter(vars, 'rx_gain', 60, validScalarPosNum);
            addParameter(vars, 'wired_ue', false, validBool);
            addParameter(vars, 'tx_freq', 3.6e9, validScalarPosNum);
            addParameter(vars, 'rx_freq', 3.6e9, validScalarPosNum);
            addParameter(vars, 'use_hub', true, validBool);
            addParameter(vars, 'is_bs', true, validBool);
            parse(vars, varargin{:});
            obj.save_inputs_to_obj(vars);
            
            obj.set_serial_numbers();
            
            % Error checking.
            assert(obj.n_antennas==numel(obj.serials), 'Number of antennas doesnt match number of serial numbers');
        end
        
        function subclass_tx(obj, data)
            for i=1:obj.n_antennas
                obj.node.sdrtx_single(data(i,:), i);  % Burn data to the BS RAM
            end
            obj.node.sdrtrigger();
        end
        
        function outputArg = subclass_rx(obj,inputArg)
            [rx_vec_iris, data0_len] = node_ue.uesdrrx(N_SAMP); % read data
        end
        
        function out = subclass_measure_noise(obj)
            
        end
        
        function report(obj)
            
            
        end
    end
    
    methods (Access=protected)
        function set_serial_numbers(obj)
            chains = ["", "" ,"" ,"" ,"" ,"" ,"","";...
                "RF3E000087", "RF3E000084", "RF3E000107", "RF3E000110", "RF3E000086", "RF3E000162", "RF3E000127", "RF3E000597"; ...];  % Chain 2
                "RF3E000346", "RF3E000543", "RF3E000594", "RF3E000404", "RF3E000616", "RF3E000622", "RF3E000601", "RF3E000602"; ...        % Chain 3
                "RF3E000146", "RF3E000122", "RF3E000150", "RF3E000128", "RF3E000168", "RF3E000136", "RF3E000213", "RF3E000142"; ...        % Chain 4
                "RF3E000356", "RF3E000546", "RF3E000620", "RF3E000609", "RF3E000604", "RF3E000612", "RF3E000640", "RF3E000551"; ...        % Chain 5
                "RF3E000208", "RF3E000636", "RF3E000632", "RF3E000568", "RF3E000558", "RF3E000633", "RF3E000566", "RF3E000635";...           % Chain 6
                "RF3E000164", ""          , ""          , ""          , ""          ,  ""         , ""          , ""]; %, "RF3D000016"]; %, "RF3E000180"];
            obj.serials = chains(obj.chain_ids, obj.node_ids);
            obj.serials = obj.serials(:);
            
            if obj.use_hub
                obj.hub_id = "FH4B000019";
            else
                obj.hub_id = [];
            end
        end
        
        function setup_board(obj)
            sdr_params = struct(...
                'id', obj.serials, ...
                'n_sdrs', obj.n_antennas, ...
                'txfreq', obj.tx_freq, ...
                'rxfreq', obj.rx_freq, ...
                'txgain', obj.tx_gain, ...
                'rxgain', obj.rx_gain, ...
                'sample_rate', obj.required_fs, ...
                'n_samp', N_SAMP, ...          % number of samples per frame time.
                'n_frame', N_FRM, ...
                'tdd_sched', bs_sched, ...     % number of zero-paddes samples
                'n_zpad_samp', N_ZPAD_PRE ...
                );
            
            obj.node = iris_py(sdr_params, obj.hub_id);
            
            if obj.is_bs
                obj.node.sdrsync();                 % Synchronize delays only for BS
                obj.node.sdr_setupbeacon();   % Burn beacon to the BS RAM
            else
                obj.node.sdr_configgainctrl();
                obj.node.sdrrxsetup();
                
                if ~obj.wired_ue
                    obj.node.sdr_setcorr();   % activate correlator
                end
                obj.node.sdr_activate_rx();   % activate reading stream
            end
        end
    end
end

