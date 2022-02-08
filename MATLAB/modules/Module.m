classdef Module < handle
    %MODULE. Main superclass for everything.
    % Enforces a common structure on all modules.
    % Use the create method as a factory method for all the modules.
    %
    % Example:
    %
    
    
    properties
        name
        required_domain
        required_fs
        index
    end
    
    methods (Abstract)
        report(obj);  % All classes should have a report method to print out some details.
    end
    
    methods
        function out_signal = use(obj, input_signal)
            input_signal.match_this(obj(1).required_domain, obj(1).required_fs);
            output_data = obj.subclass_use(input_signal.data);
            sizes = size(output_data);
            n_streams = sizes(1);
            out_signal = Signal(output_data, n_streams, obj.required_domain, ...
                obj.required_fs, input_signal.modulator);
        end
    end
    
    methods (Access=protected)
        function save_inputs_to_obj(obj, vars)
            % Will take the parsed vars and save to object in your
            % subclass.
            % Each field needs to be a valid property of the class.
            fields = fieldnames(vars.Results);
            for i = 1:numel(fields)
                obj.(fields{i}) = vars.Results.(fields{i});
            end
        end
    end
    
    methods (Static)
        function [module_dictionary] = populateModuleDictionary()
            % populate dictionary of leaf modules
            module_dictionary = containers.Map();
            
            % Modulators
            module_dictionary('OFDM') = @(varargin) OFDM(varargin{:});
            
            module_dictionary('UEs') = @(p, i) User(p, i);
            
            % Channels
            module_dictionary('RealChannel') = @(varargin) RealChannel(varargin{:});
            module_dictionary('Quadriga') = @(varargin) Quadriga(varargin{:});
            module_dictionary('LOS') = @(varargin) LOS(varargin{:});
            
            % Precoders
            module_dictionary('PrecoderBypass') = @(varargin) PrecoderBypass(varargin{:});
            module_dictionary('ZF') = @(varargin) ZF(varargin{:});
            module_dictionary('MRT') = @(varargin) MRT(varargin{:});
            
            % DPDs
            module_dictionary('DPD') = @(p, i) DPD(p, i);
            
            % Arrays
            module_dictionary('Sim_Array') = @(varargin) Sim_Array(varargin{:});
            module_dictionary('IRIS') = @(varargin) IRIS(varargin{:});
            module_dictionary('PA') = @(p, i) PA(p, i);
            module_dictionary('GMP') = @(varargin) GMP(varargin{:});
            
        end
        
        function objs = create(category_name, p, n_objs)
            if nargin == 2
                n_objs = 1;
            end
            module_dictionary = Module.populateModuleDictionary();
            moduleConstructor = module_dictionary(p.(category_name).name);
            for i = n_objs:-1:1
                objs(i) = moduleConstructor(p.(category_name), 'index', i);
                objs(i).name = p.(category_name).name;
                objs(i).required_fs = p.(category_name).required_fs;
                objs(i).required_domain = p.(category_name).required_domain;
            end
        end
    end
end