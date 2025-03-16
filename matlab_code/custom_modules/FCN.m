classdef FCN
    properties
        
    end
    methods
        function y = forward(x)
            ... # TODO
        end
    end
end

% simpler
net = importNetworkFromPyTorch("results/best_model/model_weights.pt")