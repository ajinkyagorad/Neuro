% This function is for computing the neuronal current
% Parameters 
% spike_time - the time at which the neurons are spiking
% current    - the current being carried out by output neurons at the
% time of spike

% In practise seen that the max current decays in about 200ms. So the max
% time simulation window is taken to be 200ms

function[final_current] = compute_current(weight,current, previous_time, current_time, spiking_pixels)

    % Defining the current parameters
     
    tau_leak    = 5e-3 ;
    % updating the current values of the neurons which have just spiked
    weights_to_be_concerned = sum((weight(:,spiking_pixels)),2);
    % Computing the current for a total time of 'time_simulation'     
    final_current = current.*exp(-(previous_time - current_time)/tau_leak) + weights_to_be_concerned;
    
    % Final current has been computed
    
    
    
    
   
    
   
        
        
    
    