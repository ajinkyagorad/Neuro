% Using retina sensor as input stimulus for motion sensing
% The input stimulus is a ball trajectory 
% Having a 16 X 16 pixel data as the input stimulus
clc; clear

%% ===================================== Defining the neuronal structure ===============================================
% The input stimulus matrix where the stimulus is the trajectory of the
% ball

% the output layer which will give the direction in a part of the
% trajectory therefore detecting the motion of the ball

N = 48;                            % The number of output neurons
output_layer = zeros(N,1);         % A column matrix as the output layer which has 48 neurons


% The first layer is the input AER data
% two matrices of 16 X 16 pixels due to presence of two states of the
% pixels
% The 'ON' state and the 'OFF' state
% 'ON' - Increasing intensity of the pixels
% 'OFF' - Decreasing intensity of pixels
%

AER_input_pixels = ones(16,16,2);   % A spike will be represented as a one
previous_AER_input_pixels_time = zeros(N,16*16*2);  % to store the previous spike time of pixels
previous_spiking_pixels_time = zeros(N,1); % the previous neurons which have spiked in the implementation of "compute current"

% Defining the weight matrix

w_max = 2000;
w_min = 100;

% initialising weight with mean 800 and standard deviation 20

weight = 800 + 160.*randn(16,16,2,N);

% learning rate 

alpha_plus = 100 + 20.*randn(N,16*16*2);
alpha_minus = 50 + 10.*randn(N,16*16*2);

% damping rate
% why not finite?? always add a contsant value
beta_plus = 0;
beta_minus = 0;

% Time instants for the whole simulation
time_simulation =2000e-3;
time_step = 1e-3;
times = 0:time_step:time_simulation;

% neuronal current parameters

I_threshold = 40000;
tau_leak    = 5e-3 ;
% input current for the 48 output neurons
neuronal_current = zeros(48,1);

% refractory and lateral inhibition parameters

Trefrac = 10e-3;
Tltp    = 2e-3;
Tinhibit = 1.5e-3;

% for tracking the inhibitory and refractory time periods of the neurons
inhibitory_flag = zeros(N,1);
refractory_flag = zeros(N,1);
Time_inhibit = zeros(N,1);




%% ====================================================== Simulation starts =====================================================
flag = 0;
weight = reshape(weight,16*16*2,N)';      % for convenience of using weight matrix
for i=1:length(times) - 1   
 
      % finding out which neurons spiked
      % the AER pixels
      % reshaping the AER input pixels as a column vector
      
      temp_AER_input_pixels = reshape(AER_input_pixels,16*16*2,1);      % size of matrix is 16 X 16 X 2
      spiking_pixels        = find(temp_AER_input_pixels);
      % Computing new current in the output neurons after spikes have been
      % received in the AER neurons
      
          
      lagging_neurons = find(Time_inhibit);                              % neurons which have to be inhibited
      Time_inhibit(lagging_neurons) = Time_inhibit(lagging_neurons) - 1; % decrementing by 1 in step
      
      if(~isempty(spiking_pixels))
         current_time         = i*time_step; 
         new_neuronal_current = compute_current(weight,neuronal_current,previous_spiking_pixels_time, current_time, spiking_pixels);
         indices              = find(Time_inhibit);                      % neurons which are inhibited due to lateral inhibition or refraction
         new_neuronal_current(indices,1) = neuronal_current(indices,1);               % resetting the inhibited neurons to its previous value of current
         neuronal_current     = new_neuronal_current ;
         % Replacing the current window with a new current after a spike
         % has been elicited by the AER neurons
         
         % Fiinding out the neurons which will elicit a spike after reviseing
         % the neuronal current
      
         threshold = find(neuronal_current>=I_threshold);     % checking if the current has gone above threshold 
         output_layer(threshold',1) = 1;                      % the spiking of output neurons
         neuronal_current(threshold') = 0;                    % resetting the current to zero for neurons which spiked
         previous_spiking_pixels_time(threshold') = i*time_step;
         previous_AER_input_pixels_time(threshold',spiking_pixels) = i*time_step;
         
         % Taking care of the refractory period and lateral inhibition
      
         refractory_flag(threshold',1)        = 1;                                                                  % flag for the start of refractory period
         Time_inhibit(:,1)                    = Time_inhibit(:,1) + Tinhibit/time_step;                        % Implementing lateral inhibition , adding 1 so that it can be decremented in the following code 
         Time_inhibit(threshold',1)           = Time_inhibit(threshold',1) + (Trefrac - Tinhibit)/time_step ;      %undoing the effect of time inhibit of spiking neurons
      
         inhibit                              = 1:N;
         inhibit(threshold')                  = [];
         inhibitory_flag(inhibit,1)           = 1;                                                                 % neurons which did not spike are inhibited
          
         % Weight updation rules for the synapses which have elicited a
         % spike 
         total_neurons             = 1:N;        % All indices of output neurons
         total_neurons(threshold)  = [];         % indices of neurons which have'nt spiked
         weights_to_be_incremented = (previous_AER_input_pixels_time >= ((i*time_step) - Tltp));
         weight_to_be_incremented(total_neurons,:) = 0;     % all synapses have to be depreciated for which output neurons have'nt spiked
         
         % All other neurons have to be decremented
         weights_to_be_decremented = ~(weights_to_be_incremented) ; 
         
         % Poteniating the neurons which are closely grouped together
         delta_weight_add = alpha_plus(weights_to_be_incremented).*exp(-(beta_plus)*((weight(weights_to_be_incremented) - w_min)/(w_max - w_min)));
         weight(weights_to_be_incremented) = weight(weights_to_be_incremented) + delta_weight_add; 
         
         % Depreciating the weights of neurons which are not related to the
         % current stimulus
         delta_weight_sub = alpha_minus(weights_to_be_decremented).*exp(-(beta_minus)*((weight(weights_to_be_decremented) - w_min)/(w_max - w_min)));
         weight(weights_to_be_decremented) = weight(weights_to_be_decremented) - delta_weight_sub; 
         
          flag = flag + 1;     
         
      end
   
         
         
         
        
      
     
     
      
     
      
     
      
      
     
              
end












