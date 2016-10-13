%% Single direction detection Neural Network
% Input AER data  of NxM ( 16X4 pixels ) of moving target of size less than
% 4x4 pixels.
clc;clear
% Number of Output Neurons  = 12
Nout = 12;


N = 16;
M = 16;
d = 9;
X = 1:N;
Y = 1:M;
AER_TH =20; %AER data threshold; found from plot of gradient
Data = zeros(M,N,N-d);            % twice the number  !!
DataSmooth = zeros(M,N,N-d);
RawImgInput = zeros(M,N,N-d);
for i=0:1:N-d
    x = X-i;y=Y;
    x = mod(x,N);
    y = mod(y,M);
    countourX = (heaviside(x-d)-heaviside(x)).*(x).*(x-d-1);
    countourY = (heaviside(y-d)-heaviside(y)).*(y).*(y-d-1);
    frame  = countourY'*countourX;
    if(i>0)
        frame = frame+wgn(M,N,20);
        gradient = frame-frame_;
        gradientD = sign(gradient.*floor(abs(gradient/AER_TH)));
        RawImgInput(:,:,i) = frame;
        DataSmooth(:,:,i) = gradient;
        Data(:,:,i)= gradientD;
        image(gradientD,'CDataMapping','scaled')
    end
    
    frame_ = frame;
    pause(0.01)
end

print('Data Generated');
%% Generate Test Object for STDP checking
%     countourX = (heaviside(X-3*N/4)-heaviside(X-N/4)).*(X+1).*(X-3*N/4-1);
%     countourY = (heaviside(Y-3*M/4)-heaviside(Y-M/4)).*(Y+1).*(Y-3*M/4-1);
%     frame  = countourY'*countourX;
%     image(frame,'CDataMapping','Scaled')
%% STDP Network Learn
%Network size NxM, with input size of N-d (d = target size) repeated K
%times

output_layer = zeros(Nout,1);

dataSize = size(Data);
NumFrames = dataSize(3);

AER_input = ones(N,M,NumFrames); %A spike will be represented as a (1) or (-1)
previous_AER_input_pixels_time = zeros(Nout,N*M*2);% to store the previous spike time of pixels
previous_spiking_pixels_time = zeros(Nout,1);% the previous neurons which have spiked in the implementation of "compute current"


% Defining the weight matrix

w_max = 2000;
w_min = 100;

% initialising weight with mean 800 and standard deviation 20

weight = 500 + 100.*randn(N,M,2,Nout);

% learning rate 

alpha_plus =  100+ 10.*randn(Nout,N*M*2);
alpha_minus = 20+ 5.*randn(Nout,N*M*2);

% damping rate
% why not finite?? always add a constant value
beta_plus = 0.5;
beta_minus = 0.5;

% Time instants for the whole simulation
time_simulation =2000e-3;
time_step = 1e-3;
times = 0:time_step:time_simulation;

% neuronal current parameters

I_threshold = 2E5;
tau_leak    = 50e-3 ;

% input current for the 48 output neurons
neuronal_current = zeros(Nout,1);

% refractory and lateral inhibition parameters

Trefrac = 10e-3;
Tltp    = 2e-3;
Tinhibit = 10*time_step;

% for tracking the inhibitory and refractory time periods of the neurons
inhibitory_flag = zeros(Nout,1);
refractory_flag = zeros(Nout,1);
Time_inhibit = zeros(Nout,1);

% Recording spike data
spike_dat = zeros(Nout,1);
%% ====================================================== Simulation starts =====================================================
flag = 0;
weight = reshape(weight,N*M*2,Nout)';      % for convenience of using weight matrix

for i=0:length(times) % changed index i to start from 0 rather than 1
        
        currentDataIndex = 1+mod(i,NumFrames);
        AER_input_pixels = Data(:,:,currentDataIndex);
        %For debugging purposes 
        
        %AER_input_pixels = 100*Data(:,:,1);
        
        

      % finding out which neurons spiked
      % the AER pixels
      % reshaping the AER input pixels as a column vector
      temp_AER_input_pixels = reshape(AER_input_pixels,N*M,1);      % size of matrix is N*M, 
      spiking_pixels        = find(temp_AER_input_pixels);  % get  all spikes +1 and -1
      % Computing new current in the output neurons after spikes have been
      % received in the AER neurons
      lagging_neurons = find(Time_inhibit);                              % neurons which have to be inhibited
      Time_inhibit(lagging_neurons) = Time_inhibit(lagging_neurons) - 1; % decrementing by 1 in step
      
      if(~isempty(spiking_pixels))
         current_time         = i*time_step; 
         new_neuronal_current = compute_current(weight,neuronal_current,previous_spiking_pixels_time, current_time, spiking_pixels,tau_leak);
         indices              = find(Time_inhibit>0);                      % neurons which are inhibited due to lateral inhibition or refraction
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
         spikes = zeros(Nout,1);
         spikes(threshold) = 1;
         spike_dat = [spike_dat spikes];
         % Taking care of the refractory period and lateral inhibition
         % threshold contains output neurons which spiked
         refractory_flag(threshold',1)        = 1;                                                                  % flag for the start of refractory period
         %adding  t0 to all
         Time_inhibit(:,1)                    = Time_inhibit(:,1) + Tinhibit;                        % Implementing lateral inhibition , adding 1 so that it can be decremented in the following code 
         % subtractingg t0 from those who had threshold
         Time_inhibit(threshold',1)           = Time_inhibit(threshold',1) + (Trefrac - Tinhibit) ;      %undoing the effect of time inhibit of spiking neurons
      
         % set inhibitory neurons
         inhibit                             = 1:Nout;
         inhibit(threshold')                  = [];                                         % remove those entries those spiked
         inhibitory_flag(inhibit,1)           = 1;                                          % neurons which did not spike are inhibited
          
         % Weight updation rules for the synapses which have elicited a
         % spike 
         total_neurons             = 1:Nout;        % All indices of output neurons
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
         delta_weight_sub = alpha_minus(weights_to_be_decremented).*exp(-(beta_minus)*(w_max-(weight(weights_to_be_decremented))/(w_max - w_min)));
         weight(weights_to_be_decremented) = weight(weights_to_be_decremented) - delta_weight_sub; 
         
         flag = flag + 1;
    end
      figure(1)
      %%%% Display Realtime parameters
     % if(mod(i,NumFrames)==0)
      weight_ = reshape(weight,N,M,2*Nout); % there are twice the number of weights
        subplot(2,2,1);
        surf(Data(:,:,currentDataIndex))
        title('Input Data')
        
        subplot(2,2,2)
        plot(neuronal_current);
        title(['OutputLayer(' int2str(i) ')']);
        
        subplot(2,2,3);
        image(DataSmooth(:,:,currentDataIndex),'CDataMapping','Scaled')
        title('Smooth AER data (before threshold)')
        
        subplot(2,2,4);
        image(RawImgInput(:,:,currentDataIndex),'CDataMapping','Scaled')     % the latest frame(Nth) is displayed,above from data N-1 th, Nth frame
        title('Raw Image Video Input')
        figure(2) 
        for k=1:4
            subplot(2,2,k)
            surf(weight_(:,:,k));
            %image(weight_(:,:,k)','CDataMapping','Scaled');
            title(['Weight ' int2str(k)])
        end
            
        %image(weight_(:,:,1),'CDataMapping','Scaled');colormap;
        pause(0.1);
        i % print current iteration index
      %end
              
end
