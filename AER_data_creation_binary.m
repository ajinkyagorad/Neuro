%% AERdata creation
%+1 = Spiking due to increase in intensity (first synapse fired)
%0 = no spiking
%-1  = spiking due to decrease in intensity (second synapse fired)
% uncomment 'flipud' for reverse direction
N = 16; % square grid frame length
d = 3; % size of the object/ball
X = 1:N;
slope =2;
invert = 1;
normTh = 0.23; % Normalised threshold for AER spike data generation
lastFrame = [];
for i = 0:N+d
    x = X-i;
    y = X-slope*i;            % set slope of movement /variable transformation 
    countourX = (heaviside(x-d)-heaviside(x)).*(x+1).*(x-d-1);
    countourY = (heaviside(y-d)-heaviside(y)).*(y+1).*(y-d-1);
    frame = countourY'*countourX;   % generate the ball with NxN grid
    M =max(abs(frame(:)));
    if(M==0)break;end;      % object  goes out of bound, no change, stop
    frame = frame/M;        % Normalise the values
    if(invert)frame=flipud(flipud(frame)');end;
    
    if(i>0)
        gradient = frame-lastFrame;
        AERdata = sign(gradient.*floor(abs(gradient/normTh)));     % threshold the image and generate AER data 
        image(AERdata,'CDataMapping','scaled')    % animate display
        colorbar     %display colormap
        pause(.200);
    end
    lastFrame = frame;
    %frame has the AER data; 
end