%% AERdata creation
%+1 = Spiking due to increase in intensity (first synapse fired)
%0 = no spiking
%-1  = spiking due to decrease in intensity (second synapse fired)
% uncomment 'flipud' for reverse direction
N = 16; % square grid frame length
d = 3; % size of the object/ball
X = 1:N;
slope =0;
invert = 0;
normTh = 0.23; % Normalised threshold for AER spike data generation
lastFrame = [];
for j = -N:1:N
    display(j)
    for i = 0:N+d
        x = X-i;
        y = X-slope*i;            % set slope of movement /variable transformation 
        if(j>   0)
            x=x-abs(j);
            y=y+abs(j);
        else
            x = x+abs(j);
            y = y-abs(j);
        end;
        countourX = (heaviside(x-d)-heaviside(x)).*(x+1).*(x-d-1);
        countourY = (heaviside(y-d)-heaviside(y)).*(y+1).*(y-d-1);
        frame = countourY'*countourX;   % generate the ball with NxN grid
        M =max(abs(frame(:)));
        if(M==0)if(i==0)display('out of bounds');end;break;end;      % object  goes out of bound, no change, stop
        frame = frame/M;        % Normalise the values
        if(invert)frame=flipud(flipud(frame)');end;

        if(i>0)
            gradient = frame-lastFrame;
            AERdata = sign(gradient.*floor(abs(gradient/normTh)));     % threshold the image and generate AER data 
            image(frame,'CDataMapping','scaled')    % animate display
            colorbar     %display colormap
            pause(.020);
        end
        lastFrame = frame;
        %frame has the AER data; 
    end
    
end