%% Single direction detection Neural Network
% Input AER data  of NxM ( 16X4 pixels ) of moving target of size less than
% 4x4 pixels.
% 12 output Neurons for two directions

N = 16
M = 4
d = 4
X = 1:N;
Y = 1:M;

for i=1:1:N-d
    x = X-i;y=Y;
    countourX = (heaviside(x-d)-heaviside(x)).*(x+1).*(x-d-1);
    countourY = (heaviside(y-d)-heaviside(y)).*(y+1).*(y-d-1);
    frame  = countourY'*countourX;
    image(frame,'CDataMapping','scaled')
    pause(0.2)
end
