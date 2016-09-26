%CodeStart-----------------------------------------------------------------
%Resetting MATLAB environment
    close all;
    clear all;
    clc;
%Declaring environmental variable
    r_ball=20;       %Ball's radius
%Declaring ball's initial condition
    initpos=-32 - r_ball;     %Ball's initial vertical position
    initvel= 20;      %Ball's initial vertical velocity

%Declaring animation timestep
    dt=0.0125;      %Animation timestep

 % For the mesh pixels   
    x= -32:4:32;
   
    [X,Y] = meshgrid(x);
    
%Drawing first frame
    rectangle('Position',[initpos,-r_ball/2,r_ball,r_ball],...
              'Curvature',[1,1],...
              'FaceColor','b');
    hold on; 
    plot(X,Y,'b-x',X',Y','b-x');
%Executing animation
    pos=initpos+r_ball;             %Ball's current vertical position
    vel=initvel;                    %Ball's current vertical velocity
    

    while 1
        %Declaring time counter
        t_loopstart=tic();
        %Updating ball's condition
        pos=pos+(vel*dt);           %Ball's current horizontal position
        %Adjusting ball's velocity if the ball is hitting the floow
        if pos>32
            break;      %Balls' current vertical velocity
        end
        %Clearing figure
        clf;
        %Drawing frame
        
        rectangle('Position',[pos,-r_ball/2,r_ball,r_ball],...
                  'Curvature',[1,1],...
                  'FaceColor','b');
        hold on;
         plot(X,Y,'b-x',X',Y','b-x');
        %Preserving axes
        axis([-50,50,0,2*r_ball]);
        axis('equal');
        axis('on');
        %Pausing animation
        el_time=toc(t_loopstart);
        disp(['Elapse time : ',num2str(el_time),' seconds']);
        pause(0.0125);
    
    end
%CodeEnd------------------