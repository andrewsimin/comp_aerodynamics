% Function written by Andrew Simin
% This function calculates the geometry factor of a 4 digit NACA 
% airfoil at a specified angle of attack for implementation of a 
% linear vortex panel method

%function [XB,YB,XC,YC,theta,beta,s] = panel_gen(D1,D2,D3,D4,c,alpha,nop)

% Inputs -- D1 -- First digit of NACA 4 digit airfoil
%        -- D2 -- Second digit of NACA 4 digit airfoil
%        -- D3 -- Third digit of NACA 4 digit airfoil
%        -- D4 -- Fourth digit of NACA 4 digit airfoil
%        -- c -- Chord length [m] [ft]
%        -- alpha -- Angle of attack [deg]
%        -- nop -- Number of panels for approximation
% Outputs -- XB -- X-coordinates of boundary points 
%         -- YB -- Y-coordinates of boundary points
%         -- XC -- X-coordinates of control points
%         -- YC -- Y-coordinates of control points
%         -- theta -- angle of r with respect to the x-axis
%         -- beta -- angle between normal vector and free stream velocity
%         vector
%         -- s -- length of panel
% starting from lower surface TE to upper surface TE

function [XB,YB,XC,YC,theta,beta,s] = panel_gen(D1,D2,D3,D4,c,alpha,nop)

t = str2double(append(num2str(D3),num2str(D4)))/100;                        % define thickness from symbolic imputs 
m = D1*0.01;                                                                % m coefficient
p = D2*0.1;                                                                 % p coefficient                                                                      % adjust n to accout for for loop data


%% Calculate Camber Line
x=linspace(0,2*pi,nop+1);
x=c*0.5*(cos(x)+1);

for i = 1:nop+1                                                             % loop from 1 to # of data points + 1
    if (0 <= x(i)) && (x(i) <= p*c);                                        % piecewise constraint
        y_c(i) = [(m*x(i)/p^2)*(2*p - x(i)/c)];                             % thickness equation
    elseif (p*c <= x(i)) && (x(i) <= c);                                    % piecewise constraint
        y_c(i) = [m*((c-x(i))/(1-p)^2)*(1 + (x(i)/c) - 2*p)];               % thickness equation
    end                                                                     % end conditional statements
end                                                                         % end for loop

y_c(nop+1) = 0;
x_t = x*cosd(alpha) + y_c*sind(alpha);                                      % transform x coordinates for AOA
y_c_t = -x*sind(alpha) + y_c*cosd(alpha);                                   % transform y coordinates for AOA                   

%x_t = x*cosd(alpha)                                     % transform x coordinates for AOA
%y_c_t = y_c*cosd(alpha);                                   % transform y coordinates for AOA                   


%% Calculate Theta by dyc/dx

for i = 1:nop+1                                                             % loop from 1 to # of data points + 1
    if (0 < x(i)) && (x(i) < p*c);                                          % piecewise constraint
        dy_c(i) = [(2*m/p^2)*(p-x(i)/c)];                                   % thickness equation
    elseif (p*c < x(i)) && (x(i) < c);                                      % piecewise constraint
        dy_c(i) = [((2*m)/((1-p)^2))*(p-x(i)/c)];                           % thickness equation
    end                                                                     % end conditional statement
end                                                                         % end for loop
dy_c(nop+1) = 0;
theta = atan2d(dy_c,1);                                                     % theta for airfoil coordinates
    
%% Calculate Symmetrical Airfoil Coordinates

y_t = c*5*t*[0.2969.*sqrt(x/c) - 0.1260.*x/c - 0.3516.*(x/c).^2 + 0.2843.*(x/c).^3 - (0.1036)*(x/c).^4]; % symmetrical airfoil equation

%% Calculate Camber Airfoil Coordinates
XU = (x - y_t.*sind(theta));                                                % x-coordinates for upper surface
XL = (x + y_t.*sind(theta));                                                % x-coordinates for lower surface
YU = (y_c + y_t.*cosd(theta));                                              % y-coordinates for upper surface
YL = (y_c - y_t.*cosd(theta));                                              % y-coordinates for lower surface

%% Transform Coordinates   

XU_t = XU*cosd(alpha) + YU*sind(alpha);                                     % transform XU to represent AOA                 
XL_t = XL*cosd(alpha) + YL*sind(alpha);                                     % transform XL to represent AOA
YU_t = -XU*sind(alpha) + YU*cosd(alpha);                                    % transform YU to represent AOA
YL_t = -XL*sind(alpha) + YL*cosd(alpha);                                    % transform YL to represent AOA

plot(XU)

XU_t = XU_t(1,1:nop/2);                                                     % eliminate repetition of coordinates
XL_t = XL_t(1,1:(nop/2));                                                   % eliminate repetition of coordinates
YU_t = YU_t(1,1:nop/2);                                                     % eliminate repetition of coordinates
YL_t = YL_t(1,1:(nop/2));                                                   % eliminate repetition of coordinates

%% Combine Vectors to Find Boundary Points

XB = zeros(nop+1);                                                          % populate x boundary point matrix
XB = XB(1,1:end);                                                           % turn n x n matrix to row vector
XB(1,1:(nop/2)) = (XL_t);                                                   % replace first half of XB matrix with lower surface coordinates
XB(nop+1) = 0                                                               % define origin
XB(1,(nop/2)+2:end) = flip(XU_t);                                           % replace second half of XB matrix with upper surface coordinates

YB = zeros(nop+1);                                                          % populate y boundary point matrix
YB = YB(1,1:end);                                                           % turn n x n matrix to row vector
YB(1,1:(nop/2)) = (YL_t);                                                   % replace first half of YB matrix with lower surface coordinates
YB(nop+1) = 0                                                               % define origin
YB(1,(nop/2)+2:end) = flip(YU_t);                                           % replace second half of YB matrix with upper surface coordinates

%% Control Points

for i = 1:nop                                                               % loop from 1 to # of panels
XC(i) = 0.5*(XB(i) + XB(i+1));                                              % x coordinate of control point
YC(i) = 0.5*(YB(i) + YB(i+1));                                              % y coordinate of control point
end

%% Plot

plot(XB,YB,'ro')
hold on
plot(XB,YB,'k')
hold on
plot(XC,YC,'xg')
axis equal

%% Panel Geometry

for i = 1:nop                                                               % loop from 1 to # of panels
s(i) = sqrt((XB(i+1) - XB(i))^2 + (YB(i+1) - YB(i))^2);                     % distance from boundary point to cotrol point  
thet(i) = atan2d((YB(i+1) - YB(i)),(XB(i+1) - XB(i)));                      % angle of panel to x axis                              % return positive angle
end

theta = deg2rad(thet)                                                       % deg to rad
beta = theta + pi/2                                                         % angle between normal vector and free stream velocity 
alpha = deg2rad(alpha)                                                      % deg to rad

end