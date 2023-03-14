% Function written by Andrew Simin
% This function calculates the coefficients and vortex strength for a
% linear vortex method using a NACA 4 digit airfoil

% [gamma at] = vortex_strength(XB,YB,XC,YC,theta,nop,s)

% Inputs -- XB -- X Boundary Points
%        -- YB -- Y Boundary Points
%        -- XC -- X Control Points
%        -- YC -- Y Control Points
%        -- theta -- Angle of r with respect to the x-axis
%        -- nop -- Number of panels
%        -- s -- Length of panel 
% Outputs -- gamma -- Vortex strength
%         -- at -- Tangential component of vortex strength

function [gamma at] = vortex_strength(XB,YB,XC,YC,theta,nop,s)

%% Calculate Coefficients

an=zeros(nop+1,nop+1);                                                      % populate an matrix
at = zeros(nop+1,nop+1);                                                    % populate at matrix

for i = 1:nop;                                                              % loop from 1 to number of ith panels
  for j = 1:nop;                                                            % loop from 1 to number of jth panels
    if i == j;                                                              % if i = j, return magnitude of vortex source influence coeff.
      cn1(i,j) = -1;                                                        % normal coefficient 1
      cn2(i,j) = 1;                                                         % normal coefficient 2
      ct1(i,j) = pi/2;                                                      % tangential coefficient 1
      ct2(i,j) = pi/2;                                                      % tangential coefficient 2
    else;
      A = -(XC(i) - XB(j))*cos(theta(j)) - (YC(i) - YB(j))*sin(theta(j));   % A coefficient
      B = (XC(i) - XB(j))^2 + (YC(i) - YB(j))^2;                            % B coefficient
      C = sin(theta(i) - theta(j));                                         % C coefficient
      D = cos(theta(i) - theta(j));                                         % D coefficient
      E = (XC(i) - XB(j))*sin(theta(j)) - (YC(i) - YB(j))*cos(theta(j));    % E coefficient
      F = log(1 + (s(j)^2 + 2*A*s(j))/B);                                   % F coefficient
      G = atan2((E*s(j)),(B+A*s(j)));                                       % G coefficient
      P = (XC(i) - XB(j))*sin(theta(i) - 2*theta(j)) + (YC(i) - YB(j))*cos(theta(i) - 2*theta(j)); % P coefficient
      Q = (XC(i) - XB(j))*cos(theta(i) - 2*theta(j)) - (YC(i) - YB(j))*sin(theta(i) - 2*theta(j)); % Q coefficient
      cn2(i,j) = D + 0.5*Q*F/s(j) - (A*C+D*E)*G/s(j);                       % normal coefficient 1
      cn1(i,j) = 0.5*D*F + C*G - cn2(i,j);                                  % normal coefficient 2
      ct2(i,j) = C + 0.5*P*F/s(j) + (A*D - C*E)*G/s(j);                     % tangential coefficient 1
      ct1(i,j) = 0.5*C*F - D*G - ct2(i,j);                                  % tangential coefficient 2
    end;
  end;
end;

%% Solve for Votex Strength Gamma

for i= 1:nop;                                                               % loop from 1 to number of panels
  an(i,1) = cn1(i,1);                                                       
  an(i,end) = cn2(i,nop);
  at(i,1) = ct1(i,1);
  at(i,end) = ct2(i,nop);
  for j = 2:nop;
    an(i,j) = cn1(i,j) + cn2(i,(j - 1));
    at(i,j) = ct1(i,j) + ct2(i,(j - 1));
  end;
end;

rhs = sin(theta);                                                           % Right hand side of equation

%% Kutta Condition
an(end,1) = 1;
an(end,end) = 1;
rhs(nop+1) = 0;
for j = 2:nop;
  an(end,j)=0;
end;

rhs = rhs(:);                                                               % column vector

%% Solve Linear System for Vortex Strength, Gamma

gamma = an \ rhs;
at = at
end