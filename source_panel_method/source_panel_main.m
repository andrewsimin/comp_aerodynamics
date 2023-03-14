% Script written by Andrew Simin
% Vortex Panel Method

% Program Instructions
% Enter user inputs with specified units and click "ok"
% Only 4 digit NACA airfoils (cambered or symmetric) can be used
% A plot of the panels, control points and boundary points will be
% superimposed on a plot of the NACA airfoil geometry with higher
% resolution
% A plot of the pressure coefficient at the control points will be plotted

clear all
close all

prompt = {'NACA 4 Digit Number', 'Chord Length (m)', 'Angle of Attack (deg)', 'Free Stream Velocity (m/s)', 'Number of Airfoil Data Points (For Plotting Geometric Shape Only)', 'Number of Panel Approximations (Must be Even Number and > 20)'};
title = 'User Inputs';
dims = [1 100];
defaults = {'2412','1','5','1','300','100'}
input_arr = inputdlg(prompt,title,dims,defaults);
NACA = char(input_arr(1));
D1 = str2num(NACA(1,1));
D2 = str2num(NACA(1,2));
D3 = str2num(NACA(1,3));
D4 = str2num(NACA(1,4));
str = string(input_arr)
c = str2num(str(2,1));
alpha = str2num(str(3,1));
V_inf = str2num(str(4,1));
n = str2num(str(5,1));
nop = str2num(str(6,1));

[XU_t YU_t XL_t YL_t X_chord Y_chord] = NACA_4_Digit(D1,D2,D3,D4,c,alpha,n);

[XB,YB,XC,YC,theta,beta,s] = panel_gen(D1,D2,D3,D4,c,alpha,nop);
legend('Boundary Points','Surface','Control Points')

[gamma at] = vortex_strength(XB,YB,XC,YC,theta,nop,s);

%% Velocity and Pressure Coefficients

for i=1:nop;
  smt=0;
  for j=1:nop+1;
    sm=at(i,j)*gamma(j);
    smt=smt+sm;
  end;

  v(i)=(cos(theta(i))+smt);                                                 % v is the velocity parallel to the panel at the control points
  cp(i)=1-(v(i))^2;                                                         % cp is the pressure coeff at the control points
end;

% plot the pressure coefficient at the control points
figure
hold off
plot(XC,cp);
xlabel('Control Point Along Airfoil Surface (m)')
ylabel('C_p (Pa/Pa)')
grid on
set(gca,'ydir','reverse')

% compute normal and axial force coefficients
CN = -cp.*s.*sin(beta);                                                     % Normal force coeff
CA = -cp.*s.*cos(beta);                                                     % Axial force coeff

% lift coefficient
CL = CN.*cosd(alpha) - CA.*sind(alpha);
CL = sum(CL);
