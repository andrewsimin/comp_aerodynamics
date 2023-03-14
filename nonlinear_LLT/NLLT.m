% Program Written by Andrew Simin
% This program implements a numerical lifting line theory for prediction of
% post stall aerodynamic coefficients.
% All underlines are parameters to be input by the user
% Flow properties is a section for comparison to experimental data, otherwise the
% only input needs to be velocity
% Program outputs lift and induced drag for a finite wing

clear all
close all
clc

%% Load section lift data
data = load('0015_Section_UM_Extended.txt'); % load section lift data
AOA_2D = data(1:end,1)'; % experimental angle of attack
cl_2D = data(1:end,2)'; % experimental lift coefficient
[AOA_2D unique_indices] = unique(AOA_2D) % return unique indices for interpolation
cl_2D = cl_2D(unique_indices) % return unique values of cl
clear unique_indices

%% Input data
k = 140; % number of stations
noi = 400; % number of circulation update iterations
damping_factor = 0.01; % damping factor

%% Define wing properties
b = 0.3515; % wing span
lambda = 1; % taper ratio
AR = 2.768; % aspect ratio % planform area
s = (b^2)/AR; % planform area
rc = 2*s/(b*(1+lambda)); % root chord
tc = rc*lambda; % tip chord

%% Flow properties
Re_test = 350000 % test reynolds number (optional)
nu_inf = 1.4207E-5 % kinematic viscosity [m^2/s]
char_chord_length = 0.127 % chord length [m] (optional)
V_inf = 39.1531 % free stream velocity [m/s]
%V_inf = (Re_test*nu_inf)/char_chord_length % velocity based on experimental reynolds # [m/s]

%% Loop through geometric angles of attack based on input data
for a = 1:length(AOA_2D)
AOAd_geo = AOA_2D(1,a); % geometric angle of attack [deg]
AOAr_geo = deg2rad(AOAd_geo); % geometric angle of attack [rad]

%% Anderson numerical implementation
%% 1) Discretize wing
y = 0:k; % number of y coordinates (k+1)
dy = b/k; % distance between spanwise coordinates
y = y'*dy - b/2; % y coordinates
chord = ((2*s)/((1+lambda)*b))*((1-((1-lambda)/b)*2*abs(y))) % chord length along span
%chord = sqrt(1-(2*y/(b)).^2)

%% 2. Asuume elliptical lift distribution
gamma = 1 * sqrt(1-(2*y/b).^2); % assume elliptical distribution
% allocate storage
dgamma_dy = zeros(k+1,1); % allocate storage for dgamma_dy
I = zeros(k+1,1); % allocate storage for integral
f = zeros(k+1); % allocate storage for function inside integral
Error = zeros(noi,1); % allocate storage for error between new gamma and old gamma
for i = 1:noi

%% 3. Compute induced angle of attack
for t = 1:k/2
dgamma_dy(t) = (gamma(t+1) - gamma(t)) / dy; % calculate dgamma_dy for negative y values
end
dgamma_dy(k/2+1) = 0; % define gamma = 0 at center span
for t = k/2+2 : k+1
dgamma_dy(t) = (gamma(t) - gamma(t-1)) / dy; % calculate dgamma_dy for positive y values
end
% Calculation of the integral
for n=1:k+1
for j=1:k+1
if j == n
f(n,j) = 0; % average value between j+1 and j-1 when n = j
else
f(n,j) = dgamma_dy(j)/(y(n)-y(j)); % integrand
end
end
%I(n) = trapz(y,f(n,:)); % numerically integrate
I(n) = simpsons(f(n,:),y(1,1),y(end,1),[]);
end
AOAiR = I ./ (4*pi*V_inf); % induced AOA [rad]

%% 4. Commpute alpha effective
AOAEff = AOAd_geo*pi/180 - AOAiR;

%% 5. Compute cl
%cl = 2*pi*(alphaEff - alpha0liftRad); For linear lifting line
cl = interp1(AOA_2D,cl_2D,rad2deg(AOAEff)); % interpolate section cl with effects of downwash

%% 6. New lift distribution
gammaNew = 1/2 * V_inf .* chord .* cl; % kutta joukowski circulation
gammaNew(1) = 0; % set circulation at tip to 0
gammaNew(end) = 0; % set circulation at tip to 0
if abs(gammaNew - gamma) < 0.01
break
end

%% 7. Gamma update
gamma = gamma + damping_factor*(gammaNew - gamma); % update circulation
end

%% 10. Last computations
CL = 2 / ((V_inf)*s) * trapz(y,gamma); % finite wing lift coefficient
CDi = 2 / ((V_inf)*s) * trapz(y,gamma.*AOAiR); % induced drag

%% Plots
plot(AOAd_geo,CL,'^') % plot individual data point
hold on
end
hold off
h = findobj(gca,'Type','line') % return graphics root object
AOA_NLLLT = cell2mat(get(h,'Xdata')) % retrive x data points
CL_NLLLT = cell2mat(get(h,'Ydata')) % retrieve y data points
plot(AOA_NLLLT,CL_NLLLT,'-b')
hold on
plot(AOA_2D,cl_2D,'--b') % plot input data (section lift coefficients)
data1 = load('0015EXP_Anderson.txt');
xEXP = data1(:,1);
yEXP = data1(:,2);
plot(xEXP,yEXP,'^r') % plot finite wing experimental data
data2 = load('0015LLT.txt')
AOA_LLT = data2(1:end,1)
CL_LLT = data2(1:end,2)
plot(AOA_LLT,CL_LLT,'g') % plot LLT approximation
ylabel('Lift Coefficient (C_L)')
xlabel('Angle of Attack (AOA) [deg]')
title('Lift Coefficient vs. Angle of Attack for NACA 0015 Wing')
data3 = load('Anderson_Numerical.txt')
xAnd = data3(1:end,1)
yAnd = data3(1:end,2)
plot(xAnd,yAnd,'k') % plot Andersons numerical data for comparison
legend('Current Implementation','Airfoil Section Lift Coefficient Data Re=350,000','Finite Wing Experimental Data AR=2.77, b=0.3515 [m]','LLT','Anderson Numerical')
axis([0 50 0 1.2])
[xEXP unique_indices] = unique(xEXP) % return unique indices for interpolation
yEXP = yEXP(unique_indices) % return unique values of cl
clear unique_indices
[xAnd unique_indices] = unique(xAnd) % return unique indices for interpolation
yAnd = yAnd(unique_indices) % return unique values of cl
clear unique_indices
alpha = linspace(0,50,200) % alpha for interpolation
for i = 1:200
CL_exp(i) = interp1(xEXP,yEXP,alpha(i)); % interpolate lift coefficient
CL_num(i) = interp1(AOA_NLLLT,CL_NLLLT,alpha(i)); % interpolate numerical data
CL_And(i) = interp1(xAnd,yAnd,alpha(i)) % interpolate Anderson
perror(i) = abs((CL_exp(i)-CL_num(i))/CL_exp(i))*100 % percent error present code
perrorAnd(i) = abs((CL_exp(i)-CL_And(i))/CL_exp(i))*100 % percent error Anderson
end
figure
plot(alpha,perror,alpha,perrorAnd)
title('Percent Error vs. Angle of Attack')
xlabel('Angle of Attack [Deg]')
ylabel('Percent Error From Experimental Results [%]')
legend('Present Method','Anderson')