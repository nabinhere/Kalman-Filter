%%
% Initialize simulation variables
SigmaW = 1e-5; % Process-noise covariance
SigmaV = 0.1; % Sensor-noise covariance
A = 1; B = -1e-4; % State-equation matrices
C = 0.7; D = -0.01; % Output-equation matrices
maxIter = 1000; % Number of simulation time steps
xtrue = 0.5; % Initialize true system initial state
xhat = 0.5; % Initialize Kalman filter initial estimate
SigmaX = 0; % Initialize Kalman filter covariance
u = 1; % Initial driving input, u[0]
% Reserve storage for variables we might want to plot/evaluate
xstore = zeros(maxIter+1,length(xtrue)); xstore(1,:) = xtrue;
xhatstore = zeros(maxIter,length(xhat));
SigmaXstore = zeros(maxIter,length(xhat)^2);

%%
% Main Kalman Filter steps
for k = 1:maxIter,
% KF Step 1a: State-prediction time update
xhat = A*xhat + B*u; % use prior value of "u"      
% KF Step 1b: Error-covariance time update
SigmaX = A*SigmaX*A' + SigmaW;
% [Implied operation of system in background, with
% input signal u, and output signal z]
switch k,
case 1, u = 0.5; % to match earlier example
case 2, u = 0.25; % to match earlier example
otherwise, u = randn(1); % just some interesting input
end
w = chol(SigmaW,'lower')*randn(length(xtrue)); % rand. process noise
v = chol(SigmaV,'lower')*randn(length(C*xtrue)); % rand. sensor noise
switch k,
case 1, ytrue = 3.85 - 3.5; % to match example   %output is 3.85 but 3.5 is subtracted so as to make the equation completely linear from affine
case 2, ytrue = 3.84 - 3.5; % to match example
otherwise, ytrue = C*xtrue + D*u + v; % based on present x and u
end
xtrue = A*xtrue + B*u + w; % future x is based on present u
% KF Step 1c: Estimate system output
yhat = C*xhat + D*u;
% KF Step 2a: Compute Kalman gain matrix
SigmaY = C*SigmaX*C' + SigmaV;
L = SigmaX*C'/SigmaY;
% KF Step 2b: State-estimate measurement update
xhat = xhat + L*(ytrue - yhat);
% KF Step 2c: Error-covariance measurement update
SigmaX = SigmaX - L*SigmaY*L';
% [Store information for evaluation/plotting purposes]
xstore(k+1,:) = xtrue;
xhatstore(k,:) = xhat;
SigmaXstore(k,:) = SigmaX(:);
end;

%%
figure(1); clf;
plot(0:maxIter-1,xstore(1:maxIter),'k-',0:maxIter-1,xhatstore,'b--', ...
0:maxIter-1,xhatstore+3*sqrt(SigmaXstore),'m-.',...
0:maxIter-1,xhatstore-3*sqrt(SigmaXstore),'m-.'); grid;
legend('True','Estimate','Bounds','location','northeast');
title('Kalman filter example');
xlabel('Iteration'); ylabel('State');
figure(2); clf;
plot(0:maxIter-1,xstore(1:maxIter)-xhatstore,'b-',0:maxIter-1, ...
3*sqrt(SigmaXstore),'m--',0:maxIter-1,-3*sqrt(SigmaXstore),'m--');
grid; legend('Error','Bounds','location','northeast');
title('Error with bounds');
xlabel('Iteration'); ylabel('Estimation error');

