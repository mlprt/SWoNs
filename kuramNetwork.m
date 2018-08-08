%% Run Kuramoto oscillator through network
function [theta,r, psi] = kuramNetwork (testNet,Edges,N,Lam,omega,theta,numNodes)

% This script  runs a kuramoto oscillator through a pre-generated network
% conditioned on local connections
% Numerical integration through fourth-order Runge-Kutta Method
% BC/ML/SWoNS/2018
% Code adapted from appmath.wordpress.com, courtesy
% Jeongho Kim, Mathematical Sciences, Seoul National University.
% Someone please cleanv

%% Generate network
iter = 500;
h = 0.1;
r = zeros(length(iter-1),1);
figure;
for k = 1:iter-1
    
    thetaConnect = theta(:,k)*ones(1,N)-(ones(N,1)*theta(:,k)'); %Generates phase matrix\
    indEdgeConnect = adjacency(testNet);
    thetaConnect(~indEdgeConnect) = 0; %no connection == 0
        
        f1 = kuramoto (thetaConnect,Lam,N,omega);
        f2=kuramoto(thetaConnect+0.5*h*f1(:,1),Lam,N,omega);
        f3=kuramoto(thetaConnect+0.5*h*f2(:,1),Lam,N,omega);           %4-th order Runge-Kutta method.
        f4=kuramoto(thetaConnect+h*f3(:,1),Lam,N,omega);        
        theta(:,k+1) = theta(:,k)+(h/6)*(f1(:,1))+2*f2(:,1)+2*f3(:,1)+f4(:,1);
        
        indOver = find (theta(:,k+1) > 2*pi);
        indUnder = find (theta(:,k+1) < 0);
        
        theta(indOver,k+1) = theta(indOver,k+1) - 2*pi;
        theta(indUnder,k+1) = theta(indUnder,k+1) + 2*pi;
        
        z = sum(exp(1i*theta(:,k+1)))/N;
        r(k+1) = abs (z);
        psi(k+1) = angle(z);
        
end

function f=kuramoto(x,Lam,N,omega)
 
f=omega+(Lam/N)*sum(sin(x))'; %Take out rand for noise
 
end

end