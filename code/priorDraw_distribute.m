function theta = priorDraw_distribute(num)
% function theta = priorDrawHLW(num)
%
% This function generates a draw from the priors for the Test 1 model, and the input
% parameter is restricted to do this one draw at a time anyway.
%
% Kurt Lewis
% Created: April 27. 2016
% Last-updated Time-stamp: <2018-05-23 00:33:33 (m1kfl00)>
%
    
% If they fail to give me the number 1, fill it in
if nargin < 1
    num = 1;
end
    
% Error if they try to draw more than one
if num > 1
    error('Cannot draw more than one at a time')
end

% Most of the parameters priors will be mean zero and variance 2 for now, to be drawing
% from the core of the prior for this one and only time.
m = 0;
v = 2;

% First, we draw the coefficient parameters

% Parameters that don't have to be linked to stationarity
aX = normrnd(m,v);
%rhoG = normrnd(m,v);
%rhoZ = normrnd(m,v);
muG = normrnd(m,v);
muZ = normrnd(m,v);
a1 = normrnd(m,v);
a2 = normrnd(m,v);

% b1 has to be between 0 and 1
b1 = unifrnd(0,1);

esc = 0;
while esc == 0
    rhoG = normrnd(1,0.5);
    if rhoG > 0
        esc = 1;
    else 
        esc = 0;
    end
end

esc = 0;
while esc == 0
    rhoZ = normrnd(1,0.5);
    if rhoZ > 0
        esc = 1;
    else 
        esc = 0;
    end
end


esc = 0;
while esc == 0
    ar = normrnd(m,v);
    if ar < -0.0025
        esc = 1;
    else 
        esc = 0;
    end
end

bY = normrnd(m,v);

esc1 = 0;
while esc1 == 0
    bY = normrnd(m,v);
    if bY > 0.025
        esc1 = 1;
    else 
        esc1 = 0;
    end
end

% Second we draw the shock variances from an inverse gamma.  We can draw from a gamma
% and use the relationship between gamma and inverse gamma to draw from the appropriate
% inverse gamma distribution

% To give the least informative prior, I will set both the shape and scale of the inverse
% gamma to 1.
%sh = 3;  
%sc = 0.5;

% Uniform priors on the sigmas, again favor more toward the lower part where we think
% it is as a starting point.
sigLB = sqrt(2*eps);
sigHB = 2;


% We draw from the gamma with shape set equal to the inverse gamma shape and the scale
% set to the inverse of the scale parameter in the inverse distribution.

% We need to make these obey the lambda rules as well
% Boundaries on lambdas
lambHB = 5;
lambLB = 0;

escL = 0;
while escL == 0
    s1 = unifrnd(sigLB,sigHB);
    s2 = unifrnd(sigLB,sigHB);
    s3 = unifrnd(sigLB,sigHB);
    s4 = unifrnd(sigLB,sigHB);
    s5 = unifrnd(sigLB,sigHB);
    
    %% lambda constraints 
    lambG = s5/s4;
    lambZ = (-ar*s3)/s1;
    
    % Check to see if either lambda violates either bound
    if ((lambG < lambLB) || (lambG > lambHB) || ...
        (lambZ < lambLB) || (lambZ > lambHB))
        
        % If it does violate the constraint, reject the whole draw
        escL = 0;
    else
        % If it doesn't violate the constraint, go on as usual dependent on
        % the other information from theta.
        escL = 1;
    end
    
end

    

% Assemble the output from the prior draw
theta = [a1 a2 ar b1 bY s1 s2 s3 s4 s5 aX rhoG rhoZ muG muZ];


    

    
