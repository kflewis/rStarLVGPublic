function [A,B,C,D,F,G,S] = paramMat_distribute(theta,y,u)
% function [A,B,C,D,F,G,S = paramMatHLW_levelInf(theta,y,u)
%
% This function takes in the parameter vector theta, along with the observed data from
% the model and organizes the matrices in the state space model in order to be able to
% implement the state space form.  This also generates initial estimates and variances
% for the state based on the unconditional mean and variance of the state based on the
% information present.
%
% Model structure:
%
%% We will need to provide the mapping from the model specification into of the parameters
%% into the State space, which will take the form:
%% 
%% x_t = A x_t-1 + B u_t + C w_t
%% y_t = D x_t + F u_t + G w_t
%% 
%% where w_t is distributed N(0,S).
%
% In this model, we assume that all the errors are iid normals, meaning that the model is
% fully specified by these parameters and the diagonal S matrix with its inividual sigmas
% on each diagonal element (s1, s2, s3, s4, s5, s6).
%
% This is the specification for the theta parameter vector
% theta = [a1 a2 ar b1 bY s1 s2 s3 s4 s5]
%
% Kurt Lewis
% Created for original: April 20, 2016
% Last-updated Time-stamp: <2018-05-23 01:34:54 (m1kfl00)>
%

% Just to have the order in this file
%theta = [a1 a2 ar b1 bY s1 s2 s3 s4 s5 aX rhoG rhoZ];

% Pull out the individual elements to make the matrix equations more readable
a1 = theta(1);  
a2 = theta(2);
ar = theta(3);
b1 = theta(4);   % b_pi
bY = theta(5);    
s1 = theta(6);   % sig_ytilde
s2 = theta(7);   % sig_pi    
s3 = theta(8);   % sig_z     
s4 = theta(9);   % sig_ystar 
s5 = theta(10);  % sig_g     
aX = theta(11);
rhoG = theta(12);
rhoZ = theta(13);
muG = theta(14);
muZ = theta(15);

% Get some facts about the data
[N, T] = size(u);

%% To simplify the notation within the matrix , we will create some definitions here.
%% First, we note that the coefficient on the second through fourth lags must guarantee
%% that the sum of the two coefficients on the four lags must sum to one.   
b2 = (1 - b1);


%%% Form the matrices for the state space representation:
A = [1 0 0 rhoG 0 0 0;...
     1 0 0 0 0 0 0;...
     0 1 0 0 0 0 0;...
     0 0 0 rhoG 0 0 0;...
     0 0 0 1 0 0 0;...
     0 0 0 0 0 rhoZ 0;...
     0 0 0 0 0 1 0];


B = [muG*(1-rhoG) 0 0 0 0 0 0 0;...
     0 0 0 0 0 0 0 0;...
     0 0 0 0 0 0 0 0;...
     muG*(1-rhoG) 0 0 0 0 0 0 0;...
     0 0 0 0 0 0 0 0;...
     muZ*(1-rhoZ) 0 0 0 0 0 0 0;...
     0 0 0 0 0 0 0 0];     


C = [0 0 0 s4 s5;...
     0 0 0 0 0;...
     0 0 0 0 0;...
     0 0 0 0 s5;...
     0 0 0 0 0;...
     0 0 s3 0 0;...
     0 0 0 0 0];

%% Check
D = [1 -a1 -a2 -2*ar -2*ar -(ar/2) -(ar/2);...
     0 -bY 0 0 0 0 0];

F = [0 a1 a2 (ar/2) (ar/2) 0 0 aX;...
     0 bY 0 0 0 b1 b2 0];

G = [s1 0 0 0 0;...
     0 s2 0 0 0];

S = eye(5);


