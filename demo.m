clear;clc;close all;
addpath(genpath(cd));

%% load data, generate sensing matrix and measurement
load kaist_crop256_01
x     = img;
[A,b] = gen_31(x,mask);

%% reconstruction using NLTT-TV
K       = 160;  % the # of clusters 
alpha   = 1e-2; % the weight of the LTTR term
beta    = 1e-2; % the weight of the 3DTV term
mu      = 1e-3; % the penalty parameter of ALF
maxiter = 200;  % the maximum of iterations

tic
Z  = mainsolver(A,b,K,alpha,beta,mu,maxiter,x);
toc

%% quality assessment
[psnr,ssim,sam] = quality_assessment(Z,x)
% imshow(x(:,:,10),[],'border','tight','initialmagnification','fit');  set (gcf,'Position',[0,0,512,512]);
imshow(Z(:,:,10),[],'border','tight','initialmagnification','fit');  set (gcf,'Position',[0,0,512,512]);
