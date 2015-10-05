% This script is used for reconstructing the mirror sphere by solving
% degree 2 polynomials.
% 
% copy right by Miaomiao Liu@nicta.com.au
%
%
clear all
close all
clc

%% add data, toolbox paths
addpath('data/')
addpath('BBS/')
addpath('toolbox/')

%% load data 
load('pixels.mat'); % pixels visible in the image
load('spnts.mat');  % 3D points of the visible points on the sphere. 
load('visrays.mat'); % visual rays
load('rfppnts.mat');% reflection correspondences in the camera coordinate system.
load('cameraInfo.mat') % kk: intrinsics| r: rotation matrix| t: translation vector | c: camera center


figure(10)
plot(pixels(1,:)',pixels(2,:)');
title('pixels');


figure(20)
plot3(visrays(1,:),visrays(2,:),visrays(3,:),'b.');


figure(30)
scatter3(rfppnts(1,:),rfppnts(2,:),rfppnts(3,:),'g.');
axis equal;
%% open matlabpool
isopen = matlabpool('size')>0;
if isopen==0,
   matlabpool(4); 
end


%% processing data

% define pixels_uf as points on the unit focal plane
pixels_uf = [visrays(1,:)./visrays(3,:);visrays(2,:)./visrays(3,:);ones(1,size(visrays,2))]; % pixels on the unit focal plane

figure(40)
plot3(pixels_uf(1,:),pixels_uf(2,:),pixels_uf(3,:),'k.');
% calculate reflection correspondences in the reference plane coordinate
% system. tranvec_pl: translation vector of the reference plane w.r.t the
% world coordinate system.
tranvec_pl = [0 0 2]';
tmp_spe = r'*rfppnts + repmat(c,1,size(rfppnts,2)) + repmat(tranvec_pl,1,size(rfppnts,2));

figure(50)
scatter3(tmp_spe(1,:),tmp_spe(2,:),tmp_spe(3,:),'g.');
axis equal;
% translation vector of the reference plane in the camera coordinate
% system.
tran_vec = -r*(tranvec_pl+c);
spe_corres = tmp_spe(1:2,:);

figure(60)
plot(spe_corres(1,:),spe_corres(2,:),'g.');

% Using 2D uniform cubic bspline to interpolate the 
ctrlpntsnum = 40;          % control points number
nptx_2d = ctrlpntsnum;
npty_2d = ctrlpntsnum;

margin = 0;
xmin_2d = min(pixels_uf(1,:) - margin);
xmax_2d = max(pixels_uf(1,:) + margin);
ymin_2d = min(pixels_uf(2,:) - margin);
ymax_2d = max(pixels_uf(2,:) + margin);

bbs_2d = bbs_create(xmin_2d, xmax_2d, nptx_2d, ymin_2d, ymax_2d, npty_2d, 2);

% using bbs to obtain the spline correspondences
colmat = bbs_coloc(bbs_2d, (pixels_uf(1,:))', (pixels_uf(2,:))');

% control points for spline correspondences
ctrlpnts_2d = colmat\(spe_corres(1:2,:))';

% data structure for all points
rec_3dpnts = zeros(3,size(pixels_uf,2));

% sphere information
scenter = [0 0 10]';   % sphere center in the world coordinate system
radius = 3;
scenter_cam = r*(scenter - c); % sphere center coordinates in the camera coordinate system
solutionmat = zeros(1,size(pixels,2));
twodepthind = [];

%% define strInfo for orderOneTesting
strInfo.bbs_2d = bbs_2d; strInfo.ctrlpnts_2d = ctrlpnts_2d;strInfo.kk = kk;strInfo.tran_vec = tran_vec;strInfo.r = r;

%% obtain all the points by solving polynomials\
tic
disp('start computing the depth from solving degree 2 polynomials...it will take around 41 seconds using 4 cores for 477976 points..')
parfor i = 1:size(pixels_uf,2),
    
    % set seedpnt
    seedpnt = pixels_uf(:,i);
    
    % calculate the specular correspondences in 3D (camera coordinate system)
    spe_pnt = bbs_eval(bbs_2d,ctrlpnts_2d',seedpnt(1),seedpnt(2));
    spe_pnt_3d = r*[spe_pnt;0] + tran_vec;
    
    % calculate the partial derivative of u
    spe_pnt_du = bbs_eval(bbs_2d,ctrlpnts_2d',seedpnt(1),seedpnt(2),1,0);
    
    % calculate the partial derivative of v
    spe_pnt_dv = bbs_eval(bbs_2d,ctrlpnts_2d',seedpnt(1),seedpnt(2),0,1);
    
    % calculate the partial derivative of x
    spe_pnt_dx = r*[spe_pnt_du;0];
    spe_pnt_dy = r*[spe_pnt_dv;0];
    
    % visual ray
    visray_pnt = seedpnt/norm(seedpnt);
    
    % calculate the ground truth data: visual rays intersect with sphere
    sol_gt_3dpnt = mmlLineIntersectSphereBest_MatCal(scenter_cam,radius,visray_pnt,[0 0 0]');    
    ref_depth = sol_gt_3dpnt(3); % ground truth depth
    
    % using polynomial to calculate the depth for the seed points (good)
    [sol_final,ambiindval]= mmlSolvingOrder2Polynomial(seedpnt(1),seedpnt(2),spe_pnt_dx(1),spe_pnt_dx(2),spe_pnt_dx(3),spe_pnt_dy(1),spe_pnt_dy(2),spe_pnt_dy(3),spe_pnt_3d(1),spe_pnt_3d(2),spe_pnt_3d(3),strInfo);
    
    % depth from solving degree 2 polynomials
    solutionmat(i) = sol_final;
end
toc
%%
ind = find(abs(solutionmat(:))<10e6); % valid solution

estshape = repmat(solutionmat(ind),3,1).*pixels_uf(:,ind);

% plot all the reconstructed points: (estshape: estimated 3D shape, spnts: ground truth points)
figure(1)
plot3(estshape(1,:),estshape(2,:),estshape(3,:),'r.')
hold on
plot3(spnts(1,:),spnts(2,:),spnts(3,:),'b.')
title('red: estimated 3D points; blue: ground truth points')
