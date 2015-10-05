function [spnts normals] = mmlLineIntersectSphereBest_MatCal(sc,r,visrays,cc)
%  this script is used to get the intersection of the visual rays with the
%  sphere.
%  INPUT:
%  sc: sphere center:3*1 or 1*3
%  r: radius
%  visrays: visual rays 3*n
%  cc: camera center: 3*1 or 1*3
%  spnts: sphere points: 
%  OUTPUT:
%  spnts: sphere points 3*n
%  normals: normal to the surface points 3*n
sx = sc(1);
sy = sc(2);
sz = sc(3);

cx = cc(1);
cy = cc(2);
cz = cc(3);
spnts = [];
[acr acc] = size(sc);

if acr < acc,
    sc = sc';
end

% % roots: return a solution vector
% for i = 1:size(visrays,2),  
%     i
% %%%%%% new idea%%%%%%%%%%
% vx = visrays(1,i);
% vy = visrays(2,i);
% vz = visrays(3,i);
% t = [ 
%   (sx*vx - cy*vy - cz*vz - cx*vx + sy*vy + sz*vz + (- cx^2*vy^2 - cx^2*vz^2 + 2*cx*cy*vx*vy + 2*cx*cz*vx*vz + 2*cx*sx*vy^2 + 2*cx*sx*vz^2 - 2*cx*sy*vx*vy - 2*cx*sz*vx*vz - cy^2*vx^2 - cy^2*vz^2 + 2*cy*cz*vy*vz - 2*cy*sx*vx*vy + 2*cy*sy*vx^2 + 2*cy*sy*vz^2 - 2*cy*sz*vy*vz - cz^2*vx^2 - cz^2*vy^2 - 2*cz*sx*vx*vz - 2*cz*sy*vy*vz + 2*cz*sz*vx^2 + 2*cz*sz*vy^2 + r^2*vx^2 + r^2*vy^2 + r^2*vz^2 - sx^2*vy^2 - sx^2*vz^2 + 2*sx*sy*vx*vy + 2*sx*sz*vx*vz - sy^2*vx^2 - sy^2*vz^2 + 2*sy*sz*vy*vz - sz^2*vx^2 - sz^2*vy^2)^(1/2))/(vx^2 + vy^2 + vz^2)
%  -(cx*vx + cy*vy + cz*vz - sx*vx - sy*vy - sz*vz + (- cx^2*vy^2 - cx^2*vz^2 + 2*cx*cy*vx*vy + 2*cx*cz*vx*vz + 2*cx*sx*vy^2 + 2*cx*sx*vz^2 - 2*cx*sy*vx*vy - 2*cx*sz*vx*vz - cy^2*vx^2 - cy^2*vz^2 + 2*cy*cz*vy*vz - 2*cy*sx*vx*vy + 2*cy*sy*vx^2 + 2*cy*sy*vz^2 - 2*cy*sz*vy*vz - cz^2*vx^2 - cz^2*vy^2 - 2*cz*sx*vx*vz - 2*cz*sy*vy*vz + 2*cz*sz*vx^2 + 2*cz*sz*vy^2 + r^2*vx^2 + r^2*vy^2 + r^2*vz^2 - sx^2*vy^2 - sx^2*vz^2 + 2*sx*sy*vx*vy + 2*sx*sz*vx*vz - sy^2*vx^2 - sy^2*vz^2 + 2*sy*sz*vy*vz - sz^2*vx^2 - sz^2*vy^2)^(1/2))/(vx^2 + vy^2 + vz^2)
% ];
% pnt =[cx cy cz]+ min(t(:))*[vx vy vz];
% spnts = [spnts;pnt];
% 
% end
vx = visrays(1,:);
vy = visrays(2,:);
vz = visrays(3,:);
t = [(sx*vx - cy*vy - cz*vz - cx*vx + sy*vy + sz*vz + (- cx^2*vy.^2 - cx^2*vz.^2 + 2*cx*cy*vx.*vy + 2*cx*cz*vx.*vz + 2*cx*sx*vy.^2 + 2*cx*sx*vz.^2 - 2*cx*sy*vx.*vy - 2*cx*sz*vx.*vz - cy^2*vx.^2 - cy^2*vz.^2 + 2*cy*cz*vy.*vz - 2*cy*sx*vx.*vy + 2*cy*sy*vx.^2 + 2*cy*sy*vz.^2 - 2*cy*sz*vy.*vz - cz^2*vx.^2 - cz^2*vy.^2 - 2*cz*sx*vx.*vz - 2*cz*sy*vy.*vz + 2*cz*sz*vx.^2 + 2*cz*sz*vy.^2 + r^2*vx.^2 + r^2*vy.^2 + r^2*vz.^2 - sx^2*vy.^2 - sx^2*vz.^2 + 2*sx*sy*vx.*vy + 2*sx*sz*vx.*vz - sy^2*vx.^2 - sy^2*vz.^2 + 2*sy*sz*vy.*vz - sz^2*vx.^2 - sz^2*vy.^2).^(1/2))./(vx.^2 + vy.^2 + vz.^2);
 -(cx*vx + cy*vy + cz*vz - sx*vx - sy*vy - sz*vz + (- cx^2*vy.^2 - cx^2*vz.^2 + 2*cx*cy*vx.*vy + 2*cx*cz*vx.*vz + 2*cx*sx*vy.^2 + 2*cx*sx*vz.^2 - 2*cx*sy*vx.*vy - 2*cx*sz*vx.*vz - cy^2*vx.^2 - cy^2*vz.^2 + 2*cy*cz*vy.*vz - 2*cy*sx*vx.*vy + 2*cy*sy*vx.^2 + 2*cy*sy*vz.^2 - 2*cy*sz*vy.*vz - cz^2*vx.^2 - cz^2*vy.^2 - 2*cz*sx*vx.*vz - 2*cz*sy*vy.*vz + 2*cz*sz*vx.^2 + 2*cz*sz*vy.^2 + r^2*vx.^2 + r^2*vy.^2 + r^2*vz.^2 - sx^2*vy.^2 - sx^2*vz.^2 + 2*sx*sy*vx.*vy + 2*sx*sz*vx.*vz - sy^2*vx.^2 - sy^2*vz.^2 + 2*sy*sz*vy.*vz - sz^2*vx.^2 - sz^2*vy.^2).^(1/2))./(vx.^2 + vy.^2 + vz.^2)
];

tmin = min(t);

spnts = repmat([cx;cy;cz],1,size(visrays,2)) + repmat(tmin,3,1).*visrays; % obtained sphere points

% normals = (spnts - repmat(sc',size(spnts,1),1))./repmat(sqrt(sum((spnts - repmat(sc',size(spnts,1),1)).^2,2)),1,3);
normaltmp = spnts - repmat([sx;sy;sz],1,size(spnts,2));
% normtmp = sqrt(normaltmp(1,:).^2 + normaltmp(2,:).^2 + normaltmp(3,:).^2);
normals = mmlNormalizeMatrix(normaltmp);
end



