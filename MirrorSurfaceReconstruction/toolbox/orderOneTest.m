function [err_deltax err_deltay] = orderOneTest(seedpnt,spe_pnt_3d,sol,strInfo)
% this function is used for testing whether the solution satisfies the first order partial differential equation

% pass all the global parameters in use
kk = strInfo.kk; bbs_2d = strInfo.bbs_2d; ctrlpnts_2d = strInfo.ctrlpnts_2d;tran_vec = strInfo.tran_vec; r = strInfo.r;

% examine the first order equation
tempdeltavec = kk\[1 0 0]';
scale = 0.1;
deltax = scale*tempdeltavec;

% changes in x direction
seedpnt_deltax = seedpnt + deltax;

% calculate the specular correspondences
spe_pnt_deltax = bbs_eval(bbs_2d,ctrlpnts_2d',seedpnt_deltax(1),seedpnt_deltax(2));
spe_pnt_deltax_3d = r*[spe_pnt_deltax;0] + tran_vec;

% calculate the partial derivative of u
spe_pnt_deltax_du = bbs_eval(bbs_2d,ctrlpnts_2d',seedpnt_deltax(1),seedpnt_deltax(2),1,0);

% calculate the partial derivative of v
spe_pnt_deltax_dv = bbs_eval(bbs_2d,ctrlpnts_2d',seedpnt_deltax(1),seedpnt_deltax(2),0,1);

% calculate the partial derivative of x and y
spe_pnt_deltax_dx = r*[spe_pnt_deltax_du;0];
spe_pnt_deltax_dy = r*[spe_pnt_deltax_dv;0];

% using polynomial to calculate the depth for the seed points (good)
sol_final_deltax  = mmlSolvingOrder2Polynomial_basic(seedpnt_deltax(1),seedpnt_deltax(2),spe_pnt_deltax_dx(1),spe_pnt_deltax_dx(2),spe_pnt_deltax_dx(3),spe_pnt_deltax_dy(1),spe_pnt_deltax_dy(2),spe_pnt_deltax_dy(3),spe_pnt_deltax_3d(1),spe_pnt_deltax_3d(2),spe_pnt_deltax_3d(3));

% examine \partial s/\partial x = -s(n_x/<n,v>)
% normal= (spe_pnt_3d - (sol + 1)* seedpnt)/norm(spe_pnt_3d - (sol + 1)* seedpnt);
normal_tmp= (spe_pnt_3d - sol * seedpnt)/norm(spe_pnt_3d - sol * seedpnt) - seedpnt/norm(seedpnt);
normal = normal_tmp/norm(normal_tmp);

err_deltax = [];

if max(size(sol_final_deltax))==2,
    
    tmpval = (sol_final_deltax(1) - sol)/deltax(1)+ sol * normal(1)/dot(normal,seedpnt);
    err_deltax = [err_deltax tmpval];
    
    tmpval = (sol_final_deltax(2) - sol)/deltax(1)+ sol * normal(1)/dot(normal,seedpnt);
    err_deltax = [err_deltax tmpval];    
end

if max(size(sol_final_deltax))==1,
     
    tmpval = (sol_final_deltax(1) - sol)/deltax(1)+ sol * normal(1)/dot(normal,seedpnt);
    err_deltax = [err_deltax tmpval];
end



%% changes in y direction

tempdeltayvec = kk\[0 1 0]';
scale = 0.1;
deltay = scale*tempdeltayvec;

seedpnt_deltay = seedpnt + deltay;

% calculate the specular correspondences
spe_pnt_deltay = bbs_eval(bbs_2d,ctrlpnts_2d',seedpnt_deltay(1),seedpnt_deltay(2));
spe_pnt_deltay_3d = r*[spe_pnt_deltay;0] + tran_vec;

% calculate the partial derivative of u
spe_pnt_deltay_du = bbs_eval(bbs_2d,ctrlpnts_2d',seedpnt_deltay(1),seedpnt_deltay(2),1,0);

% calculate the partial derivative of v
spe_pnt_deltay_dv = bbs_eval(bbs_2d,ctrlpnts_2d',seedpnt_deltay(1),seedpnt_deltay(2),0,1);

% calculate the partial derivative of x and y
spe_pnt_deltay_dx = r*[spe_pnt_deltay_du;0];
spe_pnt_deltay_dy = r*[spe_pnt_deltay_dv;0];

% using polynomial to calculate the depth for the seed points (good)
sol_final_deltay = mmlSolvingOrder2Polynomial_basic(seedpnt_deltay(1),seedpnt_deltay(2),spe_pnt_deltay_dx(1),spe_pnt_deltay_dx(2),spe_pnt_deltay_dx(3),spe_pnt_deltay_dy(1),spe_pnt_deltay_dy(2),spe_pnt_deltay_dy(3),spe_pnt_deltay_3d(1),spe_pnt_deltay_3d(2),spe_pnt_deltay_3d(3));

% record the first order error in y

err_deltay = [];
% test the obtained neighboring points
if max(size(sol_final_deltay))==2,
    
    tmpval = (sol_final_deltay(1) - sol)/deltay(2) + sol * normal(2)/dot(normal,seedpnt);    
    err_deltay = [err_deltay tmpval];
    
    tmpval = (sol_final_deltay(2) - sol)/deltay(2) + sol * normal(2)/dot(normal,seedpnt);
    err_deltay = [err_deltay tmpval];
end

if max(size(sol_final_deltay))==1,
    
    tmpval = (sol_final_deltay(1) - sol)/deltay(2) + sol * normal(2)/dot(normal,seedpnt);
    err_deltay = [err_deltay tmpval];
end

end



















