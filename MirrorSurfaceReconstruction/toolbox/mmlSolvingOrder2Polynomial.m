function [solution amsolind coefmat]= mmlSolvingOrder2Polynomial(x,y,w1x,w2x,w3x,w1y,w2y,w3y,w1,w2,w3,strInfo)
% This function is used to solve the polynomials and remove the ambiguity

% solving the polynomial and get the solutions with ambiguity
[sol,coefmat]= mmlSolvingOrder2Polynomial_basic(x,y,w1x,w2x,w3x,w1y,w2y,w3y,w1,w2,w3);
amsolind = 0;
if length(sol)==2,
    disp('there are two solutions')
    amsolind = 1;
end


% remove the ambiguity in the solution
% errthreshold = 20; % this threshold could be changed (20 for sphere)
errthreshold = 30; % this threshold could be changed (30 for ellipsoid)

% form the variable
seedpnt = [x y 1]'; % seed pixel
spe_pnt_3d = [w1 w2 w3]'; % seed point in 3d

solution = -10e6;
minerr_deltax = 10e6;
minerr_deltay = 10e6;
minindx = 1;
minindy = 1;

if isempty(sol)~=1,    
    for i = 1:length(sol),
        [err_deltax err_deltay] = orderOneTest(seedpnt,spe_pnt_3d,sol(i),strInfo);
        
        if isempty(err_deltax)~=1 && isempty(err_deltay)~=1,
            
            if min(abs(err_deltax)) < minerr_deltax,
                minerr_deltax = min(abs(err_deltax));
                minindx = i;
            end
            
 
            if min(abs(err_deltay)) < minerr_deltay,
                minerr_deltay = min(abs(err_deltay));
                minindy = i;
            end
        end
    end
    
    % if both x and y directions give minimum self consistency error by the
    % same depth, solution should be determined by the obtained depth from
    % the minimum error
    
    if minindx == minindy && minerr_deltax < errthreshold && minerr_deltay < errthreshold,
        solution = sol(minindx);
    end
end

end