function [u, erriter, iter, timet] = Cont_Max_Flow_Cut(im,UB_canopy,UB_background)
%   Performing the continuous max-flow algorithm to solve the 
%   continuous min-cut problem in 2D
%    
%   Usage: [u, erriter, i, timet] = CMF_Cut;
%   input:
%   -im: the image 
%   - alpha: the penalty parameter to the total-variation term.
%       For the case without incorporating image-edge weights, alpha is given
%       by the constant everywhere. For the case with image-edge weights,
%       alpha is given by the pixelwise weight function:
%
%       For example, alpha(x) = b/(1 + a*| nabla f(x)|) where b and a are positive
%       constants and |nabla f(x)| gives the strength of the local gradient.
%
%   - cc: gives the step-size of the augmented Lagrangian method.
%       The optimal range of cc is [0.3, 3].
%
%   - errbound: the error bound for convergence.
%
%   - numIter: the maximum iteration number.
%
%   - steps: the step-size for the graident-projection step to the
%       total-variation function. The optimal range of steps is [0.1,
%       0.17].
%
%   Outputs: 
%       - u: the final results u(x) in [0,1]. As the following paper,
%           the global binary result can be available by threshholding u
%           by any constant Beta in (0,1):
% 
%       - erriter: it returns the error evaluation of each iteration,
%           i.e. it shows the convergence rate. One can check the algorithm
%           performance.
%
%       - i: gives the total number of iterations, when the algorithm converges.
%
%       - timet: gives the total computation time.

%% defining parameters to Cont_Max_Flow_Cut
const_alph = 0.5; 
cc = 3;%[0.3, 3]
errbound = 1e-4;% the error bound for convergence
numIter = 300;
steps = 0.16;%step-size for the graident-projection step to the total-variation function. 
% The optimal range of steps is [0.1,0.17]
[rows,cols]=size(im);
alpha = const_alph*ones(rows,cols); 
imgSize = (rows*cols);
% build up the data terms
ulab(1) = UB_canopy;
ulab(2) = UB_background;
Cs = abs(im - ulab(1));%the source capacity
Ct = abs(im - ulab(2));%the sink capacity

% the initial value of u is set to be an initial cut, see below.
u = double((Cs-Ct) >= 0);

% the initial values of two terminal flows ps and pt are set to be the specified legal flows.
ps = min(Cs, Ct);%the initial flow from the sourc terminal to each pixel
pt = ps; %the initial flow from each pixel to the t terminal

% the initial value of the spatial flow  p = (pp1, pp2) is set to be zero.
pp1 = zeros(rows, cols+1);
pp2 = zeros(rows+1, cols);
R2L_flow = pp1(:,2:cols+1)-pp1(:,1:cols);%The pixel flow to the pixel to the left, right to left
D2U_flow = pp2(2:rows+1,:)-pp2(1:rows,:);% The flow from pixel to pixel above it down to up 
divp =R2L_flow + D2U_flow; 

erriter = zeros(numIter,1);%Initialization of the error vector to be zero

tic %measuring the time of the process

for iter = 1:numIter

    %update the total-variation function which includes in each pixel the
    %total flow from the spatial flows minus the the enters from the source
    %Plus the flow delivered to terminal t minus 
    %cc the step-size of the augmented Lagrangian method  
    pts = divp - (ps - pt  + u/cc);
    %update the spatial flow by gradient descent of delta between pts and pp and step as the step-size

    pp1(:,2:cols) = pp1(:,2:cols) + steps*(pts(:,2:cols) - pts(:,1:cols-1));%from right to left pts
    pp2(2:rows,:) = pp2(2:rows,:) + steps*(pts(2:rows,:) - pts(1:rows-1,:));%from down to up pts
    
    % creat a projection in order to stand with the constraint about the flow capacity |p(x)| <= alpha(x)
    gk = pp1(:,1:cols).^2 + pp1(:,2:cols+1).^2 + pp2(1:rows,:).^2 + pp2(2:rows+1,:).^2;%
    gk = sqrt(gk*0.5);

    gk = double(gk <= alpha) + double(~(gk <= alpha)).*(gk ./ alpha);
    gk = 1 ./ gk;

    pp1(:,2:cols) = (0.5*(gk(:,2:cols) + gk(:,1:cols-1))).*pp1(:,2:cols);%update the flows accordings to gk 
    pp2(2:rows,:) = (0.5*(gk(2:rows,:) + gk(1:rows-1,:))).*pp2(2:rows,:);


    divp = pp1(:,2:cols+1)-pp1(:,1:cols)+pp2(2:rows+1,:)-pp2(1:rows,:);%update toatal  spatial flows R2L and D2U
    
    % update the terminal source flow ps
       
    pts = divp + pt - u/cc + 1/cc;
    ps = min(pts, Cs);%according to the source capcity
   
    
    % update the sink flow pt
    pts = - divp + ps + u/cc;
    pt = min(pts, Ct);%according to the sink capcity

	% update the multiplier u
	lagra = cc*(divp - ps + pt);
	u = u - lagra;

    % evaluate the avarage error
    erriter(iter) = sum(sum(abs(lagra)))/imgSize; 
    if (erriter(iter) < errbound)
        break;
    end
end
toc
timet = toc

% figure();
% imshow(u,[]);

fprintf('number of iterations = %u. \n', iter);

end

