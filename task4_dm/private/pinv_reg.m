function pDM = pinv_reg(DM,lambda,regtype,breakPoints)

% pDM = pinv_reg(DM,lambda,regtype)
%
% returns pseudoinverse of design matrix DM with ridge regression
% using identity matrix multiplied by parameter lambda
% http://www.math.kent.edu/~reichel/publications/square.pdf

if nargin<2 % no regularisation
    lambda = 0;
end
if nargin<3
    regtype = 'ident';
end
if nargin < 4
        breakPoints = [];
end

switch regtype
    case 'ident' % aka "Zero-Order Tikhonov", aka "L0", aka "ridge"
        LtL = eye(size(DM,2));
    case 'onediff' % aka "First-Order Tikhonov"
        L = 0.5 * ([eye(size(DM,2)-1) zeros(size(DM,2)-1,1)] + ...
            [zeros(size(DM,2)-1,1) -1 * eye(size(DM,2)-1)]);
        
        %         % Hard-coding this for now. DO NOT LEAVE AS-IS
        %         L(161,162) = 0;
        %
        %         if size(L,2) > 322
        %             L(322,323) = 0;
        %         end
        for i = 1:length(breakPoints)
            L(breakPoints(i),breakPoints(i)+1) = 0;
        end
        
        LtL = L'*L;
    case 'twodiff' % aka "Second-Order Tikhonov"
        L = 0.25 * ([-1 * eye(size(DM,2)-2) zeros(size(DM,2)-2,2)] + ...
            [zeros(size(DM,2)-2,1) 2 * eye(size(DM,2)-2) zeros(size(DM,2)-2,1)] + ...
            [zeros(size(DM,2)-2,2) -1 * eye(size(DM,2)-2)]);
        for i = 1:length(breakPoints)
            L(breakPoints(i),breakPoints(i)+1) = 0;
            L(breakPoints(i),breakPoints(i)+2) = 0;
        end
        LtL = L'*L;
    case 'threediff'
        L = 0.125 * ([-1 * eye(size(DM,2)-3) zeros(size(DM,2)-3,3)] + ...
            [zeros(size(DM,2)-3,1) 3 * eye(size(DM,2)-3) zeros(size(DM,2)-3,2)] + ...
            [zeros(size(DM,2)-3,2) -3 * eye(size(DM,2)-3) zeros(size(DM,2)-3,1)] + ...
            [zeros(size(DM,2)-3,3) 1 * eye(size(DM,2)-3)]);
        LtL = L'*L;
        for i = 1:length(breakPoints)
            L(breakPoints(i),breakPoints(i)+1) = 0;
            L(breakPoints(i),breakPoints(i)+2) = 0;
            L(breakPoints(i),breakPoints(i)+3) = 0;
        end
    otherwise
        error('unrecognised entry for L matrix');
end

% pDM = inv(DM'*DM+lambda*diffmat2'*diffmat2)*DM';
pDM = inv(DM'*DM+lambda*LtL)*DM';
% pDM = pinv(DM'*DM+lambda*LtL)*DM';


