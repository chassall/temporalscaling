%% gaussian process function, for simulating data
function test_Y = try_gaussian_processes(train_X, train_Y, test_x, length_scale,doplots)

%some code to simulate data using gaussian processes, based on https://planspace.org/20181226-gaussian_processes_are_not_so_fancy/

if nargin<4
    length_scale = 1;
end

if nargin<5
    doplots = 0;
end

%% calculate euclidean distances between all elements of X

dist_XX = pdist2(train_X,train_X);
kern_XX = squared_exponential(dist_XX, length_scale);

dist_xX = pdist2(test_x,train_X);
kern_xX = squared_exponential(dist_xX, length_scale);

dist_xx = pdist2(test_x,test_x);
kern_xx = squared_exponential(dist_xx, length_scale);

dist_Xx = pdist2(train_X,test_x);
kern_Xx = squared_exponential(dist_Xx, length_scale);

test_Y = kern_xX*inv(kern_XX)*train_Y;
%can use pinv instead of inv if you have multiple observations with different Ys!

cov_x = kern_xx - kern_xX*inv(kern_XX)*kern_Xx;

var_x = diag(cov_x);

if doplots
    figure; plot(train_X,train_Y,'.','MarkerSize',20);
    hold on; plot(test_x,test_Y,'r');
    plot(test_x,test_Y+2*var_x,'k--');
    plot(test_x,test_Y-2*var_x,'k--');
end

end