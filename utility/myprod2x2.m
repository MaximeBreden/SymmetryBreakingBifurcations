function c = myprod2x2(a, b, N)

if nargin == 2
    N = size(a{1,1});
end

c = {convo(a{1,1},b{1,1},N)+convo(a{1,2},b{2,1},N), convo(a{1,1},b{1,2},N)+convo(a{1,2},b{2,2},N);
     convo(a{2,1},b{1,1},N)+convo(a{2,2},b{2,1},N), convo(a{2,1},b{1,2},N)+convo(a{2,2},b{2,2},N)};

