
function [TrueLrr, TrueK] = calculateLKfunction(points, r)

area    = (max(points(:,1))-min(points(:,1)))*(max(points(:,2))-min(points(:,2)));

DIST = calcDist(points);
L = calcLfunction(DIST,r,area);
TrueLrr=L-r';
TrueK = calcKfunction(DIST, r, area);
% figure; hold on;
% plot(r, TrueLrr,'magenta');



function [DIST] = calcDist(Points)
A = bsxfun(@minus,transpose(Points),mean(transpose(Points),2));
S = full(dot(A,A,1));
DIST = sqrt(bsxfun(@plus,S,S')-full(2*(A'*A)));
end


function L=calcLfunction(DIST,r,area)

L=sqrt(calcKfunction(DIST,r,area)/pi);

end


function K = calcKfunction(DIST, r, area)

    K = zeros(length(r),1);
    
for k=1:length(K)
    K(k) = sum(sum(DIST<r(k)));
end

K=K-size(DIST,1);
K=K/(size(DIST,1)*size(DIST,1)/(area*1.0));

end

end