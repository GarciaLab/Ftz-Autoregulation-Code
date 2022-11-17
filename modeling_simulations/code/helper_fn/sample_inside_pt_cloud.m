% adapted from https://www.mathworks.com/matlabcentral/answers/327990-generate-random-coordinates-inside-a-convex-polytope
function uvw = sample_inside_pt_cloud(x,npts)

ch = convhulln(x);

% add center of mass to list of points
x_com = mean(x);
x = [x;x_com];

% triangulate
ntri = size(ch,1);
tri = [ch,repmat(size(x,1),ntri,1)];

% figure
% % plot3(x(:,1),x(:,2),x(:,3),'bo');
% % hold on
% plot3([x(tri(:,1),1),x(tri(:,2),1),x(tri(:,3),1),x(tri(:,4),1)]', ...
%     [x(tri(:,1),2),x(tri(:,2),2),x(tri(:,3),2),x(tri(:,4),2)]', ...
%     [x(tri(:,1),3),x(tri(:,2),3),x(tri(:,3),3),x(tri(:,4),3)]','g-')

% calculate volume of each simplex
V = zeros(1,ntri);
for ii = 1:ntri
    V(ii) = abs(det(x(tri(ii,1:3),:) - x_com));
end
V = V/sum(V);

% allocate points by simplex (uniform by volume)
[~,~,simpind] = histcounts(rand(npts,1),cumsum([0,V]));

% generate points within each simplex
r1 = rand(npts,1);
uvw = x(tri(simpind,1),:).*r1 + x(tri(simpind,2),:).*(1-r1);
r2 = sqrt(rand(npts,1));
uvw = uvw.*r2 + x(tri(simpind,3),:).*(1-r2);
r3 = nthroot(rand(npts,1),3);
uvw = uvw.*r3 + x(tri(simpind,4),:).*(1-r3);

% hold on
% scatter3(uvw(:,1),uvw(:,2),uvw(:,3),'bo')
% hold off
end