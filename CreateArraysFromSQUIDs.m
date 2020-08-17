%clear
%close all

%% Set parameters
dist_to_head = 3e-3; % OPM distance to head

% these are parameters for adjusting the SQUID positions for better
% coverage
phi = 0;
scale_param = 1.05;
horizontal_transla = 0;

%% Get the geometrical information, you should modify this to your needs

%load('./geometry/bmeshes-m3.mat')
%load('./geometry/coils-csquid0.mat')

% coilp is a 102x3 matrix containing the SQUID positions
%coilp = coils.pos(3:3:306,:);
coilp = coilp;
% coiln is a 102x3 matrix containing the SQUID orientations
%coiln = coils.n(3:3:306,:);
coiln = coiln;

coilp_orig = coilp;

%clear coils;

% Next get the vertices of the skull BEM mesh to a matrix skullp
%skullp = bmeshes.meshes{2}.p;
skullp = skullp;
%vertices of the scalp BEM mesh to a matrix scallp
%scalpp = bmeshes.meshes{3}.p;
scalpp = scalpp;
%scalpnn contains the normal vectors of the BEM mesh at the vertex positions
%load('scalpnn.mat');
scalpnn=scalpnn;
%clear bmeshes;
%% This next code segment adjusts the SQUID helmet
% so that it has better coverage of the brain

cm = mean(skullp );
v = eye(3);
cms = mean(coilp);
[us,ss,vs] = svd(coilp-bsxfun(@times,cms,ones(size(coilp))),0);

cm_dist = norm(cms-cm);


I = eye(3);

v_val = v;
for i = 1:3
    dots = sum(repmat(I(:,i),1,size(v_val,2)).*v_val,1);
    idx = find(abs(dots)  == max(abs(dots)));
    v(:,i) = sign(dots(idx))*v_val(:,idx);
    v_val(:,idx) = [];
end

vs_val = vs;
for i = 1:3
    dots = sum(repmat(v(:,i),1,size(vs_val,2)).*vs_val,1);
    idx = find(abs(dots)  == max(abs(dots)));
    vs(:,i) = sign(dots(idx))*vs_val(:,idx);
    vs_val(:,idx) = [];
end

% Plot the initial position of the SQUID helmet
figure(1)
temp = bsxfun(@times,cms,ones(3,3));
quiver3(temp(:,1),temp(:,2),temp(:,3),vs(1,:)',vs(2,:)',vs(3,:)',0.1);
hold on
plot3(coilp(:,1),coilp(:,2),coilp(:,3),'xk')
quiver3(coilp(:,1),coilp(:,2),coilp(:,3),coiln(:,1),coiln(:,2),coiln(:,3),0.5);

temp = bsxfun(@times,cm,ones(3,3));
tempvec = eye(3);
quiver3(temp(:,1),temp(:,2),temp(:,3),tempvec(:,1),tempvec(:,2),tempvec(:,3),0.1);
plot3(scalpp(:,1),scalpp(:,2),scalpp(:,3),'.')
axis equal
axis tight
view(-90,0)
set(figure(1),'position',[0 200 400 400]);

R = cell(3,1);
for i = 1:2
    vec = cross(vs(:,i),v(:,i));
    vec = vec/norm(vec);
    theta = acos(dot(vs(:,i),v(:,i))/(norm(vs(:,i))*norm(v(:,i))));
    R{i} = ConstructRotationMatrix(theta,vec);
    vs = R{i}*vs;
    coilp = (R{i}*coilp')';
end
R{3} = ConstructRotationMatrix(phi*pi/180,v(:,1));
coilp = (R{3}*coilp')';
fullR = R{3}*R{2}*R{1};
vs = R{3}*vs;


cms_n = mean(coilp);
dist = norm(cm-cms);
vec = cm - cms_n;

T = [fullR vec'];
T(end+1,:) = [0 0 0 1];

coilp = T*[coilp_orig';ones(1,size(coilp_orig,1))];
coilp = coilp(1:3,:)';

transla = scale_param*cm_dist*vs(:,3);
hor_transla = horizontal_transla*vs(:,2);
coilp = coilp + repmat(transla',size(coilp,1),1) + repmat(hor_transla',size(coilp,1),1);
coils2proj = coilp; 
cms_n = mean(coilp);


% Plot the adjusted position of the SQUID helmet
figure(2)
hold on
plot3(scalpp(:,1),scalpp(:,2),scalpp(:,3),'.')
plot3(coilp(:,1),coilp(:,2),coilp(:,3),'xk')

view(-90,0)
axis equal
axis tight
view(-90,0)
set(figure(2),'position',[0 200 400 400]);


%% The following projects the SQUID positions to a distance dist_to_head from the head
% it also constructs the OPM positions to matrix OPMp and orientations
% OPMnn: normal component
% OPMtt1: #1 tangential component
% OPMtt2: #2 tangential component

rvecs = repmat(cm,size(scalpp,1),1)-scalpp;
rvecs = bsxfun(@times,rvecs,1./sqrt(sum(rvecs.^2,2)));

for i = 1:size(coils2proj,1)
    svec = cm - coils2proj(i,:);
    svec = svec/norm(svec);
    thetas = acos(sum(repmat(svec,size(scalpp,1),1).*rvecs,2));
    idx = find(thetas == min(thetas));
    OPMnn(i,:) = scalpnn(idx,:);
    OPMp(i,:) = scalpp(idx,:) + dist_to_head*OPMnn(i,:);


    z = dot(OPMnn(i,:),OPMp(i,:))/OPMnn(i,3);
    val = [0 0 z] - OPMp(i,:);
    if val(3) < 0
        val = -val;
    end
    OPMtt1(i,:) = val/norm(val);

    val = cross(OPMnn(i,:),OPMtt1(i,:));
    OPMtt2(i,:) = val/norm(val);
end

%% Visualize OPM positions
figure(3)
hold on
plot3(scalpp(:,1),scalpp(:,2),scalpp(:,3),'.')
scatter3(OPMp(:,1),OPMp(:,2),OPMp(:,3),20,'k','filled')
view(-90,0)
axis equal
axis tight
view(-90,0)
set(figure(3),'position',[0 200 400 400]);

%% Visualize OPM orientations
figure(4)
hold on
plot3(OPMp(:,1),OPMp(:,2),OPMp(:,3),'xk')
quiver3(OPMp(:,1),OPMp(:,2),OPMp(:,3),OPMnn(:,1),OPMnn(:,2),OPMnn(:,3),0.4)
quiver3(OPMp(:,1),OPMp(:,2),OPMp(:,3),OPMtt1(:,1),OPMtt1(:,2),OPMtt1(:,3),0.4)
quiver3(OPMp(:,1),OPMp(:,2),OPMp(:,3),OPMtt2(:,1),OPMtt2(:,2),OPMtt2(:,3),0.4)
view(-90,0)
axis equal
axis tight
view(-90,0)
set(figure(4),'position',[0 200 400 400]);
