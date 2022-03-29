function [Xglob] = loc_to_glob_coord_mixC_mat(X,xloc,A,no_active_particles,on_gpu)

% This is for a triangular prism
%
% X: (6x3) global coordinates of the corners [m]
% Q: (5x1) total flux on a face [m3/s]
% xloc: (3x1) local coordinates of point of observation
% 
% Q(1) face opposite nodes 1&4
% Q(2) face opposite nodes 2&5
% Q(3) face opposite nodes 3&6
% Q(4) bottom face
% Q(5) top face

% Mixed coordinate system:
%     [X(1,1) X(1,2) 0;...
%     X(2,1) X(2,2) 0; ...
%     X(3,1) X(3,2) 0; ...
%     X(4,1) X(4,2) 1;...
%     X(5,1) X(5,2) 1; ...
%     X(6,1) X(6,2) 1];


xA = reshape(X(1,1,:),[no_active_particles,1]);
yA = reshape(X(1,2,:),[no_active_particles,1]);

xB = reshape(X(2,1,:),[no_active_particles,1]);
yB = reshape(X(2,2,:),[no_active_particles,1]);

xC = reshape(X(3,1,:),[no_active_particles,1]);
yC = reshape(X(3,2,:),[no_active_particles,1]);

% Shape functions for a triangular prism in mixed coordinates:
N = [((xB.*yC - xC.*yB)./(2*A) + ((yB-yC)./(2*A)).*xloc(1,:)' + ((xC-xB)./(2*A)).*xloc(2,:)' ).*(1-xloc(3,:)'),...
    ((xC.*yA - xA.*yC)./(2*A) + ((yC-yA)./(2*A)).*xloc(1,:)' + ((xA-xC)./(2*A)).*xloc(2,:)').*(1-xloc(3,:)'),...
    ((xA.*yB - xB.*yA)./(2*A) + ((yA-yB)./(2*A)).*xloc(1,:)' + ((xB-xA)./(2*A)).*xloc(2,:)').*(1-xloc(3,:)'),...
    ((xB.*yC - xC.*yB)./(2*A) + ((yB-yC)./(2*A)).*xloc(1,:)' + ((xC-xB)./(2*A)).*xloc(2,:)').*(xloc(3,:)'),...
    ((xC.*yA - xA.*yC)./(2*A) + ((yC-yA)./(2*A)).*xloc(1,:)' + ((xA-xC)./(2*A)).*xloc(2,:)').*(xloc(3,:)'),...
    ((xA.*yB - xB.*yA)./(2*A) + ((yA-yB)./(2*A)).*xloc(1,:)' + ((xB-xA)./(2*A)).*xloc(2,:)').*(xloc(3,:)')];

if ~on_gpu
    Xglob = zeros(no_active_particles,3);
end

if ~on_gpu
for xg = 1:no_active_particles
   Xglob(xg,:) = N(xg,:)*X(:,:,xg);
end

else
     Xr_xv = X(:,1,:);
     Xr_yv = X(:,2,:);
     Xr_zv = X(:,3,:);
     
     Xr_x = reshape(Xr_xv,[6*no_active_particles,1]);
     Xr_y = reshape(Xr_yv,[6*no_active_particles,1]);
     Xr_z = reshape(Xr_zv,[6*no_active_particles,1]);
     
     Xr = [Xr_x,Xr_y,Xr_z];
     
     indices1_v = repmat([1:1:no_active_particles]',[1,6]);
     indices1 = sort(indices1_v(:))';
     indices2 = 1:1:no_active_particles*6;
     
     N_elems = reshape(N',1,[]);
     N_elems = real(N_elems);
     
     indices1 = gpuArray(indices1);
     indices2 = gpuArray(indices2);
     N_elems = gpuArray(N_elems);
     
     N_mat = sparse(indices1,indices2,N_elems);
     Xglob = N_mat*Xr;
     
end
