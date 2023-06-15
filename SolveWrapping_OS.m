function [lmt,Ld,event,Jm,wrap_points] = SolveWrapping_OS(model,Xup,V,J,midx)
% Solve the muscle wrapping
% model struct variable which contains the information of the whole
%       biomechanical system
% Xup   6D Kinematic Transforms of the system
% midx  muscle index
% u0    initial condition for the u parameter of the parametric equations
% v0    initial condition for the v parameter of the parametric equations
% lg0   initial condition for the geodesic length

widx = model.m{midx}.wgeom;

% transform all via points
for i=1:size(model.m{midx}.P,2)
    if model.m{midx}.P{1,i}.parent == model.appearance.wrapCyl{widx}.parent
        R = model.appearance.wrapCyl{widx}.R;   % Rotation from parent to wrapping geometry
        p = model.appearance.wrapCyl{widx}.tra; % Translation from parent to wrapping geometry
    else
        T = inv(pluho(Xup{model.m{midx}.P{1,i}.parent})); % Transformation from path point's parent to its parent
        R = T(1:3,1:3)'*model.appearance.wrapCyl{widx}.R;      % Rotation from parent to wrapping geometry
        p = T(1:3,1:3)'*(model.appearance.wrapCyl{widx}.tra - T(1:3,4)); % Translation from parent to wrapping geometry
    end
    wp(:,i) = R'*(model.m{midx}.P{1,i}.location - p);
end

% rotate the path points to mantain the cylinder in 0,0,0 coordinates and
% the heigth aligned with z axis, so the points will be expressed in that
% coordinate frame
p0 = (rx(pi/2)'*wp(:,model.m{midx}.wrange(1)))'*100;
q0 = (rx(pi/2)'*wp(:,model.m{midx}.wrange(2)))'*100;

% initial points on the wrap geometry (cylinder)
r = model.m{midx}.wd*model.appearance.wrapCyl{widx}.radius*100;
h = model.appearance.wrapCyl{widx}.height*100;

%% Muscle Wrapping Algorithm
out = SolveWrap(p0,q0,r);

p = out{1,1}; % p point
q = out{1,2}; % q point

p0 = p0/100; q0 = q0/100;
p = p/100; q = q/100;

if out{1,5} == 2 % if out{1,5} == 2 wrapping does occur
    
    lmt = 0; j = 1;
    for i=1:size(model.m{midx}.P,2)
        if model.m{midx}.wrange(1)==i
            j = j+1;
        elseif model.m{midx}.wrange(2)==i
%           lmt = lmt + norm(wp(:,j)-rx(-pi/2)'*q');
        else
            lmt = lmt + norm(wp(:,j+1)-wp(:,j));
            j = j+1;
        end
    end

    la = out{1,6}/100;  % curved segment length
    lmt = lmt + la;
    event = 1;
    
    % p in shoulder frame (parent frame of the wrapping geometry)
    p_s = model.appearance.wrapCyl{widx}.tra + model.appearance.wrapCyl{widx}.R*(rx(-pi/2)'*(p'));
    % q in shoulder frame (parent frame of the wrapping geometry)
    q_s = model.appearance.wrapCyl{widx}.tra + model.appearance.wrapCyl{widx}.R*(rx(-pi/2)'*(q'));

    wp = [model.m{midx}.P{1,model.m{midx}.wrange(1)}.location,...
          model.m{midx}.P{1,model.m{midx}.wrange(2)}.location,...
          p_s, q_s];
    wg = model.appearance.wrapCyl{widx}.tra;
    % wp = [ps0,q0_s,p_s,q_s]; % wrapping points are concatenated in columns, fixed points first, contact points last
    % compute wrapping points in world coordinates
    [wp_w, e, vp_w, Jm] = wrapPosVel(model,midx,widx,Xup,wp,wg,V,J);

    % lmt velocity
    lmt_d = dot(e{1},vp_w{3}) - dot(e{1},vp_w{1}) + dot(e{2},vp_w{2}) - dot(e{2},vp_w{3});

elseif out{1,5} == 0

    lmt = 0;
    for i=1:size(model.m{midx}.P,2)-1
        lmt = lmt + norm(wp(:,i+1)-wp(:,i));
    end
    
    event = 0;
    
    wp = [model.m{midx}.P{1,model.m{midx}.wrange(1)}.location,...
          model.m{midx}.P{1,model.m{midx}.wrange(2)}.location];
    
    [wp_w, e, vp_w, Jm] = nowrapPosVel(model,midx,Xup,wp,V,J);
    
    % lmt velocity
    lmt_d = dot(e{1},vp_w{2}-vp_w{1});
    
end

Ld = lmt_d;

wrap_points = out{1,3};

end

function out = SolveWrap(p0_l,q0_l,Radius)
    noWrap = 0;    
    insideRadius = 1;
    wrapped = 2;
    
    max_wrap_pts_circle_ang = (5.0/360.0)*2*pi;
    
    p = zeros(1,3);
    q = p;
    R = Radius;
    
    % Compute displacements of P and S from cylinder axis
    Px=p0_l(1); Py=p0_l(2); Pz=p0_l(3); dP=Px*Px+Py*Py; rootP=dP-R*R;
    Sx=q0_l(1); Sy=q0_l(2); Sz=q0_l(3); dS=Sx*Sx+Sy*Sy; rootS=dS-R*R;
    
    %dPx=dp0_l(1); dPy=dp0_l(2); dPz=dp0_l(3);
    %dSx=dq0_l(1); dSy=dq0_l(2); dSz=dq0_l(3);
    
    % Check P and S against cylinder, and compute x and y components of wrap points Q and T
    if( rootP<0.0 || rootS<0.0 ) % One of P or S lies within the cylinder
        out = {p,q,0,0,insideRadius};
        return
    end
    dP=R/dP; rootP=sqrt(rootP); Qx=(R*Px-rootP*Py)*dP; Qy=(R*Py+rootP*Px)*dP;
    dS=R/dS; rootS=sqrt(rootS); Tx=(R*Sx+rootS*Sy)*dS; Ty=(R*Sy-rootS*Sx)*dS;
    
    % Apply the 180-degree wrapping rule to see if contact is appropriate (i.e. wrap > 180 = no contact)
    if( R*(Qx*Ty-Qy*Tx) < 0.0 )
        path_length = norm(q0_l - p0_l);
        %d_path_length = ((q0_l - p0_l)/path_length)*(dq0_l - dp0_l)';
        out = {p,q,0,0,noWrap,path_length};
        return
    end
    
    % Compute respective wrapping segment lengths
    PQ = sqrt( (Qx-Px)*(Qx-Px) + (Qy-Py)*(Qy-Py) );
    TS = sqrt( (Tx-Sx)*(Tx-Sx) + (Ty-Sy)*(Ty-Sy) );
    QtoTang = acos( 1.0 - 0.5*( (Qx-Tx)*(Qx-Tx) + (Qy-Ty)*(Qy-Ty) )/(R*R) );
    QT = R*QtoTang;
    if(QT<0.0) 
        QT=-QT;
    end
    
    % Assign z-axis components of wrap points Q and T
    Qz = Pz + (Sz-Pz)*(PQ) / (PQ+TS+QT);
    Tz = Sz + (Pz-Sz)*(TS) / (PQ+TS+QT);
    
    PQm = sqrt( (Qx-Px)*(Qx-Px) + (Qy-Py)*(Qy-Py) + (Qz-Pz)*(Qz-Pz));
    TSm = sqrt( (Tx-Sx)*(Tx-Sx) + (Ty-Sy)*(Ty-Sy) + (Tz-Sz)*(Tz-Sz));
    QTm = sqrt(QT^2 + (Tz - Qz)^2);
    
    % Register results and return
    wrap_path_length = PQm + TSm + QTm;
    p(1) = Qx;  p(2) = Qy;  p(3) = Qz;  
    q(1) = Tx;  q(2) = Ty;  q(3) = Tz;
    
    % Generate wrap_pts sequence of points tracing out wrapping path
    Qang = atan2(Qy,Qx);               % Angle of point Q
    if QtoTang>=0 
        d = QtoTang;
    else
        d = -QtoTang;
    end
    %n = 1 + int32(d/max_wrap_pts_circle_ang); % Number of angle steps from Q to T angles
    n = 13;
    
    angDelt=(QtoTang)/double(n); % Delta angle for n steps from Q to T angles
    
    for i=1:n
        ang = Qang + (double(i-1))*angDelt; % Angle ranging from that of Q to that of T
        aPointi{i} = [R*cos(ang), R*sin(ang), Qz+(Tz-Qz)*double(i-1)/double(n)]./100;
    end
    
    out = {p,q,aPointi,n,wrapped,wrap_path_length};
    
end