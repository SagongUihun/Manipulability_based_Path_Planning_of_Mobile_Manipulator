a = 0.333;
b = 0.316;
c = 0.384;
d = 0.088;
e = 0.107;
f = 0.0825;
    
panda=loadrobot("frankaEmikaPanda", "DataFormat", "row");
removeBody(panda, "panda_rightfinger");
removeBody(panda, "panda_leftfinger");

ik = inverseKinematics('RigidBodyTree', panda, 'SolverAlgorithm', 'BFGSGradientProjection');
weights = [1,1,1,1,1,1];
initialguess = panda.homeConfiguration;

% n개의 random EE Position & Rotation

num=1;
step = 0.2;
for y = -0.855:step:0.855
    
    for x = 0.855:-step:-0.855
        position_xy(num, 1) = x;
        position_xy(num, 2) = y;
        num = num + 1;
    end

end

% x=[0.3, 0.35, 0.40, 0.45];
% y=[0.1, 0.0, 0.0, 0.0];
% num_data = size(x);

num_data= size(position_xy);

cnt = 0;

for i=1:num_data(1,1)
    EE_P(i,:) = [position_xy(i,1), position_xy(i,2), 0.625];
    EE_R1 = [0, 1, 0, pi/2];
    EE_R2 = [0, 0, 1, -pi/4];
    
    EE_R = [1,0,0,0;0,1,0,0;0,0,1,0; 0,0,0,1];
    %q(i,:) = ik('panda_link8',trvec2tform(EE_P(i,:))*axang2tform(EE_R1(1,:))*axang2tform(EE_R2(1,:)),weights,initialguess); 
    q(i,:) = ik('panda_link8',trvec2tform(EE_P(i,:))*EE_R,weights,initialguess); 
    %q(i,:) = config(i,:);
    for j=1:7
        if q(i,j) > pi
            q(i,j) = q(i,j) - 2*pi;
        end
        if q(i,j) < -pi
            q(i,j) = q(i,j) + 2*pi;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% \joint limits check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    jl_error_flag = 0;
    if abs(q(i,1) * 180/pi) > 166
        jl_error_flag = 1;
    end
    
    if abs(q(i,2) * 180/pi) > 101 
        jl_error_flag = 1;
    end
    
    if abs(q(i,3) * 180/pi) > 166 
        jl_error_flag = 1;
    end
    
    if -176 > q(i,4) * 180/pi ||  q(i,4) * 180/pi > -4
        jl_error_flag = 1;
    end
    
    if abs(q(i,5) * 180/pi) > 166 
        jl_error_flag = 1;
    end
    
    if -1 > q(i,6) * 180/pi ||  q(i,6) * 180/pi > 215
        jl_error_flag = 1;
    end
    
    if abs(q(i,7) * 180/pi) > 166
        jl_error_flag = 1;
    end
    
    
%     HAHA_T07 = make_T07(initialguess, 1);
%     nominal_xyz2 = [HAHA_T07(1,4), HAHA_T07(2,4), HAHA_T07(3,4)];
%     real_xyz2 = [d, 0, a+b+c-e];
%     error2 = abs(norm(nominal_xyz2 - real_xyz2));
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  workspace check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmp_T07 = make_T07(q,i);
    nominal_xyz = [tmp_T07(1,4), tmp_T07(2,4), tmp_T07(3,4)];
    real_xyz = [EE_P(i,1), EE_P(i,2), EE_P(i,3)];
    error = abs(norm(nominal_xyz - real_xyz));
    
    ws_error_flag = 0;
    if error > 1.0
        ws_error_flag = 1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if jl_error_flag ==1 || ws_error_flag==1
        continue;
    else
        cnt = cnt +1;
        real_q(cnt,:) = q(i,:);
        real_EE_P(cnt,:) = EE_P(i,:);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Jacobian --> Manipulability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    theta1(cnt,1) = real_q(cnt,1);
    theta2(cnt,1) = real_q(cnt,2);
    theta3(cnt,1) = real_q(cnt,3);
    theta4(cnt,1) = real_q(cnt,4);
    theta5(cnt,1) = real_q(cnt,5);
    theta6(cnt,1) = real_q(cnt,6);
    theta7(cnt,1) = real_q(cnt,7);

    DH =  [ a,  theta1(cnt,1),    0,   pi/2;
            0,  theta2(cnt,1),    0,  -pi/2;
            b,  theta3(cnt,1),    f,   pi/2;
            0,  theta4(cnt,1),    -f, -pi/2;
            c,  theta5(cnt,1),    0,   pi/2;
            0,  theta6(cnt,1),    d,   pi/2;
            e,  theta7(cnt,1),    0,   0  ];

    T01 = TF_matrix(1, DH);
    T12 = TF_matrix(2, DH);
    T23 = TF_matrix(3, DH);
    T34 = TF_matrix(4, DH);
    T45 = TF_matrix(5, DH);
    T56 = TF_matrix(6, DH);
    T67 = TF_matrix(7, DH);

    T76 = inv(T67);
    T65 = inv(T56);
    T54 = inv(T45);
    T43 = inv(T34);
    T32 = inv(T23);
    T21 = inv(T12);
    T10 = inv(T01);

    T07 = T01*T12*T23*T34*T45*T56*T67;
    T70 = inv(T07);
    
    Jacob70 = jacobian(T76, T65, T54, T43, T32, T21, T10);
    Jacob70_constraint = Jacob70(1:2,:);

    manipulability(cnt,1) = sqrt(det(Jacob70_constraint * transpose(Jacob70_constraint)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

Z = zeros(cnt,1);
plot3c(-real_EE_P(:,1), -real_EE_P(:,2), Z, manipulability(:,1), 'o', 'Manipulability');
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on
hold on 
axis equal
% show(panda);

% for i = 1:num_data(1,1)
%     plot3(real_EE_P(i,1), real_EE_P(i,2), real_EE_P(i,3));
%     hold on
%     axis equal
% end

% show(panda);
% figure(2);
% show(panda,real_q(5,:),'PreservePlot',false,'visuals','on','collision','off', 'position', [-real_EE_P(5,1), -real_EE_P(5,2),-0.625,0]);
% hold on;
% show(panda,real_q(10,:),'PreservePlot',false,'visuals','on','collision','off', 'position', [-real_EE_P(10,1), -real_EE_P(10,2),-0.625,0]);
% hold on;
% show(panda3,q(3,:),'PreservePlot',false,'visuals','on','collision','off', 'position', [-x(1,3), -y(1,3),-0.625, 0]);
% hold on;
% show(panda4,q(4,:),'PreservePlot',false,'visuals','on','collision','off', 'position', [-x(1,4), -y(1,4),-0.625,0]);
% plot(q(, y_real_out, "b");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% geoJacob = geometricJacobian(panda,q,'panda_hand');
% geoJacob = geometricJacobian(panda,randomConfiguration(panda),'panda_hand');
% panda_config=randomConfiguration(panda);
% showdetails(panda);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TF = TF_matrix(i, DH)
    d = DH(i,1);
    theta = DH(i,2);
    a = DH(i,3);
    alpha = DH(i,4);
    
    TF = [ cos(theta), -cos(alpha)*sin(theta),  sin(alpha)*sin(theta), a*cos(theta) ;
           sin(theta),  cos(alpha)*cos(theta), -sin(alpha)*cos(theta), a*sin(theta) ;
                0,            sin(alpha),              cos(alpha),          d       ;
                0,                0,                       0,               1      ];       
end

function T07 = make_T07(q, i)
    a = 0.333;
    b = 0.316;
    c = 0.384;
    d = 0.088;
    e = 0.107;
    f = 0.0825;
    DH =  [ a,  q(i,1),    0,   pi/2;
            0,  q(i,2),    0,  -pi/2;
            b,  q(i,3),    f,   pi/2;
            0,  q(i,4),    -f, -pi/2;
            c,  q(i,5),    0,   pi/2;
            0,  q(i,6),    d,   pi/2;
            e,  q(i,7),    0,   0  ];

    T01 = TF_matrix(1, DH);
    T12 = TF_matrix(2, DH);
    T23 = TF_matrix(3, DH);
    T34 = TF_matrix(4, DH);
    T45 = TF_matrix(5, DH);
    T56 = TF_matrix(6, DH);
    T67 = TF_matrix(7, DH);

    T07 = T01*T12*T23*T34*T45*T56*T67;
end

function Jacob = jacobian(T76, T65, T54, T43, T32, T21, T10)
    T75 = T76*T65;
    T74 = T75*T54;
    T73 = T74*T43;
    T72 = T73*T32;
    T71 = T72*T21;
    T70 = T71*T10;
    
    o0 = [0; 0; 0];
    o1 = T76(1:3,4);
    o2 = T75(1:3,4);
    o3 = T74(1:3,4);
    o4 = T73(1:3,4);
    o5 = T72(1:3,4);
    o6 = T71(1:3,4);
    o7 = T70(1:3,4);
    
    z0 = [0; 0; 1];
    z1 = T76(1:3,3);
    z2 = T75(1:3,3);
    z3 = T74(1:3,3);
    z4 = T73(1:3,3);
    z5 = T72(1:3,3);
    z6 = T71(1:3,3);
    z7 = T70(1:3,3);
    
    Jacob(1:3, 1) = cross(z0, o7 - o0);
    Jacob(1:3, 2) = cross(z1, o7 - o1);
    Jacob(1:3, 3) = cross(z2, o7 - o2);
    Jacob(1:3, 4) = cross(z3, o7 - o3);
    Jacob(1:3, 5) = cross(z4, o7 - o4);
    Jacob(1:3, 6) = cross(z5, o7 - o5);
    Jacob(1:3, 7) = cross(z6, o7 - o6);
    
    Jacob(4:6, 1) = z0;
    Jacob(4:6, 2) = z1;
    Jacob(4:6, 3) = z2;
    Jacob(4:6, 4) = z3;
    Jacob(4:6, 5) = z4;
    Jacob(4:6, 6) = z5;
    Jacob(4:6, 7) = z6;
end

function plot3c(x,y,z,v,marker,string)
%FUNCTION PLOT3C(X,Y,Z,V,'MARKER','TitleString') plots the values in vector v colour coded
% at the positions specified by the vectors x, y, and z in a 3-D axis
% system. A colourbar is added on the right side of the figure.
%
% The colorbar strectches from the minimum value of v to its
% maximum in 9 steps (10 values).
%
% The second last argument is optional to define the marker being used. The
% default is a point. To use a different marker (such as circles, ...) send
% its symbol to the function (which must be enclosed in '; see example):
% plot3c(X,Y,Z,V,'o')
%
% A title can be optionally added to the colorbar.:
% plot3c(X,Y,Z,V,'o','Title')
% 
% This function is an extension of PLOTC.
%
% Example:
% The seismic P-velocity (v) depends on the three parameters porosity (por) and the
% bulk moduli of the saturating fluid (kf) and the elastic frame (kd). To plot the
% velocity data as a function of these three parameters use (assuming that
% all data are given in vectors):
%
% plot3c(por,kd,kf,v,'d','Velocity')
%
% Uli Theune, University of Alberta, 2004
% utheune@phys.ualberta.ca
%

if mod(length(x)+length(y)+length(z)+length(v),4)
    disp('All vectors must be of same length')
    return
end
delete(gca)
if nargin <5
    marker='.';
end
if nargin < 6
    string=' ';
end
% Define the data range
miv=min(v);
mav=max(v);
% Get the current colormap
map=colormap;
% Plot the points
hold on
for i=1:length(x)
    in=round((v(i)-miv)*(length(map)-1)/(mav-miv));
    %--- Catch the out-of-range numbers
    if in==0;in=1;end
    if in > length(map);in=length(map);end
    plot3(x(i),y(i),z(i),marker,'color',map(in,:),'markerfacecolor',map(in,:))
end
hold off

% Re-format the colorbar
h=colorbar;

set(h,'ylim',[1 length(map)]);
yal=linspace(1,length(map),10);
set(h,'ytick',yal);
% Create the yticklabels
ytl=linspace(miv,mav,10);
s=char(10,4);
for i=1:10
    if abs(min(log10(abs(ytl)))) <= 3
        B=sprintf('%-4.3f',ytl(i));
    else
        B=sprintf('%-4.2E',ytl(i));
    end
    s(i,1:length(B))=B;
end
set(h,'yticklabel',s,'fontsize',9);
grid on
set(get(h,'title'),'string',string,'fontweight','bold')
view(3)
end