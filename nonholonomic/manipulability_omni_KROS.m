clear all;
close all;

a = 0.333;
b = 0.316;
c = 0.384;
d = 0.088;
e = 0.107;
f = 0.0825;
    
panda=loadrobot("frankaEmikaPanda", "DataFormat", "row");
removeBody(panda, "panda_rightfinger");
removeBody(panda, "panda_leftfinger");
removeBody(panda, "panda_hand");

ik = inverseKinematics('RigidBodyTree', panda, 'SolverAlgorithm', 'BFGSGradientProjection');
weights = [1,1,1,1,1,1];

homeguess = panda.homeConfiguration;
initialguess = [0,0,0,-pi/2,0,pi/3,-pi/2];
randomguess = randomConfiguration(panda);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num=1;
angle = 0;
open_angle = 75;
door_width = 0.92;
while abs(angle) <= open_angle
    x = - door_width * sin(angle*pi/180); 
    y = door_width * cos(angle*pi/180);

    handle_P(num,:) = [x, y, 1, angle];
    num = num + 1;
    angle = angle + 2.5;
end
 
mobile_P(1,:) = [-0.45, 0.92, 0.35];
mobile_next_Pose = next_P_candidate(mobile_P(1,:));

num_handle_P= size(handle_P);
num_mobile_P= size(mobile_next_Pose);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for config(1,:)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Px = handle_P(1,1)  - mobile_P(1,1);
Py = handle_P(1,2)  - mobile_P(1,2);
Pz = handle_P(1,3)  - mobile_P(1,3);

MtoH_P(1,:) = handle_P(1,1:3) - mobile_P(1,:);
MtoH_R1 = [0, 1, 0, pi/2];
MtoH_R2 = [0, 0, 1, -pi/4];
MtoH_R3 = [0, 0, 1, ( handle_P(1,4) ) *pi/180];

while 1
    randomguess = randomConfiguration(panda);
    config_q(1,:) = ik('panda_link8',trvec2tform(MtoH_P(1,:))*axang2tform(MtoH_R3(1,:))*axang2tform(MtoH_R1(1,:))*axang2tform(MtoH_R2(1,:)),weights,randomguess);
    if check_JL(config_q(1,:)) == 0
        disp("DDDDDDDD")
        break;
    end
    disp("JL")
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cnt = 0;
max_manipulability_1 = 0;
max_manipul_cnt_1 = 0;

for i=2:num_handle_P(1,1)
    
    for j=1:num_mobile_P(1,1)
        MtoH_P(j,:) = handle_P(i,1:3) - mobile_next_Pose(j,:);
        MtoH_R1 = [0, 1, 0, pi/2];
        MtoH_R2 = [0, 0, 1, -pi/4];
        MtoH_R3 = [0, 0, 1, ( handle_P(i,4) ) *pi/180];
        
        q(j,:) = ik('panda_link8',trvec2tform(MtoH_P(j,:))*axang2tform(MtoH_R3(1,:))*axang2tform(MtoH_R1(1,:))*axang2tform(MtoH_R2(1,:)),weights,config_q(i-1,:));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% joint limits check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        jl_error_flag = check_JL(q(j,:));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  workspace check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tmp_T07 = make_T07(q,j);
        nominal_xyz = [tmp_T07(1,4), tmp_T07(2,4), tmp_T07(3,4)];
        real_xyz = [MtoH_P(j,1), MtoH_P(j,2), MtoH_P(j,3)];
        error(j,1) = abs(norm(nominal_xyz - real_xyz));

        ws_error_flag = 0;
        if error(j,1) > 0.01
            ws_error_flag = 1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if  ws_error_flag==1 || jl_error_flag == 1
            continue;
        else
            cnt = cnt +1;
            real_q(cnt,:) = q(j,:);
            real_MtoH_P(cnt,:) = MtoH_P(j,:);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Jacobian --> Manipulability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            theta1(cnt,1) = real_q(cnt,1);
            theta2(cnt,1) = real_q(cnt,2);
            theta3(cnt,1) = real_q(cnt,3);
            theta4(cnt,1) = real_q(cnt,4);
            theta5(cnt,1) = real_q(cnt,5);
            theta6(cnt,1) = real_q(cnt,6);
            theta7(cnt,1) = real_q(cnt,7);

            DH =  [ a,  theta1(cnt,1),    0,  -pi/2;
                    0,  theta2(cnt,1),    0,   pi/2;
                    b,  theta3(cnt,1),    f,   pi/2;
                    0,  theta4(cnt,1),    -f, -pi/2;
                    c,  theta5(cnt,1),    0,   pi/2;
                    0,  theta6(cnt,1),    d,   pi/2;
                    e,  theta7(cnt,1) + pi/4,    0,   0  ];

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
            Jacob70_constraint = Jacob70(1,:);


            Jacob07 = jacobian(T01, T12, T23, T34, T45, T56, T67);

            Jacob07_constraint(1:2,:) = Jacob07(1:2,:);
            Jacob07_constraint(3,:) = Jacob07(6,:);
            
            
            
            manipulability_1(cnt,1)=sqrt(det(Jacob07 * transpose(Jacob07)));
%             if real_MtoH_P(cnt,1) <=0.2 ||
            if  norm(real_MtoH_P(cnt,1)) > 0.8 || norm(real_MtoH_P(cnt,1)) < 0.3
%                 manipulability_1(cnt,1) = 0;
                continue
            end
            
            if manipulability_1(cnt,1) > max_manipulability_1
                % 이전 애를 저장해둠
                max_manipulability_2 = max_manipulability_1; 
                max_manipul_cnt_2 = max_manipul_cnt_1;

                max_manipulability_1 = manipulability_1(cnt,1);
                max_manipul_cnt_1 = cnt;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end

    % 다음 자세가 선정되지 않았다는 뜻 --> 그 전에서 다른 선택지로 가야함
    if max_manipul_cnt_1 == 0
        max_manipul_cnt_1 = 1;
        break
    end
    
    mobile_P(i,:) = mobile_next_Pose(max_manipul_cnt_1, :);
    mobile_next_Pose = next_P_candidate(mobile_P(i,:));

    config_q(i,:) = real_q(max_manipul_cnt_1, :);
    selected_manipulability(i,:) = manipulability_1(max_manipul_cnt_1,1);
    
    cnt = 0;
    manipulability_1 = 0;
    max_manipulability_1 = 0;
    max_manipul_cnt_1 = 0;
end

ALL_POSITION(1:i,1:3) = mobile_P(:,1:3);
ALL_POSITION(i+1:2*i,1:3) = handle_P(:,1:3);

for order = 1 : i
    ALL_POSITION(order,4) = (order-1) * 5;
end
for order = i+1 : 2*i
    ALL_POSITION(order,4) = (order-1-i) * 5;
end

figure(1);

plot3c(ALL_POSITION(:,1), ALL_POSITION(:,2), ALL_POSITION(:,3), ALL_POSITION(:,4), 'o', 'Angle');
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on
hold on 
axis equal

% show(panda,real_q(1,:),'PreservePlot',false,'visuals','on','collision','off','position', [mobile_P(num_handle_P(1,1),1), mobile_P(num_handle_P(1,1),2),mobile_P(num_handle_P(1,1),3),0]);

function jl_error_flag = check_JL(q)
    jl_error_flag = 0;
    if abs(q(1,1)) > 2.8972
        jl_error_flag = 1;
    end

    if abs(q(1,2)) > 1.7627
        jl_error_flag = 1;
    end

    if abs(q(1,3)) > 2.8973
        jl_error_flag = 1;
    end

    if -3.0717 > q(1,4) ||  q(1,4) > -0.0699
        jl_error_flag = 1;
    end

    if abs(q(1,5)) > 2.8972
        jl_error_flag = 1;
    end

    if q(1,6)  > 3.7524 ||  q(1,6) < -0.0174
        jl_error_flag = 1;
    end

    if q(1,7) > 2.8972 - pi/4 || q(1,7) < -2.8972 - pi/4
        jl_error_flag = 1;
    end
end

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
    DH =  [ a,  q(i,1),    0,  -pi/2;
            0,  q(i,2),    0,   pi/2;
            b,  q(i,3),    f,   pi/2;
            0,  q(i,4),    -f, -pi/2;
            c,  q(i,5),    0,   pi/2;
            0,  q(i,6),    d,   pi/2;
            e,  q(i,7) + pi/4,    0,   0  ];

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
function next_P_candidate = next_P_candidate(now_P)
    step = 0.05;
    next_P_candidate(1,:) = [now_P(1,1) + step*sind(45), now_P(1,2) - step*sind(45), now_P(1,3)];
    next_P_candidate(2,:) = [now_P(1,1) + step, now_P(1,2)       , now_P(1,3)];
    next_P_candidate(3,:) = [now_P(1,1) + step*sind(45), now_P(1,2) + step*sind(45), now_P(1,3)];

    next_P_candidate(4,:) = [now_P(1,1)       , now_P(1,2) - step, now_P(1,3)];
    next_P_candidate(5,:) = [now_P(1,1)       , now_P(1,2)       , now_P(1,3)];
    next_P_candidate(6,:) = [now_P(1,1)       , now_P(1,2) + step, now_P(1,3)];

    next_P_candidate(7,:) = [now_P(1,1) - step*sind(45), now_P(1,2) - step*sind(45), now_P(1,3)];
    next_P_candidate(8,:) = [now_P(1,1) - step, now_P(1,2)       , now_P(1,3)];
    next_P_candidate(9,:) = [now_P(1,1) - step*sind(45), now_P(1,2) + step*sind(45), now_P(1,3)];

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