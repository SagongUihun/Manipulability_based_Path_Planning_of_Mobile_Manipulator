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

ik = inverseKinematics('RigidBodyTree', panda, 'SolverAlgorithm', 'LevenbergMarquardt');
weights = [1,1,1,1,1,1];
initialguess1 = [0,-pi/6,0,-pi/1.5,0,pi/2,0];
initialguess2 = [0,0,0,-0.1,0,pi/2,0];


EE_P(1,:) = [0.2,-0.755, 0.625];
EE_R1 = [0, 1, 0, pi/2];
EE_R2 = [0, 0, 1, -pi/4];

q_1(1,:) = ik('panda_link8',trvec2tform(EE_P(1,:))*axang2tform(EE_R1(1,:))*axang2tform(EE_R2(1,:)),weights,initialguess1); 



EE_P(1,:) = [0.25,-0.6, 0.625];
EE_R1 = [0, 1, 0, pi/2];
EE_R2 = [0, 0, 1, -pi/4];

q_2(1,:) = ik('panda_link8',trvec2tform(EE_P(1,:))*axang2tform(EE_R1(1,:))*axang2tform(EE_R2(1,:)),weights,initialguess1); 



EE_P(1,:) = [0.50,0, 0.625];
EE_R1 = [0, 1, 0, pi/2];
EE_R2 = [0, 0, 1, -pi/4];

q_3(1,:) = ik('panda_link8',trvec2tform(EE_P(1,:))*axang2tform(EE_R1(1,:))*axang2tform(EE_R2(1,:)),weights,initialguess1); 



EE_P(1,:) = [0.75,0, 0.625];
EE_R1 = [0, 1, 0, pi/2];
EE_R2 = [0, 0, 1, -pi/4];

EE_R = [1,0,0,0;0,1,0,0;0,0,1,0; 0,0,0,1];
q_4(1,:) = ik('panda_link8',trvec2tform(EE_P(1,:))*axang2tform(EE_R1(1,:))*axang2tform(EE_R2(1,:)),weights,initialguess1); 


%  for j=1:7
%     if initialguess1(i,j) > 2*pi
%         initialguess1(i,j) = initialguess1(i,j) - 2*pi;
%     end
%     if initialguess1(i,j) < -2*pi
%         initialguess1(i,j) = initialguess1(i,j) + 2*pi;
%     end
%  end
Jacob07_1 = make_jacobian(initialguess1,1);
Jacob07_2 = make_jacobian(initialguess2,1);
Jacob07_3 = make_jacobian(q_1,1);
Jacob07_4 = make_jacobian(q_2,1);
Jacob07_5 = make_jacobian(q_3,1);
Jacob07_6 = make_jacobian(q_4,1);

manipulability1 = sqrt(det(Jacob07_1 * transpose(Jacob07_1)));
manipulability2 = sqrt(det(Jacob07_2 * transpose(Jacob07_2)));
manipulability3 = sqrt(det(Jacob07_3 * transpose(Jacob07_3)));
manipulability4 = sqrt(det(Jacob07_4 * transpose(Jacob07_4)));
manipulability5 = sqrt(det(Jacob07_5 * transpose(Jacob07_5)));
manipulability6 = sqrt(det(Jacob07_6 * transpose(Jacob07_6)));

figure(1);
show(panda,initialguess1,'PreservePlot',false,'visuals','on','collision','off');
figure(2);
show(panda,initialguess2,'PreservePlot',false,'visuals','on','collision','off');
figure(3);
show(panda,q_1,'PreservePlot',false,'visuals','on','collision','off');
figure(4);
show(panda,q_2,'PreservePlot',false,'visuals','on','collision','off');
figure(5);
show(panda,q_3,'PreservePlot',false,'visuals','on','collision','off');
figure(6);
show(panda,q_4,'PreservePlot',false,'visuals','on','collision','off');
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

function Jacobian = jacobian(T01, T12, T23, T34, T45, T56, T67)
    T02 = T01*T12;
    T03 = T02*T23;
    T04 = T03*T34;
    T05 = T04*T45;
    T06 = T05*T56;
    T07 = T06*T67;
    
    o0 = [0; 0; 0];
    o1 = T01(1:3,4);
    o2 = T02(1:3,4);
    o3 = T03(1:3,4);
    o4 = T04(1:3,4);
    o5 = T05(1:3,4);
    o6 = T06(1:3,4);
    o7 = T07(1:3,4);
    
    z0 = [0; 0; 1];
    z1 = T01(1:3,3);
    z2 = T02(1:3,3);
    z3 = T03(1:3,3);
    z4 = T04(1:3,3);
    z5 = T05(1:3,3);
    z6 = T06(1:3,3);
    z7 = T07(1:3,3);
    
    Jacobian(1:3, 1) = cross(z0, o7 - o0);
    Jacobian(1:3, 2) = cross(z1, o7 - o1);
    Jacobian(1:3, 3) = cross(z2, o7 - o2);
    Jacobian(1:3, 4) = cross(z3, o7 - o3);
    Jacobian(1:3, 5) = cross(z4, o7 - o4);
    Jacobian(1:3, 6) = cross(z5, o7 - o5);
    Jacobian(1:3, 7) = cross(z6, o7 - o6);
    
    Jacobian(4:6, 1) = z0;
    Jacobian(4:6, 2) = z1;
    Jacobian(4:6, 3) = z2;
    Jacobian(4:6, 4) = z3;
    Jacobian(4:6, 5) = z4;
    Jacobian(4:6, 6) = z5;
    Jacobian(4:6, 7) = z6;
end
function Jacob = make_jacobian(q, i)
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
            e,  q(i,7)+pi/4,    0,   0  ];

    T01 = TF_matrix(1, DH);
    T12 = TF_matrix(2, DH);
    T23 = TF_matrix(3, DH);
    T34 = TF_matrix(4, DH);
    T45 = TF_matrix(5, DH);
    T56 = TF_matrix(6, DH);
    T67 = TF_matrix(7, DH);
    
    Jacob =jacobian(T01, T12, T23, T34, T45, T56, T67);
end

