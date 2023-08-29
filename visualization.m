% figure()
panda=loadrobot("frankaEmikaPanda", "DataFormat", "row");
removeBody(panda, "panda_rightfinger");
removeBody(panda, "panda_leftfinger");
removeBody(panda, "panda_hand");

ik = inverseKinematics('RigidBodyTree', panda, 'SolverAlgorithm', 'LevenbergMarquardt');
weights = [1,1,1,1,1,1];
initialguess2 = [0,0,0,-pi/2,0,pi,0];

W=37;
% EE_P(1,:) = [handle_P(W,1)-mobile_P(W,1),handle_P(W,2)-mobile_P(W,2), handle_P(W,3)-mobile_P(W,3)];
% EE_R1 = [0, 1, 0, pi/2];
% EE_R2 = [0, 0, 1, -pi/4];
% EE_R3 = [0, 0, 1, ( handle_P(W,4) ) *pi/180];
% EE_R3 = [0, 0, 1, ( handle_P(W,4) - mobile_P(W,4) ) *pi/180];

% q(1,:) = ik('panda_link8',trvec2tform(EE_P(1,:))*axang2tform(EE_R3(1,:))*axang2tform(EE_R1(1,:))*axang2tform(EE_R2(1,:)),weights,initialguess2);

show(panda,config_q(W,:),'PreservePlot',false,'visuals','on','collision','off','position', [mobile_P(W,1), mobile_P(W,2),mobile_P(W,3),mobile_P(W,4)*pi/180]);
% show(panda,real_q(W,:),'PreservePlot',false,'visuals','on','collision','off','position', [mobile_next_Pose(W,1), mobile_next_Pose(W,2),mobile_next_Pose(W,3),0]);
hold on