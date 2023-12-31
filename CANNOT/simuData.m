function [eRib,etib,qij,pattern,cRe_truth,cte_truth,bRp_truth,btp_truth]=simuData(n,mx,my,K,sigma_t,sigma_R,sigma_pix)

% n=20;
% mx=4;
% my=6;
% K=[912.006 0 644.224;0 911.017 355.626;0 0 1];
% sigma_pix=0.1;
% sigma_R=0.2/57.3; %转换为弧度
% sigma_t=0.0001; %单位是m
%%
cRe=angle2dcm(0.6865,0.0195,-0.7830);
cte=[0.0338;0.0570;-0.0155];
cTe=[cRe,cte;0,0,0,1];
cRe_truth=cRe;
cte_truth=cte;

bRp=angle2dcm(-1.5525,-0.0124,3.1388);
btp=[0.5362;0.1162;-0.1523];
bTp=[bRp,btp;0,0,0,1];
bRp_truth=bRp;
btp_truth=btp;

m=my*mx;
pattern=zeros(3,m);
for i=1:mx
    for j=1:my
        pattern(:,(i-1)*my+j)= [(j-1)*0.020;(i-1)*0.020;0];
    end
end


robot_roll=-1.7066+ (-0.3789+1.7066)*rand(1,n);
robot_pitch=-0.9643+ (-0.2018+0.9643)*rand(1,n);
robot_yaw=-3.0407+ (-2.0808+3.0407)*rand(1,n);
robot_x=-0.5795+(-0.1316+0.5795)*rand(1,n);
robot_y=0.3136+(0.6986-0.3136)*rand(1,n);
robot_z=-0.1054+(0.1975+0.1054)*rand(1,n);



%output
qij=zeros(2,n,m);
eRib=zeros(3,3,n);etib=zeros(3,n);
for i=1:n
    eRib(:,:,i)=angle2dcm(robot_roll(i),robot_pitch(i),robot_yaw(i));
    etib(:,i)=[robot_x(i);robot_y(i);robot_z(i)];
    eTb=[eRib(:,:,i),etib(:,i);0,0,0,1];
    for j=1:m
        qij(:,i,j)= homo(K*homo(cTe*eTb*bTp*[pattern(:,j);1]))+normrnd(0,sigma_pix,[2 1]);   
    end
      eRib(:,:,i)=angle2dcm(robot_roll(i)+normrnd(0,sigma_R),robot_pitch(i)+normrnd(0,sigma_R),robot_yaw(i)+normrnd(0,sigma_R));
      etib(:,i)=[robot_x(i)+normrnd(0,sigma_t);robot_y(i)+normrnd(0,sigma_t);robot_z(i)+normrnd(0,sigma_t)];
end

end
