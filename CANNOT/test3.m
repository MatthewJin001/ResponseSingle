

clc; clear;
currentFolder = pwd;
% 
% addpath ([currentFolder,'\dataset\'])
% addpath ([currentFolder,'\auxiliary\'])
% addpath ([currentFolder,'\method\'])
% addpath ([currentFolder,'\result\'])

n=30;mx=1;my=1;
K=[912.006 0 644.224;0 911.017 355.626;0 0 1];
fx=K(1,1);
fy=K(2,2);

sigma_pix=1;
sigma_R=1/57.3; %转换为弧度
sigma_t=0.001; %单位是m

[eRib,etib,qij,pattern,cRe_truth,cte_truth,bRp_truth,btp_truth]=simuData(n,mx,my,K,sigma_t,sigma_R,sigma_pix);

Ri=zeros(3,3,n);ti=zeros(3,n);
for i=1:n
    Ri(:,:,i)=eRib(:,:,i);
    ti(:,i)=etib(:,i);
end


step=100;

% Rk=cRe_truth;
% tk=cte_truth;
% pk=bRp_truth*pattern+btp_truth;

Rk=eye(3);
tk=zeros(3,1);
pk=zeros(3,1);

k=1;
cost=zeros(1,20);

while (step>10e-5 && k<21)

    A=zeros(9,9); B=zeros(9,1);
    for i=1:n
        g=Rk*Ri(:,:,i)*pk+Rk*ti(:,i)+tk;
        f=qij(:,i)-homo(K*g);
        J=-[fx/g(3),0,-fx*g(1)/g(3)/g(3);0,fy/g(3),-fy*g(2)/g(3)/g(3)]*...
            [-skew(g), eye(3), Rk*Ri(:,:,i)];
        A=A+J'*J;
        B=B+J'*f;
    end

    delta=-inv(A)*B;

    Rk=expm(skew(delta(1:3)))*Rk;
    tk=delta(4:6)+tk;
    pk=delta(7:9)+pk;

    step=norm(delta)
    

    cost_temp=0;
    for i=1:n
       cost_temp=cost_temp+norm(qij(:,i)-homo(K*(Rk*Ri(:,:,i)*pk+Rk*ti(:,i)+tk)))^2;
    end
    cost(k)=sqrt(cost_temp/n);
    
    k=k+1;
    
end

hold on;
plot(cost)


