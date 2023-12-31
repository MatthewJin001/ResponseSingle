

clc; clear;
currentFolder = pwd;


n=1000;mx=1;my=1;
% K=[912.006 0 644.224;0 911.017 355.626;0 0 1];
% fx=K(1,1);
% fy=K(2,2);

sigma_pos=0.001;
sigma_R=0.1/57.3; %转换为弧度
sigma_t=0.0001; %单位是m

[eRib,etib,pij,pattern,cRe_truth,cte_truth,bRp_truth,btp_truth]=simuDataPos(n,mx,my,sigma_t,sigma_R,sigma_pos);

i=3;
cRe_truth'*pij(:,i)-cRe_truth'*cte_truth
eRib(:,:,i)*(bRp_truth*pattern+btp_truth)+etib(:,i)

Ri=zeros(3,3,n);ti=zeros(3,n);
for i=1:n
    Ri(:,:,i)=eRib(:,:,i);
    ti(:,i)=etib(:,i)*1000;
end
pij=pij*1000;

for num=30:1:1000
[Rcf,tcf,pcf,Rit,tit,pit,timecf,timeit] = Alg(Ri(:,:,1:num),ti(:,1:num),pij(:,1:num))
[~,rmseit(num-29)] = RMSEfunc(Rit,tit,pit,Ri,ti,pij);
end



plot(rmseit)

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


