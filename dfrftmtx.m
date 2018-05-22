function F_a=dfrftmtx(N,a)
EYE_N=eye(N);
F_a=zeros(N,N);
for k=1:1:N
    temp=dfrft(N,a,EYE_N(k,:));
    F_a(:,k)=transpose(temp);
end
F_a