function X_k=dfrft(N,a,x_n)
%DFRFT_Eigen_decomposition_type
%fprintf('DFRFT_Eigen_decomposition_type\n');
%==========================================================================
%Settings
%N=7;
%a=1;
%==========================================================================
%Step:1 Generate matrices S and P
%fprintf('Step:1 Generate matrices S and P\n');
S=zeros(N,N);
for n=1:1:N
    C_n=2*(cos(2*pi/N*(n-1)));%-4;
	S(n,n)=C_n;
	if(n~=N)
        S(n,n+1)=1; 
    end
	if(n~=1)
        S(n,n-1)=1;
    end
end
S(1,N)=1;
S(N,1)=1;
%fprintf('matrices S:\n');
%S;

P=zeros(N,N);
if(mod(N,2)==1)
%    fprintf('N is Odd\n');
    for n=2:1:N
        P(n,N+2-n)=1;
        if (n<=((N-1)/2+1))
            P(n,n)=1;
        end
        if (n>((N-1)/2+1))
            P(n,n)=-1;
        end
    end
    P(1,1)=sqrt(2);
end
if(mod(N,2)~=1)
%    fprintf('N is Even\n');%P matrix when N is even is wrong!!
    for n=2:1:N
        P(n,N+2-n)=1;
        if (n<=(floor((N-1)/2))+1)
            P(n,n)=1;
        end
        if (n>floor((N-1)/2+1)+1)
            P(n,n)=-1;
        end
    end
    P(floor((N-1)/2)+2,floor((N-1)/2)+2)=sqrt(2);
    P(1,1)=sqrt(2);
end
P=P/sqrt(2);
%fprintf('matrices P:\n');
%P;
%==========================================================================
%Step:2 Generate the Ev and Od matrices 
%fprintf('Step:2 Generate the Ev and Od matrices\n');
P_S_INVP=P*S/P;%do not use \p 
%P_S_INVP;
Ev=P_S_INVP(1:floor(N/2)+1,1:floor(N/2)+1);
Od=P_S_INVP(floor(N/2)+2:N,floor(N/2)+2:N);
%fprintf('Ev:\n');
%Ev;
%fprintf('Od:\n');
%Od;
%==========================================================================
%Step:3 Find the eigenvectors/eigenvalues of Ev and Od 
%fprintf('Step:3 Find the eigenvectors/eigenvalues of Ev and Od \n');
[eigenVector_Ev,eigenValue_Ev]=eig(Ev);
[eigenVector_Od,eigenValue_Od]=eig(Od);
%fprintf('Ev:\n');
%eigenValue_Ev; 
%eigenVector_Ev;
%fprintf('Od:\n');
%eigenValue_Od;
%eigenVector_Od;

%Step:4 Sort the eigenvectors of Ev(Od)in the descending order of
%eigenvalues of Ev(Od) and denote the sorted eigenvectors as e_k(o_k)

[eigenValue_Ev_sort,index] = sort(diag(eigenValue_Ev),'descend');
eigenValue_Ev_sort = diag(eigenValue_Ev_sort);
eigenVector_Ev_sort = eigenVector_Ev(:,index);

[eigenValue_Od_sort,index] = sort(diag(eigenValue_Od),'descend');
eigenValue_Od_sort = diag(eigenValue_Od_sort);
eigenVector_Od_sort = eigenVector_Od(:,index);
e_k=eigenVector_Ev_sort;
o_k=eigenVector_Od_sort;
%fprintf('Ev:\n');
%eigenValue_Ev_sort ;
%e_k;
%fprintf('Od:\n')
%eigenValue_Od_sort;
%o_k;

%Step:5 Let u_2k[n]=P[e_k^T|0...0]^T
%       Let u_2k+1[n]=P[0...0|o_k^T]^T
u_2k=P*transpose([transpose(e_k),zeros(floor(N/2)+1,N-floor(N/2)-1)]);
u_2k_1=P*transpose([zeros(N-floor(N/2)-1,N-(N-floor(N/2)-1)),transpose(o_k)]);
%fprintf('u_2k[n]:\n');
%u_2k;
%fprintf('u_2k+1[n]:\n');
%u_2k_1;


u_k=zeros(N,N);
for k=0:1:floor((N-1)/2)
    u_k(:,(2*k)+1)=u_2k(:,k+1);
    if(k<=floor((N-3)/2))
        u_k(:,(2*k+1)+1)=u_2k_1(:,k+1); 
    end
end
if (mod(N,2)==0)
    u_k(:,N)=u_2k(:,floor(N/2)+1);
end

%fprintf('u_k[n]:\n');
%u_k;
%Step:6 Define F^a[m,n]=sum(on k)( u_k[m,n]*exp(-j*pi*k*a/2)*u_k[n] ) 
%            k belongs to M and k!=N-1+N_2 ; M={0,...,N-2,(N-N_2)}
F_a=zeros(N,N);


%shift the order of input signal
shift = rem((0:N-1) + fix(N/2),N)+1;
X_k=zeros(1,N);
x_n=transpose(x_n);


if (mod(N,2)==0)
    X_k(1,shift) = u_k*(exp(-j*pi/2*a*([0:N-2 N])).' .*(u_k'*x_n(shift)));
end
if (mod(N,2)==1)
     X_k(1,shift) = u_k*(exp(-j*pi/2*a*([0:N-1])).' .*(u_k'*x_n(shift)));
end
