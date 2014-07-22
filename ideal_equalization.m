function X = ideal_equalization(Y,matrix,SNR)
N=size(Y,1);
A=matrix;

Y1=A'*Y;  %匹配滤波
A1=(A'*A+1/SNR*eye(N));

%%%取3对角%%%
A3diag=zeros(N,N);
for i =1:N
    for j=1:N
        if(abs(i-j)<=1)
            A3diag(i,j)=A1(i,j);
        end
    end
end
%%%%%%%%%%%%

A3ici=A1-A3diag;




%%%%%%%迭代循环外用于保存追赶法的中间变量，以简化计算量%%
a=ones(N,1);
b=diag(A3diag);
c=ones(N-1,1);

k=2;
m=1;
for i =1:N
    for j=1:N
        if((i-j)==1)
            a(k)=A3diag(i,j);
            k=k+1;  
        end
        if((j-i)==1)
             c(m)=A3diag(i,j);
             m=m+1;
        end     
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y3diag=Y1;
%%%%%%%第一步迭代%%%%%%
%迭代次数Nitr1
Nitr1=2;
for itr1=1:Nitr1

%%%%%求解三对角方程组 A3diag^-1*X1=Y1 %%%%%%
%%%%%%%%%%  追赶法  %%%%%%%%%%%%%%%%%%%
beta=ones(N-1,1);
beta(1)=c(1)/(b(1));
for i=2:N-1
    beta(i)=c(i)/(b(i)-a(i)*beta(i-1));
end
 tmp_y=ones(N,1);
 tmp_y(1)=Y3diag(1)/(b(1));
 for i=2:N
     tmp_y(i)=(Y3diag(i)-a(i)*tmp_y(i-1))/(b(i)-a(i)*beta(i-1));
 end
 X1=zeros(N,1);
 X1(N)=tmp_y(N);
 for i=N-1:-1:1
     X1(i)=tmp_y(i)-beta(i)*X1(i+1);
 end
 %%%%%%%%%%%%% 解得X1 %%%%%%%%%%%%%%%%
 
 
%%%X1即为X的MMSE估计 消除了相邻子载波对信号的干扰%%
    %%%对X1进行硬判决  %%
X1_hat=zeros(N,1);
for i=1:N
    
 %%%16 QAM硬判决%%%
threshold_r=max(real(X1))/2;
threshold_i=max(imag(X1))/2;
if ((abs(real(X1(i)))<=threshold_r)&&(abs(imag(X1(i)))<=threshold_i))
X1_hat(i)=1*sign(real(X1(i)))+1*sign(imag(X1(i)))*1i;
elseif ((abs(real(X1(i)))<threshold_r)&&(abs(imag(X1(i)))>threshold_i))
X1_hat(i)=1*sign(real(X1(i)))+3*sign(imag(X1(i)))*1i;
elseif ((abs(real(X1(i)))>threshold_r)&&(abs(imag(X1(i)))<threshold_i))
X1_hat(i)=3*sign(real(X1(i)))+1*sign(imag(X1(i)))*1i;
elseif ((abs(real(X1(i)))>threshold_r)&&(abs(imag(X1(i)))>threshold_i))
X1_hat(i)=3*sign(real(X1(i)))+3*sign(imag(X1(i)))*1i; 
end

 
    
    %%%QPSK硬判决%%%
%     % X1_hat(i)=0.707*sign(real(X1(i)))+0.707*sign(imag(X1(i)))*1i;
%      X1_hat(i)=1*sign(real(X1(i)))+1*sign(imag(X1(i)))*1i;
%     %%%%

end
 Y3ici=A3ici*X1_hat; %ici干扰项
 Y3diag=Y1-Y3ici;
end


%%%第一步循环结束 得到X1的判决X1_hat 
Adiag=diag(diag(A1));
Aici=A1-Adiag; %Aici代表了所有干扰项


%%%%第二步迭代%%%%%%
%%%迭代次数Nitr2
Nitr2=2;

for itr2=1:Nitr2
Yici=Aici*X1_hat;
Ydiag=Y1-Yici;

for i=1:N   
  
%%%16 QAM硬判决%%%
threshold_r=max(real(X1_hat))/2;
threshold_i=max(imag(X1_hat))/2;
if ((abs(real(X1_hat(i)))<=threshold_r)&&(abs(imag(X1_hat(i)))<=threshold_i))
X1_hat(i)=1*sign(real(X1_hat(i)))+1*sign(imag(X1_hat(i)))*1i;
elseif ((abs(real(X1_hat(i)))<threshold_r)&&(abs(imag(X1_hat(i)))>threshold_i))
X1_hat(i)=1*sign(real(X1_hat(i)))+3*sign(imag(X1_hat(i)))*1i;
elseif ((abs(real(X1_hat(i)))>threshold_r)&&(abs(imag(X1_hat(i)))<threshold_i))
X1_hat(i)=3*sign(real(X1_hat(i)))+1*sign(imag(X1_hat(i)))*1i;
elseif ((abs(real(X1_hat(i)))>threshold_r)&&(abs(imag(X1_hat(i)))>threshold_i))
X1_hat(i)=3*sign(real(X1_hat(i)))+3*sign(imag(X1_hat(i)))*1i; 
end
  
%      %%%QPSK硬判决%%%
%     % X1_hat(i)=0.707*sign(real(X1(i)))+0.707*sign(imag(X1(i)))*1i;
%      X1_hat(i)=1*sign(real(X1(i)))+1*sign(imag(X1(i)))*1i;
%     %%%%

end

X1_hat= 1./((conj(diag(Adiag)).*diag(Adiag))+1/SNR).*Ydiag;

end


%%%第二步迭代结束%%%%
X=X1_hat;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function X = ideal_equalization(Y,matrix,SNR)
% X=(matrix^-1)*Y;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function X = ideal_equalization(Y,matrix,SNR)
% X=matrix'*(matrix*matrix'+1/SNR*eye(size(matrix,1)))^-1*Y;
% 
% end