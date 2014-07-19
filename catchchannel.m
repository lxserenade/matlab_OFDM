function ideal_channel_matrix=catchchannel(delaytaps,Nfft,u)

ideal_channel_matrix=zeros(Nfft,Nfft);
ideal_channel=zeros(Nfft,Nfft);

DSCchannelgain = repmat(u,Nfft,1);
ideal_channel(:,delaytaps+1)=DSCchannelgain;
for idx=1:Nfft
    for idy=1:Nfft
        ideal_channel_matrix(idx,idy)=ideal_channel(idx,mod(idx-idy,Nfft)+1);  %构造时变信道循环卷积矩阵
    end
end

w=fft(eye(Nfft))/sqrt(Nfft);
ideal_channel_matrix = w*ideal_channel_matrix*w';
end