close all;
clear all;
clc;

N=256;
N_F=1000;
M = 16;                      % Size of signal constellation
k = log2(M);                % Number of bits per symbol
n = N*N_F;                  % Number of bits to process
numSamplesPerSymbol = 1;    % Oversampling factor
rng('default')              % Use default random number generator
dataIn = randi([0 1],n,1);  % Generate vector of binary data
dataSymbolsIn = dataIn;            



%%%%%% Modulation %%%%
MOD=comm.RectangularQAMModulator;
MOD.NormalizationMethod='Average Power';
MOD.BitInput=true;
dataMod=step(MOD,dataSymbolsIn);
%dataMod = qammod(dataSymbolsIn, M);
%%%%%%%%%%%%%%%%%%%%

%%%% SC %%%%%%%%
%dataMod=ifft(dataMod);
%%%%%%%%%%%%%%



dataMod=reshape(dataMod,N,size(dataMod,1)/N);





all_data=[];
%%%%%%%  rician channel  %%%%%%
chan=ricianchan;
chan.KFactor = 7.5;
chan.ResetBeforeFiltering=0;
chan.InputSamplePeriod=(1/7.56)*1e-8;
chan.MaxDopplerShift=100;
chan.StoreHistory=0;
chan.PathDelays=[0,6.67*1e-6];
chan.AvgPathGaindB=[0,-10];

% chan.PathDelays=7.56*[0,...
%                 10*(1/7.56)*1e-6,...
%                 22*(1/7.56)*1e-6,...
%                 34*(1/7.56)*1e-6,...
%                 54*(1/7.56)*1e-6,...
%                 78*(1/7.56)*1e-6];
% chan.AvgPathGaindB=[0 -1 -9 -10 -15 -20];
    
%%%%% process  frame by frame %%%%%%
for i=1:size(dataMod,2)
   
    %%% add cp  %%%
    data=ifft(dataMod(:,i));
    receivedSignal=[data(N-1/8*N+1:end);data];  %add cp

    %%%%%%%  rician channel  %%%%%%
    receivedSignal = filter(chan,receivedSignal);
    
    %%%%%%% awgn  channel %%%%%%
    EbNo = 10;
    snr = EbNo + 10*log10(k) - 10*log10(numSamplesPerSymbol);
    receivedSignal = awgn(receivedSignal, snr, 'measured');
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    
    receivedSignal=receivedSignal(N/8+1:end);   %remove cp
    receivedSignal=fft(receivedSignal);
    
    
    %%%%%%%%%  equalization   %%%%%%%%%%%%
    matrix=catchchannel(int32(chan.PathDelays/chan.InputSamplePeriod),...
                            size(receivedSignal,1),...
                            chan.PathGains);
    receivedSignal=ideal_equalization(receivedSignal,matrix,snr);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    
    
    %%%%%%%%%% fft %%%%%%%%%%%%%%%%%%
    %%% remove cp%%%%%%

    
    
%     %%%  plot constellation %%%%%%%%%%%
%     sPlotFig = scatterplot(receivedSignal, 1, 0, 'g.');
%     hold on
%     scatterplot(dataMod(:,i), 1, 0, 'k*', sPlotFig);
%     hold on
%     scatterplot(receivedSignal, 1, 0, 'b.',sPlotFig);
%     pause;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    all_data=[all_data;receivedSignal];
end

%all_data=fft(all_data);

%%%%%  demodulation %%%%%%%%
DEM=comm.RectangularQAMDemodulator;
DEM.BitOutput=true;
%DEM.DecisionMethod='Approximate log-likelihood ratio';
DEM.NormalizationMethod='Average Power';
%DEM.OutputDataType=
dataSymbolsOut=step(DEM,all_data);
%dataSymbolsOut = qamdemod(all_data, M);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataOutMatrix=dataSymbolsOut;
%dataOutMatrix = de2bi(dataSymbolsOut,k);
dataOut = dataOutMatrix(:);              % Return data in column vector



%%%%%  calc ber  %%%%%%%%%%
[numErrors, ber] = biterr(dataIn, dataOut);
fprintf('\nThe bit error rate = %5.2e, based on %d errors\n', ...
    ber, numErrors)

