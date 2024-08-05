clc
clear
%% generation of data
Number_of_Data = 1000*128;
Number_of_Data8 = 1000*72;
Data=randi([0 1],1,Number_of_Data);
Data8=randi([0 1],1,Number_of_Data8);
%% BFSK
Eb=1;
Tb=1;
phi_1=sqrt(2/Tb);
phi_2=sqrt(2/Tb);
S1=phi_1;
S2=0+phi_2*i;
%% generating Symbols
M=2;
BPSK_Symbols = reshape(Data,[Number_of_Data/log2(M) , log2(M)]);
BFSK_Symbols = reshape(Data,[Number_of_Data/log2(M) , log2(M)]);
M = 4;
QPSK_Symbols = reshape(Data,[Number_of_Data/log2(M) , log2(M)]);
M = 8;
MPSK_Symbols = reshape(Data8,[Number_of_Data8/log2(M) , log2(M)]);
M = 16;
QAM_Symbols = reshape(Data,[Number_of_Data/log2(M) , log2(M)]);


%BPSK
M = 2;
String_BPSK = cell(Number_of_Data/log2(M),1);
String_BFSK= cell(Number_of_Data/log2(M),1);

BPSK_Dec = zeros(Number_of_Data/log2(M) ,1);
BFSK_Dec=zeros(Number_of_Data/log2(M) ,1);

Mapped_BPSK = zeros(Number_of_Data/log2(M) ,1);
Mapped_BFSK =zeros(Number_of_Data/log2(M)  ,1);

for i=1:Number_of_Data/log2(M)
    String_BPSK{i}=num2str(BPSK_Symbols(i,:));
    String_BFSK{i}=num2str(BFSK_Symbols(i,:));
    
    BPSK_Dec(i,1) = bin2dec(String_BPSK{i});
    BFSK_Dec(i,1) = bin2dec(String_BFSK{i});
    
    if(BPSK_Dec(i,1))
        Mapped_BPSK(i) = 1;
    else
        Mapped_BPSK(i) = -1;
    end
    if(BFSK_Dec(i,1))
        Mapped_BFSK(i) = S2;
    else
        Mapped_BFSK(i) = S1;
    end
    
end
%% QPSK
M = 4;
String_QPSK = cell(Number_of_Data/log2(M),1);
QPSK_Dec = zeros(Number_of_Data/log2(M) ,1);
Mapped_QPSK = zeros(Number_of_Data/log2(M) ,1);
Mapped_QPSK_notGray = zeros(Number_of_Data/log2(M),1);
for i=1:Number_of_Data/log2(M)
    String_QPSK{i}=num2str(QPSK_Symbols(i,:));
    QPSK_Dec(i,1) = bin2dec(String_QPSK{i});
    if(QPSK_Dec(i,1)==0)
        Mapped_QPSK(i) = complex(-1,-1);
        Mapped_QPSK_notGray(i)= complex(-1,-1);
    elseif(QPSK_Dec(i,1)==1)
        Mapped_QPSK(i) = complex(-1,1);
        Mapped_QPSK_notGray(i)= complex(-1,1);
    elseif(QPSK_Dec(i,1)==2)
        Mapped_QPSK(i) = complex(1,-1);
        Mapped_QPSK_notGray(i)= complex(1,1);
    else
        Mapped_QPSK(i) = complex(1,1);
        Mapped_QPSK_notGray(i)= complex(1,-1);
    end
end
%% 8PSK
M = 8;
String_MPSK = cell(Number_of_Data8/log2(M),1);
MPSK_Dec = zeros(Number_of_Data8/log2(M) ,1);
Mapped_MPSK = zeros(Number_of_Data8/log2(M) ,1);
for i=1:Number_of_Data8/log2(M)
    String_MPSK{i}=num2str(MPSK_Symbols(i,:));
    MPSK_Dec(i,1) = bin2dec(String_MPSK{i});
    if(MPSK_Dec(i,1)==0)
        Mapped_MPSK(i) = complex(1,0);
    elseif(MPSK_Dec(i,1)==1)
        Mapped_MPSK(i) = complex(1/sqrt(2),1/sqrt(2));
    elseif(MPSK_Dec(i,1)==2)
        Mapped_MPSK(i) = complex(-1/sqrt(2),1/sqrt(2));
    elseif(MPSK_Dec(i,1)==3)
        Mapped_MPSK(i) = complex(0,1);
    elseif(MPSK_Dec(i,1)==4)
        Mapped_MPSK(i) = complex(1/sqrt(2),-1/sqrt(2));
    elseif(MPSK_Dec(i,1)==5)
        Mapped_MPSK(i) = complex(0,-1);
    elseif(MPSK_Dec(i,1)==6)
        Mapped_MPSK(i) = complex(-1,0);
    else
        Mapped_MPSK(i) = complex(-1/sqrt(2),-1/sqrt(2));
    end
end
%% QAM Mapped_Symbols
QAM_Table = sqrt(0.4)*[complex(-3,-3),complex(-3,-1),complex(-3,3),complex(-3,1),complex(-1,-3),complex(-1,-1),complex(-1,3),complex(-1,1),complex(3,-3),complex(3,-1),complex(3,3),complex(3,1),complex(1,-3),complex(1,-1),complex(1,3),complex(1,1)];
M = 16;
String_QAM = cell(Number_of_Data/log2(M),1);
QAM_Dec = zeros(Number_of_Data/log2(M) ,1);
Mapped_QAM = zeros(Number_of_Data/log2(M) ,1);
SumEnergy = 0;
for i=1:Number_of_Data/log2(M)
    String_QAM{i}=num2str(QAM_Symbols(i,:));
    QAM_Dec(i,1) = bin2dec(String_QAM{i});
    if(QAM_Dec(i,1)==0)
        Mapped_QAM(i) = QAM_Table(1);
    elseif(QAM_Dec(i,1)==1)
        Mapped_QAM(i) = QAM_Table(2);
    elseif(QAM_Dec(i,1)==2)
        Mapped_QAM(i) = QAM_Table(3);
    elseif(QAM_Dec(i,1)==3)
        Mapped_QAM(i) = QAM_Table(4);
    elseif(QAM_Dec(i,1)==4)
        Mapped_QAM(i) = QAM_Table(5);
    elseif(QAM_Dec(i,1)==5)
        Mapped_QAM(i) = QAM_Table(6);
    elseif(QAM_Dec(i,1)==6)
        Mapped_QAM(i) = QAM_Table(7);
    elseif(QAM_Dec(i,1)==7)
        Mapped_QAM(i) = QAM_Table(8);
    elseif(QAM_Dec(i,1)==8)
        Mapped_QAM(i) = QAM_Table(9);
    elseif(QAM_Dec(i,1)==9)
        Mapped_QAM(i) = QAM_Table(10);
    elseif(QAM_Dec(i,1)==10)
        Mapped_QAM(i) = QAM_Table(11);
    elseif(QAM_Dec(i,1)==11)
        Mapped_QAM(i) = QAM_Table(12);
    elseif(QAM_Dec(i,1)==12)
        Mapped_QAM(i) = QAM_Table(13);
    elseif(QAM_Dec(i,1)==13)
        Mapped_QAM(i) = QAM_Table(14);
    elseif(QAM_Dec(i,1)==14)
        Mapped_QAM(i) = QAM_Table(15);
    else
        Mapped_QAM(i) = QAM_Table(16);
    end
    SumEnergy = SumEnergy+abs(Mapped_QAM(i))^2;
end
Eavg_QAM = SumEnergy/(Number_of_Data/log2(M));
%% Theoritical BER
% Define the range of Eb/No in dB
    Eb_No_dB = -4:1:14;
    % Convert Eb/No to linear scale
    Eb_No = 10.^(Eb_No_dB / 10);
    No_range=1./Eb_No;
    BER_Theoretical_BPSK=0.5*erfc(sqrt(Eb_No));
    BER_Theoretical_QPSK=0.5*erfc(sqrt(Eb_No));
    BER_Theoretical_8PSK=(1/3)*erfc(sin(pi/8)*sqrt(3*Eb_No));
    BER_Theoretical_16QAM=0.375*erfc(sqrt(Eb_No/2.5));
%% Adding Noise for BPSK
Demapped_BPSK = zeros(Number_of_Data,1);

X_I_BPSK=zeros(Number_of_Data,1);
X_Q_BPSK=zeros(Number_of_Data,1);


BPSK_Noisy_Mapped=zeros(Number_of_Data,1);

BER_BPSK_Actual=zeros(size(Eb_No,2),1);

Bit_Error_BPSK =zeros(size(Eb_No_dB,2),1);

for i=1:size(Eb_No_dB,2)
  for m=1:Number_of_Data
        X_I_BPSK(m)=real(Mapped_BPSK(m));
        X_Q_BPSK(m)=imag(Mapped_BPSK(m));
        
  end
   
        %Adding AWGN channel Noise to IN-Phase Component and Quadrature Component
        X_I_BPSK=X_I_BPSK+ sqrt((No_range(i)/2)).*randn(Number_of_Data,1);% AWGN NOISE
        X_Q_BPSK=X_Q_BPSK+ sqrt((No_range(i)/2)).*randn(Number_of_Data,1);% AWGN NOISE;
            
        BPSK_Noisy_Mapped=complex(X_I_BPSK, X_Q_BPSK);
%% BPSK Demapper
 for m=1:Number_of_Data      
    if(BPSK_Noisy_Mapped(m,1)>0)
    Demapped_BPSK(m,1) = 1;
    else
    Demapped_BPSK(m,1) = 0;
    end
 end
  
    Bit_Error_BPSK(i,1)=biterr(Demapped_BPSK,BPSK_Symbols);
    BER_BPSK_Actual(i)=Bit_Error_BPSK(i)./Number_of_Data;
    
    if(BER_BPSK_Actual(i)~=0)
        lowlim=i;
   end
end
%% Adding Noise to QPSK 
M=4;
Demapped_QPSK = zeros(Number_of_Data/2,1);
Demapped_QPSK_notGray = zeros(Number_of_Data/2,1);

Demapped_QPSK_in_Bits = zeros(Number_of_Data/2,log2(M));
Demapped_QPSK_notGray_in_Bits = zeros(Number_of_Data/2,log2(M));

X_I_QPSK=zeros(Number_of_Data/2,1);
X_Q_QPSK=zeros(Number_of_Data/2,1);

X_I_QPSK_notGray=zeros(Number_of_Data/2,1);
X_Q_QPSK_notGray=zeros(Number_of_Data/2,1);

QPSK_Noisy_Mapped=zeros(Number_of_Data/2,1);
QPSK_notGray_Noisy_Mapped=zeros(Number_of_Data/2,1);

BER_QPSK_Actual=zeros(size(Eb_No,2),1);
BER_QPSK_notGray_Actual=zeros(size(Eb_No,2),1);

Bit_Error_QPSK =zeros(size(Eb_No_dB,2),1);
Bit_Error_QPSK_notGray =zeros(size(Eb_No_dB,2),1);

for i=1:size(Eb_No_dB,2)
    for j=1:Number_of_Data/2
        X_I_QPSK(j)=real(Mapped_QPSK(j));
        X_Q_QPSK(j)=imag(Mapped_QPSK(j));
        
        X_I_QPSK_notGray(j)=real(Mapped_QPSK_notGray(j));
        X_Q_QPSK_notGray(j)=imag(Mapped_QPSK_notGray(j));
    end
    E_avg=sum(X_I_QPSK.^2+X_Q_QPSK.^2)/(Number_of_Data/2);

     X_I_QPSK=X_I_QPSK+sqrt((No_range(i))*E_avg/4).*randn(Number_of_Data/2,1);
     X_Q_QPSK=X_Q_QPSK+sqrt((No_range(i))*E_avg/4).*randn(Number_of_Data/2,1);
     
     X_I_QPSK_notGray=X_I_QPSK_notGray+sqrt((No_range(i))*E_avg/4).*randn(Number_of_Data/2,1);
     X_Q_QPSK_notGray=X_Q_QPSK_notGray+sqrt((No_range(i))*E_avg/4).*randn(Number_of_Data/2,1);
     
     QPSK_Noisy_Mapped=complex(X_I_QPSK,X_Q_QPSK);
     QPSK_notGray_Noisy_Mapped=complex(X_I_QPSK_notGray,X_Q_QPSK_notGray);
%% QPSK Demapper
for m=1:Number_of_Data/2
if(real(QPSK_Noisy_Mapped(m,1))>0 && imag(QPSK_Noisy_Mapped(m,1))>0)
    Demapped_QPSK(m,1) = 3;
elseif(real(QPSK_Noisy_Mapped(m,1))>0 && imag(QPSK_Noisy_Mapped(m,1))<0)
    Demapped_QPSK(m,1) = 2;
elseif(real(QPSK_Noisy_Mapped(m,1))<0 && imag(QPSK_Noisy_Mapped(m,1))>0)
    Demapped_QPSK(m,1) = 1;
else
    Demapped_QPSK(m,1) = 0;
end
    Demapped_QPSK_in_Bits(m,:) = dec2bin(Demapped_QPSK(m),log2(M));

if(real(QPSK_notGray_Noisy_Mapped(m,1))>0 && imag(QPSK_notGray_Noisy_Mapped(m,1))>0)
     Demapped_QPSK_notGray(m,1)= 2;
elseif(real(QPSK_notGray_Noisy_Mapped(m,1))>0 && imag(QPSK_notGray_Noisy_Mapped(m,1))<0)
    Demapped_QPSK_notGray(m,1)= 3;
elseif(real(QPSK_notGray_Noisy_Mapped(m,1))<0 && imag(QPSK_notGray_Noisy_Mapped(m,1))>0)
        Demapped_QPSK_notGray(m,1)= 1;
else
      Demapped_QPSK_notGray(m,1)= 0;
end  
    
     Demapped_QPSK_notGray_in_Bits(m,:) = dec2bin(Demapped_QPSK_notGray(m),log2(M));
end
    Demapped_QPSK_in_Bits = Demapped_QPSK_in_Bits-48;
    Demapped_QPSK_notGray_in_Bits = Demapped_QPSK_notGray_in_Bits-48;
    
    Bit_Error_QPSK(i,1)=biterr(Demapped_QPSK_in_Bits,QPSK_Symbols);
    Bit_Error_QPSK_notGray(i,1)=biterr(Demapped_QPSK_notGray_in_Bits,QPSK_Symbols);
    
    BER_QPSK_Actual(i)=Bit_Error_QPSK(i)./Number_of_Data;
    BER_QPSK_notGray_Actual(i)=Bit_Error_QPSK_notGray(i)./Number_of_Data;
    if(BER_QPSK_Actual(i)~=0)
        lowlim=i;
    end
     if(BER_QPSK_notGray_Actual(i)~=0)
        lowlim2=i;
     end
end
%% ِAdding Noise fore 8PSK 
M=8;
Demapped_MPSK = zeros(Number_of_Data8/log2(M),1);
Demapped_MPSK_in_Bits = zeros(Number_of_Data8/log2(M),log2(M));
X_I_MPSK=zeros(Number_of_Data8/log2(M),1);
X_Q_MPSK=zeros(Number_of_Data8/log2(M),1);
MPSK_Noisy_Mapped=zeros(Number_of_Data8/log2(M),1);
BER_MPSK_Actual=zeros(size(Eb_No,2),1);
Bit_Error_MPSK =zeros(size(Eb_No_dB,2),1);
for i=1:size(Eb_No_dB,2)
   for m=1:Number_of_Data8/log2(M)
        X_I_MPSK(m)=real(Mapped_MPSK(m));
        X_Q_MPSK(m)=imag(Mapped_MPSK(m));
        %Adding AWGN channel Noise to IN-Phase Component and Quadrature Component
        X_I_MPSK(m)=X_I_MPSK(m)+ sqrt((No_range(i))/(2*log2(M)))*randn;% AWGN NOISE
        X_Q_MPSK(m)=X_Q_MPSK(m)+ sqrt((No_range(i))/(2*log2(M)))*randn;% AWGN NOISE;
        MPSK_Noisy_Mapped(m)=complex(X_I_MPSK(m), X_Q_MPSK(m));
%% 8PSK Demapper
    if(angle(MPSK_Noisy_Mapped(m,1))>-(1/8)*pi && angle(MPSK_Noisy_Mapped(m,1))<(1/8)*pi)
        Demapped_MPSK(m,1) = 0;
    elseif(angle(MPSK_Noisy_Mapped(m,1))>(1/8)*pi && angle(MPSK_Noisy_Mapped(m,1))<(3/8)*pi)
        Demapped_MPSK(m,1) = 1;
    elseif(angle(MPSK_Noisy_Mapped(m,1))>(3/8)*pi && angle(MPSK_Noisy_Mapped(m,1))<(5/8)*pi)
        Demapped_MPSK(m,1) = 3;
    elseif(angle(MPSK_Noisy_Mapped(m,1))>(5/8)*pi && angle(MPSK_Noisy_Mapped(m,1))<(7/8)*pi)
        Demapped_MPSK(m,1) = 2;
    elseif(angle(MPSK_Noisy_Mapped(m,1))>(7/8)*pi || angle(MPSK_Noisy_Mapped(m,1))<-(7/8)*pi)
        Demapped_MPSK(m,1) = 6;
    elseif(angle(MPSK_Noisy_Mapped(m,1))>-(7/8)*pi && angle(MPSK_Noisy_Mapped(m,1))<-(5/8)*pi)
        Demapped_MPSK(m,1) = 7;
    elseif(angle(MPSK_Noisy_Mapped(m,1))>-(5/8)*pi && angle(MPSK_Noisy_Mapped(m,1))<-(3/8)*pi)
        Demapped_MPSK(m,1) = 5;
    else
        Demapped_MPSK(m,1) = 4;
    end
    Demapped_MPSK_in_Bits(m,:) = dec2bin(Demapped_MPSK(m),log2(M));
    end
    Demapped_MPSK_in_Bits = Demapped_MPSK_in_Bits-48;
    Bit_Error_MPSK(i,1)=biterr(Demapped_MPSK_in_Bits,MPSK_Symbols);
    BER_MPSK_Actual(i)=Bit_Error_MPSK(i)./Number_of_Data8;
    if(BER_MPSK_Actual(i)~=0)
        lowlim=i;
    end
end
%% Assing Noise to 16QAM
M=16;
X_I=zeros(Number_of_Data/4,1);
X_Q=zeros(Number_of_Data/4,1);
QAM_Noisy_Mapped=zeros(Number_of_Data,1);
Demapped_QAM = zeros(Number_of_Data/4,1);
Demapped_QAM_in_Bits = zeros(Number_of_Data/4,4);
BER_QAM_Actual=zeros(size(Eb_No,2),1);
Bit_Error_QAM =zeros(size(Eb_No_dB,2),1);
for i=1:size(Eb_No_dB,2)
 for m = 1:Number_of_Data/log2(M)
        X_I(m)=real(Mapped_QAM(m));
        X_Q(m)=imag(Mapped_QAM(m));
        %Adding AWGN channel Noise to IN-Phase Component and Quadrature Component
        X_I(m)=X_I(m)+ sqrt((No_range(i)/2)*Eavg_QAM/log2(M))*randn;% AWGN NOISE
        X_Q(m)=X_Q(m)+ sqrt((No_range(i)/2)*Eavg_QAM/log2(M))*randn;% AWGN NOISE;
        QAM_Noisy_Mapped(m)=complex(X_I(m), X_Q(m));
%% 16QAM Demapper
        [MIN,Index] = min(abs(QAM_Table-QAM_Noisy_Mapped(m,1)));
        Demapped_QAM(m) = Index-1;
        Demapped_QAM_in_Bits(m,:) = dec2bin(Demapped_QAM(m),4);
end
    Demapped_QAM_in_Bits = Demapped_QAM_in_Bits-48;
    Bit_Error_QAM(i,1)= biterr(Demapped_QAM_in_Bits,QAM_Symbols);
    BER_QAM_Actual(i)=Bit_Error_QAM(i)./Number_of_Data;
    if(BER_QAM_Actual(i)~=0)
        lowlim=i;
    end
end
%% BFSK
% Adding Noise for BFSK
BFSK_Table=[S1,S2];
Demapped_BFSK = zeros(Number_of_Data,1);

X_I_BFSK=zeros(Number_of_Data,1);
X_Q_BFSK=zeros(Number_of_Data,1);

BFSK_Noisy_Mapped=zeros(Number_of_Data,1);

BER_BFSK_Actual=zeros(size(Eb_No,2),1);

Bit_Error_BFSK =zeros(size(Eb_No_dB,2),1);

for i=1:size(Eb_No_dB,2)
  for m=1:Number_of_Data

        X_I_BFSK(m)=real(Mapped_BFSK(m));
        X_Q_BFSK(m)=imag(Mapped_BFSK(m));
  end
   
        %Adding AWGN channel Noise to IN-Phase Component and Quadrature Component
        
        X_I_BFSK=X_I_BFSK+ sqrt(No_range(i)).*randn(Number_of_Data,1);% AWGN NOISE
        X_Q_BFSK=X_Q_BFSK+ sqrt(No_range(i)).*randn(Number_of_Data,1);% AWGN NOISE;
        
        BFSK_Noisy_Mapped=complex(X_I_BFSK, X_Q_BFSK);
 for m=1:Number_of_Data      

     [MIN,Index] = min(abs(BFSK_Table-BFSK_Noisy_Mapped(m,1)));
     Demapped_BFSK(m) = Index-1;
 end

    Bit_Error_BFSK(i,1)=biterr(Demapped_BFSK,BFSK_Symbols);
    BER_BFSK_Actual(i)=Bit_Error_BFSK(i)./Number_of_Data;

    if(BER_BFSK_Actual(i)~=0)
        lowlim3=i;
    end
end
%% Figures
figure;
semilogy(Eb_No_dB, BER_Theoretical_BPSK,'--','LineWidth', 2,'color','red');
xlabel("Eb/N0 dB");
ylabel("BER");
title("BER vs EB/N0");
hold on 
semilogy(Eb_No_dB, BER_Theoretical_QPSK,'--','LineWidth', 2,'color','yellow');
hold on 
semilogy(Eb_No_dB, BER_Theoretical_8PSK,'--','LineWidth', 2,'color','green');
hold on
semilogy(Eb_No_dB, BER_Theoretical_16QAM,'--','LineWidth', 2,'color','blue');
    hold on
    semilogy(Eb_No_dB,transpose(BER_BPSK_Actual),'LineWidth', 2,'color','red');
    ylim([BER_BPSK_Actual(lowlim) 1])
    hold on
    semilogy(Eb_No_dB,transpose(BER_QPSK_Actual),'LineWidth', 2,'color','yellow');
    ylim([BER_QPSK_Actual(lowlim) 1])
    hold on
    semilogy(Eb_No_dB,transpose(BER_MPSK_Actual),'LineWidth', 2,'color','green');
    ylim([BER_MPSK_Actual(lowlim) 1])
    hold on
    semilogy(Eb_No_dB,transpose(BER_QAM_Actual),'LineWidth', 2,'color','blue');
    ylim([BER_QAM_Actual(lowlim) 1])
legend('BPSK Theoretical','QPSK Theoretical','8PSK Theoretical','16QAM Theoretical',"Actual BPSK","Actual QPSK","Actual 8PSK","Actual 16QAM");
%% QPSK Gray vs QPSK not Gray
figure;
semilogy(Eb_No_dB,transpose(BER_QPSK_Actual),'LineWidth',2,'color','red');
 ylim([BER_QPSK_Actual(lowlim) 1])
xlabel("Eb/N0 dB");
ylabel("BER");
title("BER vs EB/N0");
hold on
semilogy(Eb_No_dB,transpose(BER_QPSK_notGray_Actual),'LineWidth',2,'color','blue');
 ylim([BER_QPSK_notGray_Actual(lowlim2) 1])
hold on
semilogy(Eb_No_dB, BER_Theoretical_QPSK,'--','LineWidth', 2,'color','black');
legend("QPSK Gray","QPSK Not Gray","Theoretical QPSK");
%% BFSK BER vs EB/NO PLOT
figure;
BFSK_Theoretical=0.5*erfc(sqrt(Eb_No/2));
semilogy(Eb_No_dB, BFSK_Theoretical,'--','LineWidth', 2);
xlabel("Eb/N0 dB");
ylabel("BER");
title("BER vs EB/N0");
hold on
semilogy(Eb_No_dB,transpose(BER_BFSK_Actual),'LineWidth',2);
ylim([BER_BFSK_Actual(lowlim3) 1])
legend("BFSK Theoretical","BFSK Actual");
%% BFSK PSD
Tb=1; 
Eb=1;
delta_f=1/Tb;
Number_Realizations=1000; 
Number_Bits=100; 
SamplesPerBit=8; 
t=(0:Tb/SamplesPerBit:Tb-Tb/SamplesPerBit); 
% BaseBand equivalent Signal to 1 
s1=sqrt(2*Eb/Tb);  
% BaseBand equivalent Signal to 0 
s2=s1*exp(-2*pi*(delta_f)*t*1i); 
s2=s2'; 
% Ensemble The BFSK after mapping 
Ensemble_BFSK=zeros(Number_Realizations,Number_Bits*SamplesPerBit); 

for i=1:Number_Realizations 
data=randi([0 1],1,Number_Bits); 
Tx=repmat(data,SamplesPerBit,1); 
for j=1:Number_Bits
   if data(j)==1 
       Tx(:,j)=s1; 
    else  
       Tx(:,j)=s2; 
  end 
end
% Random initial delay 
random_start=round(SamplesPerBit*rand); 
number=randi([0,1]); 
    if number==1 
       start=repmat(s1,random_start,1); 
    else 
       start=s2(SamplesPerBit-random_start+1:end,1); 
    end 
Tx_out=reshape(Tx,size(Tx,1)* size(Tx,2),1); 
Tx_rand=[start;Tx_out;zeros(SamplesPerBit-length(start),1)]; 
Tx_rand_tr=transpose(Tx_rand);
Ensemble_BFSK(i,:)= Tx_rand_tr(1:800); 
end 
% Compute Ensemble Autocorrelation Function 
Ensemble_AutoCorr(Number_Bits*SamplesPerBit,Number_Bits*SamplesPerBit) = 0; 
for i = 1:(Number_Bits*SamplesPerBit) 
for j = 1:(Number_Bits*SamplesPerBit) 
Ensemble_AutoCorr(i,j) = sum(conj(Ensemble_BFSK(:,i)).*Ensemble_BFSK(:,j),1); 
end 
end 
EnsembleAutoCorr_sided = [conj(fliplr(Ensemble_AutoCorr(1,2:end))) Ensemble_AutoCorr(1,:)]; 
% plotting the Auto Correlation Function 
figure;
x_ACF = -(SamplesPerBit*Number_Bits-1):(SamplesPerBit*Number_Bits-1); 
plot(x_ACF,abs(EnsembleAutoCorr_sided)) 
title('Auto Correlation Function across Realizations') 
% Compute Power Spectral Density 
Ensemble_PSD = (1/length(x_ACF))*(abs(fftshift(fft(EnsembleAutoCorr_sided)))); 
% Plotting the PSD 
fs=SamplesPerBit/(Tb); 
n=(SamplesPerBit*Number_Bits+SamplesPerBit); 
f=0.5*((-n+1:n-1)*fs/n); 
f_lim=f(1,1:1599);
figure;
plot(f_lim,Ensemble_PSD);
ylim([0, 2])
xlabel("Normalized Frequency fTb");
ylabel("PSD");
title ('PSD of BFSK')