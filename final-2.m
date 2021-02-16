clear all
clc
%   step1  設定初始條件
T = 2400;
Rbar = 8.314;

% a矩陣 各元素在不同物質中之莫爾數
% 物質們 O2, CO, CO2, H2O, H, H2, OH, N2, NO
a = [0, 1, 1, 0, 0, 0, 0, 0, 0;
     2, 1, 2, 1, 0, 0, 1, 0, 1;
     0, 0, 0, 2, 1, 2, 1, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 2, 1];
 % b矩陣 各元素之莫爾數
 b = [1, 3.8, 4, 14.288];
 
 % 各物質之基本性質 其中dh=h2400-h298
h0 = [0, -110541, -393546, -241845, 217977, 0, 38985, 0, 90297];
dh = [74467, 71354, 115798, 93744, 43689, 66925, 67843, 70645, 72617];
s0 = [275.627, 265.278, 320.333, 274.225, 157.956, 194.692, 248.654, 258.598, 279.722];  

% g* = h0 + dh - T*S
for j = 1:9
    gstart(j) = h0(j) + dh(j) - T * s0(j);
end

% 先猜測平衡後各物質之莫爾數
lnn = rand(1,9);
lnN = rand(1);


%   step2  設定需要用的矩陣及方程式們，並開始迭代
for times = 1:100
    g = zeros(1,9);
    % Ｐ為5*5矩陣(4元素+1)
    P = zeros(5);
    R = zeros(5,1);
    % 將猜測之平衡莫爾數取ln
    n1 = exp(lnn);
    N = exp(lnN);
    % gibbs eq 其中P=10atm，Po=1atm故省略
    for j = 1:9
        g(j) = gstart(j) + Rbar * T * log(10*n1(j)/N);
    end
    
    
    %   step3   設定P、R矩陣之係數
    % P矩陣從n,m=1,1到4,4先設定
    for n = 1:4
        for m = 1:4
            for j = 1:9
                P(n,m) = P(n,m) + n1(j) * a(n,j) * a(m,j);
            end
        end
    end
	% P矩陣設定n=5及m=5，除了最右下角之n,m=5,5
    for i = 1:4  
        for j = 1:9
            P(i,5) = P(i,5) + n1(j) * a(i,j);
        end
        P(5,i) = P(i,5);
    end
	% P矩陣設定n,m=5,5
    for j = 1:9   
        P(5,5) = P(5,5) + n1(j);
    end
    P(5,5) = P(5,5) - N;

    % R矩陣前四項
    for i = 1:4
        for j = 1:9
            R(i) = R(i) + n1(j) * a(i,j) * (g(j)/Rbar/T - 1);
        end
        R(i) = R(i) + b(i);
    end
    % R矩陣最後一項
    for j = 1:9
        R(5) = R(5) + n1(j) * (g(j)/Rbar/T - 1);
    end
    R(5) = R(5) + N;

  
	%   step4   將P取反矩陣，並同時前乘於PQ=R，算出Q矩陣
    Q = P\R ;

    
    %   step5   求個物質之莫爾數變化量
    % 9個物質的變化量dlnn
    dlnn = zeros(1,9);
    for j = 1:9
        for i = 1:4
            % Q(i)=pi(i) (共4個元素)
            dlnn(j) = dlnn(j) + a(i,j) * Q(i);
        end
        % Q(5)=dlnN
        dlnn(j) = dlnn(j) + Q(5) - g(j)/Rbar/T;
    end

    
	%   step6   取函數迭代之收斂值
    %　取最大誤差值
    % 其中exp(lnn + dlnn)為物質當前莫爾數，n1為初始物質莫爾數 
    % 先處理1-9個物質之n
    Emax = max(abs((exp(lnn + dlnn)-n1)/n1));
    %再處理總莫爾數變化量，於n1-9及N中取一個最大的
    Emax = max(Emax,abs((exp(lnN + Q(5))-N)/N));
    
    %   函數收斂，ｃｈｅｃｋ各元素之莫爾數
    if Emax < 1E-10
        element = zeros(1,4);
        for i = 1:4
            for j = 1:9
                element(i) = element(i) + a(i,j)*n1(j);
            end
        end
        
        %   step7   秀出平衡結果
        fprintf('\t平衡後各元素之莫爾數\n')
        fprintf('\tC         O         H        N\n')
        disp(element);
        fprintf('\t平衡後各物質之莫爾數\n')
        fprintf('\tO2        CO        CO2       H2O       H         H2        OH        N2        NO\n')
        n1 = exp(lnn);
        N = exp(lnN);
        disp(n1);
        fprintf('\t平衡後物質總莫爾數\n')
        disp(N);
        break;
    else
        
        % 若誤差值太大，繼續計算
        lnn = lnn + dlnn;
        lnN = lnN + Q(5);
    end
end
