clear all
clc
%   step1  �]�w��l����
T = 2400;
Rbar = 8.314;

% a�x�} �U�����b���P���褤��������
% ����� O2, CO, CO2, H2O, H, H2, OH, N2, NO
a = [0, 1, 1, 0, 0, 0, 0, 0, 0;
     2, 1, 2, 1, 0, 0, 1, 0, 1;
     0, 0, 0, 2, 1, 2, 1, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 2, 1];
 % b�x�} �U������������
 b = [1, 3.8, 4, 14.288];
 
 % �U���褧�򥻩ʽ� �䤤dh=h2400-h298
h0 = [0, -110541, -393546, -241845, 217977, 0, 38985, 0, 90297];
dh = [74467, 71354, 115798, 93744, 43689, 66925, 67843, 70645, 72617];
s0 = [275.627, 265.278, 320.333, 274.225, 157.956, 194.692, 248.654, 258.598, 279.722];  

% g* = h0 + dh - T*S
for j = 1:9
    gstart(j) = h0(j) + dh(j) - T * s0(j);
end

% ���q�����ū�U���褧������
lnn = rand(1,9);
lnN = rand(1);


%   step2  �]�w�ݭn�Ϊ��x�}�Τ�{���̡A�ö}�l���N
for times = 1:100
    g = zeros(1,9);
    % �ެ�5*5�x�}(4����+1)
    P = zeros(5);
    R = zeros(5,1);
    % �N�q�������Ų����ƨ�ln
    n1 = exp(lnn);
    N = exp(lnN);
    % gibbs eq �䤤P=10atm�APo=1atm�G�ٲ�
    for j = 1:9
        g(j) = gstart(j) + Rbar * T * log(10*n1(j)/N);
    end
    
    
    %   step3   �]�wP�BR�x�}���Y��
    % P�x�}�qn,m=1,1��4,4���]�w
    for n = 1:4
        for m = 1:4
            for j = 1:9
                P(n,m) = P(n,m) + n1(j) * a(n,j) * a(m,j);
            end
        end
    end
	% P�x�}�]�wn=5��m=5�A���F�̥k�U����n,m=5,5
    for i = 1:4  
        for j = 1:9
            P(i,5) = P(i,5) + n1(j) * a(i,j);
        end
        P(5,i) = P(i,5);
    end
	% P�x�}�]�wn,m=5,5
    for j = 1:9   
        P(5,5) = P(5,5) + n1(j);
    end
    P(5,5) = P(5,5) - N;

    % R�x�}�e�|��
    for i = 1:4
        for j = 1:9
            R(i) = R(i) + n1(j) * a(i,j) * (g(j)/Rbar/T - 1);
        end
        R(i) = R(i) + b(i);
    end
    % R�x�}�̫�@��
    for j = 1:9
        R(5) = R(5) + n1(j) * (g(j)/Rbar/T - 1);
    end
    R(5) = R(5) + N;

  
	%   step4   �NP���ϯx�}�A�æP�ɫe����PQ=R�A��XQ�x�}
    Q = P\R ;

    
    %   step5   �D�Ӫ��褧�������ܤƶq
    % 9�Ӫ��誺�ܤƶqdlnn
    dlnn = zeros(1,9);
    for j = 1:9
        for i = 1:4
            % Q(i)=pi(i) (�@4�Ӥ���)
            dlnn(j) = dlnn(j) + a(i,j) * Q(i);
        end
        % Q(5)=dlnN
        dlnn(j) = dlnn(j) + Q(5) - g(j)/Rbar/T;
    end

    
	%   step6   ����ƭ��N�����ĭ�
    %�@���̤j�~�t��
    % �䤤exp(lnn + dlnn)�������e�����ơAn1����l��������� 
    % ���B�z1-9�Ӫ��褧n
    Emax = max(abs((exp(lnn + dlnn)-n1)/n1));
    %�A�B�z�`�������ܤƶq�A��n1-9��N�����@�ӳ̤j��
    Emax = max(Emax,abs((exp(lnN + Q(5))-N)/N));
    
    %   ��Ʀ��ġA������U������������
    if Emax < 1E-10
        element = zeros(1,4);
        for i = 1:4
            for j = 1:9
                element(i) = element(i) + a(i,j)*n1(j);
            end
        end
        
        %   step7   �q�X���ŵ��G
        fprintf('\t���ū�U������������\n')
        fprintf('\tC         O         H        N\n')
        disp(element);
        fprintf('\t���ū�U���褧������\n')
        fprintf('\tO2        CO        CO2       H2O       H         H2        OH        N2        NO\n')
        n1 = exp(lnn);
        N = exp(lnN);
        disp(n1);
        fprintf('\t���ū᪫���`������\n')
        disp(N);
        break;
    else
        
        % �Y�~�t�ȤӤj�A�~��p��
        lnn = lnn + dlnn;
        lnN = lnN + Q(5);
    end
end
