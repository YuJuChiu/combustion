clear all
clc
% step1  initialize the conditions
T = 2400;
Rbar = 8.314;

% matrix a, mole of each substance 
% O2, CO, CO2, H2O, H, H2, OH, N2, NO
a = [0, 1, 1, 0, 0, 0, 0, 0, 0;
     2, 1, 2, 1, 0, 0, 1, 0, 1;
     0, 0, 0, 2, 1, 2, 1, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 2, 1];
 % matrix b, mole of each element 
 b = [1, 3.8, 4, 14.288];
 
 % properties of each substance, dh=h2400-h298
h0 = [0, -110541, -393546, -241845, 217977, 0, 38985, 0, 90297];
dh = [74467, 71354, 115798, 93744, 43689, 66925, 67843, 70645, 72617];
s0 = [275.627, 265.278, 320.333, 274.225, 157.956, 194.692, 248.654, 258.598, 279.722];  

% g* = h0 + dh - T*S
for j = 1:9
    gstart(j) = h0(j) + dh(j) - T * s0(j);
end

% guess the mole of each substance after balancing
lnn = rand(1,9);
lnN = rand(1);


% step2  initialize the neccessary martices, start to iterate
for times = 1:100
    g = zeros(1,9);
    % Ｐ = R^(5*5), (4 elements + 1)
    P = zeros(5);
    R = zeros(5,1);
    % take logarithm of the mole
    n1 = exp(lnn);
    N = exp(lnN);
    % gibbs eq 
    for j = 1:9
        g(j) = gstart(j) + Rbar * T * log(10*n1(j)/N);
    end
    
    
    % step3   initialize P, R matrix 
    for n = 1:4
        for m = 1:4
            for j = 1:9
                P(n,m) = P(n,m) + n1(j) * a(n,j) * a(m,j);
            end
        end
    end
    
    for i = 1:4  
        for j = 1:9
            P(i,5) = P(i,5) + n1(j) * a(i,j);
        end
        P(5,i) = P(i,5);
    end
    
    for j = 1:9   
        P(5,5) = P(5,5) + n1(j);
    end
    P(5,5) = P(5,5) - N;

    % R
    for i = 1:4
        for j = 1:9
            R(i) = R(i) + n1(j) * a(i,j) * (g(j)/Rbar/T - 1);
        end
        R(i) = R(i) + b(i);
    end

    for j = 1:9
        R(5) = R(5) + n1(j) * (g(j)/Rbar/T - 1);
    end
    R(5) = R(5) + N;

  
    % step4 
    Q = P\R ;

    
    % step5   the difference of mole of each substance
    dlnn = zeros(1,9);
    for j = 1:9
        for i = 1:4
            % Q(i)=pi(i)
            dlnn(j) = dlnn(j) + a(i,j) * Q(i);
        end
        % Q(5)=dlnN
        dlnn(j) = dlnn(j) + Q(5) - g(j)/Rbar/T;
    end

    
    % step6   take the value of convergence
    Emax = max(abs((exp(lnn + dlnn)-n1)/n1));

    Emax = max(Emax,abs((exp(lnN + Q(5))-N)/N));
    
    % converge
    if Emax < 1E-10
        element = zeros(1,4);
        for i = 1:4
            for j = 1:9
                element(i) = element(i) + a(i,j)*n1(j);
            end
        end
        
        % step7   print the result
        fprintf('\tmole of each element after balancing\n')
        fprintf('\tC         O         H        N\n')
        disp(element);
        fprintf('\tmole of each substance after balancing\n')
        fprintf('\tO2        CO        CO2       H2O       H         H2        OH        N2        NO\n')
        n1 = exp(lnn);
        N = exp(lnN);
        disp(n1);
        fprintf('\ttotal mole after balancing\n')
        disp(N);
        break;
    else
        
        % if the error is too too large, keep calculating
        lnn = lnn + dlnn;
        lnN = lnN + Q(5);
    end
end
