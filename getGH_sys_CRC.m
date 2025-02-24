% Convert a generator polynomial into [G,H] binary code matrices
% Requires MATLAB toolbox

function [G,H] = getGH_sys_CRC(n,k)

        % Polynomial from https://users.ece.cmu.edu/~koopman/crc/
        if n-k==3
            hex_poly = '0x5'; % 3 bit HD3 max length 4 otherwise HD2
        elseif n-k==4
            hex_poly = '0x9'; % 4 bit HD3 max length 11, otherwise HD2
        elseif n-k==5 && k<=10
            hex_poly = '0x15'; % 5 bit HD4 max length 10
        elseif n-k==5 && k<=26
            hex_poly = '0x12'; % 5 bit HD3 max length 26
        elseif n-k==6 && k<=25
             hex_poly = '0x23'; % 6 bit HD4 max length 25
        elseif n-k==6 
             hex_poly = '0x33'; % 6 bit HD3 max length 57 or HD2 thereafter
        elseif n-k==7 && k<=4
             hex_poly = '0x72'; % 7 bit HD5 max length 4
        elseif n-k==7 && k<=56
             hex_poly = '0x5b'; % 7 bit HD4 max length 56
        elseif n-k==7 && k<=120
             hex_poly = '0x65'; % 7 bit HD4 max length 120
        elseif n-k==8 && k<=4
            hex_poly = '0x9b'; % 8 bit HD6 max length 4
        elseif n-k==8 && k<=9
            hex_poly = '0xeb'; % 8 bit HD5 max length 9
        elseif n-k==8 && k<=119
            hex_poly = '0x83'; % 8 bit HD4 max length 119
        elseif n-k==8 && k<=247
            hex_poly = '0xe7'; % 8 bit HD3 max length 247
        elseif n-k==9 && k<=8
            hex_poly = '0x13c'; % 9 bit HD6 max length 8
        elseif n-k==9 && k<=13
            hex_poly = '0x185'; % 9 bit HD5 max length 13
        elseif n-k==9 && k<=246
            hex_poly = '0x17d'; % 9 bit HD4 max length 246
        elseif n-k==10 && k<=5
            hex_poly = '0x29b'; % 10 bit HD7 max length 5
        elseif n-k==10 && k<=12
            hex_poly = '0x28e'; % 10 bit HD6 max length 12
        elseif n-k==10 && k<=26
            hex_poly = '0x2b9'; % 10 bit HD5 max length 26
        elseif n-k==10 && k<=501
            hex_poly = '0x247'; % 10 bit HD4 max length 501
        elseif n-k==11 && k<=4
            hex_poly = '0x4f5'; % 11 bit HD8 max length 4
        elseif n-k==11 && k<=12
            hex_poly = '0x571'; % 11 bit HD7 max length 12
        elseif n-k==11 && k<=22
            hex_poly = '0x532'; % 11 bit HD6 max length 22
        elseif n-k==11 && k<=26
            hex_poly = '0x5d7'; % 11 bit HD5 max length 26
        elseif n-k==11 && k<=1012
            hex_poly = '0x583'; % 11 bit HD4 max length 1012
        elseif n-k==12 && k<=11
            hex_poly = '0xa4f'; % 12 bit HD7/8 max length 11
        elseif n-k==12 && k<=27
            hex_poly = '0xb41'; % 12 bit HD6 max length 27
        elseif n-k==12 && k<=53
            hex_poly = '0xbae'; % 12 bit HD5 max length 53
        elseif n-k==12 && k<=2035
            hex_poly = '0x8f3'; % 12 bit HD4 max length 2035
        elseif n-k==13 && k<=11
            hex_poly = '0x10b7'; % 13 bit HD8 max length 11
        elseif n-k==13 && k<=12
            hex_poly = '0x12a5'; % 13 bit HD7 max length 12
        elseif n-k==13 && k<=52
            hex_poly = '0x1e97'; % 13 bit HD6 max length 52
        elseif n-k==14 && k<=11
            hex_poly = '0x2371'; % 14 bit HD8 max length 11  
        elseif n-k==14 && k<=13
            hex_poly = '0x28a9'; % 14 bit HD7 max length 13        
        elseif n-k==14 && k<=57
            hex_poly = '0x372b'; % 14 bit HD6 max length 57
        elseif n-k==14 && k<=113
            hex_poly = '0x212d'; % 14 bit HD5 max length 113
        elseif n-k==15 && k<=5
            hex_poly = '0x5a47'; % 15 bit HD9 max length 5
        elseif n-k==15 && k<=114
            hex_poly = '0x573a'; % 15 bit HD6 max length 114
        elseif n-k==15 && k<=136
            hex_poly = '0x6a8d'; % 15 bit HD5 max length 136
        elseif n-k==15 && k<=16368
            hex_poly = '0x4306'; % 15 bit HD4 max length 16368
        elseif n-k==16 && k<=5
            hex_poly = '0xed2f'; % 16 bit HD10 max length 5
        elseif n-k==16 && k<=19
            hex_poly = '0x968b'; % 16 bit HD6 max length 19
        elseif n-k==16 && k<=135
            hex_poly = '0x9eb2'; % 16 bit HD6 max length 135
        elseif n-k==16 && k<=241
            hex_poly = '0xac9a'; % 16 bit HD5 max length 241
        elseif n-k==16 && k<=32751
            hex_poly = '0xd175'; % 16 bit HD4 max length 32751
        elseif n-k==18 && k<=5
            hex_poly = '0x26a3d'; % 18 bit HD11 max length 5
        elseif n-k==18 && k<=8
            hex_poly = '0x2e7de'; % 18 bit HD10 max length 8
        elseif n-k==18 && k<=240
            hex_poly = '0x32c69'; % 18 bit HD6 max length 240
        elseif n-k==19 && k<=7
            hex_poly = '0x6d133'; % 19 bit HD11 max length 7
        elseif n-k==19 && k<=10
            hex_poly = '0x44f75'; % 19 bit HD10 max length 10
        elseif n-k==19 && k<=494
            hex_poly = '0x5685a'; % 19 bit HD6 max length 494
        elseif n-k==20 && k<=11
            hex_poly = '0x8d3cc'; % 20 bit HD11 max length 11
        elseif n-k==21 && k<=10
            hex_poly = '0x165751'; % 21 bit HD12 max length 10 
        elseif n-k==21 && k<=106
            hex_poly = '0x12faa5'; % 21 bit HD7 max length 106 
        elseif n-k==22 && k<=5
            hex_poly = '0x25d467'; % 22 bit HD13 max length 5
        elseif n-k==22 && k<=105
            hex_poly = '0x289cfe'; % 22 bit HD8 max length 105
        elseif n-k==23 && k<=105
            hex_poly = '0x469d7c'; % 23 bit HD8 max length 105
        elseif n-k==29 && k<=100
            hex_poly = '0x1e150a87'; % 29 bit HD9 max length 100
        end
        poly=koopman2matlab(hex_poly); % Convert from Koopman notation.
        % Set up the CRC check using MATLAB tool.
    GCRC = comm.CRCGenerator(poly);
    n = length(GCRC(zeros(k,1)));
    
    u=[zeros(1,k-1) 1];
    G = zeros(k,n);
    for ii=1:k
        u=circshift(u,1);
        G(ii,:) = GCRC(u')';
    end
    H =[G(:,k+1:end)' eye(n-k)];
end

% Convert a polynomial in Koopman notation into one suitable for MATLAB's
% comms package
function poly = koopman2matlab(k_poly)
    poly= dec2bin(hex2dec(k_poly))-'0';
    poly = [poly 1];
end
