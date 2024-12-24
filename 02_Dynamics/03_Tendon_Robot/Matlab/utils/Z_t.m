function displace=Z_t(t)
    zA = 0;          %  z(m)  Ramp Motor Input:
    zB = -0.025;     % zA ^ ____           ______
    t1 = 1;          %    |     \         /
    t2 = 2;          % zB |      \_______/
    t3 = 3;          %    -------------------------> t(s)
    t4 = 4;          %         t1 t2    t3 t4
    if t>t1 && t<=t2
        displace=zA + (zB-zA)*(t-t1)/(t2-t1); %Ramp lower
    elseif t>t2 && t<=t3
        displace=zB;                         %Low state
    elseif t>t3 && t<=t4
        displace=zB + (zA-zB)*(t-t3)/(t4-t3); %Ramp higher
    else
        displace=zA;  %High state
    end
end