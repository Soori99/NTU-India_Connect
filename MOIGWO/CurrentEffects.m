function [velxy, wavecurrent] = CurrentEffects(inputx,inputy,option,time,headingangle,USVvelocity)
    Bnaught = 1.2;
    epsilon = 0.3;
    omega = 0.4;
    beta = pi/2;
    k = 0.84;
    c = 0.12;
    velxy = zeros(1,2);
    if option ==1
        wavecurrent = 3;
        velxy(1,1) = USVvelocity*cosd(headingangle) + wavecurrent*cosd(headingangle);
        velxy(1,2) = USVvelocity*sind(headingangle) + wavecurrent*sind(headingangle);
    else
        B = Bnaught + epsilon*cos(omega*time + beta);
        syms X Y
        phi = 1-tanh((Y - B*cosd(k*(X-c*time)))/sqrt(1+k^2*B^2*sind(k*(X-c*time))^2)); 
        dX = matlabFunction(diff(phi,X));
        dY = matlabFunction(diff(phi,Y)); 
        temp1 = -1*dY(inputx/10,inputy/10);
        temp2 = dX(inputx/10,inputy/10);
        wavecurrent = sqrt(temp1^2+temp2^2);
        velxy(1,1) = USVvelocity*cosd(headingangle)+temp1;
        velxy(1,2) = USVvelocity*sind(headingangle)+temp2;
    end 
 end
