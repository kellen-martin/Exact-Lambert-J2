function theta = arcsin(x,y)
    theta = pi-pi/2*(1+sign(x))*(1-sign(y^2))-pi/4*(2+sign(x))*sign(y) -sign(x*y)*asin((abs(x)-abs(y))/sqrt(2*x^2+2*y^2));

end