function f=sourceTerm(x)

E=exp(-100*(x-0.5).^2);
f=400*E.*cos(2*E)-2*(-200*x+100.0).^2.*E.*cos(2*E)+4*(-200*x+100.0).^2.*E.^2.*sin(2*E);
