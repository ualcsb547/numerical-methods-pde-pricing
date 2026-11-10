function [teul, yeul, ev] = mirk4(f, intv, y0, N)

h = (intv (2) - intv (1))/N ;
y = y0 (:) ;
t = intv (1) ;
yeul = y0 (:) ;
teul = intv (1) ;
ev = 0;

for i = 1: N
    
    F1 = feval(f, t, y);
    F2 = feval(f, t + h/2, y + (h/2)*F1);
    F3 = feval(f, t + h/2, y + (h/2)*F2);
    F4 = feval(f, t + h, y + h*F3);

    phi = (1/6)*(F1 + 2.*F2 + 2.*F3 + F4);
    y1 = y + h * phi;
    t1 = t + h;
    teul = [ teul , t1 ];
    yeul = [ yeul , y1 ];
    y = y1 ;
    t =  t1 ;
    ev = ev + 4; % Es 4 porque evaluamos f 4 veces por iteración (F1...F4)
end