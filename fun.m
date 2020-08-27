function t= fun(val) 
if ((real(val) < -1) && (abs(imag(val))<1e-7)) 
t=1;
 else
t=0;
 end