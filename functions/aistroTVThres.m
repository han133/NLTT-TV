  function [z1] = aistroTVThres(x,v1,sizeD,lam,beta1,alpha1,mode)
  
  tmp  = reshape(x - v1/2*beta1, sizeD);
  lam1 = lam*alpha1/2*beta1;
   
  if lam1 ~=0
    z = anisotrTV(tmp,lam1,mode);
  else
    z = tmp;
  end
   
  z1    = z(:); 
  end