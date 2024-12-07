function EE = FromAccToEE(IAAtot,AccBufferComplete,BodyMass,EstimateEE)
persistent EEw
    if isempty(EEw)
       EEw=0;
    end
    alpha=0.104;
    beta=0.023;
    Tk=30;
    EEact=0;
    if EstimateEE==true
        if AccBufferComplete==true
          EEact= alpha + beta*IAAtot;    
        else
          EEact=0;
        end
       EEw=EEact*Tk+EEw;
    else
        EEw=0;
    end
    EE=EEw*BodyMass/4184;
    
end
    
       
    
    
       
   
