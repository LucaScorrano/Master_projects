function [IAAtot,AccBufferComplete,BufferedData] = BufferAccData(axyz)
    
    % declare persistent variables (their value is preserved between function calls)
    
    persistent axyzbuffer N datacounter;

    datacounter=datacounter+1;
    fs=50; 
    deltaT=1/fs;
    axyzbuffer(1,:)=axyz(:,:,:);
    if datacounter==1500
        IAAx=0;
        IAAy=0;
        IAAz=0;
        for i=1:length(axyzbuffer)-1
            IAAx=IAAx+deltaT/2*((axyzbuffer(i,1))+axyzbuffer(i+1,1))
            IAAy=IAAy+deltaT/2*((axyzbuffer(i,2))+axyzbuffer(i+1,2))
            IAAz=IAAz+deltaT/2*((axyzbuffer(i,3))+axyzbuffer(i+1,3))
        end
        IAAtot=IAAx+IAAy+IAAz;
        AccBufferComplete=True;
        datacounter=1;
    else
        AccBufferComplete=False;
        IAAtot=0;
    end
    BufferedData=axyzbuffer;
    
    




