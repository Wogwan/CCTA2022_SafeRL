% 
% This file generates the necessary Linear F/A-18 models for the analysis purpose 
% around the trim point specified.  
% 
% Abhijit 11/18/2009

%------------ Run the AERODYNAMIC Model
 F18_AeroData = 1; 
 f18_data; 

% ---------- Define Trim Points 

V0 = 350; 
phi0grid = [0 10 25 35]*d2r; 

%---------- Generate Models

Level = 0; 

for i1 = 1:length(phi0grid)
    
    % -------------------------------------------------
    % Coordinated  Turn 
    clear phi0
    CoordTurn = 1; 
    UnCoordTurn = 0;     
    phi0 = phi0grid(i1); 
    
    % --------------Generate Variable Name
    V1 = strcat('F18Model_Coord_Phi_',num2str(phi0*r2d));
    
    GenF18LinModel; 
    GenF18LinCLPModel;  
    
    % --------------- Save Data
    save(V1,'xtrim','utrim','Pss','CssB','CLB','CssR','CLR','linsys')
  
    % -------------------------------------------------
    % UnCoordinated  Turn 
    clear phi0 
    CoordTurn = 0; 
    UnCoordTurn = 1;     
    phi0 = phi0grid(i1); 
   
    % --------------Generate Variable Name
    V2 = strcat('F18Model_UnCoord_Phi_',num2str(phi0*r2d));
    
    GenF18LinModel; 
    GenF18LinCLPModel;  
    
    % --------------- Save Data
    save(V2,'xtrim','utrim','Pss','CssB','CLB','CssR','CLR','linsys')

end

