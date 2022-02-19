
% Load Plant and Controllers: 
% Flight Condition: 35 degree Coordinated Turn 
Bank35 = load('1_F18Model_Coord_Phi_35');

redstates = {'beta (rad)','alpha (rad)','p (rad/s)','q (rad/s)', ...
            'r (rad/s)','phi (rad)'}';
inputnames = { 'aileron (rad)','rudder (rad)','stab (rad)'}';       

outputnames = {'a_y (g)','p (rad/s)','r (rad/s)','alpha (rad)', ...
                'beta (rad)','q (rad/s)','betadot'}';        

% Plant Description 
P = ss(Bank35.PA, Bank35.PB,Bank35.PC,Bank35.PD,'statename',redstates,...
                  'inputname',inputnames, 'outputname',outputnames);

% Baseline Controller
CB = ss(Bank35.CBA, Bank35.CBB,Bank35.CBC,Bank35.CBD); 

% Revised Controller
CR = ss(Bank35.CRA, Bank35.CRB,Bank35.CRC,Bank35.CRD); 

% Form Closed Loop Plants
PclpB = feedback(P, CB);    % Baseline
PclpR = feedback(P, CR);    % Revised 
