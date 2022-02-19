function [sys, x0, str,ts] = f18_sfcnUnc(t, x ,u, flag,x0,opts)
%
%
% This is a S-function representation for the plant described in
% f18full.m .  
% This S-function is used for trimming and linearizing the aircraft.
% This S-function is being called in the f18_trim.mdl file.
%
% Abhijit 11/12/2009

switch flag
    
    case 0 % initialize
        
    
        str =[];
        ts = [0 0];
            
        s = simsizes ;
        
        s.NumContStates = 9;
        s.NumDiscStates = 0;
        s.NumOutputs = 9;
        s.NumInputs = 4;
        s.DirFeedthrough = 0;
        s.NumSampleTimes = 1; 
        
        sys = simsizes(s);
    
        
    case 1  % derivative
  
            sys = f18fullUnc(t,x,u,opts);
        
         
    case 3  % output
        
        sys = x ;
        
    case {2 4 9} 
        
        sys = [];
end