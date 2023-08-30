%==========================================================================
% Aircraft Physical Paramters
% Reference: S. B. Buttrill and P. D. Arbuckle  and K. D. Hoffler
%            Simulation model of a twin-tail, high performance airplane
%            NASA 1992, NASA TM-107601



S = 400;                % Reference Area, ft^2
b =  37.42;             % Wing Span, ft
c =  11.52;             % Aerodynamic Mean Chord, ft
rho = 1.0660e-003;      % Air Density, slugs/ft^3  --- 25C / 25000 ft
Ixx = 23000;            % Principle Moment of Intertia around X-axis, slugs*ft^2
Iyy = 151293;           % Principle Moment of Intertia around Y-axis,slugs*ft^2 
Izz = 169945;           % Principle Moment of Intertia around Z-axis,slugs*ft^2 
Ixz = - 2971;           % Principle Moment of Intertia around XZ-axis,slugs*ft^2
m = 1034.5;             % mass, slugs
g = 32.2;               % gravitational constant,ft/s^2

% Unit Conversion : Degree <--> Radian
d2r = pi/180;
r2d = 1/d2r;


%==========================================================================
% Aircraft Aero Model

% For Fast Simulation. Avoid loading MAT file 

if ~exist('F18_AeroData','var')
    F18_AeroData = 0;
end

if F18_AeroData
    F18Aero.Clb_0 = -0.0556;
    F18Aero.Clb_1 = -0.4153;
    F18Aero.Clb_2 = -0.3620;
    F18Aero.Clb_3 = 2.3843;
    F18Aero.Clb_4 = -1.6196;

    F18Aero.Cldr_0= 0.0129;
    F18Aero.Cldr_1= 0.0014;
    F18Aero.Cldr_2= 0.0083;
    F18Aero.Cldr_3= -0.0274;

    F18Aero.Clda_0= 0.1424;
    F18Aero.Clda_1= -0.0516;
    F18Aero.Clda_2= -0.2646;
    F18Aero.Clda_3= 0.1989;

    F18Aero.Clp_0= -0.3540;
    F18Aero.Clp_1= 0.2377;

    F18Aero.Clr_0= 0.1983;
    F18Aero.Clr_1= 0.7804;
    F18Aero.Clr_2= -1.0871;

    F18Aero.Cnb_0= 0.0885;
    F18Aero.Cnb_1= 0.0329;
    F18Aero.Cnb_2= -0.3816;

    F18Aero.Cndr_0= -0.0780;
    F18Aero.Cndr_1= -0.0176;
    F18Aero.Cndr_2= 0.5564;
    F18Aero.Cndr_3= -0.8980;
    F18Aero.Cndr_4= 0.3899;

    F18Aero.Cnda_0= 0.0104;
    F18Aero.Cnda_1= 0.0584;
    F18Aero.Cnda_2= -0.3413;
    F18Aero.Cnda_3= 0.2694;

    F18Aero.Cnr_0= -0.4326;
    F18Aero.Cnr_1= -0.1307;

    F18Aero.Cnp_0= 0.0792;
    F18Aero.Cnp_1= -0.0881;

    F18Aero.Cyb_0= -0.7344;
    F18Aero.Cyb_1= 0.2654;
    F18Aero.Cyb_2= -0.1926;

    F18Aero.Cyda_0= -0.1656;
    F18Aero.Cyda_1= -0.2403;
    F18Aero.Cyda_2= 1.5317;
    F18Aero.Cyda_3= -0.8500;

    F18Aero.Cydr_0= 0.2054;
    F18Aero.Cydr_1= 0.4082;
    F18Aero.Cydr_2= -1.6921;
    F18Aero.Cydr_3= 0.9351;

    F18Aero.Cma_0= -0.0866;
    F18Aero.Cma_1= 0.5110;
    F18Aero.Cma_2= -1.2897;

    F18Aero.Cmds_0= -0.9051;
    F18Aero.Cmds_1= -0.3245;
    F18Aero.Cmds_2= 0.9338;

    F18Aero.Cmq_0= -4.1186;
    F18Aero.Cmq_1= 10.9921;
    F18Aero.Cmq_2= -68.5641;
    F18Aero.Cmq_3= 64.7190;

    F18Aero.CLds_0= 0.5725;
    F18Aero.CLds_1= 0.4055;
    F18Aero.CLds_2= -2.6975;
    F18Aero.CLds_3= 2.1852;

    F18Aero.Cdds_0= 0.0366;
    F18Aero.Cdds_1= -0.2739;
    F18Aero.Cdds_2= 4.2360;
    F18Aero.Cdds_3= -3.8578;
end

