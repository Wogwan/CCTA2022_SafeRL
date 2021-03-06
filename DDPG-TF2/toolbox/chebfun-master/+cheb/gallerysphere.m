function varargout = gallerysphere(name)
%CHEB.GALLERYSPHERE   Spherefun example functions.
%   F = CHEB.GALLERYSPHERE(NAME) returns a spherefun corresponding to NAME.
%   See the listing below for available names.
%
%   For example,  plot(cheb.gallerysphere('football')) plots the classic
%   icosahedral pattern of the Addias Telestar football (or soccer ball for
%   Americans).  For details of how each function is constructed try 
%   'type +cheb/gallerysphere' or 'edit cheb.gallerysphere'.
%
%   [F,FA] = CHEB.GALLERYSPHERE(NAME) also returns the anonymous function
%   FA used to define the function. Some gallery functions are generated by
%   operations beyond the usual Spherefun constructor, so FA in those cases
%   is equal to F.
%
%   CHEB.GALLERYSPHERE with no input argument returns a function chosen at
%   random from the gallery.
%
%   CHEB.GALLERYSPHERE with no output argument creates a plot of the selected
%   function.
%
%   Available names:
%
%   deathstar   A function resembling the surface of the Death Star.
%   gaussian    Gaussian function on the sphere centered at Gauss's birth place.
%   geomag      Radial component of the International Geomagnetic Reference
%               field from the IGRF-12 model for 2015.
%   football    Icosahedral pattern found on a traditional soccer ball.
%   jet         A zonal jet stream over the mid-latitudes of the northern hemisphere.
%   moire       Moire pattern from waves generated at two point sources.
%   neamtu      Function created by Mike Neamtu for testing various spline
%               interpolation methods on the sphere (see Alfeld, Neamtu,
%               Schumaker, J. Comput. Appl. Math. 1996).
%   peaks       Peaks like function on the sphere taken from the geopeaks
%               function in the MATLAB mapping toolbox.
%   randn       Random linear combination of all real spherical harmonics 
%               of exact degree 40.  The coefficients are generated from a 
%               i.i.d. Gaussian (normal) distribution with std=1.
%   reprodkern  Reproducing kernel for spherical harmonics of degree 10
%               centered at (x,y,z) = (-1/sqrt(3),-1/sqrt(3),1/sqrt(3)).
%   soccerball  Same as football, but for the American users.
%   stripes     Alternating striped pattern.
%   vortices    Two antipodal vortices taken from  Nair, Cote, and 
%               Staniforth (1999).
%   :)          A function to make you happy
%
%   Gallery functions are subject to change in future releases of Chebfun.
%
% See also CHEB.GALLERY, CHEB.GALLERYTRIG, CHEB.GALLERY2, CHEB.GALLERY3, CHEB.GALLERYDISK.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% If the user did not supply an input, return a function chosen at random
% from the gallery.
if ( nargin == 0 )
    names = {'football','jet','soccerball','deathstar', 'vortices'...
        'gaussian','reprodkern','geomag','peaks','neamtu','randn',...
        'moire','stripes',':)'};
    name = names{randi(length(names))};
end

type=1;                  % Default plotting type
addEarthPlot = 0;        % Flag on whether or not to plot the earth (1=yes)
viewAngle = [-37.5 30];  % Default viewing angle
clrmap = parula(64);     % Default colormap
clraxis = [];            % Color axis for the colormap; empty means use default

% The main switch statement.
switch lower(name)    
    % The classic football (or soccerball) pattern
    case {'football','soccerball'}
        f = spherefun.sphharm(6,0) + sqrt(14/11)*spherefun.sphharm(6,5);
        fa = f;
        cntrlvl = -[0.25 0.25];
        type = 3;
        
    % A function resembling the surface of the Death Star
    case 'deathstar'
        fa = @(x,y,z) -(exp(-30*((y+sqrt(3)/2).^2 + x.^2 + (z-1/2).^2)) + exp(-100*z.^2));
        f = spherefun(fa);
        type = 1;
        clrmap = flipud(jet);
        viewAngle = [-35 8];
        
    % Two diametrically oposed vortices taken from
    % R. D. Nair, J. Cote, A. Staniforth, Cascade interpolation for 
    % semi-Lagrangian advection over the sphere, Quart. J. Roy. Meteor. 
    % Soc. 125 (1999) 1445-1468
    case 'vortices'
        rho = @(th) 3*(sin(th)); w = @(th) (3*sqrt(2)/2*sech(rho(th)).^2.*tanh(rho(th)))./(rho(th)+eps);
        fa = @(lam,th) -tanh(rho(th)/5.*sin(lam-20*w(th)));
        f = spherefun(fa);
        cntrlvl = [0 0];
        type = 3;
        clrmap = jet;
        
    % A simple Gaussian
    case 'gaussian'
        % Coordinates of Braunschweig, where Gauss was born;
        coords = [10.516667 52.266667]/180*pi;
        [xc,yc,zc] = sph2cart(coords(1),coords(2),1);
        fa = @(x,y,z,xc,yc,zc) exp(-20*((x-xc).^2 + (y-yc).^2 + (z-zc).^2));
        f = spherefun(@(x,y,z) fa(x,y,z,xc,yc,zc));
        addEarthPlot = 1;
        viewAngle = [80 15];
        
    % The reproducing kernel for all spherical harmonics of degree 20
    case 'reprodkern'
        [lam0,th0] = cart2sph(-1/sqrt(3),-1/sqrt(3),1/sqrt(3));
        f = spherefun(@(lam,th) sphRPK(lam,th,lam0,pi/2-th0,20));
        fa = f;
        viewAngle = [-50 5];        
        
    % The IGRF-12 geomagnetic field
    case 'geomag'
        % Spherical harmonic coefficients for the IGRF-12 2015 Geomagnetic
        % Model; see http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
        c = [-29442.0, -1501.0, 4797.1, -2445.1, 3012.9, -2845.6, 1676.7,...
            -641.9, 1350.7, -2352.3, -115.3, 1225.6, 244.9, 582.0, -538.4,...
            907.6, 813.7, 283.3, 120.4, -188.7, -334.9, 180.9, 70.4, -329.5,...
            -232.6, 360.1, 47.3, 192.4, 197.0, -140.9, -119.3, -157.5, 16.0,...
            4.1, 100.2, 70.0, 67.7, -20.8, 72.7, 33.2, -129.9, 58.9, -28.9,...
            -66.7, 13.2, 7.3, -70.9, 62.6, 81.6, -76.1, -54.1, -6.8, -19.5,...
            51.8, 5.7, 15.0, 24.4, 9.4, 3.4, -2.8, -27.4, 6.8, -2.2, 24.2, 8.8,...
            10.1, -16.9, -18.3, -3.2, 13.3, -20.6, -14.6, 13.4, 16.2, 11.7,...
            5.7, -15.9, -9.1, -2.0, 2.1, 5.4, 8.8, -21.6, 3.1, 10.8, -3.3, 11.8,...
            0.7, -6.8, -13.3, -6.9, -0.1, 7.8, 8.7, 1.0, -9.1, -4.0, -10.5, 8.4,...
            -1.9, -6.3, 3.2, 0.1, -0.4, 0.5, 4.6, -0.5, 4.4, 1.8, -7.9, -0.7,...
            -0.6, 2.1, -4.2, 2.4, -2.8, -1.8, -1.2, -3.6, -8.7, 3.1, -1.5, -0.1,...
            -2.3, 2.0, 2.0, -0.7, -0.8, -1.1, 0.6, 0.8, -0.7, -0.2, 0.2, -2.2,...
            1.7, -1.4, -0.2, -2.5, 0.4, -2.0, 3.5, -2.4, -1.9, -0.2, -1.1, 0.4,...
            0.4, 1.2, 1.9, -0.8, -2.2, 0.9, 0.3, 0.1, 0.7, 0.5, -0.1, -0.3, 0.3,...
            -0.4, 0.2, 0.2, -0.9, -0.9, -0.1, 0.0, 0.7, 0.0, -0.9, -0.9, 0.4,...
            0.4, 0.5, 1.6, -0.5, -0.5, 1.0, -1.2, -0.2, -0.1, 0.8, 0.4, -0.1,...
            -0.1, 0.3, 0.4, 0.1, 0.5, 0.5, -0.3, -0.4, -0.4, -0.3, -0.8];
        % Compute the vertical component of the magnetic field
        f = 0*spherefun.sphharm(0,0);
        k = 1;
        for l=1:13
            f = f + -((l+1))*sqrt(4*pi/(2*l+1))*c(k)*spherefun.sphharm(l,0);
            k = k + 1;
            for m=1:l
                f = f - (-1)^m*((l+1))*sqrt(4*pi/(2*l+1))*...
                        (c(k)*spherefun.sphharm(l,m) + ...
                        c(k+1)*spherefun.sphharm(l,-m));        
                k = k + 2;
            end
        end
        fa = f;
        type = 3;
        addEarthPlot = 1;
        cntrlvl = -60000:5000:60000;
        
    % A "peaks-like" function for the sphere
    case 'peaks'
        fa = @(x,y,z) 8*(1-x).^2.*exp(-4*(x - 0.059).^2 - 2*(y + 0.337).^2 - 2*(z + 0.940).^2) - ...
            30*(z/10 - x.^3 - y.^5) .* exp(-3*(x - 0.250).^2 - 2*(y - 0.433).^2 - 3*(z - 0.866).^2) + ...
            (20*y - 8*z.^3) .* exp(-2*(x + 0.696).^2 - 3*(y + 0.123).^2 - 2*(z - 0.707).^2) + ...
            (7*y - 10*x + 10*z.^3) .* exp(-3*(x - 0.296).^2 - 3*(y + 0.814).^2 - 3*(z + 0.5).^2);
        f = spherefun(fa);
        
    % A function developed by Mike Neamtu of Vanderbilt University to test
    % approximation schemes on the sphere.
    case 'neamtu'
        fa = @(x,y,z) 1 + x.^8 + exp(2*y.^3) + exp(2*z.^2) + 10*x.*y.*z;
        f = spherefun(fa);
        
    % A Guassian random function on the sphere suggested by Dmitry Belyaev
    % at Oxford.
    case 'randn'
        % Compute a random spherefun
        f = randnfunsphere(0.079,'monochromatic');
        fa = f;
        type = 4;
        clrmap = gray(2);
        % Set the color axis so that the transition from white to black is
        % exactly at the value zero.
        minmax = minandmax2est(f);
        clraxis = [-1 1]*max(abs(minmax));

    % An interesting Moire pattern generated by two sources on the sphere
    case 'moire'
        % Centers of the beacons
        boise = [-116.237651 43.613739]*pi/180;
        oxford = [-1.257778 51.751944]*pi/180;
        % ithaca = [-76.5 42.443333]*pi/180;
        % stellenbosh = [18.86 -33.92]*pi/180;
        [xb,yb,zb] = sph2cart(boise(1),boise(2),1);
        [xo,yo,zo] = sph2cart(oxford(1),oxford(2),1);
        % Pick the number of oscillations and make each of the "waves" 
        % vanish at the anti-podal points from their centers.
        omega = besselroots(0,30); omega = omega(end)/2;
        % Use a combination of the J0 bessel functions centered at Boise
        % and Oxford to generate the Moire pattern.
        fa = @(x,y,z,omega) 2 + besselj(0,omega*sqrt((x-xb).^2+(y-yb).^2+(z-zb).^2)) + ...
            2 + besselj(0,omega*sqrt((x-xo).^2+(y-yo).^2+(z-zo).^2));
        f = spherefun(@(x,y,z) fa(x,y,z,omega));
        type = 1;
        addEarthPlot = 1;
        viewAngle = [32 8];
        
    % A zipper-like stripe pattern for the sphere    
    case 'stripes'
        fa = @(x,y,z) (1 + cos(10*pi*x)).*exp(-exp(-20*z)) + (1 - cos(10*pi*x)).*exp(-exp(20*z));
        f = spherefun(fa);
        type = 3;
        cntrlvl = [1 1];
        viewAngle = [-36 8];
        clrmap = hot;
        
    % A zonal jet stream in the mid-latitudes of the northern hemisphere.
    case 'jet'
        fa = spherefun(@(lam,th) sin(th-(pi/4+0.01*cos(12*lam))).*exp(-300*(1-cos(th-(pi/4+0.01*cos(12*lam))))));
        f = spherefun(fa);
        type = 1;
        clrmap = jet;
        viewAngle = [32 8];
        addEarthPlot = 1;
        
    % A smiley face emoticon
    case ':)'
        fa = @(x,y,z) exp(-20*(4*(x+cos(pi/6)*sin(pi/4)).^2 + 4*(y+sin(pi/6)*sin(pi/4)).^2 + (z-cos(pi/4)).^2)) + ...
                      exp(-20*(4*(x+cos(-pi/6)*sin(pi/4)).^2 + 4*(y+sin(-pi/6)*sin(pi/4)).^2 + (z-cos(pi/4)).^2)) + ...
                   9*exp(-100*((x + 1).^2 - x.*z + (z - 1/2).^2 - 1/10).^2);
        f = spherefun(fa);
        viewAngle = [-75 20];
        
    otherwise
        error('CHEB:GALLERYSPHERE:unknown:unknownFunction', ...
            'Unknown function.')
end

% Only return something if there is an output argument.
if ( nargout > 0 )
    varargout = {f, fa};
    return;
end

ptitle = [name ', rank = ' num2str(length(f))];
% Determine the type of plot to make
if type==1
    % Plot the function with a grid
    surf(f,'grid','k-')
    axis off, title(ptitle)
elseif type==2
    contour(f)
    axis off, title(ptitle)
elseif type==3
    plot(f), hold on
    contour(f,cntrlvl,'k-'), hold off
    axis off, title(ptitle)
else
    plot(f)
    axis off, title(ptitle)
end

view(viewAngle);
colormap(clrmap);

% Set the color axis if required.
if ( ~isempty(clraxis) )
    caxis(clraxis);
end

if ( addEarthPlot )
    hold on, spherefun.plotEarth('w-'), hold off
end

end

function f = sphRPK(lam,th,lam0,th0,deg)
%SPHRPK  Reproducing kernel for the spherical harmonics of a given degree.
%   F = SPHRPK(LAM,TH,LAM0,TH0,DEG) is the reproducing kernel centered at
%   (LAM0,TH0) for the space of spherical harmonics of degree DEG.

t = cos(lam).*sin(th).*cos(lam0).*sin(th0) + sin(lam).*sin(th).*sin(lam0).*sin(th0) + cos(th).*cos(th0);
pl = legpoly(0:deg,[-1,1]);
[m,n] = size(t);
t = t(:);
c = ones(length(t),1)*((2*(0:deg)+1)/4/pi);
f = reshape(sum(c.*pl(t),2),m,n);

end
