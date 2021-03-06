% Test file for @chebmatrix/waterfall.m.

function pass = test_waterfall(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

hfig = figure('Visible', 'off');

Q = chebmatrix({chebfun(@(x) sin(pi*x), pref), ...
                chebfun(@(x) sin(pi*(x - 0.3)), pref), ...
                chebfun(@(x) sin(pi*(x - 0.6)), pref)});
pass(1) = doesNotCrash(@() waterfall(Q));
pass(2) = doesNotCrash(@() waterfall(Q, [0 0.3 0.6]));
pass(3) = doesNotCrash(@() waterfall(Q, 'LineWidth', 2, 'FaceAlpha', 0.5, ...
    'FaceColor', 'r'));

close(hfig);

end

function pass = doesNotCrash(fn)

try
    fn();
    pass = true;
catch ME;
    pass = false;
end

end
