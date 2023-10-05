function Isp = specificImpulse(I, mp, g0)
    Isp = I/(mp*g0);
end

function T = thrust(mdot, ue, Ae, pe, pa)
    T = mdot*ue + Ae*(pe-pa);
end

function M = machNumber(v, a)
    M = v / a;
end

function a = speedOfSound(adiabaticConstant, R, T)
    a = sqrt(adiabaticConstant * R * T);
end

function stagnationTemp = stagnationTemp(adiabaticConstant, M)
    stagnationTemp = 1 + (adiabaticConstant - 1)/2 * M ^ 2;
end

function stagnationPressure = stagnationPressure(adiabaticConstant, M)
    stagnationPressure = (1+(y-1)/2*M^2)^(y/(y-1));
end

function stagnationDensity = stagnationDensity(adiabaticConstant, M)
    stagnationDensity = (1+(y-1)/2*M^2)^(1/(y-1));
end
