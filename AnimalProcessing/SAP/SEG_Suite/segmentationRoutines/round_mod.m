function ra=round_mod(a,low,high)
% GOYA Yûki
% v0.1
    ra=round(a);
    if ra<low, ra=low; end
    if ra>high, ra=high; end
end