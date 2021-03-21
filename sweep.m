%bl = [0.2 0.4 0.6 0.8 1] * 1e-7;
bl = 40e-9;
bw = linspace(0,0.45,10) * 1e-7;

I = [];

for n = 1 : length(bw)
    I(n) = MC_Current(bw(n),bl);
end

plot(bw,I);