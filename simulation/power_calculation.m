function N = power_calculation(alpha,power_level,p,delta)
if p>0.5
    p = 1-p;
end
talpha2 = norminv(1-alpha/2);
tbeta = norminv(power_level);
sd1 = sqrt(2*p*(1-p));
sd2 = sqrt(p*(1-p) + (p + delta) * (1-p-delta));
N = (talpha2*sd1 + tbeta*sd2) * (talpha2*sd1 + tbeta*sd2) / (delta*delta);
N = round(N);
end