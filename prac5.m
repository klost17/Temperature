% No-slip boundary conditions
Ra_ns_even = 6000; ns_Kcr_even = []; ns_Racr_even = [];
for K = 1:0.001:6
    [solRa,it,res]=prac5newton(@prac5funEven,Ra_ns_even,K,1e-6,100);
    Ra_ns_even = solRa;
    ns_Kcr_even = [ns_Kcr_even K];
    ns_Racr_even= [ns_Racr_even solRa];
end

% Stress-free boundary conditions
n = 1;
sf_Kcr = 0.01:0.001:6;
sf_Racr = ((n*pi)^2+sf_Kcr.^2).^3./sf_Kcr.^2;
sfKcr = n*pi/sqrt(2);
sfRacr = 27*(n*pi)^4/4;

% Critical value of K for which the Ra curve is minimum.
global Ra0
Ra0 = 1500;
nsKcr_even = fminbnd(@prac5criticalReEven,3,3.5);
nsRacr_even = prac5criticalReEven(nsKcr_even);

% Plot everything
figure(1)
plot(ns_Kcr_even,ns_Racr_even,sf_Kcr,sf_Racr,nsKcr_even,nsRacr_even,'ob',sfKcr,sfRacr,'or'); axis([0 6 0 3000])
xlabel('K'); ylabel('Ra_{cr}(K)'); grid; title('Neutral stability curve')
legend('No-slip boundary conditions','Stress-free boundary conditions')