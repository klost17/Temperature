function Racr = prac5criticalReEven(K0)

global Ra0

[solRa,it,res]=prac5newton(@prac5funEven,Ra0,K0,1e-6,100);

if it<100 && res<1e-5
    Racr = solRa;

end