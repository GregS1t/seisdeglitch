L = length(keepstat_glitch)

sumMonoGlitches = zeros(L);
sumPolyGlitches   = zeros(L);
for i=1:L
    sumMonoGlitches(i) = sum([keepstat_glitch(i).U, keepstat_glitch(i).V, keepstat_glitch(i).W ]);
    sumPolyGlitches(i)   = sum([keepstat_glitch(i).Tilt, keepstat_glitch(i).ECStrain]);  
end

disp("toto")
disp(length(sumMonoGlitches))
disp(length(sumPolyGlitches))
disp(length([keepstat_glitch.sol]))

figure(200)
bar([keepstat_glitch.sol], sumMonoGlitches)
grid on

figure(201)
bar([keepstat_glitch.sol], sumPolyGlitches)
grid on

