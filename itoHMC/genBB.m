#!/Applications/Mathematica.app/Contents/MacOS/MathematicaScript -script
argv=$ScriptCommandLine;

eps =  ToExpression[argv[[2]]];
seed = ToExpression[argv[[3]]];

(* Gaussian random numbers *)
grand := Sqrt[ -2*Log[1. - Random[]]]*Cos[ 2 Pi Random[]];



SeedRandom[seed]
Nu = 20000;
du = .005;
U = Nu d;
nb = Nu + 1;

x = Table[0, {Nu + 1}];

Do[ x[[i + 1]] = x[[i]] + Sqrt[2.*eps*du] grand;, {i, 1, Nu}]

Do[ x[[i + 1]] -=  (i/Nu) x[[Nu + 1]];
 , {i, 1, Nu}]

Export["randBBlist.dat",x,"Data"]
