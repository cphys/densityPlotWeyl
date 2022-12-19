ClearAll["Global`*"]

tau[dim_, numb_] := 
 KroneckerProduct[PauliMatrix[numb], IdentityMatrix[dim]]

sigma[dim_, numb_] := 
 KroneckerProduct[IdentityMatrix[dim], PauliMatrix[numb]]

k = {0, ky, kz};

Hdeg = vf*Sum[k[[i]]*sigma[2, i], {i, 1, 3}] . tau[2, 3];

HD[vf_, g_, Bz_] := 
  vf*Sum[k[[i]]*sigma[2, i], {i, 1, 3}] . tau[2, 3] + 
   g*Bz*(sigma[2, 3] . tau[2, 0]);

Eigs[Bz_, g_, vf_, ky_, kz_] :=
 Eigs[Bz, g, vf, ky, kz] = {
   -Sqrt[Bz^2 g^2 - 2 Bz g kz vf + (ky^2 + kz^2) vf^2],
   Sqrt[Bz^2 g^2 - 2 Bz g kz vf + (ky^2 + kz^2) vf^2],
   -Sqrt[Bz^2 g^2 + 2 Bz g kz vf + (ky^2 + kz^2) vf^2],
   Sqrt[Bz^2 g^2 + 2 Bz g kz vf + (ky^2 + kz^2) vf^2]}

tabEig[Bz_, i_] := 
  tabEig[Bz, i] = 
   N[ParallelTable[
     N[Eigs[Bz, 1., 1., ky, kz][[i]]], {ky, -3.14, 
      3.14, .1}, {kz, -3.14, 3.14, .1}]];

graphicsTable = Table[
   ListPlot3D[
    {tabEig[Bz, 1], tabEig[Bz, 2], tabEig[Bz, 3], tabEig[Bz, 4]},
    BoxRatios -> {1, 1, 1},
    PlotTheme -> "Detailed",
    DataRange -> {{-\[Pi], \[Pi]}, {-\[Pi], \[Pi]}},
    ViewPoint -> {Pi/2, -Pi, .95},
    AxesLabel -> {"\!\(\*SubscriptBox[\(k\), \(y\)]\)", 
      "\!\(\*SubscriptBox[\(k\), \(z\)]\)", "E"},
    Ticks -> {{-\[Pi], 0, \[Pi]}, {\[Pi], 0, \[Pi]}, {-5, -3, -1, 0, 
       1, 3, 5}},
    PlotRange -> {-5, 5}],
   {Bz, 0, 3, .01}];

Export[
 FileNameJoin[{NotebookDirectory[], "testWeylgif2.gif"}],
 graphicsTable,
 "DisplayDurations" -> 0.04]
