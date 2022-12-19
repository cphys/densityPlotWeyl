ClearAll["Global`*"]

(* define the pauli matricies for particle hole and spin spaces *)
tau[dim_, numb_] := 
 tau[dim, numb] = 
  KroneckerProduct[PauliMatrix[numb], IdentityMatrix[dim]]
sigma[dim_, numb_] := 
 sigma[dim, numb] = 
  KroneckerProduct[IdentityMatrix[dim], PauliMatrix[numb]]
  
(* The Hamiltonian for the system *)
H4Lat[kx_, ky_, kz_, t_, tz_, b_, \[Lambda]_, 
  Bz_] := (t*(Cos[kx] + Cos[ky] - 2) + tz*(Cos[kz] - b))*
   sigma[2, 0] . tau[2, 3] + \[Lambda]*Sin[kx]*
   sigma[2, 1] . tau[2, 1] + \[Lambda]*Sin[ky]*
   sigma[2, 2] . tau[2, 1] + Bz*(sigma[2, 3] . tau[2, 0])
   
(* Find the Eigenvalues of the system *)
eigs[kx_, ky_, kz_, i_, Bz_] := 
 eigs[kx, ky, kz, i, Bz] = 
  N[Eigenvalues[H4Lat[kx, ky, kz, 1., 1., 1., 1., Bz]]][[i]]
  
tabEig[i_, Bz_] := 
 tabEig[i, Bz] = 
  N[ParallelTable[
    eigs[kx, ky, kz, i, Bz], {kx, -3.14, 3.14, .05}, {ky, -3.14, 
     3.14, .05}, {kz, -3.14, 3.14, .05}]]
     
totalE[Bz_] := Sum[(-1)^i*Abs[tabEig[i, Bz]], {i, 1, 4}]

densPlot[numbs_, Bz_] :=
 ListDensityPlot3D[
  (-1)^numbs*Abs[tabEig[numbs, Bz]],
  ImageSize -> 400,
  ImagePadding -> {{50, 50}, {30, 10}},
  BoxRatios -> {1, 1, 1},
  PlotTheme -> "Detailed",
  AxesLabel -> {Style["\!\(\*SubscriptBox[\(k\), \(x\)]\)", 
     FontWeight -> Bold, FontSize -> 24],
    Style["\!\(\*SubscriptBox[\(k\), \(y\)]\)", FontWeight -> Bold, 
     FontSize -> 24],
    Style["\!\(\*SubscriptBox[\(k\), \(z\)]\)", FontWeight -> Bold, 
     FontSize -> 24]},
  Ticks -> {{-\[Pi], 0, \[Pi]}, {-\[Pi], 0, \[Pi]}, {-\[Pi], 
     0, \[Pi]}},
  TicksStyle -> Directive[14],
  PlotLegends -> Placed[BarLegend[{"SunsetColors", {0, 3}}], Below],
  DataRange -> {{-3.2, 3.2}, {-3.2, 3.2}, {-3.2, 3.2}},
  OpacityFunction -> Interval[{-.5, .5}],
  OpacityFunctionScaling -> False ,
  ColorFunction -> "SunsetColors"]
  
contPlot[numbs_, Bz_] :=
 ListDensityPlot3D[
  (-1)^numbs*Abs[tabEig[numbs, Bz]],
  ImageSize -> 400,
  ImagePadding -> {{50, 50}, {30, 10}},
  PlotTheme -> "Detailed",
  AxesLabel -> {Style["\!\(\*SubscriptBox[\(k\), \(x\)]\)", 
     FontWeight -> Bold, FontSize -> 24],
    Style["\!\(\*SubscriptBox[\(k\), \(y\)]\)", FontWeight -> Bold, 
     FontSize -> 24],
    Style["\!\(\*SubscriptBox[\(k\), \(z\)]\)", FontWeight -> Bold, 
     FontSize -> 24]},
  Ticks -> {{-\[Pi], 0, \[Pi]}, {-\[Pi], 0, \[Pi]}, {-\[Pi], 
     0, \[Pi]}},
  TicksStyle -> Directive[14],
  PlotLegends -> Placed[BarLegend[{"SunsetColors", {0, 7}}], Below],
  DataRange -> {{-3.2, 3.2}, {-3.2, 3.2}, {-3.2, 3.2}},
  PlotRange -> {{-3.2, 3.2}, {0, 3.2}, {-3.2, 3.2}} ,
  BoxRatios -> Automatic,
  ColorFunction -> "SunsetColors",
  OpacityFunction -> None,
  OpacityFunctionScaling -> False]
  
densPlotTotalE[Bz_] :=
 ListDensityPlot3D[
  totalE[Bz],
  BoxRatios -> {1, 1, 1},
  PlotTheme -> "Detailed",
  ViewPoint -> {Pi/2, -Pi, .65},
  AxesLabel -> {"k_x", "\!\(\*SubscriptBox[\(k\), \(y\)]\)", 
    "\!\(\*SubscriptBox[\(k\), \(z\)]\)"},
  Ticks -> {{-\[Pi], 0, \[Pi]}, {\[Pi], 0, \[Pi]}, {\[Pi], 
     0, \[Pi]}},
  DataRange -> {{-3.2, 3.2}, {-3.2, 3.2}, {-3.2, 3.2}}]
  
graphicsTable = 
  Table[Row[{contPlot[2, B], contPlot[4, B], densPlot[2, B], 
     densPlot[4, B]}], {B, 0, 2, .25}];
     
Export[
  FileNameJoin[{NotebookDirectory[], "weylDensity.gif"}],
  graphicsTable,
  "DisplayDurations" -> 0.4];
