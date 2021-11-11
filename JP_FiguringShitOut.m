clear,clc,clf;

inputs = GetGoogleSpreadsheet('1mX9oFI3Zd5SyJR2twwYRWJ417U7gUv60qco82aeOLdE');
inputs = str2double(inputs);
inputs = num2cell(inputs);
%Import Vars from Spreadsheet --->
[density, temp, dynamicViscosity, wingSpan, velocity, wingChord, wingXOverC, wingTOverC, wingSweepAngle, Qwing, fuselageLength] = inputs{2,:};
[LVT, CVT, sht, LHT, cht, lht] = inputs{6,2:7};
[quarterChordVertStab, chordVertStab, vertTOverC, vertSweepAngle, quarterChordHorzStab, chordHorzStab, horzToverC, horzSweepAngle] = inputs{10,1:8};
[frontStruts_Sfront, frontStruts_dOverq, wheel_Sfront, wheel_dOverq, backStrut_Sfront, backStrut_dOverq, E_density, E_battery] = inputs{14,1:8};
[alat0, anot, Cl_max, W_e, W_p, n_prop, n_motor] = inputs{18,1:7};
[Power_Max, g, R, Tsea, density_sea, a, m, altitude_max] = inputs{22,1:8};
[n_Strut_pos, n_Strut_neg] = inputs{25,1:2};

