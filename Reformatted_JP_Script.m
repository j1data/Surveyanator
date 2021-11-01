clear,clc,clf
%AerE 261
%Jr.JPL Blake, Ellie, Jeremy, Justin, Nicole
%Surveyanator
%scaling factor 0.28:1 to the bonanza
%based of Beechcraft Bonanza and Cessna 172


inputs = GetGoogleSpreadsheet('1mX9oFI3Zd5SyJR2twwYRWJ417U7gUv60qco82aeOLdE');
inputs = str2double(inputs);
inputs = num2cell(inputs);
%Import Vars from Spreadsheet --->
[density, temp, dynamicViscosity, wingSpan, velocity, wingChord, wingXOverC, wingTOverC, wingSweepAngle, Qwing, fuselageLength] = inputs{2,:};
[LVT, CVT, SHT, LHT, cht, lht] = inputs{6,2:7};
[quarterChordVertStab, chordVertStab, vertTOverC, vertSweepAngle, quarterChordHorzStab, chordHorzStab, horzToverC, horzSweepAngle] = inputs{10,1:8};
[frontStruts_Sfront, frontStruts_dOverq, wheel_Sfront, wheel_dOverq, backStrut_Sfront, backStrut_dOverq, E_density, E_battery] = inputs{14,1:8};
[alat0, anot, Cl_max, W_e, W_p, n_prop, n_motor] = inputs{18,1:7};
[Power_Max, g, R, Tsea, density_sea, a, m] = inputs{22,1:8};

%General Calculations ---->
Mach_val = Mach(velocity, temp);

%Main Wing Calcs -->
Swing = wingSpan*wingChord;
ReWing = Reynolds(density,velocity,dynamicViscosity,wingChord);
cfWing = FrictionCoefficient(ReWing,Mach_val);
FFWing = FormFactor(wingXOverC,wingTOverC,wingSweepAngle,Mach_val);
wingSWetted = (1.977+0.52*wingTOverC)*Swing; %(m^2)
cd_o_Wing = cfWing*FFWing*(wingSWetted/Swing)*Qwing;

%Fuselage Calcs -->
f = fuselageLength/(sqrt(4/pi)*pi);
ffFuse = 0.9 + (5/(f^1.5))+(f/400);
ReFuse = Reynolds(density,velocity,dynamicViscosity,fuselageLength);
sWettedFuse = pi*0.8*fuselageLength; %(m^2)
cfFuse = FrictionCoefficient(ReFuse,Mach_val);
cd_o_Fuse = cfFuse*ffFuse*1*(sWettedFuse/2.15);

%Vert tail sizing calcs -->
[svt] = TailVertCoefficient(CVT,LVT,wingSpan,Swing);

%Horz tail sizing calcs -->
Cavg = Swing/wingSpan;
[sht] = TailHorizCoefficient(cht,lht,Cavg,Swing);

%Vertical Stabilizer calcs -->
ReVert = Reynolds(density,velocity,dynamicViscosity,chordVertStab);
cfVert = FrictionCoefficient(ReVert,Mach_val);
ffVert = FormFactor(quarterChordVertStab,vertTOverC,vertSweepAngle,Mach_val);
sWettedVert = CVT/Swing; %(m^2)
cd_o_Vert = cfVert*ffVert*(sWettedVert/2.15)*1.05;

%Horizontal Stabilizer calcs -->
ReHoriz = Reynolds(density,velocity,dynamicViscosity,chordHorzStab);

cfHoriz = FrictionCoefficient(ReHoriz,Mach_val);
ffHoriz = FormFactor(chordHorzStab,horzToverC,horzSweepAngle,Mach_val);
sWettedHoriz = cht/Swing; %meters^2
cd_o_Horz = cfHoriz*ffHoriz*(sWettedHoriz/2.15)*1.05;

%Landing gear struts (front 2) calcs -->
cd_o_landF = frontStruts_dOverq*(frontStruts_Sfront/Swing);

%Landing gear struts (back 1)  calcs -->
cd_o_landB = backStrut_dOverq*(backStrut_Sfront/Swing);

%Landing Gear Wheels (3 wheels) calcs -->
cd_o_wheels = wheel_dOverq*(wheel_Sfront/Swing);

%Drag buildup calcs -->
dragBuildUp = cd_o_Wing + cd_o_Fuse + cd_o_Vert + cd_o_Horz + cd_o_landB + cd_o_landB + cd_o_wheels;

%Lift calculations -->
AR = wingSpan*wingChord;
e_o = 1.78*(1-0.045*(AR^0.68))-0.64; %Unitless
K = 1/(pi*e_o*AR);
a3D = anot/(1+((57.3*anot)/(pi*0.7*AR)));

%Battery Calculations -->
batteryWeight = (E_battery/E_density)*9.81; %(N)

%Cl & CL calcs -->
W_total = W_e+W_p+batteryWeight; %(Newtons)
CLift = (W_total)/(0.5*density*(velocity^2)*Swing);

%Weight Fraction Calcs -->
frac_W_e = W_e/W_total;
frac_W_p = W_p/W_total;
frac_W_f = batteryWeight/W_total;

%AoA @SLF
alpha3D_SLF = (CLift/a3D)+alat0;

%Range calculations -->
CLoCD_max = 1/(4*K*dragBuildUp);
CLoCD_3half_max = (((3*dragBuildUp)/K)^1.5)/(4*dragBuildUp);

endurance = (E_battery*n_prop*n_motor*((density*Swing)^0.5)*CLoCD_3half_max)/((2^(0.5))*(W_total^(1.5)))/60;
range = (((E_battery*n_motor*n_prop)/W_total)*CLoCD_max)/1000; %km
V_maxrange =((2/density)*(W_total/Swing)*(K/(3*dragBuildUp)^(0.5)))^(.5); %m/s
V_stall = ((2*W_total)/(density*Swing)*((K/(3*dragBuildUp))^(.5)))^(.5); %m/s

%Part E

%Input to spreadsheet


altitude= 0:1:10000; %meters

Temp_alt = Tsea+a.*(altitude); %Kelvin

density_alt = density_sea.*(Temp_alt/Tsea).^((-g/(a*R))-1); %kg/m^3

v_infin_SLF = (((2*W_total)./(density_alt.*Swing).*(K./(3*dragBuildUp)).^(.5)).^(0.5)); %m/s
Power_req = .5.*density_alt.*((v_infin_SLF).^(3)).*Swing*dragBuildUp+((2*K*((W_total)^(2)))./(density_alt.*v_infin_SLF*Swing)); %watts
Power_avail = Power_Max.*((density_alt./density).^(m)); %watts
Power_excess = Power_avail-Power_req; %watts

plot(altitude,Power_req,'color','r')
hold on
plot(altitude,Power_avail,'color','g')
hold on
plot(altitude,Power_excess,'color','b')
hold off
legend('Power Required','Power Available','Power Excess')
xlabel('Altitude (m)')
ylabel('Power (Watt)')
title('Power vs. Altitude')

ROC = Power_excess./W_total; %m/s
VMax_ROC = (((2*W_total)./(density_alt.*Swing)).*((K/(3*dragBuildUp))^(.5))).^(0.5); %m/s
ROC_Max = ((n_prop.*Power_avail)./(W_total))-VMax_ROC.*((1.155)./(CLoCD_max)); %m/s
Service_ceiling = 100/(60*3.28084); %m/s

figure()
plot(altitude,ROC,'color','r')
hold on
plot(altitude,ROC_Max,'color','k')
hold on
yline(Service_ceiling,'color','b')
hold off
xlabel('Altitude (m)')
ylabel('Rate of Climb (m/s)')
legend('ROC','ROC Max','Service Ceiling')
title('Rate of Climb vs Altitude')

% Finding the altitude where ROC max equals the service ceiling
Intersections=find(abs(ROC_Max-Service_ceiling)<=(0.0001));

SC=altitude(Intersections); %Service ceiling in meters


%Formatted output:
%Values of Interest:
%F = figure(1);
%T = table()

%Displaying values of interest
fprintf('CL = %g*(alpha-(%g))\n',a3D,alpha3D_SLF) %3Dlift equation
fprintf('Aspect Ratio = %g\n',AR)
fprintf('Planform Area = %g m^2\n',Swing)
fprintf('Mach_val = %g\n',Mach_val)
fprintf('K = %g\n',K)
fprintf('SVT is %g\n',svt)
fprintf('SHT is %g\n',sht)
fprintf('CD0 for the wing is %g\n',cd_o_Wing)
fprintf('CD0 for the fuse is %g\n',cd_o_Fuse)
fprintf('CD0 for the vertical stabilizer is %g\n',cd_o_Vert)
fprintf('CD0 for the horizontal stabilizer is %g\n',cd_o_Horz)
fprintf('CD0 for the front 2 landing gear struts is %g\n',cd_o_landF)
fprintf('CD0 for the back landing gear is %g\n',cd_o_landB)
fprintf('CD0 for the landing gear wheels is %g\n',cd_o_wheels)
fprintf('CD0 for the whole plane is %g\n\n',dragBuildUp)

%Part D  
fprintf('The empty weight of our aircraft is %g Newtons \n',W_e) %Part D
fprintf('The payload weight of our aircraft is %g Newtons \n',W_p) %Part D
fprintf('The battery weight of our aircraft is %g Newtons \n',batteryWeight) %Part D
fprintf('The fractional empty weight is %g \n',frac_W_e) %Part D
fprintf('The fractional payload weight is %g \n',frac_W_p) %Part D
fprintf('The fractional battery weight is %g \n',frac_W_f) %Part D
fprintf('The Cl max is %g \n',Cl_max) %Part D
fprintf('The CL value for our aircraft is %g at steady level flight. \n',CLift) %Part D
fprintf('Alpha at steady level flight is %g degrees at a CL of %g \n\n',alpha3D_SLF,CLift) %Part D

%Part E
fprintf('The endurance is %g minutes \n',endurance)
fprintf('The range is %g kilometers \n',range)
fprintf('The velocity to achieve max range is %g m/s \n',V_maxrange)
fprintf('The stall velocity is %g m/s \n\n',V_stall)

%Part F
fprintf('The service ceiling is %g meters \n', SC)

%Functions used in the program ---->

function [Re] = Reynolds(density,velocity,mu,length)
    %Reynolds Number (density, velocity, dynamic viscosity, and length in %metric)
    Re=(density*velocity*length)/mu;
end

function [M] = Mach(velocity,Temp)
    %Gives the mach number from the velocity (m/s) and Temp (K)
    gamma = 1.4; %unitless
    R = 287; %J/(mol*K)
    a = sqrt(gamma*R*Temp); %speed of sound in m/s
    M = velocity/a; %Mach number
end

function [Cf] = FrictionCoefficient(Re,Mach)
    %calculation for the friction coefficient
    %   Reynolds, velocity, and temp K in SI
    Cf = (0.455)/(((log10(Re))^2.58)*(1+0.144*Mach^2)^0.65);
end

function [FF] = FormFactor(xovercmax,toverc,sweep,Mach)
    %Function for Form Factor
    %   (x/c)max, t/c, sweep angle(degrees) 
    FF=(1+((0.6)/xovercmax)*(toverc)+100*toverc^4)*((1.34*Mach^.18)*(cosd(sweep)^0.28));
end

function [svt] = TailVertCoefficient(cvt,lvt,b,Swing)
    %Calculation for helping to determine the size of the vertical tail
    %   Enter a coefficient, length between cg and qc, span, and planform area to
    %   determine the exposed side of the vertical tail wing
    svt=(cvt*b*Swing)/(lvt);
end


function [sht] = TailHorizCoefficient(cht,lht,c,Swing)
    %Calculation for helping to determine the size of the horizontal tail
    %   Enter a coefficient, length between cg and qc, span, and planform area to
    %   determine the exposed side of the vertical tail wing
    sht=(cht*c*Swing)/(lht);
end

