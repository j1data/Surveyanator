clear,clc
%AerE 261
%Jr.JPL Blake, Ellie, Jeremy, Justin, Nicole
%Surveyanator
%scaling factor 0.28:1 to the bonanza
%based of Beechcraft Bonanza and Cessna 172

%values used
density = 1.0687; %density in Kg/m^3
Temp = 278.41; %Temp Kelvin
mu = 17.4*10^-6; %dynamic viscosity in pa*s
b = 6; %wing span meters
velocity = 50; %velocity in m/s

%Mach Calculation
Mach_val = Mach(velocity,Temp);


%Wing Calculations
cwing = 0.6; %chord length in meters
Swing = b*cwing; %Wing planform area meters^2
wingxovercmax=0.4;
wingtoverc=0.152;
wingsweepangle=0;
Qwing = 1.25; % For a high mounted wing that doesnt fillet with the fuselage

Rewing= Reynolds(density,velocity,mu,cwing);
cfwing = FrictionCoefficient(Rewing,Mach_val);
FFwing = FormFactor(wingxovercmax,wingtoverc,wingsweepangle,Mach_val);
wingswetted = (1.977+0.52*wingtoverc)*Swing; %meters^2 Estimation 
cdowing = cfwing*FFwing*(wingswetted/Swing)*Qwing;


%Fuselage
lengthfuse = 2.36; %length of fuselage in meters

f = lengthfuse/(sqrt((4/pi)*pi));

fffuse = .9+(5/f^1.5)+(f/400);
Refuse=Reynolds(density,velocity,mu,lengthfuse);
swettedfuse = pi*0.8*lengthfuse; %in meters^2
cffuse = FrictionCoefficient(Refuse,Mach_val);
cdofuse = cffuse*fffuse*1*(swettedfuse/2.15);


%Sizing vertical tail
% SVT=0.2; %meters^2
lvt=2; %meters
% CVT=(SVT*LVT)/(b*Swing);
cvt = 0.032;

[svt] = TailVertCoefficient(cvt,lvt,b,Swing);

%Sizing horizontal tail
SHT=.5; %meters^2
LHT=3; %meters
Cavg=Swing/b;
%CHT=(SHT*LHT)/(cwing*Swing);
cht = 0.70;
lht = 2;

[sht] = TailHorizCoefficient(cht,lht,Cavg,Swing);

%Vertical Stabilizer
vertxovercmax = 0.3; % quarter chord length in percent 
cvert = 1; % chord for vert in meters
verttoverc = 0.12;
sweepvert=0;
Revert = Reynolds(density,velocity,mu,cvert);

cfvert = FrictionCoefficient(Revert,Mach_val);
ffvert = FormFactor(vertxovercmax,verttoverc,sweepvert,Mach_val);
swettedvert = cvt/Swing; %meters^2 %Estimation
cdovert = cfvert*ffvert*(swettedvert/2.15)*1.05;

%Horizontal Stabilizer
horizxovercmax = 0.3; % quarter chord length in percent 
choriz = 1; % chord for horz in meters
horiztoverc = 0.12;
sweephoriz = 0;
Rehoriz = Reynolds(density,velocity,mu,choriz);

cfhoriz = FrictionCoefficient(Rehoriz,Mach_val);
ffhoriz = FormFactor(horizxovercmax,horiztoverc,sweephoriz,Mach_val);
swettedhoriz = cht/Swing; %meters^2 %Estimation 
cdohorz = cfhoriz*ffhoriz*(swettedhoriz/2.15)*1.05;

%Landing Gear Struts Main Part D
Sfrontalgear = 2*0.04445*0.2286 ;%2 struts*width(m)*length(m) = meters^2
Doverq = 0.05;
cdolandm = Doverq*(Sfrontalgear/Swing);

%Gear Wheels Part D
Doverqwheel = 0.13;
sfrontalwheel = 3*0.127*0.050038; %multiply sfrontal by 3 for each wheel %meters^2
cdowheels = Doverqwheel*(sfrontalwheel/Swing);

%Landing Gear Strut Back Part D
sfrontalgearback = 0.04445*0.2286; %width(m)*length(m) = meters^2
Doverqbackgear = 0.25;
cdogearback = Doverqbackgear*(sfrontalgearback/Swing);

%Drag buildup lets get it!!!
dragbuildup = cdowing+cdofuse+cdovert+cdohorz+cdolandm+cdogearback; %JP-Should wheels be included in this??

%Lift Equation Calculations
AR=b*cwing;
enot = 1.78*(1-0.045*(AR^0.68))-0.64; %unitless
K = 1/(pi*enot*AR); %unitless
alat0 = -5.41; %degs %angle of attack at lift equals 0
anot = 0.114; %lift slope
a3D = anot/(1+((57.3*anot)/(pi*0.7*AR)));

%Battery Info Part E
E_density = 1044000;%J/kg
Ebattery = 500000000; %Joules
weightbattery = Ebattery/E_density; %kg



%Cl and CL calculations Part D
% L=W
weightempty=300; %Newtons
weightpayload = 200; %Newtons
weight_total = weightempty+weightpayload+weightbattery;

CLift = (weight_total)/(.5*density*(velocity^2)*Swing);

%Fractional weight Calculation %Part D
frac_weightempty = weightempty/weight_total;
frac_weightpayload = weightpayload/weight_total;
frac_weightbattery = weightbattery/weight_total;

%Cl max from XFLR5 Part D
Clmax=1.5382;

%Alpha at steady level flight %Part D
alpha3D_at_SLF = (CLift/a3D)+alat0;



%Part E

%Range 
n_prop = .7; 
n_emotor = .85;
CLoverCD_max = 1/(4*K*dragbuildup);
CLoverCD3halfs_max = (((3*dragbuildup)/K)^(1.5))/(4*dragbuildup);

endurance = (Ebattery*n_prop*n_emotor*((density*Swing)^(.5))*CLoverCD3halfs_max)/((2^(0.5))*(weight_total^(1.5)))/60; %minutes
range = (((Ebattery*n_emotor*n_prop)/(weight_total))*CLoverCD_max)/1000; %km
V_maxrange =((2/density)*(weight_total/Swing)*(K/(3*dragbuildup)^(0.5)))^(.5); %m/s
V_stall = ((2*weight_total)/(density*Swing)*((K/(3*dragbuildup))^(.5)))^(.5); %m/s




%Displaying values of interest
fprintf('CL = %g*(alpha-(%g))\n',a3D,alpha3D_at_SLF) %3Dlift equation
fprintf('Aspect Ratio = %g\n',AR)
fprintf('Planform Area = %g m^2\n',Swing)
fprintf('Mach_val = %g\n',Mach_val)
fprintf('K = %g\n',K)
fprintf('SVT is %g\n',svt)
fprintf('SHT is %g\n',sht)
fprintf('CD0 for the wing is %g\n',cdowing)
fprintf('CD0 for the fuse is %g\n',cdofuse)
fprintf('CD0 for the vertical stabilizer is %g\n',cdovert)
fprintf('CD0 for the horizontal stabilizer is %g\n',cdohorz)
fprintf('CD0 for the front 2 landing gear struts is %g\n',cdolandm)
fprintf('CD0 for the back landing gear is %g\n',cdogearback)
fprintf('CD0 for the landing gear wheels is %g\n',cdowheels)
fprintf('CD0 for the whole plane is %g\n\n',dragbuildup)

%Part D  
fprintf('The empty weight of our aircraft is %g Newtons \n',weightempty) %Part D
fprintf('The payload weight of our aircraft is %g Newtons \n',weightpayload) %Part D
fprintf('The battery weight of our aircraft is %g Newtons \n',weightbattery) %Part D
fprintf('The fractional empty weight is %g \n',frac_weightempty) %Part D
fprintf('The fractional payload weight is %g \n',frac_weightpayload) %Part D
fprintf('The fractional battery weight is %g \n',frac_weightbattery) %Part D
fprintf('The Cl max is %g \n',Clmax) %Part D
fprintf('The CL value for our aircraft is %g at steady level flight. \n',CLift) %Part D
fprintf('Alpha at steady level flight is %g degrees at a CL of %g \n\n',alpha3D_at_SLF,CLift) %Part D

%Part E
fprintf('The endurance is %g minutes \n',endurance)
fprintf('The range is %g kilometers \n',range)
fprintf('The velocity to achieve max range is %g m/s \n',V_maxrange)
fprintf('The stall velocity is %g m/s \n',V_stall)

function [Re] = Reynolds(density,velocity,mu,length)
%Reynolds Number
%density, velocity, dynamic viscosity, and length in metric

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