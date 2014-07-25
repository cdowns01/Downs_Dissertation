clear all
close all

%Load data from stability test runs for n-type material
load ('stability test/rand1/1000terms/iteration1.mat')
%Cpjs contains data points for p-type carrier density
%Cnjs contains data points for n-type carrier density

%p- and n-type carriers in n-type material
nCpjs = Cpjs;
nCnjs = Cnjs;
%Calculate for n-type side

%%%%%%%%%%%%%%%%%%
%CHECK UNITS
%CHECK UNITS
%%%%%%%%%%%%%%%%%%
nDopingLevel = 10^18;%n-type doping level (in 1/cm^3)
intrinsicLevel = 2*10^6; %intrinsic carier concentration (in 1/cm^3)


nCnjs_final = nCnjs+nDopingLevel+intrinsicLevel;
nCpjs_final = nCpjs+intrinsicLevel;

Ei = (Ec+Ev)/2 + 3/4*keV*T*log(meffh/meffe);
Ef1 = Ei+ksi*T/q*log(Cnjs./Cpjs);

ndiff1 = ksi*T/q*log(nCnjs./nCpjs);
ndiff2 = ksi*T/q*log(nCnjs_final./nCpjs_final);

nEc1 = Ec-ndiff1;
nEv1 = Ev-ndiff1;

nEc2 = Ec-ndiff2;
nEv2 = Ev-ndiff2;

% plot(indexn,nEc1,indexn,nEv1,indexn,Ec,indexn,Ev,indexn,Ef1,indexn,Ei),title('n-type (no doping)')
% figure
% plot(indexn,nEc2,indexn,nEv2,indexn,Ec,indexn,Ev,indexn,Ei),title('n-type (doping)')
% figure


%Load data from stability test runs for p-type material
load ('stability test/rand1/1000terms/iteration1.mat')

%p- and n-type carriers in p-type material
pCpjs = Cpjs;
pCnjs = Cnjs;

%Calculate for p-type side

%%%%%%%%%%%%%%%%%%
%CHECK UNITS
%CHECK UNITS
%%%%%%%%%%%%%%%%%%
pDopingLevel = 10^18; %n-type doping level (in 1/cm^3)
intrinsicLevel = 2*10^6; %intrinsic carier concentration (in 1/cm^3)

indexn2 = indexn+Wn;

pCnjs_final = pCnjs+intrinsicLevel;
pCpjs_final = pCpjs+pDopingLevel+intrinsicLevel;

Ei = (Ec+Ev)/2 + 3/4*keV*T*log(meffh/meffe);

pdiff1 = -ksi*T/q*log(pCpjs./pCnjs);
pdiff2 = -ksi*T/q*log(pCpjs_final./pCnjs_final);

pEc1 = Ec-pdiff1;
pEv1 = Ev-pdiff1;

pEc2 = Ec-pdiff2;
pEv2 = Ev-pdiff2;


% plot(indexn2,pEc1,indexn2,pEv1,indexn2,Ec,indexn2,Ev,indexn2,Ei),title('p-type (no doping)')
% figure
% plot(indexn2,pEc2,indexn2,pEv2,indexn2,Ec,indexn2,Ev,indexn2,Ei),title('p-type (Doping)')
% figure

totalIndex = [indexn,indexn2];
EcTotal1 = [nEc1, pEc1];
EcTotal2 = [nEc2, pEc2];

EvTotal1 = [nEv1, pEv1];
EvTotal2 = [nEv2, pEv2];



% plot(totalIndex,EcTotal2,totalIndex,EvTotal2,totalIndex,Ec,totalIndex,Ev,totalIndex,Ei),title('combined'),xlim([0,2*Wn])


%%%%%%%%%%%%%%%%%%%
%Calculate ratio of widths on each side of the depletion region, if it
%exists.
%Dependent on the doping level of the materials.
%%%%%%%%%%%%%%%%%%

%Define the wn/wp for the depletion region
widthRatio = pDopingLevel/nDopingLevel;

%nIndexCount iterates through each index point on the n-side
%pIndexCount holds the corresponding index point (per the doping levels) on
%the p-side
nIndexCount = linspace(1,length(indexn),length(indexn));
pIndexCount = round(nIndexCount*widthRatio);


index = 1;
flag=true;

while (flag == true)
   
    %First, test the maximum built in voltage at the index points selected
    
    %Calculate voltage drop on n-type side for depletion region using
    %formula V(x) - V(xn) = -q*Nd/(2*epsilon)*(x-xn)^2
    %1/100 term is to convert length in cm to m
    Vjn = -q*nDopingLevel/2/(e0*er)*(nIndexCount(index)*intervalLength/100)^2;
    
    %Calculate voltage drop on p-type side for depletion region using
    %formula V(xp) - V(x) = -q*Na/(2*epsilon)*(xp-x)^2
    %1/100 term is to convert length in cm to m
    Vjp = -q*pDopingLevel/2/(e0*er)*(pIndexCount(index)*intervalLength/100)^2;
    
    Vtotal_doping = abs(Vjp + Vjn);
    
    %Find the actual Vbi at the index points tested
    Vtotal_bands = pEc2(pIndexCount(index)) - nEc2(length(nEc2)-nIndexCount(index)+1);
    
    %Check to see if the value of the dopant Vbi is greater than the actual
    %Vbi at those points.  This corresponds to the equality point between
    %the terms.
    if(Vtotal_doping >Vtotal_bands)
        Vout = Vtotal_doping;
        flag = false;

    end
    
    %Iterate to the next index value
    index = index +1;
    
    if(index == length(nEc2))
        Vout = 0;
        flag = false;
    end
    

end
        
Vtotal_doping;
Vtotal_bands;

%Undo last iteration
index = index-1;

%Fix edges of depletion region
nIndex = index;
pIndex = index;
%Set initial conditions before adding in depletion region
Ec_out = EcTotal2;


%Endpoints of depletion region
startPoint = (length(nEc2)-nIndex)*intervalLength;
endPoint = (length(nEc2)+pIndex-1)*intervalLength;

%Total width of depletion region
totalWidth = endPoint-startPoint;
%Vtotalbands is total voltage across depletion region

%Calculate width on each side of junction using wn/wp = Na/Nd
nWidth = totalWidth/(1/widthRatio + 1);
pWidth = totalWidth - nWidth;

%Calculate voltage on each side of junction using Vj,n/Vj,p = Na/Nd
nVoltage = Vtotal_bands/(1/widthRatio+1);
pVoltage = Vtotal_bands - nVoltage;

%Energy level at beginning of depletion region
cnStartVoltage = Ec_out(length(nEc2)-nIndex+1);
cpStartVoltage = Ec_out(length(nEc2)+pIndex);


%Define conduction band curve in depletion region on n-side using
%V(x)=-q Nd/(2 epsilon) (x-xn)^2 + V(xn)

for loopIndex = length(nEc2)-nIndex:1:length(nEc2)
   Ec_out(loopIndex) = q*nDopingLevel/2/(e0*er)*(((length(nEc2)-nIndex)-loopIndex)*intervalLength/100)^2 +cnStartVoltage;
end

for loopIndex = length(nEc2):1:length(nEc2)+pIndex-1
  Ec_out(loopIndex+1) = cpStartVoltage-q*nDopingLevel/2/(e0*er)*((pIndex-(loopIndex-length(nEc2)))*intervalLength/100)^2;
end



Ev_out = EvTotal2;

vnStartVoltage = Ev_out(length(nEc2)-nIndex+1);
vpStartVoltage = Ev_out(length(nEc2)+pIndex);

Vjvn = zeros(length(Ev_out),1);
Vjvp = zeros(length(Ev_out),1);

for loopIndex = length(nEv2)-nIndex:1:length(nEv2)
   Ev_out(loopIndex) = vnStartVoltage+q*nDopingLevel/2/(e0*er)*(((length(nEv2)-nIndex)-nIndexCount(loopIndex))*intervalLength/100)^2;
   Vjvn(loopIndex) = q*nDopingLevel/2/(e0*er)*(((length(nEv2)-nIndex)-nIndexCount(loopIndex))*intervalLength/100)^2;
end

for loopIndex = length(nEv2):1:length(nEv2)+pIndex
  Ev_out(loopIndex+1) = vpStartVoltage-q*nDopingLevel/2/(e0*er)*((pIndex-(loopIndex-length(nEv2)))*intervalLength/100)^2;
  Vjvp(loopIndex+1) =q*nDopingLevel/2/(e0*er)*((pIndex-(loopIndex-length(nEv2)))*intervalLength/100)^2;
end


figure
plot(totalIndex,Ec_out,totalIndex,Ev_out,totalIndex,EvTotal2,totalIndex,EcTotal2),xlim([0,2*Wn]),title('p-n junction band structure')
















