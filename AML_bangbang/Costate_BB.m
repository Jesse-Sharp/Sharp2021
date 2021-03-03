% Author: Jesse Sharp; Last Update: 03/03/2021
% 
function f = Costate_BB(Lambda,y,U,parameters)

ps = parameters.ps; 
pa = parameters.pa; 
pl = parameters.pl; 
gs = parameters.gs; 
ga = parameters.ga; 
gl = parameters.gl; 
ux = parameters.ux; 
ut = parameters.ut; 
Alpha = parameters.Alpha; 
Gamma = parameters.Gamma; 
a2 = parameters.a2; 

f = [-(-2*y(1)*Lambda(1)*ps-gs*Lambda(1)+gs*Lambda(2)+Lambda(1)*ps),-(-2*y(2)*Lambda(2)*pa-y(3)*Lambda(2)*pa-y(3)*Lambda(3)*pl-ga*Lambda(2)+Lambda(2)*pa),-(a2-pa*y(2)*Lambda(2)-Lambda(3)*pl*y(2)-2*pl*y(3)*Lambda(3)+Lambda(3)*pl-Lambda(3)*gl-Lambda(3)*Alpha/(Gamma+y(3))+Lambda(3)*Alpha*y(3)/(Gamma+y(3))^2-U*Lambda(3))];