%{***************************************************************************}
%{* Nom     : function f2bis pour la resolution d'une equation differentielle
%{* Descri?: via un scha Runge-Kutta
%{* Entree  :
%{* Sortie  : -
%{* Utilise :
%{* A faire : -
%{* Portee  : private
%{* Version :
%{***************************************************************************}
function [varf]=f2bis(sigma,caracPhase,Hz,HzPt,tauOx_1,r1_1,r2,vitesseT)

%constantes adimensionns
thetaPt=(HzPt(2)-HzPt(1))/(Hz(3)-Hz(2))-((Hz(2)-Hz(1))*(HzPt(3)-HzPt(2))/((Hz(3)-Hz(2))^2));
theta=(Hz(2)-Hz(1))/(Hz(3)-Hz(2));
part_therm=(caracPhase(1).alphaT-caracPhase(2).alphaT)*vitesseT; %terme du la dilatation thermique

%calcul de la diff
% varf=-(+tauOx_1*(abs(sigma)^caracPhase(2).N)*sign(sigma)-r1_1*(thetaPt/(theta)^2)*sigma+caracPhase(2).CoeffD*(HzPt(3)-HzPt(2))-tauOx_1*r2*((caracPhase(2).E/(1-caracPhase(2).nu))^(caracPhase(1).N-caracPhase(2).N))*(abs(-sigma/theta)^caracPhase(1).N)*sign(-sigma/theta)-part_therm)/(1+r1_1/theta);%version benoit
varf=-(+tauOx_1*(abs(sigma)^caracPhase(2).N)*sign(sigma)/2-r1_1*(thetaPt/(theta)^2)*sigma+caracPhase(2).CoeffD*(HzPt(3)-HzPt(2))-tauOx_1*r2*((caracPhase(2).E/(1-caracPhase(2).nu))^(caracPhase(1).N-caracPhase(2).N))*(abs(-sigma/theta)^caracPhase(1).N)*sign(-sigma/theta)/2-part_therm)/(1+r1_1/theta);%version zhaojun
