function [ut,utlmst,u]=extractanderase(iau,ibu,ut,utlmst,u)
tempt=ut;
temptlmst=utlmst;
tempu=u;
clear ut utlmst u
ut=tempt(iau:ibu);
utlmst=temptlmst(iau:ibu);
u=tempu(iau:ibu);
end
