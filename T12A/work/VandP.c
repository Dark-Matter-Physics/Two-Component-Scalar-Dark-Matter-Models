#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../CalcHEP_src/include/extern.h"
#include "../../CalcHEP_src/include/VandP.h"
#include "autoprot.h"
extern int  FError;
/*  Special model functions  */

int nModelParticles=25;
static ModelPrtclsStr ModelPrtcls_[25]=
{
  {"~Phip","~Phim",0, 2004, "MPhip","WPhip",0,1,1,3}
, {"h","h",1, 25, "Mh","Wh",0,1,1,0}
, {"~A0","~A0",1, 2003, "MA0","WA0",0,1,1,0}
, {"~P01","~P01",1, 2001, "MP01","WP01",0,1,1,0}
, {"~P02","~P02",1, 2002, "MP02","WP02",0,1,1,0}
, {"g","g",1, 21, "0","0",2,8,16,0}
, {"A","A",1, 22, "0","0",2,1,2,0}
, {"Z","Z",1, 23, "MZ","WZ",2,1,3,0}
, {"Wp","Wm",0, 24, "MWp","WWp",2,1,3,3}
, {"~FXm","~fXm",0, 3004, "MFXm","WFXm",1,1,2,-3}
, {"d1","D1",0, 1, "Md1","Wd1",1,3,6,-1}
, {"d2","D2",0, 3, "Md2","Wd2",1,3,6,-1}
, {"d3","D3",0, 5, "Md3","Wd3",1,3,6,-1}
, {"u1","U1",0, 2, "Mu1","Wu1",1,3,6,2}
, {"u2","U2",0, 4, "Mu2","Wu2",1,3,6,2}
, {"u3","U3",0, 6, "Mu3","Wu3",1,3,6,2}
, {"e1","E1",0, 11, "Me1","We1",1,1,2,-3}
, {"e2","E2",0, 13, "Me2","We2",1,1,2,-3}
, {"e3","E3",0, 15, "Me3","We3",1,1,2,-3}
, {"nu1","nu1",1, 12, "Mnu1","Wnu1",1,1,2,0}
, {"nu2","nu2",1, 14, "Mnu2","Wnu2",1,1,2,0}
, {"nu3","nu3",1, 16, "Mnu3","Wnu3",1,1,2,0}
, {"~FX01","~FX01",1, 3001, "MFX01","WFX01",1,1,2,0}
, {"~FX02","~FX02",1, 3002, "MFX02","WFX02",1,1,2,0}
, {"~FX03","~FX03",1, 3003, "MFX03","WFX03",1,1,2,0}
};
ModelPrtclsStr *ModelPrtcls=ModelPrtcls_; 
int nModelVars=22;
int nModelFunc=281;
static int nCurrentVars=21;
int*currentVarPtr=&nCurrentVars;
static char*varNames_[303]={
 "WZ","WWp","Wd1","Wd2","Wd3","Wu1","Wu2","Wu3","We1","We2"
,"We3","Wnu1","Wnu2","Wnu3","Maux","sqrt2","Pi","Q","alfSMZ","aS"
,"aEWinv","Gf","rd","MPhip","Mh","MA0","MP01","MP02","MZ","MFXm"
,"Md1","Md2","Md3","Mu1","Mu2","Mu3","Me1","Me2","Me3","Mnu1"
,"Mnu2","Mnu3","MFX01","MFX02","MFX03","Rl4S","Il4S","RlSH","IlSH","RlSP"
,"IlSP","RlPp","IlPp","RlPhi","IlPhi","Rl4P","Il4P","RlPpp","IlPpp","RcTri"
,"IcTri","Ryk1","Iyk1","Ryk2","Iyk2","Rcg21","Icg21","Rcg22","Icg22","Rcg23"
,"Icg23","Rcg31","Icg31","Rcg32","Icg32","Rcg33","Icg33","Rcg11","Icg11","Rcg12"
,"Icg12","Rcg13","Icg13","RZDL11","IZDL11","RZDL12","IZDL12","RZDL13","IZDL13","RZDL21"
,"IZDL21","RZDL22","IZDL22","RZDL23","IZDL23","RZDL31","IZDL31","RZDL32","IZDL32","RZDL33"
,"IZDL33","RZDR11","IZDR11","RZDR12","IZDR12","RZDR13","IZDR13","RZDR21","IZDR21","RZDR22"
,"IZDR22","RZDR23","IZDR23","RZDR31","IZDR31","RZDR32","IZDR32","RZDR33","IZDR33","RZUL11"
,"IZUL11","RZUL12","IZUL12","RZUL13","IZUL13","RZUL21","IZUL21","RZUL22","IZUL22","RZUL23"
,"IZUL23","RZUL31","IZUL31","RZUL32","IZUL32","RZUL33","IZUL33","RZUR11","IZUR11","RZUR12"
,"IZUR12","RZUR13","IZUR13","RZUR21","IZUR21","RZUR22","IZUR22","RZUR23","IZUR23","RZUR31"
,"IZUR31","RZUR32","IZUR32","RZUR33","IZUR33","RZEL11","IZEL11","RZEL12","IZEL12","RZEL13"
,"IZEL13","RZEL21","IZEL21","RZEL22","IZEL22","RZEL23","IZEL23","RZEL31","IZEL31","RZEL32"
,"IZEL32","RZEL33","IZEL33","RZER11","IZER11","RZER12","IZER12","RZER13","IZER13","RZER21"
,"IZER21","RZER22","IZER22","RZER23","IZER23","RZER31","IZER31","RZER32","IZER32","RZER33"
,"IZER33","RUV11","IUV11","RUV12","IUV12","RUV13","IUV13","RUV21","IUV21","RUV22"
,"IUV22","RUV23","IUV23","RUV31","IUV31","RUV32","IUV32","RUV33","IUV33","RZS11"
,"IZS11","RZS12","IZS12","RZS21","IZS21","RZS22","IZS22","RZX11","IZX11","RZX12"
,"IZX12","RZX13","IZX13","RZX21","IZX21","RZX22","IZX22","RZX23","IZX23","RZX31"
,"IZX31","RZX32","IZX32","RZX33","IZX33","HPP1","HGG1","QCDok","g3","el"
,"MWp","TW","STW","CTW","TTW","g1","g2","vvSM","RLam","ILam"
,"RYd11","IYd11","RYd12","IYd12","RYd13","IYd13","RYd21","IYd21","RYd22","IYd22"
,"RYd23","IYd23","RYd31","IYd31","RYd32","IYd32","RYd33","IYd33","RYe11","IYe11"
,"RYe12","IYe12","RYe13","IYe13","RYe21","IYe21","RYe22","IYe22","RYe23","IYe23"
,"RYe31","IYe31","RYe32","IYe32","RYe33","IYe33","RYu11","IYu11","RYu12","IYu12"
,"RYu13","IYu13","RYu21","IYu21","RYu22","IYu22","RYu23","IYu23","RYu31","IYu31"
,"RYu32","IYu32","RYu33"};
char**varNames=varNames_;
static REAL varValues_[303]={
   2.495200E+00,  2.141000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  1.510000E+00,  0.000000E+00,  0.000000E+00
,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  1.000000E+00,  1.414214E+00,  3.141593E+00,  1.000000E+02,  1.172000E-01,  1.190000E-01
,  1.370360E+02,  1.166390E-05};
REAL*varValues=varValues_;
int calcMainFunc(void)
{
   int i;
   static REAL * VV=NULL;
   static int iQ=-1;
   static int cErr=1;
   REAL *V=varValues;
   FError=0;
   if(VV && cErr==0)
   { for(i=0;i<nModelVars;i++) if(i!=iQ && VV[i]!=V[i]) break;
     if(i==nModelVars)      {if(iQ>=0 && VV[iQ]!=V[iQ]) goto FirstQ; else return 0;} 
   }
  cErr=1;
   nCurrentVars=22;
   V[22]=slhaRead("SPheno.spc.T12AMaud",0);
   if(!isfinite(V[22]) || FError) return 22;
 FirstQ:
 cErr=1;
   nCurrentVars=23;
   V[23]=slhaVal("MASS",V[17],1,2004);
   if(!isfinite(V[23]) || FError) return 23;
   nCurrentVars=24;
   V[24]=slhaVal("MASS",V[17],1,25);
   if(!isfinite(V[24]) || FError) return 24;
   nCurrentVars=25;
   V[25]=slhaVal("MASS",V[17],1,2003);
   if(!isfinite(V[25]) || FError) return 25;
   nCurrentVars=26;
   V[26]=slhaVal("MASS",V[17],1,2001);
   if(!isfinite(V[26]) || FError) return 26;
   nCurrentVars=27;
   V[27]=slhaVal("MASS",V[17],1,2002);
   if(!isfinite(V[27]) || FError) return 27;
   nCurrentVars=28;
   V[28]=slhaVal("MASS",V[17],1,23);
   if(!isfinite(V[28]) || FError) return 28;
   nCurrentVars=29;
   V[29]=slhaVal("MASS",V[17],1,3004);
   if(!isfinite(V[29]) || FError) return 29;
   nCurrentVars=30;
   V[30]=slhaVal("MASS",V[17],1,1);
   if(!isfinite(V[30]) || FError) return 30;
   nCurrentVars=31;
   V[31]=slhaVal("MASS",V[17],1,3);
   if(!isfinite(V[31]) || FError) return 31;
   nCurrentVars=32;
   V[32]=slhaVal("MASS",V[17],1,5);
   if(!isfinite(V[32]) || FError) return 32;
   nCurrentVars=33;
   V[33]=slhaVal("MASS",V[17],1,2);
   if(!isfinite(V[33]) || FError) return 33;
   nCurrentVars=34;
   V[34]=slhaVal("MASS",V[17],1,4);
   if(!isfinite(V[34]) || FError) return 34;
   nCurrentVars=35;
   V[35]=slhaVal("MASS",V[17],1,6);
   if(!isfinite(V[35]) || FError) return 35;
   nCurrentVars=36;
   V[36]=slhaVal("MASS",V[17],1,11);
   if(!isfinite(V[36]) || FError) return 36;
   nCurrentVars=37;
   V[37]=slhaVal("MASS",V[17],1,13);
   if(!isfinite(V[37]) || FError) return 37;
   nCurrentVars=38;
   V[38]=slhaVal("MASS",V[17],1,15);
   if(!isfinite(V[38]) || FError) return 38;
   nCurrentVars=39;
   V[39]=slhaVal("MASS",V[17],1,12);
   if(!isfinite(V[39]) || FError) return 39;
   nCurrentVars=40;
   V[40]=slhaVal("MASS",V[17],1,14);
   if(!isfinite(V[40]) || FError) return 40;
   nCurrentVars=41;
   V[41]=slhaVal("MASS",V[17],1,16);
   if(!isfinite(V[41]) || FError) return 41;
   nCurrentVars=42;
   V[42]=slhaVal("MASS",V[17],1,3001);
   if(!isfinite(V[42]) || FError) return 42;
   nCurrentVars=43;
   V[43]=slhaVal("MASS",V[17],1,3002);
   if(!isfinite(V[43]) || FError) return 43;
   nCurrentVars=44;
   V[44]=slhaVal("MASS",V[17],1,3003);
   if(!isfinite(V[44]) || FError) return 44;
   nCurrentVars=45;
   V[45]=slhaVal("HDM",V[17],1,3);
   if(!isfinite(V[45]) || FError) return 45;
   nCurrentVars=46;
   V[46]=slhaVal("IMHDM",V[17],1,3);
   if(!isfinite(V[46]) || FError) return 46;
   nCurrentVars=47;
   V[47]=slhaVal("HDM",V[17],1,2);
   if(!isfinite(V[47]) || FError) return 47;
   nCurrentVars=48;
   V[48]=slhaVal("IMHDM",V[17],1,2);
   if(!isfinite(V[48]) || FError) return 48;
   nCurrentVars=49;
   V[49]=slhaVal("HDM",V[17],1,15);
   if(!isfinite(V[49]) || FError) return 49;
   nCurrentVars=50;
   V[50]=slhaVal("IMHDM",V[17],1,15);
   if(!isfinite(V[50]) || FError) return 50;
   nCurrentVars=51;
   V[51]=slhaVal("HDM",V[17],1,8);
   if(!isfinite(V[51]) || FError) return 51;
   nCurrentVars=52;
   V[52]=slhaVal("IMHDM",V[17],1,8);
   if(!isfinite(V[52]) || FError) return 52;
   nCurrentVars=53;
   V[53]=slhaVal("HDM",V[17],1,7);
   if(!isfinite(V[53]) || FError) return 53;
   nCurrentVars=54;
   V[54]=slhaVal("IMHDM",V[17],1,7);
   if(!isfinite(V[54]) || FError) return 54;
   nCurrentVars=55;
   V[55]=slhaVal("HDM",V[17],1,6);
   if(!isfinite(V[55]) || FError) return 55;
   nCurrentVars=56;
   V[56]=slhaVal("IMHDM",V[17],1,6);
   if(!isfinite(V[56]) || FError) return 56;
   nCurrentVars=57;
   V[57]=slhaVal("HDM",V[17],1,9);
   if(!isfinite(V[57]) || FError) return 57;
   nCurrentVars=58;
   V[58]=slhaVal("IMHDM",V[17],1,9);
   if(!isfinite(V[58]) || FError) return 58;
   nCurrentVars=59;
   V[59]=slhaVal("HDM",V[17],1,10);
   if(!isfinite(V[59]) || FError) return 59;
   nCurrentVars=60;
   V[60]=slhaVal("IMHDM",V[17],1,10);
   if(!isfinite(V[60]) || FError) return 60;
   nCurrentVars=61;
   V[61]=slhaVal("HDM",V[17],1,13);
   if(!isfinite(V[61]) || FError) return 61;
   nCurrentVars=62;
   V[62]=slhaVal("IMHDM",V[17],1,13);
   if(!isfinite(V[62]) || FError) return 62;
   nCurrentVars=63;
   V[63]=slhaVal("HDM",V[17],1,14);
   if(!isfinite(V[63]) || FError) return 63;
   nCurrentVars=64;
   V[64]=slhaVal("IMHDM",V[17],1,14);
   if(!isfinite(V[64]) || FError) return 64;
   nCurrentVars=65;
   V[65]=slhaVal("COUPLING2",V[17],1,1);
   if(!isfinite(V[65]) || FError) return 65;
   nCurrentVars=66;
   V[66]=slhaVal("IMCOUPLING2",V[17],1,1);
   if(!isfinite(V[66]) || FError) return 66;
   nCurrentVars=67;
   V[67]=slhaVal("COUPLING2",V[17],1,2);
   if(!isfinite(V[67]) || FError) return 67;
   nCurrentVars=68;
   V[68]=slhaVal("IMCOUPLING2",V[17],1,2);
   if(!isfinite(V[68]) || FError) return 68;
   nCurrentVars=69;
   V[69]=slhaVal("COUPLING2",V[17],1,3);
   if(!isfinite(V[69]) || FError) return 69;
   nCurrentVars=70;
   V[70]=slhaVal("IMCOUPLING2",V[17],1,3);
   if(!isfinite(V[70]) || FError) return 70;
   nCurrentVars=71;
   V[71]=slhaVal("COUPLING3",V[17],1,1);
   if(!isfinite(V[71]) || FError) return 71;
   nCurrentVars=72;
   V[72]=slhaVal("IMCOUPLING3",V[17],1,1);
   if(!isfinite(V[72]) || FError) return 72;
   nCurrentVars=73;
   V[73]=slhaVal("COUPLING3",V[17],1,2);
   if(!isfinite(V[73]) || FError) return 73;
   nCurrentVars=74;
   V[74]=slhaVal("IMCOUPLING3",V[17],1,2);
   if(!isfinite(V[74]) || FError) return 74;
   nCurrentVars=75;
   V[75]=slhaVal("COUPLING3",V[17],1,3);
   if(!isfinite(V[75]) || FError) return 75;
   nCurrentVars=76;
   V[76]=slhaVal("IMCOUPLING3",V[17],1,3);
   if(!isfinite(V[76]) || FError) return 76;
   nCurrentVars=77;
   V[77]=slhaVal("COUPLING1",V[17],1,1);
   if(!isfinite(V[77]) || FError) return 77;
   nCurrentVars=78;
   V[78]=slhaVal("IMCOUPLING1",V[17],1,1);
   if(!isfinite(V[78]) || FError) return 78;
   nCurrentVars=79;
   V[79]=slhaVal("COUPLING1",V[17],1,2);
   if(!isfinite(V[79]) || FError) return 79;
   nCurrentVars=80;
   V[80]=slhaVal("IMCOUPLING1",V[17],1,2);
   if(!isfinite(V[80]) || FError) return 80;
   nCurrentVars=81;
   V[81]=slhaVal("COUPLING1",V[17],1,3);
   if(!isfinite(V[81]) || FError) return 81;
   nCurrentVars=82;
   V[82]=slhaVal("IMCOUPLING1",V[17],1,3);
   if(!isfinite(V[82]) || FError) return 82;
   nCurrentVars=83;
   V[83]=slhaVal("UDLMIX",V[17],2,1,1);
   if(!isfinite(V[83]) || FError) return 83;
   nCurrentVars=84;
   V[84]=slhaVal("IMUDLMIX",V[17],2,1,1);
   if(!isfinite(V[84]) || FError) return 84;
   nCurrentVars=85;
   V[85]=slhaVal("UDLMIX",V[17],2,1,2);
   if(!isfinite(V[85]) || FError) return 85;
   nCurrentVars=86;
   V[86]=slhaVal("IMUDLMIX",V[17],2,1,2);
   if(!isfinite(V[86]) || FError) return 86;
   nCurrentVars=87;
   V[87]=slhaVal("UDLMIX",V[17],2,1,3);
   if(!isfinite(V[87]) || FError) return 87;
   nCurrentVars=88;
   V[88]=slhaVal("IMUDLMIX",V[17],2,1,3);
   if(!isfinite(V[88]) || FError) return 88;
   nCurrentVars=89;
   V[89]=slhaVal("UDLMIX",V[17],2,2,1);
   if(!isfinite(V[89]) || FError) return 89;
   nCurrentVars=90;
   V[90]=slhaVal("IMUDLMIX",V[17],2,2,1);
   if(!isfinite(V[90]) || FError) return 90;
   nCurrentVars=91;
   V[91]=slhaVal("UDLMIX",V[17],2,2,2);
   if(!isfinite(V[91]) || FError) return 91;
   nCurrentVars=92;
   V[92]=slhaVal("IMUDLMIX",V[17],2,2,2);
   if(!isfinite(V[92]) || FError) return 92;
   nCurrentVars=93;
   V[93]=slhaVal("UDLMIX",V[17],2,2,3);
   if(!isfinite(V[93]) || FError) return 93;
   nCurrentVars=94;
   V[94]=slhaVal("IMUDLMIX",V[17],2,2,3);
   if(!isfinite(V[94]) || FError) return 94;
   nCurrentVars=95;
   V[95]=slhaVal("UDLMIX",V[17],2,3,1);
   if(!isfinite(V[95]) || FError) return 95;
   nCurrentVars=96;
   V[96]=slhaVal("IMUDLMIX",V[17],2,3,1);
   if(!isfinite(V[96]) || FError) return 96;
   nCurrentVars=97;
   V[97]=slhaVal("UDLMIX",V[17],2,3,2);
   if(!isfinite(V[97]) || FError) return 97;
   nCurrentVars=98;
   V[98]=slhaVal("IMUDLMIX",V[17],2,3,2);
   if(!isfinite(V[98]) || FError) return 98;
   nCurrentVars=99;
   V[99]=slhaVal("UDLMIX",V[17],2,3,3);
   if(!isfinite(V[99]) || FError) return 99;
   nCurrentVars=100;
   V[100]=slhaVal("IMUDLMIX",V[17],2,3,3);
   if(!isfinite(V[100]) || FError) return 100;
   nCurrentVars=101;
   V[101]=slhaVal("UDRMIX",V[17],2,1,1);
   if(!isfinite(V[101]) || FError) return 101;
   nCurrentVars=102;
   V[102]=slhaVal("IMUDRMIX",V[17],2,1,1);
   if(!isfinite(V[102]) || FError) return 102;
   nCurrentVars=103;
   V[103]=slhaVal("UDRMIX",V[17],2,1,2);
   if(!isfinite(V[103]) || FError) return 103;
   nCurrentVars=104;
   V[104]=slhaVal("IMUDRMIX",V[17],2,1,2);
   if(!isfinite(V[104]) || FError) return 104;
   nCurrentVars=105;
   V[105]=slhaVal("UDRMIX",V[17],2,1,3);
   if(!isfinite(V[105]) || FError) return 105;
   nCurrentVars=106;
   V[106]=slhaVal("IMUDRMIX",V[17],2,1,3);
   if(!isfinite(V[106]) || FError) return 106;
   nCurrentVars=107;
   V[107]=slhaVal("UDRMIX",V[17],2,2,1);
   if(!isfinite(V[107]) || FError) return 107;
   nCurrentVars=108;
   V[108]=slhaVal("IMUDRMIX",V[17],2,2,1);
   if(!isfinite(V[108]) || FError) return 108;
   nCurrentVars=109;
   V[109]=slhaVal("UDRMIX",V[17],2,2,2);
   if(!isfinite(V[109]) || FError) return 109;
   nCurrentVars=110;
   V[110]=slhaVal("IMUDRMIX",V[17],2,2,2);
   if(!isfinite(V[110]) || FError) return 110;
   nCurrentVars=111;
   V[111]=slhaVal("UDRMIX",V[17],2,2,3);
   if(!isfinite(V[111]) || FError) return 111;
   nCurrentVars=112;
   V[112]=slhaVal("IMUDRMIX",V[17],2,2,3);
   if(!isfinite(V[112]) || FError) return 112;
   nCurrentVars=113;
   V[113]=slhaVal("UDRMIX",V[17],2,3,1);
   if(!isfinite(V[113]) || FError) return 113;
   nCurrentVars=114;
   V[114]=slhaVal("IMUDRMIX",V[17],2,3,1);
   if(!isfinite(V[114]) || FError) return 114;
   nCurrentVars=115;
   V[115]=slhaVal("UDRMIX",V[17],2,3,2);
   if(!isfinite(V[115]) || FError) return 115;
   nCurrentVars=116;
   V[116]=slhaVal("IMUDRMIX",V[17],2,3,2);
   if(!isfinite(V[116]) || FError) return 116;
   nCurrentVars=117;
   V[117]=slhaVal("UDRMIX",V[17],2,3,3);
   if(!isfinite(V[117]) || FError) return 117;
   nCurrentVars=118;
   V[118]=slhaVal("IMUDRMIX",V[17],2,3,3);
   if(!isfinite(V[118]) || FError) return 118;
   nCurrentVars=119;
   V[119]=slhaVal("UULMIX",V[17],2,1,1);
   if(!isfinite(V[119]) || FError) return 119;
   nCurrentVars=120;
   V[120]=slhaVal("IMUULMIX",V[17],2,1,1);
   if(!isfinite(V[120]) || FError) return 120;
   nCurrentVars=121;
   V[121]=slhaVal("UULMIX",V[17],2,1,2);
   if(!isfinite(V[121]) || FError) return 121;
   nCurrentVars=122;
   V[122]=slhaVal("IMUULMIX",V[17],2,1,2);
   if(!isfinite(V[122]) || FError) return 122;
   nCurrentVars=123;
   V[123]=slhaVal("UULMIX",V[17],2,1,3);
   if(!isfinite(V[123]) || FError) return 123;
   nCurrentVars=124;
   V[124]=slhaVal("IMUULMIX",V[17],2,1,3);
   if(!isfinite(V[124]) || FError) return 124;
   nCurrentVars=125;
   V[125]=slhaVal("UULMIX",V[17],2,2,1);
   if(!isfinite(V[125]) || FError) return 125;
   nCurrentVars=126;
   V[126]=slhaVal("IMUULMIX",V[17],2,2,1);
   if(!isfinite(V[126]) || FError) return 126;
   nCurrentVars=127;
   V[127]=slhaVal("UULMIX",V[17],2,2,2);
   if(!isfinite(V[127]) || FError) return 127;
   nCurrentVars=128;
   V[128]=slhaVal("IMUULMIX",V[17],2,2,2);
   if(!isfinite(V[128]) || FError) return 128;
   nCurrentVars=129;
   V[129]=slhaVal("UULMIX",V[17],2,2,3);
   if(!isfinite(V[129]) || FError) return 129;
   nCurrentVars=130;
   V[130]=slhaVal("IMUULMIX",V[17],2,2,3);
   if(!isfinite(V[130]) || FError) return 130;
   nCurrentVars=131;
   V[131]=slhaVal("UULMIX",V[17],2,3,1);
   if(!isfinite(V[131]) || FError) return 131;
   nCurrentVars=132;
   V[132]=slhaVal("IMUULMIX",V[17],2,3,1);
   if(!isfinite(V[132]) || FError) return 132;
   nCurrentVars=133;
   V[133]=slhaVal("UULMIX",V[17],2,3,2);
   if(!isfinite(V[133]) || FError) return 133;
   nCurrentVars=134;
   V[134]=slhaVal("IMUULMIX",V[17],2,3,2);
   if(!isfinite(V[134]) || FError) return 134;
   nCurrentVars=135;
   V[135]=slhaVal("UULMIX",V[17],2,3,3);
   if(!isfinite(V[135]) || FError) return 135;
   nCurrentVars=136;
   V[136]=slhaVal("IMUULMIX",V[17],2,3,3);
   if(!isfinite(V[136]) || FError) return 136;
   nCurrentVars=137;
   V[137]=slhaVal("UURMIX",V[17],2,1,1);
   if(!isfinite(V[137]) || FError) return 137;
   nCurrentVars=138;
   V[138]=slhaVal("IMUURMIX",V[17],2,1,1);
   if(!isfinite(V[138]) || FError) return 138;
   nCurrentVars=139;
   V[139]=slhaVal("UURMIX",V[17],2,1,2);
   if(!isfinite(V[139]) || FError) return 139;
   nCurrentVars=140;
   V[140]=slhaVal("IMUURMIX",V[17],2,1,2);
   if(!isfinite(V[140]) || FError) return 140;
   nCurrentVars=141;
   V[141]=slhaVal("UURMIX",V[17],2,1,3);
   if(!isfinite(V[141]) || FError) return 141;
   nCurrentVars=142;
   V[142]=slhaVal("IMUURMIX",V[17],2,1,3);
   if(!isfinite(V[142]) || FError) return 142;
   nCurrentVars=143;
   V[143]=slhaVal("UURMIX",V[17],2,2,1);
   if(!isfinite(V[143]) || FError) return 143;
   nCurrentVars=144;
   V[144]=slhaVal("IMUURMIX",V[17],2,2,1);
   if(!isfinite(V[144]) || FError) return 144;
   nCurrentVars=145;
   V[145]=slhaVal("UURMIX",V[17],2,2,2);
   if(!isfinite(V[145]) || FError) return 145;
   nCurrentVars=146;
   V[146]=slhaVal("IMUURMIX",V[17],2,2,2);
   if(!isfinite(V[146]) || FError) return 146;
   nCurrentVars=147;
   V[147]=slhaVal("UURMIX",V[17],2,2,3);
   if(!isfinite(V[147]) || FError) return 147;
   nCurrentVars=148;
   V[148]=slhaVal("IMUURMIX",V[17],2,2,3);
   if(!isfinite(V[148]) || FError) return 148;
   nCurrentVars=149;
   V[149]=slhaVal("UURMIX",V[17],2,3,1);
   if(!isfinite(V[149]) || FError) return 149;
   nCurrentVars=150;
   V[150]=slhaVal("IMUURMIX",V[17],2,3,1);
   if(!isfinite(V[150]) || FError) return 150;
   nCurrentVars=151;
   V[151]=slhaVal("UURMIX",V[17],2,3,2);
   if(!isfinite(V[151]) || FError) return 151;
   nCurrentVars=152;
   V[152]=slhaVal("IMUURMIX",V[17],2,3,2);
   if(!isfinite(V[152]) || FError) return 152;
   nCurrentVars=153;
   V[153]=slhaVal("UURMIX",V[17],2,3,3);
   if(!isfinite(V[153]) || FError) return 153;
   nCurrentVars=154;
   V[154]=slhaVal("IMUURMIX",V[17],2,3,3);
   if(!isfinite(V[154]) || FError) return 154;
   nCurrentVars=155;
   V[155]=slhaVal("UELMIX",V[17],2,1,1);
   if(!isfinite(V[155]) || FError) return 155;
   nCurrentVars=156;
   V[156]=slhaVal("IMUELMIX",V[17],2,1,1);
   if(!isfinite(V[156]) || FError) return 156;
   nCurrentVars=157;
   V[157]=slhaVal("UELMIX",V[17],2,1,2);
   if(!isfinite(V[157]) || FError) return 157;
   nCurrentVars=158;
   V[158]=slhaVal("IMUELMIX",V[17],2,1,2);
   if(!isfinite(V[158]) || FError) return 158;
   nCurrentVars=159;
   V[159]=slhaVal("UELMIX",V[17],2,1,3);
   if(!isfinite(V[159]) || FError) return 159;
   nCurrentVars=160;
   V[160]=slhaVal("IMUELMIX",V[17],2,1,3);
   if(!isfinite(V[160]) || FError) return 160;
   nCurrentVars=161;
   V[161]=slhaVal("UELMIX",V[17],2,2,1);
   if(!isfinite(V[161]) || FError) return 161;
   nCurrentVars=162;
   V[162]=slhaVal("IMUELMIX",V[17],2,2,1);
   if(!isfinite(V[162]) || FError) return 162;
   nCurrentVars=163;
   V[163]=slhaVal("UELMIX",V[17],2,2,2);
   if(!isfinite(V[163]) || FError) return 163;
   nCurrentVars=164;
   V[164]=slhaVal("IMUELMIX",V[17],2,2,2);
   if(!isfinite(V[164]) || FError) return 164;
   nCurrentVars=165;
   V[165]=slhaVal("UELMIX",V[17],2,2,3);
   if(!isfinite(V[165]) || FError) return 165;
   nCurrentVars=166;
   V[166]=slhaVal("IMUELMIX",V[17],2,2,3);
   if(!isfinite(V[166]) || FError) return 166;
   nCurrentVars=167;
   V[167]=slhaVal("UELMIX",V[17],2,3,1);
   if(!isfinite(V[167]) || FError) return 167;
   nCurrentVars=168;
   V[168]=slhaVal("IMUELMIX",V[17],2,3,1);
   if(!isfinite(V[168]) || FError) return 168;
   nCurrentVars=169;
   V[169]=slhaVal("UELMIX",V[17],2,3,2);
   if(!isfinite(V[169]) || FError) return 169;
   nCurrentVars=170;
   V[170]=slhaVal("IMUELMIX",V[17],2,3,2);
   if(!isfinite(V[170]) || FError) return 170;
   nCurrentVars=171;
   V[171]=slhaVal("UELMIX",V[17],2,3,3);
   if(!isfinite(V[171]) || FError) return 171;
   nCurrentVars=172;
   V[172]=slhaVal("IMUELMIX",V[17],2,3,3);
   if(!isfinite(V[172]) || FError) return 172;
   nCurrentVars=173;
   V[173]=slhaVal("UERMIX",V[17],2,1,1);
   if(!isfinite(V[173]) || FError) return 173;
   nCurrentVars=174;
   V[174]=slhaVal("IMUERMIX",V[17],2,1,1);
   if(!isfinite(V[174]) || FError) return 174;
   nCurrentVars=175;
   V[175]=slhaVal("UERMIX",V[17],2,1,2);
   if(!isfinite(V[175]) || FError) return 175;
   nCurrentVars=176;
   V[176]=slhaVal("IMUERMIX",V[17],2,1,2);
   if(!isfinite(V[176]) || FError) return 176;
   nCurrentVars=177;
   V[177]=slhaVal("UERMIX",V[17],2,1,3);
   if(!isfinite(V[177]) || FError) return 177;
   nCurrentVars=178;
   V[178]=slhaVal("IMUERMIX",V[17],2,1,3);
   if(!isfinite(V[178]) || FError) return 178;
   nCurrentVars=179;
   V[179]=slhaVal("UERMIX",V[17],2,2,1);
   if(!isfinite(V[179]) || FError) return 179;
   nCurrentVars=180;
   V[180]=slhaVal("IMUERMIX",V[17],2,2,1);
   if(!isfinite(V[180]) || FError) return 180;
   nCurrentVars=181;
   V[181]=slhaVal("UERMIX",V[17],2,2,2);
   if(!isfinite(V[181]) || FError) return 181;
   nCurrentVars=182;
   V[182]=slhaVal("IMUERMIX",V[17],2,2,2);
   if(!isfinite(V[182]) || FError) return 182;
   nCurrentVars=183;
   V[183]=slhaVal("UERMIX",V[17],2,2,3);
   if(!isfinite(V[183]) || FError) return 183;
   nCurrentVars=184;
   V[184]=slhaVal("IMUERMIX",V[17],2,2,3);
   if(!isfinite(V[184]) || FError) return 184;
   nCurrentVars=185;
   V[185]=slhaVal("UERMIX",V[17],2,3,1);
   if(!isfinite(V[185]) || FError) return 185;
   nCurrentVars=186;
   V[186]=slhaVal("IMUERMIX",V[17],2,3,1);
   if(!isfinite(V[186]) || FError) return 186;
   nCurrentVars=187;
   V[187]=slhaVal("UERMIX",V[17],2,3,2);
   if(!isfinite(V[187]) || FError) return 187;
   nCurrentVars=188;
   V[188]=slhaVal("IMUERMIX",V[17],2,3,2);
   if(!isfinite(V[188]) || FError) return 188;
   nCurrentVars=189;
   V[189]=slhaVal("UERMIX",V[17],2,3,3);
   if(!isfinite(V[189]) || FError) return 189;
   nCurrentVars=190;
   V[190]=slhaVal("IMUERMIX",V[17],2,3,3);
   if(!isfinite(V[190]) || FError) return 190;
   nCurrentVars=191;
   V[191]=slhaVal("UVMIX",V[17],2,1,1);
   if(!isfinite(V[191]) || FError) return 191;
   nCurrentVars=192;
   V[192]=slhaVal("IMUVMIX",V[17],2,1,1);
   if(!isfinite(V[192]) || FError) return 192;
   nCurrentVars=193;
   V[193]=slhaVal("UVMIX",V[17],2,1,2);
   if(!isfinite(V[193]) || FError) return 193;
   nCurrentVars=194;
   V[194]=slhaVal("IMUVMIX",V[17],2,1,2);
   if(!isfinite(V[194]) || FError) return 194;
   nCurrentVars=195;
   V[195]=slhaVal("UVMIX",V[17],2,1,3);
   if(!isfinite(V[195]) || FError) return 195;
   nCurrentVars=196;
   V[196]=slhaVal("IMUVMIX",V[17],2,1,3);
   if(!isfinite(V[196]) || FError) return 196;
   nCurrentVars=197;
   V[197]=slhaVal("UVMIX",V[17],2,2,1);
   if(!isfinite(V[197]) || FError) return 197;
   nCurrentVars=198;
   V[198]=slhaVal("IMUVMIX",V[17],2,2,1);
   if(!isfinite(V[198]) || FError) return 198;
   nCurrentVars=199;
   V[199]=slhaVal("UVMIX",V[17],2,2,2);
   if(!isfinite(V[199]) || FError) return 199;
   nCurrentVars=200;
   V[200]=slhaVal("IMUVMIX",V[17],2,2,2);
   if(!isfinite(V[200]) || FError) return 200;
   nCurrentVars=201;
   V[201]=slhaVal("UVMIX",V[17],2,2,3);
   if(!isfinite(V[201]) || FError) return 201;
   nCurrentVars=202;
   V[202]=slhaVal("IMUVMIX",V[17],2,2,3);
   if(!isfinite(V[202]) || FError) return 202;
   nCurrentVars=203;
   V[203]=slhaVal("UVMIX",V[17],2,3,1);
   if(!isfinite(V[203]) || FError) return 203;
   nCurrentVars=204;
   V[204]=slhaVal("IMUVMIX",V[17],2,3,1);
   if(!isfinite(V[204]) || FError) return 204;
   nCurrentVars=205;
   V[205]=slhaVal("UVMIX",V[17],2,3,2);
   if(!isfinite(V[205]) || FError) return 205;
   nCurrentVars=206;
   V[206]=slhaVal("IMUVMIX",V[17],2,3,2);
   if(!isfinite(V[206]) || FError) return 206;
   nCurrentVars=207;
   V[207]=slhaVal("UVMIX",V[17],2,3,3);
   if(!isfinite(V[207]) || FError) return 207;
   nCurrentVars=208;
   V[208]=slhaVal("IMUVMIX",V[17],2,3,3);
   if(!isfinite(V[208]) || FError) return 208;
   nCurrentVars=209;
   V[209]=slhaVal("ZSCALAR",V[17],2,1,1);
   if(!isfinite(V[209]) || FError) return 209;
   nCurrentVars=210;
   V[210]=slhaVal("IMZSCALAR",V[17],2,1,1);
   if(!isfinite(V[210]) || FError) return 210;
   nCurrentVars=211;
   V[211]=slhaVal("ZSCALAR",V[17],2,1,2);
   if(!isfinite(V[211]) || FError) return 211;
   nCurrentVars=212;
   V[212]=slhaVal("IMZSCALAR",V[17],2,1,2);
   if(!isfinite(V[212]) || FError) return 212;
   nCurrentVars=213;
   V[213]=slhaVal("ZSCALAR",V[17],2,2,1);
   if(!isfinite(V[213]) || FError) return 213;
   nCurrentVars=214;
   V[214]=slhaVal("IMZSCALAR",V[17],2,2,1);
   if(!isfinite(V[214]) || FError) return 214;
   nCurrentVars=215;
   V[215]=slhaVal("ZSCALAR",V[17],2,2,2);
   if(!isfinite(V[215]) || FError) return 215;
   nCurrentVars=216;
   V[216]=slhaVal("IMZSCALAR",V[17],2,2,2);
   if(!isfinite(V[216]) || FError) return 216;
   nCurrentVars=217;
   V[217]=slhaVal("ZX",V[17],2,1,1);
   if(!isfinite(V[217]) || FError) return 217;
   nCurrentVars=218;
   V[218]=slhaVal("IMZX",V[17],2,1,1);
   if(!isfinite(V[218]) || FError) return 218;
   nCurrentVars=219;
   V[219]=slhaVal("ZX",V[17],2,1,2);
   if(!isfinite(V[219]) || FError) return 219;
   nCurrentVars=220;
   V[220]=slhaVal("IMZX",V[17],2,1,2);
   if(!isfinite(V[220]) || FError) return 220;
   nCurrentVars=221;
   V[221]=slhaVal("ZX",V[17],2,1,3);
   if(!isfinite(V[221]) || FError) return 221;
   nCurrentVars=222;
   V[222]=slhaVal("IMZX",V[17],2,1,3);
   if(!isfinite(V[222]) || FError) return 222;
   nCurrentVars=223;
   V[223]=slhaVal("ZX",V[17],2,2,1);
   if(!isfinite(V[223]) || FError) return 223;
   nCurrentVars=224;
   V[224]=slhaVal("IMZX",V[17],2,2,1);
   if(!isfinite(V[224]) || FError) return 224;
   nCurrentVars=225;
   V[225]=slhaVal("ZX",V[17],2,2,2);
   if(!isfinite(V[225]) || FError) return 225;
   nCurrentVars=226;
   V[226]=slhaVal("IMZX",V[17],2,2,2);
   if(!isfinite(V[226]) || FError) return 226;
   nCurrentVars=227;
   V[227]=slhaVal("ZX",V[17],2,2,3);
   if(!isfinite(V[227]) || FError) return 227;
   nCurrentVars=228;
   V[228]=slhaVal("IMZX",V[17],2,2,3);
   if(!isfinite(V[228]) || FError) return 228;
   nCurrentVars=229;
   V[229]=slhaVal("ZX",V[17],2,3,1);
   if(!isfinite(V[229]) || FError) return 229;
   nCurrentVars=230;
   V[230]=slhaVal("IMZX",V[17],2,3,1);
   if(!isfinite(V[230]) || FError) return 230;
   nCurrentVars=231;
   V[231]=slhaVal("ZX",V[17],2,3,2);
   if(!isfinite(V[231]) || FError) return 231;
   nCurrentVars=232;
   V[232]=slhaVal("IMZX",V[17],2,3,2);
   if(!isfinite(V[232]) || FError) return 232;
   nCurrentVars=233;
   V[233]=slhaVal("ZX",V[17],2,3,3);
   if(!isfinite(V[233]) || FError) return 233;
   nCurrentVars=234;
   V[234]=slhaVal("IMZX",V[17],2,3,3);
   if(!isfinite(V[234]) || FError) return 234;
   nCurrentVars=235;
   V[235]=slhaVal("EFFHIGGSCOUPLINGS",V[17],3,25,22,22);
   if(!isfinite(V[235]) || FError) return 235;
   nCurrentVars=236;
   V[236]=slhaVal("EFFHIGGSCOUPLINGS",V[17],3,25,21,21);
   if(!isfinite(V[236]) || FError) return 236;
   nCurrentVars=237;
   V[237]=initQCD(V[18],V[34],V[32],V[35]);
   if(!isfinite(V[237]) || FError) return 237;
   nCurrentVars=238;
   V[238]=Sqrt(alphaQCD(V[17])*4*3.1415927)*1;
   if(!isfinite(V[238]) || FError) return 238;
   nCurrentVars=239;
   V[239]=2*Sqrt(1/(V[20]))*Sqrt(V[16]);
   if(!isfinite(V[239]) || FError) return 239;
   nCurrentVars=240;
   V[240]=Sqrt(Pow(V[28],2)/(2.)+Sqrt(Pow(V[28],4)/(4.)-Pow(V[28],2)*V[16]/(V[15]*V[20]*V[21])));
   if(!isfinite(V[240]) || FError) return 240;
   nCurrentVars=241;
   V[241]=Asin(Sqrt(1-Pow(V[240],2)/(Pow(V[28],2))));
   if(!isfinite(V[241]) || FError) return 241;
   nCurrentVars=242;
   V[242]=Sin(V[241]);
   if(!isfinite(V[242]) || FError) return 242;
   nCurrentVars=243;
   V[243]=Cos(V[241]);
   if(!isfinite(V[243]) || FError) return 243;
   nCurrentVars=244;
   V[244]=Tan(V[241]);
   if(!isfinite(V[244]) || FError) return 244;
   nCurrentVars=245;
   V[245]=V[239]*1/(Cos(V[241]));
   if(!isfinite(V[245]) || FError) return 245;
   nCurrentVars=246;
   V[246]=V[239]*1/(Sin(V[241]));
   if(!isfinite(V[246]) || FError) return 246;
   nCurrentVars=247;
   V[247]=2*Sqrt(Pow(V[240],2)/(Pow(V[246],2)));
   if(!isfinite(V[247]) || FError) return 247;
   nCurrentVars=248;
   V[248]=Pow(V[24],2)/(Pow(V[247],2));
   if(!isfinite(V[248]) || FError) return 248;
   nCurrentVars=249;
   V[249]=0.;

   nCurrentVars=250;
   V[250]=V[15]*V[30]/(V[247]);
   if(!isfinite(V[250]) || FError) return 250;
   nCurrentVars=251;
   V[251]=0.;

   nCurrentVars=252;
   V[252]=0;

   nCurrentVars=253;
   V[253]=0.;

   nCurrentVars=254;
   V[254]=0;

   nCurrentVars=255;
   V[255]=0.;

   nCurrentVars=256;
   V[256]=0;

   nCurrentVars=257;
   V[257]=0.;

   nCurrentVars=258;
   V[258]=V[15]*V[31]/(V[247]);
   if(!isfinite(V[258]) || FError) return 258;
   nCurrentVars=259;
   V[259]=0.;

   nCurrentVars=260;
   V[260]=0;

   nCurrentVars=261;
   V[261]=0.;

   nCurrentVars=262;
   V[262]=0;

   nCurrentVars=263;
   V[263]=0.;

   nCurrentVars=264;
   V[264]=0;

   nCurrentVars=265;
   V[265]=0.;

   nCurrentVars=266;
   V[266]=V[15]*MbEff(V[17])*1/(V[247]);
   if(!isfinite(V[266]) || FError) return 266;
   nCurrentVars=267;
   V[267]=0.;

   nCurrentVars=268;
   V[268]=V[15]*V[36]/(V[247]);
   if(!isfinite(V[268]) || FError) return 268;
   nCurrentVars=269;
   V[269]=0.;

   nCurrentVars=270;
   V[270]=0;

   nCurrentVars=271;
   V[271]=0.;

   nCurrentVars=272;
   V[272]=0;

   nCurrentVars=273;
   V[273]=0.;

   nCurrentVars=274;
   V[274]=0;

   nCurrentVars=275;
   V[275]=0.;

   nCurrentVars=276;
   V[276]=V[15]*V[37]/(V[247]);
   if(!isfinite(V[276]) || FError) return 276;
   nCurrentVars=277;
   V[277]=0.;

   nCurrentVars=278;
   V[278]=0;

   nCurrentVars=279;
   V[279]=0.;

   nCurrentVars=280;
   V[280]=0;

   nCurrentVars=281;
   V[281]=0.;

   nCurrentVars=282;
   V[282]=0;

   nCurrentVars=283;
   V[283]=0.;

   nCurrentVars=284;
   V[284]=V[15]*V[38]/(V[247]);
   if(!isfinite(V[284]) || FError) return 284;
   nCurrentVars=285;
   V[285]=0.;

   nCurrentVars=286;
   V[286]=V[15]*V[33]/(V[247]);
   if(!isfinite(V[286]) || FError) return 286;
   nCurrentVars=287;
   V[287]=0.;

   nCurrentVars=288;
   V[288]=0;

   nCurrentVars=289;
   V[289]=0.;

   nCurrentVars=290;
   V[290]=0;

   nCurrentVars=291;
   V[291]=0.;

   nCurrentVars=292;
   V[292]=0;

   nCurrentVars=293;
   V[293]=0.;

   nCurrentVars=294;
   V[294]=V[15]*McEff(V[17])*1/(V[247]);
   if(!isfinite(V[294]) || FError) return 294;
   nCurrentVars=295;
   V[295]=0.;

   nCurrentVars=296;
   V[296]=0;

   nCurrentVars=297;
   V[297]=0.;

   nCurrentVars=298;
   V[298]=0;

   nCurrentVars=299;
   V[299]=0.;

   nCurrentVars=300;
   V[300]=0;

   nCurrentVars=301;
   V[301]=0.;

   nCurrentVars=302;
   V[302]=V[15]*MtEff(V[17])*1/(V[247]);
   if(!isfinite(V[302]) || FError) return 302;
   if(VV==NULL) 
   {  VV=malloc(sizeof(REAL)*nModelVars);
      for(i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0) iQ=i;
   }
   for(i=0;i<nModelVars;i++) VV[i]=V[i];
   cErr=0;
   nCurrentVars++;
   return 0;
}
