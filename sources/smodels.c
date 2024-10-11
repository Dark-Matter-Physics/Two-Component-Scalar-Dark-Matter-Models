#include"micromegas.h"
#include"micromegas_aux.h"
#include "../CalcHEP_src/include/version.h"

static int smH=-1; 

static int newPDG(int pdg)
{ 
   if(smH>=nModelParticles) {if(pdg==25) return 25252525; else return pdg;}
   if(ModelPrtcls[smH].NPDG==pdg) return 25; 
   if(pdg==25) return 11111111;
   return pdg;   
}   

static int newSlhaDecayPrint(char * name, int dVirt, FILE*f)
{
   txtList all;
   int i,dim; 
   int PDG;
   char N[5][P_NAME_SIZE];
   int id[5];
           
   PDG=newPDG(pNum(name));
   if(!PDG) return 0;
   fprintf(f,"DECAY %d  %E  # %s\n",PDG,pWidth(name,&all),name);
   for(;all;all=all->next)
   {  
      char pn[20], *chB,*chE;
      double br;
      
      sscanf(all->txt,"%s",pn);
      sscanf(pn,"%lf",&br);
      chB=strstr(all->txt,"->");
      chB+=2;
      for(dim=0,chE=chB ; chE;dim++, chE=strchr(chE+1,',')) 
      {  sscanf(chE+1,"%[^,]",N[dim]); 
         trim(N[dim]); 
         id[dim]=newPDG(pNum(N[dim]));
      }
      if(dVirt && dim==2 && pMass(name)<= pMass(N[0])+pMass(N[1]))  
      {  
         int v[2],k; 
         for(k=0;k<2;k++) v[k] = (id[k]==23 || abs(id[k])==24);
          
         if(v[0]||v[1])
         { 
            txtList LV;
            if(id[0]!=id[1]) br/=2;   
            for(k=0;k<2;k++) if(v[k])
            { 
               pWidth(N[k],&LV);
               for(;LV;LV=LV->next)
               {  double brV;
                  char* chD=strstr(LV->txt,"->")+2;
                  char name1[20],name2[20];
                  sscanf(chD,"%[^,],%s", name1,name2);
                  trim(name1);
                  sscanf(LV->txt,"%lf",&brV);
                  fprintf(f,"   %e  3  %d  %d  %d # %s,%s->%s\n",br*brV , id[1-k] , newPDG(pNum(name1)),newPDG(pNum(name2)),N[1-k],N[k],chD);
               } 
               if(id[0]==id[1]) break;
            }
            continue;   
         }   
      }     
      
      fprintf(f,"   %s   %d  ",pn,dim);
      for(i=0;i<dim;i++) fprintf(f," %d", id[i] ); 
      chB=strstr(all->txt,"->");
      fprintf(f,"  # %s \n",chB+2);
   } 
   fprintf(f,"\n");
   return PDG;
} 


   

int smodels(int Run, int nf,double csMinFb,  char*fileName, char*version, int wrt) 
{ 
   int SMP[18]={1,2,3,4,5,6, 11,12,13,14,15,16, 21,22,23,24, 111, 211};
   int i,j;
   int nChan=0;
   double PcmMax=6500; // LHC pcm
   double PcmMin=4000;
   if((Run & (LHC8|LHC13)) == 0) 
   { printf("SMODELS: The third parameter has to be either  LHC8 or  LHC13 or  LHC8+LHC13 \n");      
     return 1;
   }  
   if((Run & LHC8) == 0)  PcmMin=6500;
   if((Run & LHC13) == 0) PcmMax=4000;
   
   FILE*f=fopen(fileName,"w");  
   int np=0;
   char**plist=NULL;
   char* gluname=NULL;
   char* phname=NULL;
   char* bname=NULL;
   char* Bname=NULL;
   char* lname=NULL;
   char* Lname=NULL;
   char* wname=NULL;
   char* Wname=NULL;
   char* zname=NULL;
   int  VZdecay_=VZdecay,VWdecay_=VWdecay;
   
   fprintf(f,"####################################################################\n");   
   fprintf(f,"# SLHA input file for SModelS, automatically created by micrOMEGAs #\n");   
   fprintf(f,"####################################################################\n\n");   
    
   for(i=0;i<nModelParticles;i++)
   {
      if(ModelPrtcls[i].NPDG== 21)   gluname=ModelPrtcls[i].name;
      if(ModelPrtcls[i].NPDG== 22)   phname=ModelPrtcls[i].name;
      if(ModelPrtcls[i].NPDG==  5) { bname=ModelPrtcls[i].name;  Bname=ModelPrtcls[i].aname;}
      if(ModelPrtcls[i].NPDG== -5) { bname=ModelPrtcls[i].aname; Bname=ModelPrtcls[i].name; }  
      if(ModelPrtcls[i].NPDG== 15) { lname=ModelPrtcls[i].name;  Lname=ModelPrtcls[i].aname;}
      if(ModelPrtcls[i].NPDG==-15) { Lname=ModelPrtcls[i].aname; lname=ModelPrtcls[i].name; }
      if(ModelPrtcls[i].NPDG== 23) { zname=ModelPrtcls[i].name;  }
      if(ModelPrtcls[i].NPDG== 24) { wname=ModelPrtcls[i].name;  Wname=ModelPrtcls[i].aname;} 
      if(ModelPrtcls[i].NPDG==-24) { wname=ModelPrtcls[i].aname; Wname=ModelPrtcls[i].name; }   
   }

   
   if(!VZdecay || !VWdecay) { VZdecay=1; VWdecay=1; cleanDecayTable();} 

// look for SM Higgs 

   if(gluname && bname && lname && wname && zname)
   for(smH=0;smH<nModelParticles;smH++) if( ModelPrtcls[smH].spin2==0 && ModelPrtcls[smH].cdim==1 
   && ModelPrtcls[smH].name[0]!='~'  && strcmp(ModelPrtcls[smH].name,ModelPrtcls[smH].aname)==0  )
   {  double w,ggBr,bbBr,llBr,zzBr,wwBr; 
      txtList L;
      double ren;
//                       gg     bb      ll     zz     ww       
      double Br128SM[5]={0.077, 0.54,  0.059, 0.034, 0.264};    double  wSM_128=4.71E-3;
      double Br125SM[5]={0.080, 0.59,  0.065, 0.027, 0.213};    double  wSM_125=4.24E-3;
      double Br123SM[5]={0.081, 0.62,  0.068, 0.022, 0.182};    double  wSM_123=3.99E-3;
      double prec=0.9;  
      char chan[50];
      double hMass=pMass(ModelPrtcls[smH].name);
      if(hMass<123 || hMass>128)  continue;
      double a123= (hMass-125)*(hMass-128)/(2*5), a125=-(hMass-123)*(hMass-128)/(2*3), a128=(hMass-123)*(hMass-125)/(3*5);
      double  ggBrSM=a123*Br123SM[0]+a125*Br125SM[0]+a128*Br128SM[0],
              bbBrSM=a123*Br123SM[1]+a125*Br125SM[1]+a128*Br128SM[1],
              llBrSM=a123*Br123SM[2]+a125*Br125SM[2]+a128*Br128SM[2],
              zzBrSM=a123*Br123SM[3]+a125*Br125SM[3]+a128*Br128SM[3],
              wwBrSM=a123*Br123SM[4]+a125*Br125SM[4]+a128*Br128SM[4],
              wSM   =a123*wSM_123+a125*wSM_125+a128*wSM_128;

      w=pWidth(ModelPrtcls[smH].name,&L);
//printf("h2 width=%e\n",w);      
      sprintf(chan,"%s,%s",gluname,gluname);
      ggBr=findBr(L, chan);
      sprintf(chan,"%s,%s",lname,Lname);
      llBr=findBr(L, chan);
      sprintf(chan,"%s,%s",bname,Bname);     
      bbBr=findBr(L, chan);
      sprintf(chan,"%s,%s",wname,Wname);
      wwBr=findBr(L, chan);
      if(wwBr==0) { sprintf(chan,"%s *",wname); wwBr=2*findBr(L, chan);}
         
      sprintf(chan,"%s,%s",zname,zname);      
      zzBr=findBr(L, chan);                
      if(zzBr==0) {sprintf(chan,"%s *",zname); zzBr=findBr(L, chan);}    
      if(ggBr==0)
      { ren=w/(w+ ggBrSM*wSM);
        llBr*=ren;
        bbBr*=ren;
        wwBr*=ren;
        zzBr*=ren;
      }     
//printf("     ggBr=%.2e   llBr=%.2e bbBr=%e wwBr=%.2e zzBr=%.2e \n", ggBr,  llBr, bbBr, wwBr, zzBr);  
//printf(" SM  ggBr=%.2e   llBr=%.2e bbBr=%e wwBr=%.2e zzBr=%.2e \n", ggBrSM,llBrSM, bbBrSM, wwBrSM, zzBrSM);   
      if(    bbBrSM*prec< bbBr && bbBr<bbBrSM*(2-prec) 
          && llBrSM*prec< llBr && llBr<llBrSM*(2-prec)
//          && zzBrSM*prec< zzBr && zzBr<llBrSM*(2-prec)
          && zzBr <0.05
          && wwBrSM*prec< wwBr && wwBr<wwBrSM*(2-prec)
       ) break;       
   }
   
   if(smH<nModelParticles) printf("found SM-like Higgs = %s\n",ModelPrtcls[smH].name);
   else  printf("warning: no SM-like Higgs\n");

   printf("writing mass block and decay tables ... \n");

//  QNUMBERS
   for(i=0;i<nModelParticles;i++)
   {  
      // check whether in SM list 
      int isBSMparticle=1;
      for(j=0;j<18;j++) {
         if(abs(ModelPrtcls[i].NPDG)==SMP[j]) {
           isBSMparticle=0;
           continue;
         }  
      } 
      if(abs(newPDG(ModelPrtcls[i].NPDG))==25) isBSMparticle=0; // Higgs is SM particle

      // write QNUMBERS for BSM particles
      if(isBSMparticle==1) 
      {
        fprintf(f,"BLOCK QNUMBERS %d  # %s\n", newPDG(ModelPrtcls[i].NPDG), ModelPrtcls[i].name);   
        fprintf(f," 1  %d # 3*el.charge\n 2  %d # 2*spin+1\n 3  %d # color dim\n 4  %d # 0={ self-conjugated}\n",
             ModelPrtcls[i].q3, 
             ModelPrtcls[i].spin2+1, 
             ModelPrtcls[i].cdim, 
             strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname)? 1:0   );
        if(ModelPrtcls[i].name[0]=='~') fprintf(f,"11  1 # Z2-like symmetry factor (0=even,1=odd)\n\n"); 
        else                            fprintf(f,"11  0 # Z2-like symmetry factor (0=even,1=odd)\n\n"); 

      } 
   }   

   fprintf(f,"BLOCK MASS\n");
   for(i=0;i<nModelParticles;i++) if(pMass(ModelPrtcls[i].name) <PcmMax)
   { 
     for(j=0;j<16;j++) if(abs(ModelPrtcls[i].NPDG)==SMP[j]) break; 
     if(j==16 )
     { 
        np++; 
        plist=realloc(plist,np*sizeof(char*));
        plist[np-1]=ModelPrtcls[i].name;
        if(strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname))
        { np++;
          plist=realloc(plist,np*sizeof(char*));
          plist[np-1]=ModelPrtcls[i].aname;
        }    
        fprintf(f,"  %d  %E  # %s  \n",newPDG(ModelPrtcls[i].NPDG),findValW(ModelPrtcls[i].mass),ModelPrtcls[i].name);   
     }
   }
   fprintf(f,"\n");

   for(i=0;i<nModelParticles;i++) 
   {  for(j=0;j<16;j++) if(ModelPrtcls[i].NPDG==SMP[j]) break;
      
      if(j==16) newSlhaDecayPrint(ModelPrtcls[i].name,0,f); 
   }

   printf("computing LHC cross sections ... \n");

   for(i=0;i<np;i++) for(j=i;j<np;j++) if(pMass(plist[i])+pMass(plist[j])<PcmMax)
    if(plist[i][0]=='~' && plist[j][0]=='~') 
    {  int q31,q32,q3,c1,c2;

       qNumbers(plist[i], NULL, &q31,&c1);
       qNumbers(plist[j], NULL, &q32,&c2);
       
       q3=q31+q32;

       if(q3<0) { q3*=-1; if(abs(c1)==3) c1*=-1; if(abs(c2)==3)  c2*=-1;}
       if(c1>c2){ int c=c1; c1=c2;c2=c;}
       if(c1==8) c1=1;
       if(c2==8) c2=1;

       int ok=0;

       switch(q3)
       {  case 0:  if( (c1==-3 && c2==3 ) || (c1==1  && c2==1) ) ok=1;
                   break;
          case 1:  if( (c1==3  && c2==3 ) || (c1==-3 && c2==1) ) ok=1;
                   break;
          case 2:  if( (c1==-3 && c2==-3) || (c1==1  && c2==3) ) ok=1;
                   break;
          case 3:  if( (c1==-3 && c2==-3) || (c1==1  && c2==1) ) ok=1;
                   break;
          case 4:  if( (c1==3  && c2==3) ) ok=1;
                   break;
       }
   
       if(ok==0) continue;
       //------------------------- MSSM EWinos only -----------------------------------
       //if( abs(pNum(plist[i]))<1000022 || abs(pNum(plist[i]))>1000037 || 
       //    abs(pNum(plist[j]))<1000022 || abs(pNum(plist[j]))>1000037) continue;
       //------------------------------------------------------------------------------
       {  double dcs;
          double Qf=0.5*(pMass(plist[i])+pMass(plist[j]));
          for(double Pcm=PcmMin; ; ) 
          {    
             dcs=hCollider(Pcm,1,nf,Qf,Qf,plist[i],plist[j],0,wrt);             
             if(dcs>csMinFb*0.001)
             {
               fprintf(f,"XSECTION  %E   2212  2212  2  %d  %d\n",2*Pcm, newPDG(pNum(plist[i])),newPDG(pNum(plist[j]))); 
/*pb*/         fprintf(f," 0  0  0  0  0  0 %E # CalcHEP %s\n\n", dcs,VERSION);
               nChan++;
             }
             if(Pcm==PcmMax) break; else Pcm=PcmMax;
          }  
       }
    }    

  fclose(f);
  free(plist);

  printf("SLHA input file done.\n\n");
  

  if(VZdecay!= VZdecay_ || VWdecay!=VWdecay_ ) { VZdecay=VZdecay_; VWdecay=VWdecay_; cleanDecayTable();}  
  if(nChan) return 0; else return 1;
}

int  smodels_(int*LHCrun, int *nf,double *csMinFb, char*fileName, char*version, int *wrt,int len1,int len2)
{
  char cName1[30];
  fName2c(fileName,cName1, len1);
  char cName2[30];
  fName2c(version,cName2, len2);
     
  return smodels(*LHCrun,*nf, *csMinFb, cName1,cName2,*wrt);
}
