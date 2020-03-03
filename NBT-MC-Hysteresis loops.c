#include<math.h>
#include<stdio.h>
#include<stdlib.h>

int main(){
	FILE *fichierx;
	fichierx=fopen("25Dy-80k.dat","w");
	FILE *fichiery;
	fichiery=fopen("poly.dat","w");
	FILE *fichierz;
	fichierz=fopen("polz.dat","w");
	

	int L=4,N=L*L*L,   // L=taille du réseau selon les axes x et y et N= nombre de sites
	MCSTEPS=100000,  // Etapes Monte Carlo à l'équilibre
	 MCDISC=10000;   // Etapes MC à éliminer (avant équilibre)
	double kboltz=0.000086173303,dipx[L][L][L],dipy[L][L][L],dipz[L][L][L], polx,comp, poly, polz, polx_moy,pol_moy4,pol_moy2,polx_moy2, polx_moy4, poly_moy,poly_moy2,poly_moy4, polz_moy,polz_moy2,polz_moy4, susceptx, suscepty, susceptz, suscept, En, Enx, Eny, Enz, En_moy, Enx_moy, Eny_moy, Enz_moy, En_moy2, Enx_moy2, Eny_moy2,Enz_moy2,chal_sp,chal_spx, chal_spy,chal_spz,Cum, Cumx, Cumy, Cumz,theta[L][L][L],phi[L][L][L],theta_f,phi_f, susceptibi,T=80.,champmin=-0.2,champmax=0.2,champstep=0.01, theta_c=0., phi_c=1.,Pb=0.15,Pn=0.1,Pd=0.005, cd=0.25,Th=1.;

	double polxn,polyn,polzn,polx_moysn,poly_moysn,polz_moysn, polx_moy2sn,poly_moy2sn , polz_moy2sn, polx_moy4sn, poly_moy4sn, polz_moy4sn, susceptxsn, susceptysn, susceptzsn, susceptsn,Cumsn,Cumxsn,Cumysn, Cumzsn;
	double polxb,polyb,polzb,polx_moysb,poly_moysb,polz_moysb, polx_moy2sb,poly_moy2sb , polz_moy2sb, polx_moy4sb, poly_moy4sb, polz_moy4sb, susceptxsb, susceptysb, susceptzsb, susceptsb,Cumsb,Cumxsb,Cumysb, Cumzsb;
	double theta_initial_N=Th*M_PI/4., phi_initial_N=Th*M_PI/4.,theta_initial_B=Th*M_PI/4., phi_initial_B=Th*M_PI/4.;
        double delta=0.,champ,champx=0,champy=0,champz=0;
	int i,j,k,u,v, sig=1., CC=-1., site, compt;
	double factbb=0.066*3, factnn=0.29*3, factbn=0.61*3;
	double Jbbx=factbb,Jbby=factbb,Jbbz=factbb, Jnnx=factnn,Jnny=factnn,Jnnz=factnn,Jbnx=factbn, Jbny=factbn, Jbnz=factbn;
	double polxn_moy2, polxn_moy, polxn_moy4, polyn_moy2, polyn_moy, polyn_moy4, polzn_moy2, polzn_moy, polzn_moy4, polxb_moy2, polxb_moy, polxb_moy4, polyb_moy2, polyb_moy, polyb_moy4, polzb_moy2, polzb_moy, polzb_moy4, poln_moy, polb_moy;
	//***************  Etape d'initialisation à la config. Ferro *****************

		// Bi 
	for(i=0;i<=L-1;i=i+1){
		for(j=0;j<=L-1;j=j+1){
			for(k=0;k<=L-1;k=k+2){
				dipx[i][j][k]=sig*Pb*1./2.;
				dipy[i][j][k]=sig*Pb*1./2.;
				dipz[i][j][k]=sig*Pb*sqrt(2.)/2.;	
				theta[i][j][k]=theta_initial_B;
				phi[i][j][k]=phi_initial_B;
			
	}
		}
	}
			
			

			compt=0;
		for(i=0;i<=L-1;i=i+1){
					for(j=0;j<=L-1;j=j+1){
						for(k=0;k<=L-1;k=k+2){
			double x= (double)rand()/RAND_MAX;
	
				if (x<=cd) {

		
				 dipx[i][j][k]=(Pd/Pb)*dipx[i][j][k];
				 dipy[i][j][k]=(Pd/Pb)*dipy[i][j][k];
				 dipz[i][j][k]=(Pd/Pb)*dipz[i][j][k];
				theta[i][j][k]=theta_initial_B;
				phi[i][j][k]=phi_initial_B;


}
}
	}

}
		
			
	
		//Na
	for(i=0;i<=L-1;i=i+1){
		for(j=0;j<=L-1;j=j+1){
			for(k=1;k<=L-1;k=k+2){
				dipx[i][j][k]=sig*Pn*1./2.;
				dipy[i][j][k]=sig*Pn*1./2.;
				dipz[i][j][k]=sig*Pn*sqrt(2.)/2.;	
				theta[i][j][k]=theta_initial_N;
				phi[i][j][k]=phi_initial_N;
		
				
		}
		}
	}


	polxn=sig*1./2.*Pn*N*0.5;    
	 polyn=sig*1./2.*Pn*N*0.5;
	 polzn=sig*sqrt(2)/2.*Pn*N*0.5; 
	

	polxb=sig*1./2.*N*0.5*(Pb*(1-cd)+Pd*cd);    
	 polyb=sig*1./2.*N*0.5*(Pb*(1-cd)+Pd*cd);
	 polzb=sig*sqrt(2)/2.*N*0.5*(Pb*(1-cd)+Pd*cd); 

				for(champ=champmin;champ<=champmax;champ=champ+champstep){  //Electric field configuration

                champx=CC*champ*sin(theta_c)*cos(phi_c);	
		champy=CC*champ*sin(theta_c)*sin(phi_c);
		champz=CC*champ*cos(theta_c);	

	En=-(2./N)*(Jbbx*(polxb*polxb)+Jbby*(polyb*polyb)+Jbbz*(polzb*polzb)+Jnnx*(polxn*polxn)+Jnny*(polyn*polyn)+Jnnz*(polzn*polzn))-(1./N)*(Jbnx*(polxb*polxn)+Jbny*(polyb*polyn)+Jbnz*(polzb*polzn))-(1./2.)*(champx*(polxn+polxb)+champy*(polyn+polyb)+champz*(polzn+polzb)); 	


      
	 

	//****************************************************************************
	int mc;

		En_moy=0.;
		Enx_moy=0.;
		Eny_moy=0.;
		Enz_moy=0.;

		En_moy2=0.;
		Enx_moy2=0.;
		Eny_moy2=0.;
		Enz_moy2=0.;

		// sous réseaux Na 
		polxn_moy=0.;
		polyn_moy=0.;
		polzn_moy=0.;

		polxn_moy2=0.;
		polyn_moy2=0.;
		polzn_moy2=0.;

		polxn_moy4=0.;
		polyn_moy4=0.;
		polzn_moy4=0.;
		poln_moy=0.;
		polb_moy=0.;		

		// sous réseaux Bi
		polxb_moy=0.;
		polyb_moy=0.;
		polzb_moy=0.;

		polxb_moy2=0.;
		polyb_moy2=0.;
		polzb_moy2=0.;

		polxb_moy4=0.;
		polyb_moy4=0.;
		polzb_moy4=0.;		


		//---------------------------------------------------------------------------
		for(mc=1;mc<=MCSTEPS+MCDISC;mc++){  // Boucle sur les configurations
				compt=0;
				
				for(site=1;site<=N;site++){   // boucle sur N sites par configuration
			//choix aléatoire d'un site de coordonnés (i,j,k) tels que 0<=i<=L-1 et 0<=j<=L-1
				int i=rand()%L;
				int j=rand()%L;
				int k=rand()%L;
				int e=i+j+k;
					

				double dipx_initial=dipx[i][j][k];
				double dipy_initial=dipy[i][j][k];
				double dipz_initial=dipz[i][j][k];
				double theta_initial=theta[i][j][k];
				double phi_initial=phi[i][j][k];
				
                                // les 6 sites voisins du site (i,j,k) choisi en tenant 
				// comptes conditions aux bords périodiques
				int isuiv=i+1; if(i==L-1)isuiv=0;
				int iprec=i-1; if(i==0)iprec=L-1;
				int jsuiv=j+1; if(j==L-1)jsuiv=0;
				int jprec=j-1; if(j==0)jprec=L-1;
				int ksuiv=k+1; if(k==L-1)ksuiv=0;
				int kprec=k-1; if(k==0)kprec=L-1;
	

			double som1x=dipx[isuiv][j][k]+dipx[iprec][j][k];
			double som1y=dipy[isuiv][j][k]+dipy[iprec][j][k];
			double som1z=dipz[isuiv][j][k]+dipz[iprec][j][k];

				double som2x=dipx[i][jsuiv][k]+dipx[i][jprec][k];
				double som2y=dipy[i][jsuiv][k]+dipy[i][jprec][k];
				double som2z=dipz[i][jsuiv][k]+dipz[i][jprec][k];

			double som3x=dipx[i][j][ksuiv]+dipx[i][j][kprec];
			double som3y=dipy[i][j][ksuiv]+dipy[i][j][kprec];
			double som3z=dipz[i][j][ksuiv]+dipz[i][j][kprec];

		

                    

                        //domaine 
			double rt=(double)rand()/RAND_MAX; // nombre aléatoire 0<=r<1
			double rp=(double)rand()/RAND_MAX; // nombre aléatoire 0<=r<1

			double dtheta=(2.*rt-1.)*(M_PI/12.); 
			double dphi=(2.*rp-1.)*(M_PI/12.); 
		
			double theta_f=theta_initial+dtheta; 
			double phi_f=phi_initial+dphi; 



	if (k%2 == 0){

			double dipx_final=dipx_initial*sin(theta_f)*cos(phi_f)/(sin(theta_initial)*cos(phi_initial));
			double dipy_final=dipy_initial*sin(theta_f)*sin(phi_f)/(sin(theta_initial)*sin(phi_initial));
			double dipz_final=dipz_initial*cos(theta_f)/(cos(theta_initial)); 
	
//--------------------------------------------------------------------------------------------------
double dEB=-(dipx_final-dipx_initial)*(Jbbx*(som1x+som2x)+Jbnx*(som3x))-(dipy_final-dipy_initial)*(Jbby*(som1y+som2y)+Jbny*(som3y))-(dipz_final-dipz_initial)*(Jbbz*(som1z+som2z)+Jbnz*(som3z))-champx*(dipx_final-dipx_initial)-champy*(dipy_final-dipy_initial)-champz*(dipz_final-dipz_initial);
	
			 if(dEB<=0){ //dans ce cas le dipole (i,j,k) sera flipé, on calcule la polarisation
            				// de la nouvelle configuration
					polxb=polxb+dipx_final-dipx_initial;
					polyb=polyb+dipy_final-dipy_initial;
					polzb=polzb+dipz_final-dipz_initial;
                         
					En=En+dEB;
					dipx[i][j][k]=dipx_final;	 //flip du dipole (i,j,k)
					dipy[i][j][k]=dipy_final;
					dipz[i][j][k]=dipz_final;
					theta[i][j][k]=theta_f;
					phi[i][j][k]=phi_f;
	
				}else{
					double p=exp(-dEB/(kboltz*T));   // probabilité de flip
					double r=(double)rand()/RAND_MAX; // nombre aléatoire 0<=r<1
					if(r<=p){ //dans ce cas aussi le dipole (i,j,k) sera flipé, on calcule la polarisation
                  // de la nouvelle configuration
					polxb=polxb+dipx_final-dipx_initial;
					polyb=polyb+dipy_final-dipy_initial;
					polzb=polzb+dipz_final-dipz_initial;
					
					En=En+dEB;
					dipx[i][j][k]=dipx_final;	 //flip du dipole (i,j,k)
					dipy[i][j][k]=dipy_final;
					dipz[i][j][k]=dipz_final;
					theta[i][j][k]=theta_f;
					phi[i][j][k]=phi_f;
	
		
					
		}		
}
	}

//Impair
				else{


	    	double dipx_final=Pn*sin(theta_f)*cos(phi_f);
			double dipy_final=Pn*sin(theta_f)*sin(phi_f);
			double dipz_final=Pn*cos(theta_f);


          double dEN=-(dipx_final-dipx_initial)*(Jnnx*(som1x+som2x)+Jbnx*(som3x))-(dipy_final-dipy_initial)*(Jnny*(som1y+som2y)+Jbny*(som3y))-(dipz_final-dipz_initial)*(Jnnz*(som1z+som2z)+Jbnz*(som3z))-champx*(dipx_final-dipx_initial)-champy*(dipy_final-dipy_initial)-champz*(dipz_final-dipz_initial);

	                if(dEN<=0){ 
                                
					polxn=polxn+dipx_final-dipx_initial;
					polyn=polyn+dipy_final-dipy_initial;
					polzn=polzn+dipz_final-dipz_initial;
                          

					En=En+dEN;
					dipx[i][j][k]=dipx_final;	 //flip du dipole (i,j,k)
					dipy[i][j][k]=dipy_final;
					dipz[i][j][k]=dipz_final;
					theta[i][j][k]=theta_f;
					phi[i][j][k]=phi_f;

				}else{
					double p=exp(-dEN/(kboltz*T));   // probabilité de flip
					double r=(double)rand()/RAND_MAX; // nombre aléatoire 0<=r<1
					if(r<=p){ //dans ce cas aussi le dipole (i,j,k) sera flipé, on calcule la polarisation
  // de la nouvelle configuration
					polxn=polxn+dipx_final-dipx_initial;
					polyn=polyn+dipy_final-dipy_initial;
					polzn=polzn+dipz_final-dipz_initial;
           
					
					En=En+dEN;
					dipx[i][j][k]=dipx_final;	 //flip du dipole (i,j,k)
					dipy[i][j][k]=dipy_final;
					dipz[i][j][k]=dipz_final;
					theta[i][j][k]=theta_f;
					phi[i][j][k]=phi_f;
		

					}
				}
				}

}

			if(mc>MCDISC){  // accumulation des sommes pour le calcul des valeurs moyenne
					// cette accumulation se fait uniquement pour les configurations d'équilibre
				polxn_moy=polxn_moy+polxn;
				polyn_moy=polyn_moy+polyn;
				polzn_moy=polzn_moy+polzn;

				polxn_moy2=polxn_moy2+polxn*polxn;
				polyn_moy2=polyn_moy2+polyn*polyn;
				polzn_moy2=polzn_moy2+polzn*polzn;

				polx_moy4=polxn_moy4+polxn*polxn*polxn*polxn;
				poly_moy4=polyn_moy4+polyn*polyn*polyn*polyn;
				polz_moy4=polzn_moy4+polzn*polzn*polzn*polzn;
//
				polxb_moy=polxb_moy+polxb;
				polyb_moy=polyb_moy+polyb;
				polzb_moy=polzb_moy+polzb;

				polxb_moy2=polxb_moy2+polxb*polxb;
				polyb_moy2=polyb_moy2+polyb*polyb;
				polzb_moy2=polzb_moy2+polzb*polzb;

				polxb_moy4=polxb_moy4+polxb*polxb*polxb*polxb;
				polyb_moy4=polyb_moy4+polyb*polyb*polyb*polyb;
				polzb_moy4=polzb_moy4+polzb*polzb*polzb*polzb;
			
				En_moy=En_moy+En;
				En_moy2=En_moy2+En*En;
			}
		}
		// A la fin des étapes MC, les valeurs moyennes sont calculés en divisant les sommes
		// calculées par le nombre de configurations générées à l'équilibre
		polxn_moy=polxn_moy/MCSTEPS;
		polxn_moy2=polxn_moy2/MCSTEPS;
		polxn_moy4=polxn_moy4/MCSTEPS;
		polyn_moy=polyn_moy/MCSTEPS;
		polyn_moy2=polyn_moy2/MCSTEPS;
		polyn_moy4=polyn_moy4/MCSTEPS;
		polzn_moy=polzn_moy/MCSTEPS;
		polzn_moy2=polzn_moy2/MCSTEPS;
		polzn_moy4=polzn_moy4/MCSTEPS;
		
		polxb_moy=polxb_moy/MCSTEPS;
		polxb_moy2=polxb_moy2/MCSTEPS;
		polxb_moy4=polxb_moy4/MCSTEPS;
		polyb_moy=polyb_moy/MCSTEPS;
		polyb_moy2=polyb_moy2/MCSTEPS;
		polyb_moy4=polyb_moy4/MCSTEPS;
		polzb_moy=polzb_moy/MCSTEPS;
		polzb_moy2=polzb_moy2/MCSTEPS;
		polzb_moy4=polzb_moy4/MCSTEPS;

		En_moy=En_moy/MCSTEPS;
		En_moy2=En_moy2/MCSTEPS;

		
		//---------------
	

///////////////////////////////////////////////////////////////////////////////////////////////////////////		
poln_moy=2.*sqrt((polxn_moy/N)*(polxn_moy/N)+(polyn_moy/N)*(polyn_moy/N)+(polzn_moy/N)*(polzn_moy/N));
polb_moy=2.*sqrt((polxb_moy/N)*(polxb_moy/N)+(polyb_moy/N)*(polyb_moy/N)+(polzb_moy/N)*(polzb_moy/N));



printf("%.6lf \t %.6lf \t %.6lf \t %.6lf \t %.6lf \t %.6lf \n",-champ,2.*(polxn_moy/N)+2.*(polxb_moy/N),2.*(polyn_moy/N)+2.*(polyb_moy/N),2.*(polzn_moy/N)+2.*(polzb_moy/N),poln_moy+polb_moy,En_moy);
	fprintf(fichierx,"%.6lf \t %.6lf \t %.6lf \t %.6lf \t %.6lf \t %.6lf\n",-champ,2.*(polxn_moy/N)+2.*(polxb_moy/N),2.*(polyn_moy/N)+2.*(polyb_moy/N),2.*(polzn_moy/N)+2.*(polzb_moy/N),poln_moy+polb_moy,En_moy);
	

		printf("%.6lf \t %.6lf \t %.6lf \t %.6lf \t %.6lf \t %.6lf \n",champ,2.*polxb_moy/N,2.*polyb_moy/N,2.*polzb_moy/N,polb_moy,En_moy);
	fprintf(fichiery,"%.6lf \t %.6lf \t %.6lf \t %.6lf \t %.6lf \t %.6lf \n",champ,2.*polxb_moy/N,2.*polyb_moy/N,2.*polzb_moy/N,polb_moy,En_moy);

//printf("%.6lf \n",2./N);
/**	printf("%.6lf \t %.6lf \t %.6lf \t %.6lf \t %.6lf\n",T,polx_moy/N,susceptibi ,chal_sp,Cumx);
       fprintf(fichierx,"%.6lf \t %.6lf \t %.6lf \t %.6lf \t %.6lf\n",T,polx_moy/N,susceptibi,chal_sp,Cumx);
	
	printf("%.6lf \t %.6lf \t %.6lf \t %.6lf \t %.6lf\n",T,poly_moy/N,susceptibi,chal_sp,Cumy);
       fprintf(fichiery,"%.6lf \t %.6lf \t %.6lf \t %.6lf \t %.6lf\n",T,poly_moy/N,susceptibi,chal_sp,Cumy);
		
	printf("%.6lf \t %.6lf \t %.6lf \t %.6lf \t %.6lf\n",T,polz_moy/N,susceptibi,chal_sp,Cumz);
       fprintf(fichierz,"%.6lf \t %.6lf \t %.6lf \t %.6lf \t %.6lf \t %.6lf \n",T,polz_moy/N,susceptibi,chal_sp,Cumz,En);
**/

}
	fclose(fichierx);
	fclose(fichiery);
	fclose(fichierz);
	
	getchar();
	return 0;
}


