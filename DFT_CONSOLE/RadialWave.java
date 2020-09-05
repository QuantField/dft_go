import javax.swing.JFrame;

import java.util.*;


public class RadialWave 
{
	    int      N ;
	    double   h;
	    double   PI = ath.PI;
	    double   Rmax;
	    double[] Sol;
	    double[] r   ;
	    double[] rho ;
	    double[] wave;
	    double[] Vcoulomb;
	    double[] Vcentrifugal;
	    double[] Veff;  // Total Potential
	    double[] Veff0; // used to store the Coulomb potential or any given one
     	double[] K;
     	double[][] EnergyBrackets;
     	int        EnergyLevels;
     	double E;
     	int Z;
	   
	    /*
	     * details about the bisection method.
	     * The two values we use for the secant method it we wish to use it
	     * bisec_E0
	     * bisec_E1  
	     */
	    double bisec_E0 , bisec_E1;
	     
	    public RadialWave(int N, int Z0, double Rmax)
	    {
	    	this.N = N;
	        Z = Z0;
	        this.Rmax = Rmax;
	        h = Rmax/(double)N;
	   
	        r   = new double[N+1];
	        rho = new double[N+1];
	        Sol = new double[N+1];
		    Vcoulomb     = new double[N+1];
		    Vcentrifugal = new double[N+1];
		    Veff         = new double[N+1]  ; 
		    Veff0        = new double[N+1]  ;
	        
	        K    = new double[N+1];
	        // maximum of 10 levels of energies
	        EnergyBrackets = new double[10][2];
	        for(int i=0; i<=N; i++)
	        {
	        	r[i]        = i*h;
	        	Vcoulomb[i] = -Z/r[i];	        	
	        }	        
	        
	        // this is the constructor for the computation of energy levels of hydrogenoid atoms
	        // therefore the non centrifugal potential is the usual Coulomb one, need to set Veff0 to it:
	        setVeff0( Vcoulomb);
	        
	        // if however we want to use a different potential, in the DFT case for example, setVeff0 need to be
	        // set to the new potential.
	    }
	    
	    void setVcentrifugal(int j)
	    {
	        for(int i=0; i<=N; i++)
	        {
	        	Vcentrifugal[i] = 0.5*j*(j+1)/(r[i]*r[i]);	        	
	        }   	
	    }
	    
	    void setVeff0(double[] Pot)
	    {
	    	for(int i=0; i<=N; i++)	
	    	Veff0[i] = Pot[i];
	    }

	    
	    void setVeff()
	    {
	    	for(int i=0; i<=N; i++)	
	    	Veff[i] = Veff0[i]+Vcentrifugal[i];
	    }
        
	    public ArrayList<WaveFunction> calculatesStates( int Nmax)
	    {
	    	// Nmax is the highest principal quantum number  max principal number.
	    	
	    	// Orb contains all the wavefunctions with associated data (psi, n,l,E)
	    	ArrayList<WaveFunction> Orb = new ArrayList<WaveFunction>();
	    	
	    	//System.out.println("Finding Energy Brackets :");
	    	
	    	for (int l=0; l<=Nmax-1; l++)
	    	{
	    		//System.out.println("l = " + l);
	          	setVcentrifugal(l);
		    	setVeff();
		    	// -Z^2*Ry is the lowest bracket.
		    	double Emin = -0.5*Z*Z/((l+1)*(l+1))-3.; // for Ry energy units 0.5 becomes 1 
		    	findEigenEnergies( Emin, 0.05, 0,Nmax-l, false);// if true displays on console
		    	
		    	for(int i=0; i<EnergyLevels;i++)
		    	{
		    		double Eig = getEnergy(EnergyBrackets[i][0], EnergyBrackets[i][1]);
		    		Numerov(Eig);
		    		int n=i+1+l;
		    		//System.out.printf("n =%d    E =%14.10f  \n", n, E);
		    		Orb.add(makeState(n,l));		    		
		    	}
	    	}
 	    	// Sorting the Eigen States.
	    	  
	    	Collections.sort(Orb);
	    	
	    	/*
	    	System.out.println("===========================");
	    	System.out.println("  Sorted Quantum States:   ");
	    	System.out.println("===========================");
	    	for (int i=0; i<Orb.size(); i++)
	    	{
	    		System.out.printf(" n=%d \t l=%d \t E=%10.6f \n", 
	    				Orb.get(i).n, Orb.get(i).l, Orb.get(i).Energy);
	    	}
	    	*/
	    	return Orb;	    	
	  	}
	    
	    public void test1()
	    {
	    	RadialWave H = new RadialWave(4000, 1, Rmax);// Z = 1;
	    	
	    	H.setVcentrifugal(0);
	    	H.setVeff();
	    	
	    	System.out.println("Finding Energy Brackets :");
	    	
	    	H.findEigenEnergies( -10, 0.1, 0, 3,true);
	    	
	    	
	    	
	    	WaveFunction[] Spectrum = new WaveFunction[H.EnergyLevels];
	    	
	    	System.out.println("-------------------------------------------------------");
	    	System.out.println("Energy Levels(Ry) :");
	    	int l=0;
	    	for(int i=0; i<H.EnergyLevels;i++)
	    	{
	    		//H.secant(H.EnergyBrackets[i][0]);
	    		double Eig = H.getEnergy(H.EnergyBrackets[i][0], H.EnergyBrackets[i][1]);
	    		H.Numerov(Eig);
	    		System.out.printf(" E =%14.10f  \n", H.E);
	    			    		
	    		Spectrum[i] = H.makeState(i+1,l);// the lowest state we call it 1 not 0;
	    		int n=i+1;
	    		//H.plot(H.r, Spectrum[i].Density,"Probability Density",  "n = " + n);
	    	}
	    }
	    
	    // secant method after few iterations of the bisection method.
	    double getEnergy(double Ea,double Eb)
	    {
	    	bisection(Ea, Eb, 1e-4);	    		  
	  	    double E  = bisec_E0;
	  	    double E1 = bisec_E1;
	  	    double E2 = 0;
	  	    int iter = 0;
	  	    while ( Math.abs((E2-E1)/E2)>1e-8 || iter>30)
	  	    {   iter++;
	  	    	E2 = E1 - (E1-E)*F(E1)/(F(E1)-F(E));
	  	    	E  = E1;
	  	    	E1 = E2;
	  	   }
	  	   //System.out.printf(" iter =%3d   E=%14.10f\n",iter,E2);	  	    
	  	   return E2;
	  	 }

	    
	 /**
	  *    
	  * @param Emin
	  * @param dE
	  * @param Emax
	  * 
	  * We find the bracket that contain the Eigen energy E : 
	  * Ea = EnergyBrackets[i][0] and Eb= EnergyBrackets[i][0]
	  * i.e  Ea<E<Eb
	  * we can then either use bisection(Ea,Eb) or secant(Ea).
	  * n is the maximum principal number.
	  */
	 void findEigenEnergies(double Emin, double dE, double Emax, int n, boolean show)
	 {
		 int M = (int)((Emax-Emin)/dE);
		 double[] SolN = new double[M+1];
		 for(int i=0; i<=M; i++) SolN[i] = F(Emin+i*dE);
		 
		 EnergyLevels = 0;
		 for(int i=0; i<M; i++)	 {
			 if (SolN[i]*SolN[i+1]<=0)  {
				 EnergyBrackets[EnergyLevels][0] = Emin+i*dE;
				 EnergyBrackets[EnergyLevels][1] = Emin+(i+1)*dE;
				 EnergyLevels++;
				 if (EnergyLevels>=n) break; 
			 }
		 }
		 if (show) {
			for(int i=0; i<EnergyLevels; i++) {
				System.out.printf("Ea = %16.10f  \t  Eb = %16.10f \n",EnergyBrackets[i][0],EnergyBrackets[i][1]);
			}
		 }
		 
	 }
	    
	       
	 void Numerov(double E0)
	    {
	        this.E = E0;
	        
	        for (int i=1;i<=N;i++)
	        {
	          K[i] =  2*(E0-Veff[i]); // in Hartree  or E0-2*Veff[i] for Ryleigs;          
	        }
		    double  C = h*h/12.;
		    
		    double alph = Math.sqrt(-2*E0);
			Sol[N]   = r[N]*Math.exp(-alph*r[N]);//0; // r[N]*Math.exp(-Z*r[N]) ;  // Sol(xb) = 0;
			Sol[N-1] = r[N-1]*Math.exp(-alph*r[N-1]);//1e-8; // r[N-1]*Math.exp(-Z*r[N-1]) ; //1e-6; //    
			//Sol[N]   =  r[N]*Math.exp(-Z*r[N]) ;  
			//Sol[N-1] =  r[N-1]*Math.exp(-Z*r[N-1]) ;    
			
			for(int i=N-1; i>=1;i--) 
			{    	    		     
				Sol[i-1] = ( 2*(1. - 5.*C*K[i])*Sol[i] - 
						(1. + C*K[i+1] )*Sol[i+1] )/(1.+C*K[i-1]);
						
			}
	 }
	 
	 double F(double E0)
	  {
	    Numerov( E0);
	    return Sol[0]; //Sol[N]; this for Numerov starting from the Left
	  }
	 
	 
	  void bisection(double Ea0, double Eb0, double epsilon)
	  {
	      double Ea = Ea0 ;
	      double Eb = Eb0;
	      double s  = 1;
	      double s0 = 0.5;
	      int iter = 0;
	      if (F(Ea)*F(Eb)<0)
	      {  while ( Math.abs((s-s0)/s) >epsilon )
	          {   iter++;
	           s0 = s;
	           s = (Ea+Eb)/2;
	           bisec_E0 = s0; // nothing of use
	 	       bisec_E1 = s;  // here, is just to speed up the secant method
	           if (F(s)*F(Ea)<0) Eb=s;
	           if (F(s)*F(Eb)<0) Ea=s;
	           //System.out.printf(" iter = %3d E = %14.10f\n",iter,s);
	          } 
	          //System.out.printf(" iter = %3d E = %14.10f\n",iter,s);       
	      }
	      else System.out.println("No solutions");
	   }
	 
	
		
	public  double integrate(double[] x, double[] y)
	{
		int N = x.length-1;
		if (N%2 !=0) 
		{ 
		   System.out.println("method integrate error:N must be even."); 
		   System.exit(1);
		}
		if (x.length!=y.length) 
		{
		   System.out.println("method integrate error: x and y must have same length."); 
		   System.exit(1);}
		double h = x[1]-x[0];
		// Composite Simpson integration
 		double  S = y[0]+ y[N];
		for(int i=1;i<=N;i=i+2)   S+= 4*y[i];
		for(int i=2;i<=N-1;i=i+2) S+= 2*y[i];
		return  h*S/3;		    	
	}

	
	public void normilise()
	{
		double[] p = new double[N+1];
		for(int i=0; i<=N; i++) p[i]=Sol[i]*Sol[i];
		double Norm = integrate(r, p);
		Norm =Math.sqrt(Norm);
		for(int i=0; i<=N; i++) Sol[i]=Sol[i]/Norm;		
	}
	
	
	 WaveFunction makeState(int k, int l)
	 {
		 normilise();
		 double[] denVol = new double[N+1];
		 wave   = new double[N+1];
		 rho[0]    = Sol[0]*Sol[0];// because the loop start from i=1
		 wave[0]   = 0;
		 
		 for(int i=1;i<=N;i++)
		 {
			 rho[i]    = Sol[i]*Sol[i];
			 wave[i]   = Sol[i]/r[i];
			 denVol[i] = rho[i]/(4*PI*r[i]*r[i]); // normalised per volume				 
		 }
		 denVol[0] = denVol[1] ;
		 //wave[0]=wave[1];// discontinuity at r=0 ;
		 
		 WaveFunction Wve = new WaveFunction(this.wave, this.rho, denVol,  this.E, k, l);
		 return Wve;
	 } 
}	 

	 






