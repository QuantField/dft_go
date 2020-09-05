import java.util.ArrayList;


public class LDA 
{
	 int      N ;
	 int      Z ;
	 double   h;
	 double   PI = Math.PI;
	 double   Rmax;
	 double[] Sol;
	 double[] r   ;
	 double[] rho ;
	 double[] newrho;
	 double[] Veff0;
	 double[] Veff;
	 double[] Vxc;
	 double[] Uhartree;
	 
	 double Eb    ;   // sum Eigen energies
     double Etot  ;
     double Ekin  ;
     double Ecoul ;
     double Excor ;
     double Enuclei;
     
     int    N_nl; // number of orbitals (n,l)
     
     ExchangeCorrelation XC;
	 
	 
	 public LDA( int z, double Rmax, int N )
	 {
			this.Z     = z; // Number of electrons in the atom: Z=1,2,10,18,36
			this.Rmax = Rmax;
			this.N    = N;
			h = Rmax/(double)N;
		
			r     = new double[N+1];
			rho   = new double[N+1];
			newrho   = new double[N+1];
			Sol   = new double[N+1];
			Veff0 = new double[N+1];
			Veff  = new double[N+1];
			Uhartree = new double[N+1];
			Vxc      = new double[N+1];
			
			// setting the grid 
			for(int i=0; i<=N; i++) r[i]=i*h;
			
			// Initial guess of rho
			// rho normalised to Z , i.e. integ [0 Infinity] of 4*PI*r^2*rho = Z
			double Nrm =Z*Z*Z*Z/(64*PI);
			for(int i=0; i<=N; i++)
			{
				rho[i]=Nrm*Math.exp(-Z/2.*r[i]);
				Uhartree[i] = 0;
			}
			
			XC = new ExchangeCorrelation(3);
			// 3 corresponds to S.H.Vosko, L.Wilk, and M.Nusair, Can.J.Phys.58, 1200 (1980) 
			
		}
	 
	 
	   
	 
	   public void mainDemo()
	   {
		   
		RadialWave W = new RadialWave( N, Z, Rmax);
		
		double Etot0=10;
		Etot =20;
		int iter=0;
		
		int n =2;// this should be calculated automatically, this is just for Z<=10 
		ArrayList<WaveFunction> spectrum = null; // orbitals are stored here.
		
		double[] rhoOld = new double[N+1];
		
		while(distance(rho,newrho)>1e-8  && iter<=50) //this is more accurate but lengthier
		{	
		   iter++;
		   solvePoisson(); // we get Uhartree, from the rho, run first time with the initial guessed density.
		   calulate_Vxc(); // result in Vxc.. exchange correlation term potential
		   
		   // calculating the Vks term without the centrifugal potential.
		   for (int i=0; i<=N; i++) 
				  Veff0[i] = (-Z+Uhartree[i])/r[i]+Vxc[i];
		   Veff0[0]=Veff0[1];
		   	   
		   W.setVeff0(Veff0);
		   
		   spectrum = W.calculatesStates(n);
		   
		   BuildNewRho(Z, spectrum); // build new rho and sums the eigen energies of the occupied states (Eb);
		   
		   Etot0 = Etot;
		   calculateTotalEnergy(); // Etot is calculated
		   
		   System.out.printf("iter = %d   Etot = %14.10f\n", iter, Etot);
		   
		   double admix = 0.9;
		   for (int k=0; k<rho.length; k++)
		   {
		       //rhoOld[k] = rho[k];
			   rho[k] = (1.-admix)*rho[k]+admix*newrho[k]; // linear mixing
		       //rho[k] = newrho[k];
		   }
		}
		//while(Math.abs((Etot-Etot0)/Etot)>1e-8  && iter<=50);
		
		System.out.printf("\n**********************************\n");
		System.out.printf("Etot    =%10.6f\n", Etot);
		System.out.printf("Ekin    =%10.6f\n", Ekin);
		System.out.printf("Ecoul   =%10.6f\n", Ecoul);
		System.out.printf("Excor   =%10.6f\n", Excor);
		System.out.printf("Enuclei =%10.6f\n", Enuclei);
		System.out.printf("**********************************\n");
		
		for(int i = 0; i<=N_nl; i++)
		{
		  System.out.printf("n=%d \t l=%d \t E=%10.6f  \n",
				  spectrum.get(i).n, spectrum.get(i).l, spectrum.get(i).Energy );
		}
		System.out.printf("**********************************\n");
		
	   }
	 
	 
	 
	    public void firstTest()
	    {
	    	RadialWave W = new RadialWave( N, Z, Rmax);
	    	
	    	// get states up to principal quantum number n
	    	int n = 3;
		    ArrayList<WaveFunction> spectrum = W.calculatesStates(n);
		    
		    // experiencing with Z=10
		    int Z = 2;
		    // filling the states n,l with Z electrons, and summing up their Eigen Energies
			BuildNewRho(Z, spectrum);
			
			// just testing the density is correct
			double[] rhotest = new double[newrho.length];
			for(int i=0; i<newrho.length;i++)  rhotest[i] = 4*Math.PI*r[i]*r[i]*newrho[i];
			 
			//plot(r, newrho, "new dens");
			System.out.printf("Testing with an atom with %d electrons...\n",Z); 
			System.out.printf("Integration of the density over the whole volulme : %14.10f\n", integrate(r,rhotest));
			System.out.printf("The result should be %d\n",Z);		    		    
	    }    
	    
	    
		public  static  double integrate(double[] x, double[] y)
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
		
        // another technique that might be more accurate than the above
		public double integrate2(double[] r, double[] F)
		{ // Forth order integration routine
		  double dr= r[1]-r[0];
		  int Np = F.length;
		  double coeff[] = {17./48., 59./48.,43./48.,49./48.};
		  double sum = coeff[0]*F[0]+coeff[1]*F[1]+coeff[2]*F[2]+coeff[3]*F[3];
		  for (int i=4; i<F.length-4; i++) sum += F[i];
		  sum += coeff[0]*F[Np-1]+coeff[1]*F[Np-2]+coeff[2]*F[Np-3]+coeff[3]*F[Np-4];
		  return sum*dr;
		}

		
		public  double  calulate_rs(double rho)
		{ 
			return Math.pow(3./(4*PI*rho),1./3.) ;
		}
		
		
		// Calculating exchange-correlation part
		public void calulate_Vxc()
		{
			for (int i=0; i<=N; i++)
	        {
	           double rs = calulate_rs(rho[i]);
	           Vxc[i] = XC.Vx(rs) + XC.Vc(rs);
	        }
	   }
		
		//Calculates Uhartree Potential Vhartree= Uhartree/r
		public   void solvePoisson()
		{
			double[] k    = new double[N+1];
			for (int i=0;i<=N;i++) 	k[i] = -4*PI*r[i]*rho[i];
			
			double    C = h*h/12.;
			Uhartree[0] = 0 ;
			Uhartree[1] = Uhartree[0]+h;     
			for(int i=1; i<=N-1;i++) 
			{	    	    		     
				Uhartree[i+1] = -Uhartree[i-1] + 2.*Uhartree[i] 
						        + C*(k[i+1]+10*k[i]+k[i-1]);				
			}			
			// setting the boundary conditin U(infinity)=Z
			double alpha = (Z -Uhartree[N])/r[N];
			for (int i=0;i<=N;i++) Uhartree[i]=Uhartree[i]+alpha*r[i];			
		}
		
		
		public void BuildNewRho(int Z, ArrayList<WaveFunction> states)
		{// Knowing the energies of eigenstates, finds chemical potential and new charge density
		  		  
		  for (int k=0; k<newrho.length; k++) newrho[k]=0;
		  
		  Eb =0;        // Sum of eigenvalues
		  int    Nt =0;        // number of electrons added
		  for (int k=0; k<states.size(); k++)
		  {
		    int l = states.get(k).l;
		    int dN = 2*(2*l+1);  // degeneracy of each radial wave level
		    double deg = (double)dN;
		    double ferm = Nt+dN<=Z  ? 1. : (Z-Nt)/(2.*(2.*l+1)); // if shell is not fully-filled, take only part of charge
		  
		    for (int i=0; i<newrho.length; i++) 
		    	newrho[i] += ferm*deg*states.get(k).DensityPerVolume[i];
		    
		    Eb += deg*ferm*states.get(k).Energy; // Sum of eigenvalues times degeneracy
		    
		    Nt += dN;
		    
		    System.out.printf("adding .. n=%d  l=%d  E=%10.6f  OccupPerc=%2.2f == sumE=%10.6f\n", 
		    		            states.get(k).n, l,states.get(k).Energy, ferm, Eb);
		    
		    if (Nt>=Z){ N_nl = k; break;} // Finish when enough electrons added
		  }
		  
		}
		
		
		public void calculateTotalEnergy()
		{
			// preparing for integration dV = S*dr , where S =4*pi*r^2
			double[] S = new double[N+1];
			for (int i=0; i<=N; i++) S[i]=4.*Math.PI*r[i]*r[i];
			// to use temporarily
			double[] tmp = new double[N+1];
			
			
			// E_hartree = 1/2*V_hartree*n(r)			
			for (int i=0; i<=N; i++) tmp[i] = 0.5*S[i]*Uhartree[i]*rho[i]/r[i]; //new why newrho?
			tmp[0]=tmp[1];
	        double Ehartree = integrate(r,tmp);
	    
	    
	        // Adding exchange-correlation energy part
	        //ExchangeCorrelation XC = new ExchangeCorrelation(3);
	        // here not to confuse EcVc with Vc and ExVx with Vx
			tmp[0]=0;
	        for (int i=1; i<=N; i++)
	        {  double rs = calulate_rs(rho[i]); 
	           tmp[i]    = S[i]*rho[i]*(XC.ExVx(rs)+XC.EcVc(rs));  //new why newrho?
	         }
	        double Excor0 = integrate(r,tmp); // need to know about the difference between this and the other Excor
			
			
			tmp[0]=0;
	        for (int i=1; i<=N; i++)
	        {   
	           tmp[i]    = S[i]*rho[i]*Vxc[i];  //new why newrho?
	         }
	        double Excor2 = integrate(r,tmp);
			Excor = Excor0 + Excor2;
	    
	        // Enuclei = <-Z/r>
	        for (int i=0; i<=N; i++) tmp[i] = -Z*rho[i]*S[i]/r[i] ; 
	        tmp[0]=tmp[1]; 
	        Enuclei = integrate(r,tmp);
	    
	        // Epotential = <V_KS>
	        for (int i=0; i<=N; i++) tmp[i]  = S[i]*rho[i]*Veff0[i]; 
	        tmp[0]=tmp[1]; 
	        double Epotential = integrate(r,tmp);

	        // Total energy on output density
	        Etot  = Eb-Ehartree+Excor0 ;
	        Ekin  = Eb-Epotential;
	        Ecoul = Ehartree;	       	        
		}
		
		
		
		public double distance(double[] x, double[] y)
		{
			double dist =0;
			double denom = 0;
			if(x.length != y.length) System.out.println("method disantce : x and y must have same length");
			for(int i=0;i<x.length;i++)
			{
				dist += (x[i]-y[i])*(x[i]-y[i]);
				denom += x[i]*x[i];
			}
			//System.out.println(Math.sqrt(dist));
			return Math.sqrt(dist)/Math.sqrt(denom);
		}
		
		
/*		
		public void plot(double[] x, double[] y, String Title)
		{
			 // create your PlotPanel (you can use it as a JPanel)
		    Plot2DPanel plot = new Plot2DPanel();
		   
		    // add a line plot to the PlotPanel
		    plot.addLinePlot(Title, x, y);
		    //plot.addBaseLabel("text", Color.BLUE, 0) ;

		    // put the PlotPanel in a JFrame, as a JPanel
		    JFrame frame = new JFrame(Title);
		    frame.setSize(500, 500);
		    frame.setContentPane(plot);
		    frame.setVisible(true);
		}
*/		
		
		
}












