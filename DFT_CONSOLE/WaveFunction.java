
	 // Psi_nl
	 // This class implements Comparable..for the sorting purposes
	 public class WaveFunction implements Comparable<WaveFunction>
	 {
		 int       n;
		 int       l;
		 double    Energy;
		 double[]  Psi;
		 double[]  Density;
		 double[]  DensityPerVolume;
		 
		 WaveFunction(double[] fun, double[] den, double[] denVol, double E, int k, int l)
		 {
			 int M  = fun.length;
			 n  = k;
			 Energy = E;
			 this.l = l;
			 Psi    = new double[M] ;
			 Density = new double[M];
			 DensityPerVolume = new double[M];

			 for(int i=0; i<M; i++)
			 {
				 Psi[i]     = fun[i];
				 Density[i] = den[i];
				 DensityPerVolume[i] = denVol[i]; 
			 }
		 }
		 
		 // it is better to use @override,..the compiler will stop you if
		 // we are not overriding the named method..
		 @Override public int compareTo(WaveFunction o2)
		 { 
		     int comp = 0;
		     if (this.n > o2.n) comp = 1;
		     if (this.n < o2.n) comp = -1;
		     return comp;	 
		 }  

	 }
	
	
