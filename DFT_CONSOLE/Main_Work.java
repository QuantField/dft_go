
import java.util.ArrayList;
// https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations/atomic-reference-data-electronic-7
//import RadialWave.WaveFunction;

public class Main_Work 
{	
	 
	 public static void main(String[] args)
	 {
		 /*
		 LDA P = new LDA(2,10,1000);
		 P.plot(P.r, P.rho, "Density");
		 P.calulate_Vxc();
		 P.plot(P.r, P.Vxc, "Vxc");
		 P.solvePoisson();
		 P.plot(P.r, P.Uhartree, "Uhartree");	 
		 */
		 
		  
		 int    Z    = 2; 
		 double Rmax = 20; 
		 int    N    = 20000; // with 100 thousands is better 
		 
		 LDA P = new LDA(Z,Rmax,N);
		 System.out.println(" Z = "+Z);
		 P.mainDemo();
		 System.out.println("\n-------------- Finished ------------\n"); 
		 
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
 
	
}
