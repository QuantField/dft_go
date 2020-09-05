/*
 * Author Kamel Saadi
 * Date : 26 Sep 2013
 * 
 * This class is taken as is(almost) from Kristjan Haule.  
 */


class ExchangeCorrelation
{
//******************************************************************************/
//  Calculates Exchange&Correlation Energy and Potential                       */ 
//  type=0 - due to U.von.Barth and L.Hedin, J.Phys.C5, 1629 (1972)            */
//  type=1 - O.E.Gunnarsson and S.Lundqvist,  Phys.Rev.B                       */
//  type=2 - V.L.Moruzzi, J.F.Janak, and A.R.Williams, Calculated              */
//           Electronic Properties of Metals (New York, Pergamon Press, 1978)  */
//  type=3 - S.H.Vosko, L.Wilk, and M.Nusair, Can.J.Phys.58, 1200 (1980)       */
//  type=4 - Correlation of Perdew and Wang 1991                               */
//******************************************************************************/
  int type;
  double A, C;
  final  double alphax = 0.610887057710857;//(3/(2 Pi))^(2/3)
  final double Aw     = 0.0311;
  final double Bw     = -0.048;
  final double Cw     = 0.002;
  final double D      = -0.0116;
  final double gamma  = -0.1423;
  final double beta1  =  1.0529;
  final double beta2  =  0.3334;
  final double Ap     =  0.0621814;
  final double xp0    = -0.10498;
  final double bp     =  3.72744;
  final double cp     =  12.9352;
  final double Qp     =  6.1519908;
  final double cp1    =  1.2117833;
  final double cp2    =  1.1435257;
  final double cp3    = -0.031167608;

  public ExchangeCorrelation(int type)
  {
    this.type = type;
	switch(type)
    {
      case 0: C = 0.0504; A = 30; break;
      case 1: C = 0.0666; A = 11.4; break;
      case 2: C = 0.045;  A = 21; break;
    }
  };
  
  double sqr(double x){return x*x;}
  
  double Vx(double rs){return -alphax/rs;}
  
  double ExVx(double rs){return 0.25*alphax/rs;}
  
  double Ex(double rs){return -0.75*alphax/rs;}
  
  double Vc(double rs)
  {
    if (type<3)
    {
      double x = rs/A;
      return -0.5*C*Math.log(1+1/x);
    }else if(type<4){// type=3 WVN
      double x    =Math.sqrt(rs);
      double xpx  =x*x+bp*x+cp;
      double atnp =Math.atan2(Qp,(2*x+bp));
      double ecp = 0.5*Ap*(Math.log(x*x/xpx)+cp1*atnp-cp3*(Math.log(sqr(x-xp0)/xpx)+cp2*atnp));
      return ecp - Ap/6.*(cp*(x-xp0)-bp*x*xp0)/((x-xp0)*xpx);
    }else{
      if (rs>1) return gamma/(1+beta1*Math.sqrt(rs)+beta2*rs)*(1+7/6.*beta1*Math.sqrt(rs)+
    		           beta2*rs)/(1+beta1*Math.sqrt(rs)+beta2*rs);
      else return Aw*Math.log(rs)+Bw-Aw/3.+2/3.*Cw*rs*Math.log(rs)+(2*D-Cw)*rs/3.;
    }
  }
  double EcVc(double rs)
  {
    if (type<3){
      double x = rs/A;
      double epsilon = -0.5*C*((1+x*x*x)*Math.log(1+1/x)+0.5*x-x*x-1/3.);
      return epsilon-Vc(rs);
    } else if (type<4){ // type=3 WVN
      double x=Math.sqrt(rs);
      return Ap/6.*(cp*(x-xp0)-bp*x*xp0)/((x-xp0)*(x*x+bp*x+cp));
    }else{
      if (rs>1) return 2*gamma/(1+beta1*Math.sqrt(rs)+beta2*rs)-Vc(rs);
      else return Aw*Math.log(rs)+Bw+Cw*rs*Math.log(rs)+D*rs-Vc(rs);
    }
  }
};
