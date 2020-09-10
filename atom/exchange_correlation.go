/*
 * Author Kamel Saadi
 * Date :  Sep 2020
 *
 * This class is taken as is(almost) from Kristjan Haule.
 */

//******************************************************************************/
//  Calculates Exchange&Correlation Energy and Potential                       */
//  type=0 - due to U.von.Barth and L.Hedin, J.Phys.C5, 1629 (1972)            */
//  type=1 - O.E.Gunnarsson and S.Lundqvist,  Phys.Rev.B                       */
//  type=2 - V.L.Moruzzi, J.F.Janak, and A.R.Williams, Calculated              */
//           Electronic Properties of Metals (New York, Pergamon Press, 1978)  */
//  type=3 - S.H.Vosko, L.Wilk, and M.Nusair, Can.J.Phys.58, 1200 (1980)       */
//  type=4 - Correlation of Perdew and Wang 1991                               */
//******************************************************************************/

package atom

import (
	"log"
	"math"
	"strconv"
)

const (
	alphax = 0.610887057710857 //(3/(2 Pi))^(2/3)
	Aw     = 0.0311
	Bw     = -0.048
	Cw     = 0.002
	D      = -0.0116
	gamma  = -0.1423
	beta1  = 1.0529
	beta2  = 0.3334
	Ap     = 0.0621814
	xp0    = -0.10498
	bp     = 3.72744
	cp     = 12.9352
	Qp     = 6.1519908
	cp1    = 1.2117833
	cp2    = 1.1435257
	cp3    = -0.031167608
)

// ExchangeCorrelation 
type ExchangeCorrelation struct {
	C, A  float64
	etype int
}

// NewExchangeCorrelation
func NewExchangeCorrelation(typ int) *ExchangeCorrelation {
	res := []ExchangeCorrelation{
		{0.0504, 30, typ},
		{0.0666, 11.4, typ},
		{0.045, 21, typ},
		{0, 0, typ},
		{0, 0, typ},
	}
	if typ > len(res)-1 {
		log.Printf(":[Error] etype should be in [0,4] found " + strconv.Itoa(typ))
		return nil
	} else {
		return &res[typ]
	}
}

// Vx
func (e *ExchangeCorrelation) Vx(rs float64) float64 {
	return -alphax / rs
}

// ExVx
func (e *ExchangeCorrelation) ExVx(rs float64) float64 {
	return 0.25 * alphax / rs
}

// Ex
func (e *ExchangeCorrelation) Ex(rs float64) float64 {
	return -0.75 * alphax / rs
}

// Vc
func (e *ExchangeCorrelation) Vc(rs float64) float64 {

	if e.etype < 3 {
		x := rs / e.A
		return -0.5 * e.C * math.Log(1+1/x)
	} else if e.etype < 4 { // type=3 WVN
		x := math.Sqrt(rs)
		xpx := x*x + bp*x + cp
		atnp := math.Atan2(Qp, (2*x + bp))
		ecp := 0.5 * Ap * (math.Log(x*x/xpx) +
			cp1*atnp -
			cp3*(math.Log((x-xp0)*(x-xp0)/xpx)+cp2*atnp))
		return ecp - Ap/6.*(cp*(x-xp0)-bp*x*xp0)/((x-xp0)*xpx)
	} else {
		if rs > 1 {
			return gamma / (1 + beta1*math.Sqrt(rs) + beta2*rs) *
				(1 + 7/6.*beta1*math.Sqrt(rs) + beta2*rs) /
				(1 + beta1*math.Sqrt(rs) + beta2*rs)
		} else {
			return Aw*math.Log(rs) + Bw - Aw/3. + 2/3.*Cw*rs*math.Log(rs) + (2*D-Cw)*rs/3.
		}
	}
}

// EcVc
func (e *ExchangeCorrelation) EcVc(rs float64) float64 {
	if e.etype < 3 {
		x := rs / e.A
		epsilon := -0.5 * e.C * ((1+x*x*x)*math.Log(1+1/x) + 0.5*x - x*x - 1/3.)
		return epsilon - e.Vc(rs)
	} else {
		if e.etype < 4 { // type=3 WVN
			x := math.Sqrt(rs)
			return Ap / 6. * (cp*(x-xp0) - bp*x*xp0) / ((x - xp0) * (x*x + bp*x + cp))
		} else {
			if rs > 1 {
				return 2*gamma/(1+beta1*math.Sqrt(rs)+beta2*rs) - e.Vc(rs)
			} else {
				return Aw*math.Log(rs) + Bw + Cw*rs*math.Log(rs) + D*rs - e.Vc(rs)
			}
		}
	}
}

/*
func main() {

	for _, i := range []int{0, 1, 2, 3, 4, 5} {
		e := NewExchangeCorrelation(i)
		if e != nil {
			fmt.Println(e)
			fmt.Println(e.EcVc(1), e.Ex(1), e.ExVx(1), e.Vc(1), e.Vx(1))
		}
	}
}
*/
