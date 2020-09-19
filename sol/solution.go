/*------------------------------------------------------------------------------
 Using Numerov Algorithm
 https://en.wikipedia.org/wiki/Numerov%27s_method#:~:text=Numerov's%20method%20(
	 also%20called%20Cowell's,the%20differential%20equation%20is%20linear.

to solve Shrodinger Equation for Central Radial Potential
Using variable subsitution  u(r) = r*Psi(r)

{-alp*d^2/dr^2 + [V(r) + alp*l(l+1)/r^2]}u(r) = Eu(r)

where
alp= -h_bar^2/(2*m)

or

{-alp*d^2/dr^2 + Veff(r) }u(r) = Eu(r)

where  Veff(r) = V(r) + alp*l(l+1)/r^2

In a more practical form for calculation

(d^2/dr^2)u(r) = (2*m/h_bar^2)(E - Veff(r))u(r)

in E units (Hart) :

   u"(r) = 2(E-Veff)u(r)  --- (1)

Now if we use logarithmic grid:

	 x = log(Zr) ,   log is base 10
	 and y = log(u(r))/sqrt(r)
	 
(1) becomes 

  y"(x)  + [2*r^2*(E-V(r)) - (l+1/2)^2]*y(x) = 0

  Here V(r) excludes centrifugel potential


------------------------------------------------------------------------------*/

package sol

import (
	"encoding/csv"
	"fmt"
	"math"
	"os"	
)

type interval struct {
	Lower float64
	Upper float64
}

// Solution contains the necessary data to use Numerov algorithm
type Solution struct {	
	Grid            *LogGrid
	GridStep        float64  // h = dr
	R               []float64
	V               []float64
	Vcent           float64
	E               float64
	K               []float64
	Sol             []float64
	sol0, sol1      float64
}

// NewSolution create a Solution object
// Numerov is iterated backward. The goad is to tune E in order to make
// the solution =0 at r=0, i.e. u(r[boundaryElementIndex])=0
// where boundaryElementIndex =0
func NewSolution( grid *LogGrid, V []float64) *Solution {
	N := len(grid.Val)
	s := Solution{
		Grid: grid, 
		GridStep: grid.Step,
		R:        grid.Val,
		V:        V,		
		K:        make([]float64, N),
		Sol:      make([]float64, N),		
		sol0:  0,
		sol1:  0,
	}
	return &s
}

// InitialConditions 
func (s *Solution) InitialConditions(E float64) {	
	N := len(s.R)
	Rsqr:= s.Grid.ValSqr
	for i := 1; i < N; i++ {
		s.K[i] = 2*Rsqr[i]*(E - s.V[i]) - s.Vcent// in Hartree  or E0-2*V[i] for Ryleigs;
	}	
	alph := math.Sqrt(-2 * E)
	s.sol0 = math.Sqrt(s.R[N-1]) * math.Exp(-alph*s.R[N-1])
	s.sol1 = math.Sqrt(s.R[N-2]) * math.Exp(-alph*s.R[N-2])
}


// Numerov  Method for solving y'' = A*y
// if we start iterating fomr r[0] to r[n] ( n index of last element)
// then sol0 = Sol(r[0]), sol1 = Sol(r[1])
// if we start iterating from r[n] down to r[0] then
// sol0 = r[n], sol1 = r[n-1]
func (s *Solution) Numerov(E float64) {

	s.InitialConditions(E)
	
	h := s.Grid.Step
	N := len(s.R)
	C := h*h / 12	
	
	s.Sol[N-1], s.Sol[N-2] = s.sol0, s.sol1
	for i := N - 2; i >= 1; i-- {
		s.Sol[i-1] = (2*(1.-5.*C*s.K[i])*s.Sol[i] -
			(1.+C*s.K[i+1])*s.Sol[i+1]) / (1. + C*s.K[i-1])
	}
}

// GetOrigSol transforma back to u(r) by mutlipying by sqrt(r)
func (s *Solution) GetOrigSol() []float64{
	u := make([]float64, len(s.Sol))
	for i:=0; i<len(u); i++ {
		u[i] = s.Grid.ValSqrRoot[i]*s.Sol[i]
	}
	return u
}



// BoundaryValueResult Returns the last computed value
// This value  should be =0 at r = 0
func (s *Solution) BoundaryValueResult(E float64) float64 {		
	s.Numerov(E)
	return s.Sol[0]
}

// FindSolutionIntervals checks for solutions between Emin an Emax 
// search is done by steps of length dE
func (s *Solution) FindSolutionIntervals(Emin, Emax, dE float64) map[int]interval {
	solIntervals := make(map[int]interval)
	n := 1
	for val := Emin; val < Emax; val += dE {
		// fmt.Println("------ Searching in ", interval{val, val + dE}, "..." )
		if s.BoundaryValueResult(val)*s.BoundaryValueResult(val+dE) < 0 {
			solIntervals[n] = interval{val, val + dE}
			n++
		}
	}
	return solIntervals
}

// FindEigenState given an interval [Elow, Eupp] that contains a solution 
// (f(E)=0) this function finds the exact value of E
func (s *Solution) FindEigenState(Elow, Eupp float64) (float64, []float64) {
	psi := make([]float64, len(s.Sol))
	E, iter := Brent(s.BoundaryValueResult, Elow, Eupp)
	if iter != -1 {
		copy(psi, s.Sol)
		return E, psi
	} 
	return 0, nil
	
}

// SaveData save on the form r, f(r)
func SaveData(x, y []float64, filename string) {	
	var S [][]string	
	for i := 0; i < len(x); i++ {
		q := fmt.Sprintf("%f \t %f", x[i], y[i])		
		S = append(S, []string{q})
	}
	file, _ := os.Create(filename)
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	err := writer.WriteAll(S)
	if err != nil {
		fmt.Println("error")
	}
}
