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

------------------------------------------------------------------------------*/

package main

import (
	"encoding/csv"
	"fmt"
	"math"
	"os"
)

type interval struct {
	lower float64
	upper float64
}

type Solution struct {
	h                      float64
	r                      []float64
	V                      []float64
	E                      float64
	K                      []float64
	Sol                    []float64
	boundary_element_index int // either r[0] or r[len(r)-1]
}

// Numerov is iterated backward. The goad is to tune E in order to make
// the solution =0 at r=0, i.e. u(r[boundary_element_index])=0
// where boundary_element_index =0
func NewSolution(radius *Radius, V []float64) *Solution {
	N := len(radius.array)
	s := Solution{
		h:                      radius.h,
		r:                      radius.array,
		V:                      V,
		K:                      make([]float64, N),
		Sol:                    make([]float64, N),
		boundary_element_index: 0,
	}
	return &s
}

// if E==0, s0 and s1 are used
// if not they are calculated based on some logic
func (s *Solution) SetInitialValues(E, s0, s1 float64) (float64, float64) {
	N := len(s.r)
	sol0, sol1 := s0*s.r[N-1], s1*s.r[N-2]
	if E != 0 || s0*s1 == 0 {
		alph := math.Sqrt(-2 * E)
		sol0 = s.r[N-1] * math.Exp(-alph*s.r[N-1])
		sol1 = s.r[N-2] * math.Exp(-alph*s.r[N-2])
	}
	return sol0, sol1
}

// if we start iterating fomr r[0] to r[n] ( n index of last element)
// then sol0 = Sol(r[0]), sol1 = Sol(r[1])
// if we start iterating from r[n] down to r[0] then
// sol0 = r[n], sol1 = r[n-1]
func (s *Solution) Numerov(E, sol0, sol1 float64) {

	N := len(s.r)
	C := s.h * s.h / 12

	for i := 1; i < N; i++ {
		s.K[i] = 2 * (E - s.V[i]) // in Hartree  or E0-2*V[i] for Ryleigs;
	}

	s.Sol[N-1], s.Sol[N-2] = sol0, sol1
	for i := N - 2; i >= 1; i-- {
		s.Sol[i-1] = (2*(1.-5.*C*s.K[i])*s.Sol[i] -
			(1.+C*s.K[i+1])*s.Sol[i+1]) / (1. + C*s.K[i-1])
	}
}

func (s *Solution) BoundaryValueResult(E float64) float64 {	
	sol0, sol1 := s.SetInitialValues(E, 0, 0)
	//sol0, sol1 := s.SetInitialValues(0, 7e-06, 7.5e-06) // not working
	// fmt.Println("E:", E , "sol0: ", sol0, " ", "sol1: ", sol1)
	s.Numerov(E, sol0, sol1)
	return s.Sol[s.boundary_element_index]
}

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

func (s *Solution) FindEigenState(Elow, Eupp float64) (float64, []float64) {
	psi := make([]float64, len(s.Sol))
	E, iter := Brent(s.BoundaryValueResult, Elow, Eupp)
	if iter != -1 {
		copy(psi, s.Sol)
		return E, psi
	} else {
		return 0, nil
	}
}

func SaveData(x, y []float64, filename string) {
	//S := make([][]string, len(x), len(x)+10)
	var S [][]string
	// S = append(S, []string{"# Solution for the Hydrogen atome n=1 l=0"})
	// S = append(S, []string{"#   r        r*psi"})
	// S = append(S, []string{"# ---------------------------------------"})
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
