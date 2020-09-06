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
	h     float64
	r     []float64
	V     []float64
	E     float64
	K     []float64
	Sol   []float64
	index int
}

func NewSolution(radius *Radius, V []float64) *Solution {
	N := len(radius.array)
	s := Solution{
		h:     radius.h,
		r:     radius.array,
		V:     V,
		K:     make([]float64, N),
		Sol:   make([]float64, N),
		index: 0,
	}
	return &s
}

// if E==0 s0 and s1 are used
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

// if we start iterating fomr r[0] to r[n] :(index of last element)
// sol0 = Sol(r[0]), sol1 = Sol(r[1])
// if we start iterating from r[n] down to r[0] :
// sol0 = r[n], sol1 = r[n-1]
func (s *Solution) Numerov(E, sol0, sol1 float64) {

	N := len(s.r)
	C := s.h * s.h / 12

	for i := 1; i < N; i++ {
		s.K[i] = 2 * (E - s.V[i]) // in Hartree  or E0-2*V[i] for Ryleigs;
	}

	s.index = 0 // because we are sweeping from the high end, we are starting
	// from r=Rmax and coming down to r =0
	// if we are starting from r=0 then s.index = N-1, in which case the for loop
	// must be reversed as well.
	// s.index is used in BoundaryValueResult

	s.Sol[N-1], s.Sol[N-2] = sol0, sol1
	for i := N - 2; i >= 1; i-- {
		s.Sol[i-1] = (2*(1.-5.*C*s.K[i])*s.Sol[i] -
			(1.+C*s.K[i+1])*s.Sol[i+1]) / (1. + C*s.K[i-1])
	}
}

func (s *Solution) BoundaryValueResult(E float64) float64 {
	// in our case index = 0
	sol0, sol1 := s.SetInitialValues(E, 0, 0)
	//sol0, sol1 := s.SetInitialValues(0, 7e-06, 7.5e-06) // not working
	// fmt.Println("E:", E , "sol0: ", sol0, " ", "sol1: ", sol1)
	s.Numerov(E, sol0, sol1)
	return s.Sol[s.index]
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
	E, iter := brent(s.BoundaryValueResult, Elow, Eupp)
	if iter != -1 {
		copy(psi, s.Sol)
		return E, psi
	} else {
		return 0, nil
	}
}

func GenerateTextToSave(x, y []float64, filename string) {
	//S := make([][]string, len(x), len(x)+10)
	var S [][]string
	// S = append(S, []string{"# Solution for the Hydrogen atome n=1 l=0"})
	// S = append(S, []string{"#   r        r*psi"})
	// S = append(S, []string{"# ---------------------------------------"})
	for i := 0; i < len(x); i++ {
		q := fmt.Sprintf("%f \t %f", x[i], y[i])
		//S[i] = []string{q}
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

func main() {

	// defining the ratius with max length of 10 units
	// grid of N+1 points
	Rmax := 10.0
	N := 1000
	radius := NewRadius(Rmax, N)

	// Atom Z = 1
	atom := NewAtom(1, radius)
	atom.SetVcoulomb()
	atom.SetVcentrifugal(0) // l=0
	atom.SetVeff(nil)

	sol := NewSolution(radius, atom.Veff)

	fmt.Println("Solution: ")

	q := sol.FindSolutionIntervals(-10, 0, 1)
	fmt.Println(q)
	E, psi := sol.FindEigenState(-1, 0)
	if psi != nil {
		fmt.Println("E = ", E)
		fmt.Println("Psi = ", psi[0:5])
		fmt.Println("Psi = ", psi[len(psi)-5:])
	}

	GenerateTextToSave(sol.r, sol.Sol, "psi.dat")
	GenerateTextToSave(sol.r, sol.V, "pot.dat")

}
