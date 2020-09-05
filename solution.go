package main

import (
	"fmt"
	"math"
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

func (s *Solution) Numerov(E float64) {

	N := len(s.r)
	C := s.h * s.h / 12

	for i := 1; i < N; i++ {
		s.K[i] = 2 * (E - s.V[i]) // in Hartree  or E0-2*V[i] for Ryleigs;
	}

	s.index = 0 // because we are sweeping from the high end, we are starting
	// from r=Rmax and coming down to r =0
	// if we are starting from r=0 then s.index = N-1
	// s.index is used in BoundaryValueResult
	alph := math.Sqrt(-2 * E)
	s.Sol[N-1] = s.r[N-1] * math.Exp(-alph*s.r[N-1])
	s.Sol[N-2] = s.r[N-2] * math.Exp(-alph*s.r[N-2])

	for i := N - 2; i >= 1; i-- {
		s.Sol[i-1] = (2*(1.-5.*C*s.K[i])*s.Sol[i] -
			(1.+C*s.K[i+1])*s.Sol[i+1]) / (1. + C*s.K[i-1])
	}
}

func (s *Solution) BoundaryValueResult(E float64) float64 {
	// in our case index = 0
	s.Numerov(E)
	return s.Sol[s.index]
}

func (s *Solution) FindSolutionIntervals(Emin, Emax, dE float64) map[int]interval {

	solIntervals := make(map[int]interval)
	n := 1
	for val := Emin; val < Emax; val += dE {
		if s.BoundaryValueResult(val)*s.BoundaryValueResult(val+dE) < 0 {
			solIntervals[n] = interval{val, val + dE}
			n++
		}
	}
	return solIntervals
}

func (s *Solution) FindEigenState(Elow, EUpp float64) (float64, []float64) {
	psi := make([]float64, len(s.Sol))
	E, iter := brent(s.BoundaryValueResult, Elow, EUpp)
	if iter != -1 {
		copy(psi, s.Sol)
		return E, psi
	} else {
		return 0, nil
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

	sol.Numerov(-0.5)

	fmt.Println("Solution: ")

	q := sol.FindSolutionIntervals(-10, 0, 1)
	fmt.Println(q)
	E, psi := sol.FindEigenState(-1, 0)
	if psi != nil {
		fmt.Println("E = ", E)
		fmt.Println("Psi = ", psi[0:5])
		fmt.Println("Psi = ", psi[len(psi)-5:])
	}
}
