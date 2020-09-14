package sol

import "math"

// EigenState : Hamiltonian eigen solution u(r) = r*psi(r)
type EigenState struct {
	Nnumber, Lnumber  int
	E     float64
	Value []float64
	h     float64 // r discretization step
}

// NewEigenEigenState : EigenState constructor
func NewEigenEigenState(E float64, n, l int, sol []float64, h float64) *EigenState {
	psi := Normalize(sol, h)
	state := EigenState{
		E:     E,
		Nnumber:     n,
		Lnumber:     l,
		Value: psi,
		h:     h,
	}
	return &state
}

//ProbabilityDensity : generates radial probabilty density
func (e *EigenState) ProbabilityDensity() []float64 {
	p := make([]float64, len(e.Value))
	for i := 0; i < len(e.Value); i++ {
		p[i] = e.Value[i] * e.Value[i]
	}
	return p
}

//Density : generates volume probabilty density
func (e *EigenState) Density(r []float64 ) []float64 {
	rho := make([]float64, len(e.Value))
	denom:= 4*math.Pi
	rho[0]=0
	for i := 1; i < len(e.Value); i++ {
		psi:= e.Value[i]/r[i]
		rho[i] = psi*psi/denom
	}
	return rho
}


