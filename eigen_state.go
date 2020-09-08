package main

// EigenState : Hamiltonian eigen solution
type EigenState struct {
	n, l  int
	E     float64
	value []float64
	h     float64 // r discretization step
}

// NewEigenEigenState : EigenState constructor
func NewEigenEigenState(E float64, n, l int, sol []float64, h float64) *EigenState {
	psi := Normalize(sol, h)
	state := EigenState{
		E:     E,
		n:     n,
		l:     l,
		value: psi,
		h:     h,
	}
	return &state
}

//ProbabilityDensity : generate probabilty density
func (e *EigenState) ProbabilityDensity() []float64 {
	p := make([]float64, len(e.value))
	for i := 0; i < len(e.value); i++ {
		p[i] = e.value[i] * e.value[i]
	}
	return p
}
