package sol

// EigenState : Hamiltonian eigen solution
type EigenState struct {
	n, l  int
	E     float64
	Value []float64
	h     float64 // r discretization step
}

// NewEigenEigenState : EigenState constructor
func NewEigenEigenState(E float64, n, l int, sol []float64, h float64) *EigenState {
	psi := Normalize(sol, h)
	state := EigenState{
		E:     E,
		n:     n,
		l:     l,
		Value: psi,
		h:     h,
	}
	return &state
}

// these n , l functions no necessary just need to find 
// a good name starting with Caps

// QuantumNumberN  n 
func (e *EigenState) QuantumNumberN() int {
	return e.n 
}

// QuantumNumberL l
func (e *EigenState) QuantumNumberL() int {
	return e.n 
}



//ProbabilityDensity : generate probabilty density
func (e *EigenState) ProbabilityDensity() []float64 {
	p := make([]float64, len(e.Value))
	for i := 0; i < len(e.Value); i++ {
		p[i] = e.Value[i] * e.Value[i]
	}
	return p
}
