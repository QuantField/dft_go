package sol

// EigenState : Hamiltonian eigen solution
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

//ProbabilityDensity : generates probabilty density
func (e *EigenState) ProbabilityDensity() []float64 {
	p := make([]float64, len(e.Value))
	for i := 0; i < len(e.Value); i++ {
		p[i] = e.Value[i] * e.Value[i]
	}
	return p
}
