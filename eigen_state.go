package main

type EigenState struct {
	n, l  int
	E     float64
	value []float64
	h     float64 // r discretization step
}

func NewEigenEigenState(E float64, n, l int, sol []float64, h float64) *EigenState {
	psi := make([]float64, len(sol))
	copy(psi, sol)
	state := EigenState{
		E:     E,
		n:     n,
		l:     l,
		value: psi,
		h: h,
	}
	return &state
}

func (e *EigenState) normalize() {
	e.value = normalize(e.value, e.h)
}
