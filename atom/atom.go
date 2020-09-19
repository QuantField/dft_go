package atom

import (
	
	sol "quantum/sol"

)

//import "fmt"



//============================= Atom =============================

// Atom stuctue
type Atom struct {
	Z            float64	
	Vcoulomb     []float64
	//Vcentrifugal []float64
	// Total Potential, Veff = Veff0 + Vcentrifugal
	Veff0 []float64 // Veff0 place holer for varying  Vks
	Veff  []float64
	grid *sol.LogGrid
}

// NewAtom constructor
func NewAtom(Z float64, grid *sol.LogGrid) *Atom {
	// r0 := make([]float64, len(r))
	// copy(r0, r) //deep copy
	
	a := Atom{
		Z: Z,
		grid: grid,
	}
	n := len(grid.Val)
	a.Vcoulomb = make([]float64, n)
	//a.Vcentrifugal = make([]float64, n)
	a.Veff0 = make([]float64, n)
	a.Veff = make([]float64, n)
	return &a
}

// SetVcoulomb  Sets the Coulomb Potential -Z/r
func (a *Atom) SetVcoulomb() {
	r := a.grid.Val 
	for i := 1; i < len(r); i++ {
		a.Vcoulomb[i] = -a.Z /r[i]
	}
	a.Vcoulomb[0] = 0 //INFINITY
}

// SetVcentrifugal Sets the centrifugal  Potential 0.5l(l+1)/r^2
/* Won't be used as it is included in Soltion Initial Values
func (a *Atom) SetVcentrifugal(j int) {
	l := float64(j)
	r2 := a.grid.ValSqr
	for i := 1; i < len(r2); i++ {
		a.Vcentrifugal[i] = 0.5*l*(l + 1)/r2[i]
	}
	a.Vcentrifugal[0] = 0 //INFINITY
}
*/

// SetVcentrifugal Constant in this case
func (a *Atom) SetVcentrifugal(j int) float64{
	l := float64(j)
	return (l+0.5)*(l+0.5)
}


// SetVeff Sets effective potential Veff = Veff0 + 0.5l(l+1)/r^2
// if Veff0 is nil then Veff0 becomes the Coulomb potential 
// already calculated and stored in Vcoulom. 
func (a *Atom) SetVeff(Veff0 []float64) {
	var V []float64
	V = Veff0	
	if Veff0 == nil {
		V = a.Vcoulomb
	}
	a.Veff = V 

	// for i := 0; i < len(r); i++ {
	// 	a.Veff[i] = V[i] + a.Vcentrifugal[i]
	// }
}


/*
func main() {
	//r := NewRadius(10, 100)
	Rmax:= 10.0
	N := 100

	atom:= NewAtom(1, NewRadius(Rmax, N).Array)
	//fmt.Println(atom.r)
	print_arr(atom.r)
	print_arr(nil)
}
*/
