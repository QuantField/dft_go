package atom

//import "fmt"

//============================= Radius =============================

// INFINITY represents infinity at the center r=0
// This is set in order to work with Potentials at r=0
// All infinite values are set to 0 at r=0, it doesn't affect the solution
// as we work with r*psi(r) which = 0 anyway
const INFINITY = 0.0 

// Radius reprsent r values from 0 to Rmax with step h and N+1 points
type Radius struct {
	Rmax  float64
	N     int
	GridStep   float64
	Array []float64
}

// NewRadius Radius constructor
func NewRadius(Rmax float64, N int) *Radius {

	dr := Rmax / float64(N)
	rad := make([]float64, N+1)
	for i := 0; i < len(rad); i++ {
		rad[i] = float64(i) * dr
	}

	r := Radius{
		Rmax:  Rmax,
		N:     N,
		GridStep:     dr,
		Array: rad,
	}
	return &r
}


//============================= Atom =============================

// Atom stuctue
type Atom struct {
	Z            float64
	r            []float64
	Vcoulomb     []float64
	Vcentrifugal []float64
	// Total Potential, Veff = Veff0 + Vcentrifugal
	Veff0 []float64 // Veff0 place holer for varying  Vks
	Veff  []float64
}

// NewAtom constructor
func NewAtom(Z float64, radius *Radius) *Atom {
	// r0 := make([]float64, len(r))
	// copy(r0, r) //deep copy
	r := radius.Array
	a := Atom{
		Z: Z,
		r: r,
	}
	a.Vcoulomb = make([]float64, len(r))
	a.Vcentrifugal = make([]float64, len(r))
	a.Veff0 = make([]float64, len(r))
	a.Veff = make([]float64, len(r))
	return &a
}

// SetVcoulomb  Sets the Coulomb Potential -Z/r
func (a *Atom) SetVcoulomb() {
	for i := 1; i < len(a.r); i++ {
		a.Vcoulomb[i] = -a.Z / a.r[i]
	}
	a.Vcoulomb[0] = INFINITY
}

// SetVcentrifugal Sets the centrifugal  Potential 0.5l(l+1)/r^2
func (a *Atom) SetVcentrifugal(j int) {
	l := float64(j)
	for i := 1; i < len(a.r); i++ {
		a.Vcentrifugal[i] = 0.5 * l * (l + 1) / (a.r[i] * a.r[i])
	}
	a.Vcentrifugal[0] = INFINITY
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
	for i := 0; i < len(a.r); i++ {
		a.Veff[i] = V[i] + a.Vcentrifugal[i]
	}
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
