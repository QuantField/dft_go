package main

//import "fmt"

//============================= Radius =============================

const INFINITY = 0.0 // this is to use when deviding by 0
// doesn't matter as we work with r*f(r)
// solutions never use V(r=0)

type Radius struct {
	Rmax  float64
	N     int
	h     float64
	array []float64
}

func NewRadius(Rmax float64, N int) *Radius {

	dr := Rmax / float64(N)
	rad := make([]float64, N+1)
	for i := 0; i < len(rad); i++ {
		rad[i] = float64(i) * dr
	}

	r := Radius{
		Rmax:  Rmax,
		N:     N,
		h:     dr,
		array: rad,
	}
	return &r
}

//============================= Atom =============================

type Atom struct {
	Z            float64
	r            []float64
	Vcoulomb     []float64
	Vcentrifugal []float64
	// Total Potential, Veff = Veff0 + Vcentrifugal
	Veff0 []float64 // Veff0 place holer for varying  Vks
	Veff  []float64
}

func NewAtom(Z float64, radius *Radius) *Atom {
	// r0 := make([]float64, len(r))
	// copy(r0, r) //deep copy
	r := radius.array
	a := Atom{
		Z: Z,
		r: r,
	}

	a.Vcoulomb = make([]float64, len(r))
	a.Vcentrifugal = make([]float64, len(r))
	a.Veff0 = make([]float64, len(r))
	a.Veff = make([]float64, len(r))
	//a.SetVcoulomb()
	//a.SetVeff(nil) // start initially with Coulomb Potential
	return &a
}

func (a *Atom) SetVcoulomb() {
	for i := 1; i < len(a.r); i++ {
		a.Vcoulomb[i] = -a.Z / a.r[i]
	}
	a.Vcoulomb[0] = INFINITY
}

func (a *Atom) SetVcentrifugal(j int) {
	l := float64(j)
	for i := 1; i < len(a.r); i++ {
		a.Vcentrifugal[i] = 0.5 * l * (l + 1) / (a.r[i] * a.r[i])
	}
	a.Vcentrifugal[0] = INFINITY
}

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

// func main(){
//    const N = 20
//    r:= make([]float64, N)
//    for i:=0; i<N; i++ {
// 	   r[i] = float64(i)*0.01 + 0.01
//    }

//    a:= NewAtom(1, r)
// //    fmt.Println(a.r)
// //    fmt.Println(a.Vcoulomb)
// //    fmt.Println(a.Veff0)
//    r[10] = 11111
//    q:=r
//    q[10] = 22222
//    //fmt.Println(&a.r,"   ", &r, "  ", &q)
//    fmt.Println(a.r[10], "  ", r[10])

// }

// func print_arr(x []float64){
// 	if (x == nil) {
// 		fmt.Println("can't print")
// 	}else 	{
// 	fmt.Println(x)
// 	}
// }

/*
func main() {
	//r := NewRadius(10, 100)
	Rmax:= 10.0
	N := 100

	atom:= NewAtom(1, NewRadius(Rmax, N).array)
	//fmt.Println(atom.r)
	print_arr(atom.r)
	print_arr(nil)



}
*/
