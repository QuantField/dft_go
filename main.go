package main

import "fmt"

func main() {

	// defining the ratius with max length of 10 units
	// grid of N+1 points
	Rmax := 10.0
	N := 4000
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
		fmt.Println("Psi = ", psi[0:2])
		fmt.Println("Psi = ", psi[len(psi)-2:])

		u := NewEigenEigenState(E, 1 ,0, psi, sol.h)
		
		SaveData(sol.r, u.value, "misc/psi.dat")
		SaveData(sol.r, u.ProbabilityDensity(), "misc/dens.dat")
		SaveData(sol.r, sol.V, "misc/pot.dat")
		fmt.Println(Integrate(u.ProbabilityDensity(), sol.h ) )
	}

}
