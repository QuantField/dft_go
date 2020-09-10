package main

import (
	"fmt"
	 atom "quantum/atom"
	 solEng  "quantum/sol"
)

func main() {

	// defining the ratius with max length of 10 units
	// grid of N+1 points
	Rmax := 10.0
	N := 4000
	radius := atom.NewRadius(Rmax, N)

	// Atom Z = 1
	H := atom.NewAtom(1, radius)
	H.SetVcoulomb()
	H.SetVcentrifugal(0) // l=0
	H.SetVeff(nil)

	sol := solEng.NewSolution(radius, H.Veff)

	fmt.Println("Solution: ")

	q := sol.FindSolutionIntervals(-10, 0, 1)
	fmt.Println(q)
	E, psi := sol.FindEigenState(-1, 0)
	if psi != nil {
		fmt.Println("E = ", E)
		fmt.Println("Psi = ", psi[0:2])
		fmt.Println("Psi = ", psi[len(psi)-2:])

		u := solEng.NewEigenEigenState(E, 1 ,0, psi,  sol.GridStep)
		
		solEng.SaveData(sol.R, u.Value, "misc/psi.dat")
		solEng.SaveData(sol.R, u.ProbabilityDensity(), "misc/dens.dat")
		solEng.SaveData(sol.R, sol.V, "misc/pot.dat")
		fmt.Println(solEng.Integrate(u.ProbabilityDensity(), sol.GridStep ) )
	}

}
