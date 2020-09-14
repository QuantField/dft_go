package main

import (
	"fmt"
	"math"
	atom "quantum/atom"
	sol "quantum/sol"
)




func main() {

	N := 1000
	// 10^-5..10^1
	grid := sol.NewLogGrid(-5,1, N)
	// fmt.Println(grid.Val)
	// fmt.Println(grid.LogVal)
	// fmt.Println(grid.ValSqr )
	// fmt.Println(grid.ValSqrRoot)


	
	

	// Atom Z = 1
	H := atom.NewAtom(1, grid)
	H.SetVcoulomb()
	H.SetVcentrifugal(0) // l=0
	H.SetVeff(nil)

	sol := sol.NewSolution(radius, H.Veff)

	fmt.Println("Solution: ")

	q := sol.FindSolutionIntervals(-10, -0.0001, .01)

	for key, interval := range q {
		E, _ := sol.FindEigenState(interval.Lower, interval.Upper)
		fmt.Println("n = ", key, " Solution interval", interval, "E =", E)
	}
	//fmt.Println(q)
	/*
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
	*/
 
}
