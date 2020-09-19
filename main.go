package main

import (
	"fmt"
	atom "quantum/atom"
	lab "quantum/sol"
)


func main() {

	N := 200
	// 10^-5..10^1
	grid := lab.NewLogGrid(-5, 1.2, N)
	// fmt.Println(grid.Val)
	// fmt.Println(grid.LogVal)
	// fmt.Println(grid.ValSqr )
	// fmt.Println(grid.ValSqrRoot)

	// Atom Z = 1
	H := atom.NewAtom(1, grid)
	H.SetVcoulomb()	
	H.SetVeff(nil)
	// i:=10
	// fmt.Println(H.Vcoulomb[i], H.Veff[i], grid.Val[i])
   

	sols := lab.NewSolution(grid, H.Veff)
	sols.Vcent = H.SetVcentrifugal(0)
	
	
	
	E := -0.5 // -0.5, -0.25
	//TODO  the graph looks fine, but cant find E via searching with Brent method
	sols.Numerov(E)
	u:= sols.GetOrigSol()
	u1:= lab.Normalize2(grid.Val, u)
	u1 = lab.Square(u1)
    SaveData(sols.Grid.Val, u1, "misc/psi_1.dat")


	E = -1.0/4.0 // -0.5, -0.25	
	sols.Numerov(E)
	u= sols.GetOrigSol()
	u1= lab.Normalize2(grid.Val, u)
	u1 = lab.Square(u1)
	SaveData(sols.Grid.Val, u1, "misc/psi_2.dat")
	
	E = -1.0/9.0 // -0.5, -0.25	
	sols.Numerov(E)
	u= sols.GetOrigSol()
	u1= lab.Normalize2(grid.Val, u)
	u1 = lab.Square(u1)
    SaveData(sols.Grid.Val, u1, "misc/psi_3.dat")







	// fmt.Println("Solution: ")

	// q := sol.FindSolutionIntervals(-0.6, -0.4, 0.01)

	// for key, interval := range q {
	// 	E, _ := sol.FindEigenState(interval.Lower, interval.Upper)
	// 	fmt.Println("n = ", key, " Solution interval", interval, "E =", E)
	// }
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
	fmt.Println("end")

}
