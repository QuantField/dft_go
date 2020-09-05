## Numerov Method for Shordinger Equation with Radial Potential 
* Testing solution.go
```

func main() {

	// defining the ratius with max length of 10 units
	// grid of N+1 points
	Rmax := 10.0
	N := 1000
	radius := NewRadius(Rmax, N)

	// Atom Z = 1
	atom := NewAtom(1, radius)
	atom.SetVcoulomb()
	atom.SetVcentrifugal(0) // l=0
	atom.SetVeff(nil)

	sol := NewSolution(radius, atom.Veff)

	sol.Numerov(-0.5)

	fmt.Println("Solution: ")

	q := sol.FindSolutionIntervals(-10, 0, 1)
	fmt.Println(q)
	E, psi := sol.FindEigenState(-1, 0)
	if psi != nil {
		fmt.Println("E = ", E)
		fmt.Println("Psi = ", psi[0:5])
		fmt.Println("Psi = ", psi[len(psi)-5:])
	}
}

-------------------------------------------------------------------------------------------------------------------
go build solution.go atom.go go_numerics.go

Solution:
map[1:{-1 0}]
E =  -0.49996700459597737
Psi =  [1.2587153541687712e-14 0.009881949702694893 0.019583866984099516 0.029091909927566852 0.038408914952632656]
Psi =  [0.00047079194465100165 0.0004665756207883756 0.0004623965920504494 0.0004582545336864157 0.000454149123715689]

```
