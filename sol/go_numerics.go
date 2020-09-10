package sol

import (
	"fmt"
	"log"
	"math"
)

type realValuedFunction func(float64) float64

// Secant : Brent method for finding zeros of a function
func Secant(f realValuedFunction, a, b float64) (float64, int) {
	epsilon, MaxIter := 1.0E-6, 50
	fa, fb := f(a), f(b)
	iter := 2
	x2, f2 := b, 0.0
	x1, f1 := a, 0.0
	if fa*fb > 0 {
		fmt.Println("No solution")
		return 0.0, -1
	}
	for math.Abs((x1-x2)/x2) > epsilon && iter <= MaxIter {
		if x2 < a || x2 > b || iter == 2 {
			x1, f1 = a, fa
			xm := (a + b) / 2
			iter++
			fxm := f(xm)
			if fa*fxm < 0 {
				b, fb = xm, fxm
			} else {
				a, fa = xm, fxm
			}
			x2, f2 = xm, fxm
		}
		iter++
		xv, fv := x2, f2
		x2 += -(x2 - x1) * f2 / (f2 - f1)
		f2 = f(x2)
		x1, f1 = xv, fv
		fmt.Printf("x = %4.8f  f(x) = %4.8f \n", x2, f2)
	}
	return x2, iter

}

// Brent : Brent method for finding zeros of a function
func Brent(f realValuedFunction, a, b float64) (float64, int) {
	epsilon, MaxIter := 1.0E-6, 50
	fa, fb := f(a), f(b)
	iter := 2

	x1, f1, x2, f2 := a, fa, b, fb
	x3 := x2 + 10
	f3 := 0.0
	if fa*fb > 0 {
		//print("No solution")
		return 0.0, -1
	}
	for math.Abs((x2-x3)/x2) > epsilon && iter <= MaxIter {
		if x3 < a || x3 > b || iter == 2 {
			x1, f1, x2, f2 = a, fa, b, fb
			xm := (a + b) / 2
			iter++
			fxm := f(xm)
			if fa*fxm < 0 {
				b, fb = xm, fxm
			} else {
				a, fa = xm, fxm
			}
			x3 = (a + b) / 2.
			iter++
			f3 = f(x3)
		}
		iter++
		xv, fv := x3, f3
		x3 = f1*f2*x3/((f3-f2)*(f3-f1)) +
			f2*f3*x1/((f1-f2)*(f1-f3)) +
			f3*f1*x2/((f2-f1)*(f2-f3))
		f3 = f(x3)
		x1, f1 = x2, f2
		x2, f2 = xv, fv
		//fmt.Printf("x = %4.8f  f(x) = %4.8f \n", x3, f3)
	}
	return x3, iter
}

// Integrate : Composite Simpson integration
func Integrate(y []float64, h float64) float64 {
	N := len(y) - 1
	if N%2 != 0 {
		log.Panic("method integrate error: length of array  must be odd.")
	}
	S := y[0] + y[N]
	for i := 1; i <= N; i = i + 2 {
		S += 4 * y[i]
	}
	for i := 2; i <= N-1; i = i + 2 {
		S += 2 * y[i]
	}
	return h * S / 3
}

// Normalize : integ psi*psi = 1
func Normalize(y []float64, h float64) []float64 {
	N := len(y) - 1
	yNorm := make([]float64, len(y))
	if N%2 != 0 {
		log.Panic("method integrate error: length of array  must be odd.")
	}
	S := y[0]*y[0] + y[N]*y[N]
	for i := 1; i <= N; i = i + 2 {
		S += 4 * y[i] * y[i]
	}
	for i := 2; i <= N-1; i = i + 2 {
		S += 2 * y[i] * y[i]
	}
	A := h * S / 3
	A = math.Sqrt(A)
	for i := 0; i <= N; i++ {
		yNorm[i] = y[i] / A
	}
	return yNorm
}

/* func main() {

	cubic := func(x float64) float64 {
		return x*x*x - 2
	}
	s, d := -1.0, 3.0

	secant(cubic, s, d)
	fmt.Println()
	brent(cubic, s, d)

	// testing simpson integration


	y := []float64{0, 1, 4, 9, 16}
	fmt.Println()
	fmt.Println(integrate(y, 1)) // 4^3/3

}
*/
