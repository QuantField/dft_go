package sol

import (
	"math"
)

// LogGrid Creates necessary data for the radial
type LogGrid struct {
	LogRmax, LogRmin        float64
	Rmax, Rmin              float64
	N                       int
	Step                    float64
	Val, ValSqr, ValSqrRoot []float64
	LogVal                  []float64
}

// NewLogGrid Radius constructor
func NewLogGrid(logRmin, logRmax float64, N int) *LogGrid {
	dr := (logRmax - logRmin) / float64(N)
	Rmin := math.Pow(10, logRmin)
	Rmax := math.Pow(10, logRmax)
	r := make([]float64, N+1)
	rLog := make([]float64, N+1)
	rSqr := make([]float64, N+1)
	rSqrRoot := make([]float64, N+1)
	v := logRmin
	for i := 0; i < len(r); i++ {
		r[i] = math.Pow(10, v)
		rLog[i] = v
		rSqr[i] = r[i] * r[i]
		rSqrRoot[i] = math.Sqrt(r[i])
		v += dr
	}
	gr := LogGrid{
		LogRmax:    logRmax,
		LogRmin:    logRmin,
		Rmax:       Rmax,
		Rmin:       Rmin,
		N:          N,
		Step:       dr,
		Val:        r,
		ValSqr:     rSqr,
		ValSqrRoot: rSqrRoot,
		LogVal:     rLog,
	}
	return &gr
}

/*
func main() {

	grid := NewLogGrid(-5, 2, 20)
	for i := 0; i < len(grid.Val); i++ {
		fmt.Printf("%4.6f  %4.6f  %4.6f  %4.6f \n ", grid.Val[i], grid.LogVal[i], grid.ValSqr[i], grid.ValSqrRoot[i])
	}
}
*/
