package main

import (
	"fmt"
	"math"
)

const PLUTO_NUM_STATES int = 51
const PLUTO_TIME_STEP int = 29200
const PLUTO_DT int = 146

const PLUTO_NSTEPS int = 201

type terse_vector_t struct {
	x float64
	y float64
	z float64
}

type body_grav_calc_t struct {
	tt float64        /* J2000 terrestrial time [days] */
	r  terse_vector_t /* position [au] */
	v  terse_vector_t /* velocity [au/day] */
	a  terse_vector_t /* acceleration [au/day^2] */
}

//name of precess_dir_t
const (
	FROM_2000 = 0
	INTO_2000 = 0
)

func main() {
	fmt.Println("going to start astronomy.go")
	fmt.Println(math.Pi)
	fmt.Println(PLUTO_NSTEPS)
}