package astronomy

import (
	"runtime"
	"testing"
)

func gravSimSimulate(t *testing.T) StateVector {
	t0 := TimeFromCalendar(2020, 1, 1, 0, 0, 0.0)
	body := StateVector{X: 1.5, Y: 0.0, Z: 0.0, Vx: 0.0, Vy: 0.017, Vz: 0.0, T: t0}
	sim, err := NewGravSim(Ssb, t0, 1, []StateVector{body})
	if err != nil {
		t.Fatal(err)
	}
	out := make([]StateVector, 1)
	t1 := TimeFromCalendar(2020, 1, 2, 0, 0, 0.0)
	if err := sim.Update(t1, 1, out); err != nil {
		t.Fatal(err)
	}
	t2 := TimeFromCalendar(2020, 1, 11, 0, 0, 0.0)
	if err := sim.Update(t2, 1, out); err != nil {
		t.Fatal(err)
	}
	return out[0]
}

func TestGravSim(t *testing.T) {
	// Reference from the C implementation.
	state := gravSimSimulate(t)
	if diff := absDiff(state.X, 1.4934812762399952e+00); diff > 1.0e-12 {
		t.Errorf("GravSim x error = %g", diff)
	}
	if diff := absDiff(state.Y, 1.6969407525209854e-01); diff > 1.0e-12 {
		t.Errorf("GravSim y error = %g", diff)
	}
	if diff := absDiff(state.Vx, -1.3032131698856176e-03); diff > 1.0e-15 {
		t.Errorf("GravSim vx error = %g", diff)
	}
	if diff := absDiff(state.Vy, 1.6932189338763851e-02); diff > 1.0e-15 {
		t.Errorf("GravSim vy error = %g", diff)
	}
}

func TestGravSimParallelDeterministic(t *testing.T) {
	// The parallel acceleration computation must produce bit-identical results
	// regardless of the number of goroutines used.
	serial := gravSimSimulate(t)
	old := runtime.GOMAXPROCS(1)
	one := gravSimSimulate(t)
	runtime.GOMAXPROCS(old)
	// Compare the numeric fields exactly (the AstroTime cache fields are NaN
	// and never compare equal).
	if serial.X != one.X || serial.Y != one.Y || serial.Z != one.Z ||
		serial.Vx != one.Vx || serial.Vy != one.Vy || serial.Vz != one.Vz {
		t.Errorf("GravSim result differs between GOMAXPROCS=1 and default: serial=%+v gomax1=%+v", serial, one)
	}
}

func TestGravSimMultiBody(t *testing.T) {
	// Simulate two interacting bodies.
	t0 := TimeFromCalendar(2020, 1, 1, 0, 0, 0.0)
	states := []StateVector{
		{X: 1.0, Y: 0.0, Z: 0.0, Vx: 0.0, Vy: 0.012, Vz: 0.0, T: t0},
		{X: -1.0, Y: 0.0, Z: 0.0, Vx: 0.0, Vy: -0.012, Vz: 0.0, T: t0},
	}
	sim, err := NewGravSim(Ssb, t0, 2, states)
	if err != nil {
		t.Fatal(err)
	}
	out := make([]StateVector, 2)
	t1 := TimeFromCalendar(2020, 1, 15, 0, 0, 0.0)
	if err := sim.Update(t1, 2, out); err != nil {
		t.Fatal(err)
	}
	if out[0].X == 0.0 || out[1].X == 0.0 {
		t.Error("multi-body simulation produced invalid positions")
	}
	if !isfinite(out[0].X) || !isfinite(out[1].Y) {
		t.Error("multi-body simulation produced non-finite values")
	}
	if sim.NumBodies() != 2 {
		t.Errorf("NumBodies = %d, expected 2", sim.NumBodies())
	}
	if sim.Origin() != Ssb {
		t.Errorf("Origin = %d, expected SSB", sim.Origin())
	}
	if diff := absDiff(sim.Time().Ut, t1.Ut); diff > 1.0e-12 {
		t.Error("simulation time does not match")
	}
	if _, err := sim.BodyState(Mars); err != nil {
		t.Errorf("BodyState(Mars) error: %v", err)
	}
}
