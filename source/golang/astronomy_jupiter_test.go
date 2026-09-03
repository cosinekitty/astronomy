package astronomy

import "testing"

func TestJupiterMoons(t *testing.T) {
	// Reference values from the Python implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	jm := JupiterMoons(time)
	cases := []struct {
		name                string
		state               StateVector
		x, y, z, vx, vy, vz float64
	}{
		{"Io", jm.Io, -9.7373664896067871e-04, 2.3908777996268209e-03, 1.1263197947485458e-03, -9.3896892165541354e-03, -3.1033537170850585e-03, -1.6308039869167847e-03},
		{"Europa", jm.Europa, 3.7695222468582096e-03, -2.1635951163459979e-03, -9.6055253067530745e-04, 4.2090017912289994e-03, 6.0845834293992143e-03, 3.0377092940355411e-03},
		{"Ganymede", jm.Ganymede, 6.9752734841741519e-03, -1.4329482664352322e-03, -5.8382645489647126e-04, 1.3652697756182637e-03, 5.5317593202497762e-03, 2.6712567332974543e-03},
		{"Callisto", jm.Callisto, -1.2673764002035187e-02, 2.0127360021468456e-04, -9.2260852936812231e-05, -5.7092254624222187e-05, -4.2537975113751316e-03, -2.0070565947053242e-03},
	}
	for _, c := range cases {
		if diff := absDiff(c.state.X, c.x); diff > 1.0e-16 {
			t.Errorf("%s: x error = %g", c.name, diff)
		}
		if diff := absDiff(c.state.Y, c.y); diff > 1.0e-16 {
			t.Errorf("%s: y error = %g", c.name, diff)
		}
		if diff := absDiff(c.state.Z, c.z); diff > 1.0e-16 {
			t.Errorf("%s: z error = %g", c.name, diff)
		}
		if diff := absDiff(c.state.Vx, c.vx); diff > 1.0e-18 {
			t.Errorf("%s: vx error = %g", c.name, diff)
		}
		if diff := absDiff(c.state.Vy, c.vy); diff > 1.0e-18 {
			t.Errorf("%s: vy error = %g", c.name, diff)
		}
		if diff := absDiff(c.state.Vz, c.vz); diff > 1.0e-18 {
			t.Errorf("%s: vz error = %g", c.name, diff)
		}
	}
}

func TestCalcSolarSystemBarycenter(t *testing.T) {
	// Reference values from the Python implementation (HelioVector SSB).
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	ssb := calcSolarSystemBarycenter(time)
	if diff := absDiff(ssb.X, 2.1939404379624889e-03); diff > 1.0e-16 {
		t.Errorf("ssb x error = %g", diff)
	}
	if diff := absDiff(ssb.Y, -6.9993151829519885e-03); diff > 1.0e-16 {
		t.Errorf("ssb y error = %g", diff)
	}
	if diff := absDiff(ssb.Z, -3.0124878121858311e-03); diff > 1.0e-16 {
		t.Errorf("ssb z error = %g", diff)
	}
}

func TestCalcEarthVsop(t *testing.T) {
	// Reference: Python HelioVector(Earth, 2019-06-24T15:45:37Z).
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	earth := calcEarth(time)
	if diff := absDiff(earth.X, 4.6002645875635940e-02); diff > 1.0e-13 {
		t.Errorf("earth x error = %g", diff)
	}
	if diff := absDiff(earth.Y, -9.3162259987610496e-01); diff > 1.0e-13 {
		t.Errorf("earth y error = %g", diff)
	}
	if diff := absDiff(earth.Z, -4.0385968694770696e-01); diff > 1.0e-13 {
		t.Errorf("earth z error = %g", diff)
	}
}
