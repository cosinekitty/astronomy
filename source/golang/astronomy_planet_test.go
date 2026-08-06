package astronomy

import "testing"

func TestCalcVsopPosVel(t *testing.T) {
	// Reference values from the Python implementation (HelioState), which uses
	// identical VSOP87 data and math.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	cases := []struct {
		body                Body
		x, y, z, vx, vy, vz float64
	}{
		{Mars, -1.0302875769127506e+00, 1.1571636339310609e+00, 5.5857129023057051e-01, -1.0386931758615763e-02, -6.9799137521443918e-03, -2.9212564623814236e-03},
		{Venus, 3.7081766155378248e-01, 5.7321176836530141e-01, 2.3444589155700651e-01, -1.7416504345155941e-02, 8.9967248164195546e-03, 5.1505550411433711e-03},
		{Jupiter, -8.8896878776369170e-01, -4.8053424284456838e+00, -2.0380229598703803e+00, 7.3553365993393929e-03, -7.7287444431664522e-04, -5.1020083489960772e-04},
	}
	for _, c := range cases {
		state := calcVsopPosVel(&vsop[c.body], time.Tt)
		if diff := absDiff(state.r.X, c.x); diff > 1.0e-13 {
			t.Errorf("%s: x error = %g", BodyName(c.body), diff)
		}
		if diff := absDiff(state.r.Y, c.y); diff > 1.0e-13 {
			t.Errorf("%s: y error = %g", BodyName(c.body), diff)
		}
		if diff := absDiff(state.r.Z, c.z); diff > 1.0e-13 {
			t.Errorf("%s: z error = %g", BodyName(c.body), diff)
		}
		if diff := absDiff(state.v.X, c.vx); diff > 1.0e-16 {
			t.Errorf("%s: vx error = %g", BodyName(c.body), diff)
		}
		if diff := absDiff(state.v.Y, c.vy); diff > 1.0e-16 {
			t.Errorf("%s: vy error = %g", BodyName(c.body), diff)
		}
		if diff := absDiff(state.v.Z, c.vz); diff > 1.0e-16 {
			t.Errorf("%s: vz error = %g", BodyName(c.body), diff)
		}
	}
}

func TestSynodicPeriod(t *testing.T) {
	cases := []struct {
		body     Body
		expected float64
	}{
		{Mars, 779.934250724223},
		{Venus, 583.923648792288},
		{Jupiter, 398.883614706403},
	}
	for _, c := range cases {
		got, err := synodicPeriod(c.body)
		if err != nil {
			t.Fatal(err)
		}
		if diff := absDiff(got, c.expected); diff > 1.0e-9 {
			t.Errorf("synodicPeriod(%s) = %g, expected %g", BodyName(c.body), got, c.expected)
		}
	}
	moon, err := synodicPeriod(Moon)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(moon, 29.530588); diff > 1.0e-12 {
		t.Errorf("synodicPeriod(Moon) = %g", moon)
	}
	if _, err := synodicPeriod(Earth); err == nil {
		t.Error("Earth should not have a synodic period")
	}
}

func TestCorrectLightTravel(t *testing.T) {
	// A stationary body 1 light-day away must have its observed time backdated
	// by exactly one day.
	time := TimeFromCalendar(2020, 1, 1, 0, 0, 0.0)
	pos, err := correctLightTravel(func(t AstroTime) (AstroVector, error) {
		return AstroVector{X: SpeedOfLightAuPerDay, Y: 0, Z: 0, T: t}, nil
	}, time)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(pos.X, SpeedOfLightAuPerDay); diff > 1.0e-9 {
		t.Errorf("light-travel corrected position error = %g", diff)
	}
	// Too far away must fail.
	_, err = correctLightTravel(func(t AstroTime) (AstroVector, error) {
		return AstroVector{X: 2 * SpeedOfLightAuPerDay, Y: 0, Z: 0, T: t}, nil
	}, time)
	if err == nil {
		t.Error("expected error for distance > 1 light-day")
	}
}
