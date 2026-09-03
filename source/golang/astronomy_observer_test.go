package astronomy

import (
	"testing"
)

func TestCalcPluto(t *testing.T) {
	// Reference from the C implementation: HelioVector(Pluto) and HelioState.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	b := calcPluto(time, true)
	if diff := absDiff(b.r.X, 1.2407758350790576e+01); diff > 1.0e-12 {
		t.Errorf("Pluto x error = %g", diff)
	}
	if diff := absDiff(b.r.Y, -2.8782898678669827e+01); diff > 1.0e-12 {
		t.Errorf("Pluto y error = %g", diff)
	}
	if diff := absDiff(b.r.Z, -1.2720848478603108e+01); diff > 1.0e-12 {
		t.Errorf("Pluto z error = %g", diff)
	}
	if diff := absDiff(b.v.X, 2.9988119448312298e-03); diff > 1.0e-15 {
		t.Errorf("Pluto vx error = %g", diff)
	}
	if diff := absDiff(b.v.Y, 8.1065069906927458e-04); diff > 1.0e-15 {
		t.Errorf("Pluto vy error = %g", diff)
	}
	if diff := absDiff(b.v.Z, -6.4828321027241942e-04); diff > 1.0e-15 {
		t.Errorf("Pluto vz error = %g", diff)
	}
}

func TestObserverVector(t *testing.T) {
	// Reference from the C implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	observer := Observer{Latitude: 29.0, Longitude: -81.0, Height: 10.0}
	vec, err := ObserverVector(&time, observer, EquatorJ2000)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(vec.X, 1.4220220212736747e-05); diff > 1.0e-15 {
		t.Errorf("observer x error = %g", diff)
	}
	if diff := absDiff(vec.Y, 3.4518863661798197e-05); diff > 1.0e-15 {
		t.Errorf("observer y error = %g", diff)
	}
	if diff := absDiff(vec.Z, 2.0522054090147405e-05); diff > 1.0e-15 {
		t.Errorf("observer z error = %g", diff)
	}
}

func TestObserverState(t *testing.T) {
	// Reference from the C implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	observer := Observer{Latitude: 29.0, Longitude: -81.0, Height: 10.0}
	state, err := ObserverState(&time, observer, EquatorJ2000)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(state.Vx, -2.1748442314656669e-04); diff > 1.0e-15 {
		t.Errorf("observer vx error = %g", diff)
	}
	if diff := absDiff(state.Vy, 8.9352122805111330e-05); diff > 1.0e-15 {
		t.Errorf("observer vy error = %g", diff)
	}
	if diff := absDiff(state.Vz, 4.0652095280833167e-07); diff > 1.0e-15 {
		t.Errorf("observer vz error = %g", diff)
	}
}

func TestVectorObserverRoundTrip(t *testing.T) {
	// ObserverVector followed by VectorObserver must recover the location.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	observer := Observer{Latitude: 29.0, Longitude: -81.0, Height: 10.0}
	vec, err := ObserverVector(&time, observer, EquatorJ2000)
	if err != nil {
		t.Fatal(err)
	}
	obs2 := VectorObserver(&vec, EquatorJ2000)
	if diff := absDiff(obs2.Latitude, observer.Latitude); diff > 1.0e-6 {
		t.Errorf("latitude round-trip error = %g", diff)
	}
	if diff := absDiff(obs2.Longitude, observer.Longitude); diff > 1.0e-6 {
		t.Errorf("longitude round-trip error = %g", diff)
	}
}

func TestSearchMoonNode(t *testing.T) {
	// Reference from the C implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	node, err := SearchMoonNode(time)
	if err != nil {
		t.Fatal(err)
	}
	if node.Kind != AscendingNode {
		t.Errorf("expected ascending node, got %d", node.Kind)
	}
	want := TimeFromCalendar(2019, 7, 3, 6, 53, 17.241)
	if diff := absDiff(node.Time.Ut, want.Ut); diff > 0.01 {
		t.Errorf("moon node time error = %g days", diff)
	}

	next, err := NextMoonNode(node)
	if err != nil {
		t.Fatal(err)
	}
	if next.Kind != DescendingNode {
		t.Errorf("expected descending node, got %d", next.Kind)
	}
	wantNext := TimeFromCalendar(2019, 7, 16, 9, 5, 36.525)
	if diff := absDiff(next.Time.Ut, wantNext.Ut); diff > 0.01 {
		t.Errorf("next moon node time error = %g days", diff)
	}
}

func TestRotation_ObserverFrames(t *testing.T) {
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	observer := Observer{Latitude: 29.0, Longitude: -81.0, Height: 10.0}
	ident := IdentityMatrix()

	pairs := []struct {
		a, b RotationMatrix
	}{
		{Rotation_EQJ_HOR(time, observer), Rotation_HOR_EQJ(time, observer)},
		{Rotation_ECL_EQD(time), Rotation_EQD_ECL(time)},
		{Rotation_ECL_HOR(time, observer), InverseRotation(Rotation_ECL_HOR(time, observer))},
	}
	for _, p := range pairs {
		rt := CombineRotation(p.a, p.b)
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				if diff := absDiff(rt.Rot[i][j], ident.Rot[i][j]); diff > 1.0e-14 {
					t.Errorf("rotation roundtrip element [%d][%d] diff = %g", i, j, diff)
				}
			}
		}
	}
}
