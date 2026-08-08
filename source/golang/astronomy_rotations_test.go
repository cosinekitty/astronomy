package astronomy

import (
	"testing"
)

func TestRotation_EQD_EQJ(t *testing.T) {
	// EQJ->EQD and EQD->EQJ must be inverses.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	eqjEqd := Rotation_EQJ_EQD(time)
	eqdEqj := Rotation_EQD_EQJ(time)
	roundTrip := CombineRotation(eqjEqd, eqdEqj)
	ident := IdentityMatrix()
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if diff := absDiff(roundTrip.Rot[i][j], ident.Rot[i][j]); diff > 1.0e-15 {
				t.Errorf("EQJ<->EQD roundtrip element [%d][%d] diff = %g", i, j, diff)
			}
		}
	}
}

func TestRotation_ECT_EQD(t *testing.T) {
	// ECT->EQD and EQD->ECT must be inverses.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	ectEqd := Rotation_ECT_EQD(time)
	eqdEct := Rotation_EQD_ECT(time)
	roundTrip := CombineRotation(ectEqd, eqdEct)
	ident := IdentityMatrix()
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if diff := absDiff(roundTrip.Rot[i][j], ident.Rot[i][j]); diff > 1.0e-15 {
				t.Errorf("ECT<->EQD roundtrip element [%d][%d] diff = %g", i, j, diff)
			}
		}
	}
}

func TestRotation_HOR_EQD(t *testing.T) {
	observer := Observer{Latitude: 29.0, Longitude: -81.0, Height: 10.0}
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	horEqd := Rotation_HOR_EQD(time, observer)
	eqdHor := Rotation_EQD_HOR(time, observer)
	roundTrip := CombineRotation(horEqd, eqdHor)
	ident := IdentityMatrix()
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if diff := absDiff(roundTrip.Rot[i][j], ident.Rot[i][j]); diff > 1.0e-15 {
				t.Errorf("HOR<->EQD roundtrip element [%d][%d] diff = %g", i, j, diff)
			}
		}
	}
}

func TestSearchLunarApsis(t *testing.T) {
	// Reference values from the Python implementation.
	time := TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
	apsis, err := SearchLunarApsis(time)
	if err != nil {
		t.Fatal(err)
	}
	if apsis.Kind != ApsisKindPericenter {
		t.Errorf("expected perigee, got kind %d", apsis.Kind)
	}
	// Reference time: 2019-07-05T04:58:49.287Z (UT).
	expected := TimeFromCalendar(2019, 7, 5, 4, 58, 49.287)
	if diff := absDiff(apsis.Time.Ut, expected.Ut); diff > 0.01 {
		t.Errorf("perigee time error = %g days", diff)
	}
	if diff := absDiff(apsis.DistAu, 0.002431300424); diff > 1.0e-12 {
		t.Errorf("perigee dist_au error = %g", diff)
	}
	if diff := absDiff(apsis.DistKm, 363717.366435353761); diff > 1.0e-6 {
		t.Errorf("perigee dist_km error = %g", diff)
	}

	next, err := NextLunarApsis(apsis)
	if err != nil {
		t.Fatal(err)
	}
	if next.Kind != ApsisKindApocenter {
		t.Errorf("expected apogee after perigee, got kind %d", next.Kind)
	}
	// Reference time: 2019-07-20T23:58:55.163Z.
	expectedNext := TimeFromCalendar(2019, 7, 20, 23, 58, 55.163)
	if diff := absDiff(next.Time.Ut, expectedNext.Ut); diff > 0.01 {
		t.Errorf("apogee time error = %g days", diff)
	}
	if diff := absDiff(next.DistAu, 0.002710395712); diff > 1.0e-12 {
		t.Errorf("apogee dist_au error = %g", diff)
	}
}
