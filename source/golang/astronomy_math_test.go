package astronomy

import (
	"math"
	"testing"
)

func absDiff(x, y float64) float64 {
	return math.Abs(x - y)
}

func compareVectors(label string, v AstroVector, ex, ey, ez float64, tol float64, t *testing.T) {
	dx := ex - v.X
	dy := ey - v.Y
	dz := ez - v.Z
	diff := math.Sqrt(dx*dx + dy*dy + dz*dz)
	if diff > tol {
		t.Errorf("%s: excessive error = %g; got (%g, %g, %g)", label, diff, v.X, v.Y, v.Z)
	}
}

func TestInverseRotation(t *testing.T) {
	a := RotationMatrix{}
	a.Rot[0] = [3]float64{1.0, 2.0, 3.0}
	a.Rot[1] = [3]float64{4.0, 5.0, 6.0}
	a.Rot[2] = [3]float64{7.0, 8.0, 9.0}
	b := InverseRotation(a)
	// b must be the transpose of a.
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if diff := absDiff(b.Rot[i][j], a.Rot[j][i]); diff != 0.0 {
				t.Errorf("Rotation_MatrixInverse: element [%d][%d] diff = %g", i, j, diff)
			}
		}
	}
}

func TestCombineRotation(t *testing.T) {
	// Same test data as ctest.c Rotation_MatrixMultiply, converted to Go's
	// column layout Rot[c][r] (column c, row r).
	a := RotationMatrix{}
	a.Rot[0] = [3]float64{1.0, 4.0, 7.0}
	a.Rot[1] = [3]float64{2.0, 5.0, 8.0}
	a.Rot[2] = [3]float64{3.0, 6.0, 9.0}

	b := RotationMatrix{}
	b.Rot[0] = [3]float64{10.0, 13.0, 16.0}
	b.Rot[1] = [3]float64{11.0, 14.0, 17.0}
	b.Rot[2] = [3]float64{12.0, 15.0, 18.0}

	v := RotationMatrix{}
	v.Rot[0] = [3]float64{84.0, 201.0, 318.0}
	v.Rot[1] = [3]float64{90.0, 216.0, 342.0}
	v.Rot[2] = [3]float64{96.0, 231.0, 366.0}

	// CombineRotation(b, a) with a=1..9 and b=10..18 must give v.
	c := CombineRotation(b, a)
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if diff := absDiff(c.Rot[i][j], v.Rot[i][j]); diff > 1.0e-15 {
				t.Errorf("Rotation_MatrixMultiply: element [%d][%d] diff = %g", i, j, diff)
			}
		}
	}
}

func TestPivot(t *testing.T) {
	// Test #1: pivot identity 90 degrees counterclockwise around the z-axis.
	r, err := Pivot(IdentityMatrix(), 2, +90.0)
	if err != nil {
		t.Fatal(err)
	}
	a := RotationMatrix{}
	a.Rot[0] = [3]float64{0.0, 1.0, 0.0}
	a.Rot[1] = [3]float64{-1.0, 0.0, 0.0}
	a.Rot[2] = [3]float64{0.0, 0.0, 1.0}
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if diff := absDiff(r.Rot[i][j], a.Rot[i][j]); diff > 1.0e-15 {
				t.Errorf("Rotation_Pivot #1: element [%d][%d] diff = %g", i, j, diff)
			}
		}
	}

	// Test #2: pivot -30 deg around x-axis, then +180 deg around y-axis; rotate a vector.
	r, err = Pivot(r, 0, -30.0)
	if err != nil {
		t.Fatal(err)
	}
	r, err = Pivot(r, 1, +180.0)
	if err != nil {
		t.Fatal(err)
	}
	time := TimeFromCalendar(2000, 1, 1, 0, 0, 0.0)
	v1 := AstroVector{X: 1.0, Y: 2.0, Z: 3.0, T: time}
	v2 := RotateVector(r, v1)
	compareVectors("Rotation_Pivot #2", v2, 2.0, 2.3660254037844390, -2.0980762113533156, 1.0e-15, t)
}

func rotateX(angle float64) RotationMatrix {
	ca := math.Cos(RadiansFromDegrees(angle))
	sa := math.Sin(RadiansFromDegrees(angle))
	var m RotationMatrix
	m.Rot[0] = [3]float64{1.0, 0.0, 0.0}
	m.Rot[1] = [3]float64{0.0, ca, sa}
	m.Rot[2] = [3]float64{0.0, -sa, ca}
	return m
}

func rotateY(angle float64) RotationMatrix {
	ca := math.Cos(RadiansFromDegrees(angle))
	sa := math.Sin(RadiansFromDegrees(angle))
	var m RotationMatrix
	m.Rot[0] = [3]float64{ca, 0.0, -sa}
	m.Rot[1] = [3]float64{0.0, 1.0, 0.0}
	m.Rot[2] = [3]float64{sa, 0.0, ca}
	return m
}

func rotateZ(angle float64) RotationMatrix {
	ca := math.Cos(RadiansFromDegrees(angle))
	sa := math.Sin(RadiansFromDegrees(angle))
	var m RotationMatrix
	m.Rot[0] = [3]float64{ca, sa, 0.0}
	m.Rot[1] = [3]float64{-sa, ca, 0.0}
	m.Rot[2] = [3]float64{0.0, 0.0, 1.0}
	return m
}

func TestSpin(t *testing.T) {
	time := TimeFromCalendar(2019, 5, 5, 12, 30, 0.0)
	cases := []struct {
		xrot, yrot, zrot float64
		sx, sy, sz       float64
		tx, ty, tz       float64
	}{
		{0.0, 0.0, 90.0, 1, 2, 3, -2, 1, 3},
		{0.0, 0.0, 0.0, 1.0, 2.0, -3.0, 1.0, 2.0, -3.0},
		{90.0, 0.0, 0.0, 1, 2, 3, 1, -3, 2},
		{0.0, 90.0, 0.0, 1, 2, 3, 3, 2, -1},
		{180.0, 180.0, 180.0, 1, 2, 3, 1, 2, 3},
		{0.0, 0.0, -45.0, 1, 0, 0, 0.7071067811865476, -0.7071067811865476, 0},
	}
	for _, c := range cases {
		m := CombineRotation(rotateX(c.xrot), rotateY(c.yrot))
		m = CombineRotation(m, rotateZ(c.zrot))
		sv := AstroVector{X: c.sx, Y: c.sy, Z: c.sz, T: time}
		tv := RotateVector(m, sv)
		compareVectors("TestSpin", tv, c.tx, c.ty, c.tz, 1.0e-15, t)
	}
}

func TestVectorFromAngles(t *testing.T) {
	time := TimeFromCalendar(2015, 3, 5, 12, 0, 0.0)
	cases := []struct {
		lat, lon float64
		x, y, z  float64
	}{
		{0.0, 0.0, 1.0, 0.0, 0.0},
		{0.0, 90.0, 0.0, 1.0, 0.0},
		{0.0, 180.0, -1.0, 0.0, 0.0},
		{0.0, 270.0, 0.0, -1.0, 0.0},
		{+90.0, 0.0, 0.0, 0.0, 1.0},
		{-90.0, 0.0, 0.0, 0.0, -1.0},
		{-30.0, +60.0, 0.43301270189221946, 0.75, -0.5},
	}
	for _, c := range cases {
		v := VectorFromSphere(Spherical{Lat: c.lat, Lon: c.lon, Dist: 1.0}, time)
		compareVectors("TestVectorFromAngles", v, c.x, c.y, c.z, 2.0e-16, t)
	}
}

func TestAnglesFromVector(t *testing.T) {
	time := TimeFromCalendar(2015, 3, 5, 12, 0, 0.0)
	cases := []struct {
		lat, lon float64
		x, y, z  float64
	}{
		{0.0, 0.0, 1.0, 0.0, 0.0},
		{0.0, 90.0, 0.0, 1.0, 0.0},
		{0.0, 180.0, -1.0, 0.0, 0.0},
		{0.0, 270.0, 0.0, -1.0, 0.0},
		{+90.0, 0.0, 0.0, 0.0, 1.0},
		{-90.0, 0.0, 0.0, 0.0, -1.0},
		{-30.0, +60.0, 0.43301270189221946, 0.75, -0.5},
	}
	for _, c := range cases {
		sphere, err := SphereFromVector(AstroVector{X: c.x, Y: c.y, Z: c.z, T: time})
		if err != nil {
			t.Fatal(err)
		}
		if latdiff := absDiff(sphere.Lat, c.lat); latdiff > 8.0e-15 {
			t.Errorf("TestAnglesFromVector: lat diff = %g for (%g, %g, %g)", latdiff, c.x, c.y, c.z)
		}
		if londiff := absDiff(sphere.Lon, c.lon); londiff > 8.0e-15 {
			t.Errorf("TestAnglesFromVector: lon diff = %g for (%g, %g, %g)", londiff, c.x, c.y, c.z)
		}
	}
}

func TestEquatorFromVector(t *testing.T) {
	time := TimeFromCalendar(2015, 3, 5, 12, 0, 0.0)
	// Vector on the celestial equator: RA should be 3 sidereal hours (45 degrees / 15).
	v := AstroVector{X: math.Sqrt(2.0) / 2.0, Y: math.Sqrt(2.0) / 2.0, Z: 0.0, T: time}
	equ, err := EquatorFromVector(v)
	if err != nil {
		t.Fatal(err)
	}
	if diff := absDiff(equ.Ra, 3.0); diff > 1.0e-14 {
		t.Errorf("RA error = %g", diff)
	}
	if diff := absDiff(equ.Dec, 0.0); diff > 1.0e-14 {
		t.Errorf("Dec error = %g", diff)
	}
	if diff := absDiff(equ.Dist, 1.0); diff > 1.0e-15 {
		t.Errorf("Dist error = %g", diff)
	}
}

func TestRotateState(t *testing.T) {
	time := TimeFromCalendar(2000, 6, 1, 12, 0, 0.0)
	// Rotate a state vector 90 degrees around the z-axis.
	r := rotateZ(90.0)
	state := StateVector{X: 1.0, Y: 2.0, Z: 3.0, Vx: 4.0, Vy: 5.0, Vz: 6.0, T: time}
	rot := RotateState(r, state)
	// 90-degree z-rotation: (x,y) -> (-y,x).
	if diff := absDiff(rot.X, -2.0); diff > 1.0e-15 {
		t.Errorf("X error = %g", diff)
	}
	if diff := absDiff(rot.Y, 1.0); diff > 1.0e-15 {
		t.Errorf("Y error = %g", diff)
	}
	if diff := absDiff(rot.Z, 3.0); diff > 1.0e-15 {
		t.Errorf("Z error = %g", diff)
	}
	if diff := absDiff(rot.Vx, -5.0); diff > 1.0e-15 {
		t.Errorf("Vx error = %g", diff)
	}
	if diff := absDiff(rot.Vy, 4.0); diff > 1.0e-15 {
		t.Errorf("Vy error = %g", diff)
	}
	if diff := absDiff(rot.Vz, 6.0); diff > 1.0e-15 {
		t.Errorf("Vz error = %g", diff)
	}
}

func TestQuadInterp(t *testing.T) {
	// Parabola f(x) = (x-10)^2 - 4, sampled at x = -1, 0, +1 with tm=0, dt=1.
	// fa=117, fm=96, fb=77. Both roots (x=8 and x=12) are outside [-1, +1].
	if _, _, ok := quadInterp(0.0, 1.0, 117.0, 96.0, 77.0); ok {
		t.Error("expected no unique in-bounds root")
	}

	// Linear function f(x)=x-0.25 sampled at x=0,1,2 -> fa=-0.25, fm=0.75, fb=1.75.
	// Root at x=0.25 (normalized x=-0.75 relative to tm=1).
	troot, dfdt, ok := quadInterp(1.0, 1.0, -0.25, 0.75, 1.75)
	if !ok {
		t.Error("expected a valid root")
	}
	if diff := absDiff(troot, 0.25); diff > 1.0e-14 {
		t.Errorf("quadInterp root error = %g", diff)
	}
	if diff := absDiff(dfdt, 1.0); diff > 1.0e-14 {
		t.Errorf("quadInterp derivative error = %g", diff)
	}
}

func TestLongitudeHelpers(t *testing.T) {
	if got := normalizeLongitude(350.0); got != 350.0 {
		t.Errorf("normalizeLongitude(350) = %g", got)
	}
	if got := normalizeLongitude(370.0); got != 10.0 {
		t.Errorf("normalizeLongitude(370) = %g", got)
	}
	if got := normalizeLongitude(-10.0); got != 350.0 {
		t.Errorf("normalizeLongitude(-10) = %g", got)
	}
	if got := normalizeLongitude(-720.0); got != 0.0 {
		t.Errorf("normalizeLongitude(-720) = %g", got)
	}
	if got := longitudeOffset(190.0); got != -170.0 {
		t.Errorf("longitudeOffset(190) = %g", got)
	}
	if got := longitudeOffset(-190.0); got != 170.0 {
		t.Errorf("longitudeOffset(-190) = %g", got)
	}
	if got := longitudeOffset(45.0); got != 45.0 {
		t.Errorf("longitudeOffset(45) = %g", got)
	}
}

func TestPivotErrors(t *testing.T) {
	if _, err := Pivot(IdentityMatrix(), 3, 45.0); err == nil {
		t.Error("expected error for invalid axis")
	}
	if _, err := Pivot(IdentityMatrix(), 0, math.NaN()); err == nil {
		t.Error("expected error for NaN angle")
	}
}

func TestSphereFromVectorZero(t *testing.T) {
	time := TimeFromUniversalDays(0.0)
	if _, err := SphereFromVector(AstroVector{X: 0, Y: 0, Z: 0, T: time}); err == nil {
		t.Error("expected error for zero-length vector")
	}
}
