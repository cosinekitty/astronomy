package astronomy

import (
	"bufio"
	"math"
	"os"
	"strconv"
	"strings"
	"testing"
)

func atof(s string) float64 {
	x, err := strconv.ParseFloat(s, 64)
	if err != nil {
		panic("Cannot parse float: " + s)
	}
	return x
}

func TestBodyFunctions(t *testing.T) {
	// BodyName / BodyCode round trips.
	names := []string{"Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto", "Sun", "Moon", "EMB", "SSB"}
	for _, name := range names {
		body := BodyCode(name)
		if BodyName(body) != name {
			t.Errorf("BodyName(BodyCode(%s)) = %s", name, BodyName(body))
		}
	}
	if BodyCode("not a body") != InvalidBody {
		t.Error("expected InvalidBody for unknown name")
	}
	if BodyName(InvalidBody) != "" {
		t.Error("expected empty name for invalid body")
	}

	// Orbital periods (from ctest.c / astronomy.c).
	periods := map[Body]float64{
		Mercury: 87.969, Venus: 224.701, Earth: 365.256, Mars: 686.980,
		Jupiter: 4332.589, Saturn: 10759.22, Uranus: 30685.4, Neptune: 60189.0, Pluto: 90560.0,
	}
	for body, expected := range periods {
		if diff := absDiff(PlanetOrbitalPeriod(body), expected); diff > 1.0e-12 {
			t.Errorf("PlanetOrbitalPeriod(%s) = %g, expected %g", BodyName(body), PlanetOrbitalPeriod(body), expected)
		}
	}
	if PlanetOrbitalPeriod(InvalidBody) != 0.0 {
		t.Error("expected 0 for invalid body orbital period")
	}

	// Mass products (from astronomy.c).
	masses := map[Body]float64{
		Sun:     0.2959122082855911e-03,
		Mercury: 0.4912547451450812e-10,
		Venus:   0.7243452486162703e-09,
		Earth:   0.8887692390113509e-09,
		Mars:    0.9549535105779258e-10,
		Jupiter: 0.2825345909524226e-06,
		Saturn:  0.8459715185680659e-07,
		Uranus:  0.1292024916781969e-07,
		Neptune: 0.1524358900784276e-07,
		Pluto:   0.2188699765425970e-11,
	}
	for body, expected := range masses {
		if diff := absDiff(MassProduct(body), expected); diff > 1.0e-22 {
			t.Errorf("MassProduct(%s) = %g, expected %g (diff %g)", BodyName(body), MassProduct(body), expected, diff)
		}
	}
	if diff := absDiff(MassProduct(Moon), masses[Earth]/81.30056); diff > 1.0e-22 {
		t.Errorf("Moon GM mismatch: %g", diff)
	}
	if diff := absDiff(MassProduct(Emb), masses[Earth]+masses[Earth]/81.30056); diff > 1.0e-22 {
		t.Errorf("EMB GM mismatch: %g", diff)
	}

	if isSuperiorPlanet(Mars) != true || isSuperiorPlanet(Earth) != false {
		t.Error("isSuperiorPlanet logic error")
	}
}

func TestDefineStar(t *testing.T) {
	// Star must be invalid before definition.
	if userDefinedStar(Star1) != nil {
		t.Error("Star1 should be undefined before DefineStar")
	}
	if err := DefineStar(Star1, 3.5, 45.0, 100.0); err != nil {
		t.Fatal(err)
	}
	star := userDefinedStar(Star1)
	if star == nil {
		t.Fatal("Star1 should be defined now")
	}
	if diff := absDiff(star.ra, 3.5); diff > 1.0e-15 {
		t.Errorf("star.ra = %g", star.ra)
	}
	if diff := absDiff(star.dec, 45.0); diff > 1.0e-15 {
		t.Errorf("star.dec = %g", star.dec)
	}
	if diff := absDiff(star.dist, 100.0*AuPerLightYear); diff > 1.0e-9 {
		t.Errorf("star.dist = %g", star.dist)
	}

	// Invalid arguments must be rejected without changing the definition.
	if err := DefineStar(InvalidBody, 1.0, 1.0, 1.0); err == nil {
		t.Error("expected error for invalid body")
	}
	if err := DefineStar(Star2, -1.0, 0.0, 10.0); err == nil {
		t.Error("expected error for negative RA")
	}
	if err := DefineStar(Star2, 24.0, 0.0, 10.0); err == nil {
		t.Error("expected error for RA >= 24")
	}
	if err := DefineStar(Star2, 0.0, 91.0, 10.0); err == nil {
		t.Error("expected error for Dec > 90")
	}
	if err := DefineStar(Star2, 0.0, 0.0, 0.5); err == nil {
		t.Error("expected error for distance < 1 ly")
	}
	if userDefinedStar(Star2) != nil {
		t.Error("Star2 must remain undefined after rejected definitions")
	}
}

func TestFrameRotations(t *testing.T) {
	// EQJ <-> ECL must be inverses of each other.
	eqjEcl := Rotation_EQJ_ECL()
	eclEqj := Rotation_ECL_EQJ()
	roundTrip := CombineRotation(eqjEcl, eclEqj)
	ident := IdentityMatrix()
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if diff := absDiff(roundTrip.Rot[i][j], ident.Rot[i][j]); diff > 1.0e-15 {
				t.Errorf("EQJ->ECL->EQJ roundtrip element [%d][%d] diff = %g", i, j, diff)
			}
		}
	}

	// EQJ <-> GAL must be inverses of each other.
	eqjGal := Rotation_EQJ_GAL()
	galEqj := Rotation_GAL_EQJ()
	roundTrip = CombineRotation(eqjGal, galEqj)
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if diff := absDiff(roundTrip.Rot[i][j], ident.Rot[i][j]); diff > 1.0e-15 {
				t.Errorf("EQJ->GAL->EQJ roundtrip element [%d][%d] diff = %g", i, j, diff)
			}
		}
	}

	// Rotating the equatorial north pole (0,0,1) from EQJ to ECL should yield
	// a point on the ecliptic with latitude 90 - obliquity (the ecliptic pole
	// is tilted from the equatorial pole by the obliquity).
	time := TimeFromCalendar(2000, 1, 1, 0, 0, 0.0)
	v := RotateVector(eqjEcl, AstroVector{X: 0.0, Y: 0.0, Z: 1.0, T: time})
	sphere, err := SphereFromVector(v)
	if err != nil {
		t.Fatal(err)
	}
	obliquity := DegreesFromRadians(math.Atan2(sinObliq2000, cosObliq2000))
	expectedLat := 90.0 - obliquity
	if diff := absDiff(sphere.Lat, expectedLat); diff > 1.0e-9 {
		t.Errorf("ecliptic pole latitude = %g, expected %g", sphere.Lat, expectedLat)
	}
}

func TestEqjGalJpl(t *testing.T) {
	// Compare EQJ -> GAL rotation against JPL Horizons reference data,
	// exactly as ctest.c Test_EQJ_GAL_JPL does.
	const filename = "../../generate/galactic/mars.txt"
	const thresholdArcsec = 23.0
	file, err := os.Open(filename)
	if err != nil {
		t.Fatalf("Cannot open file: %s", filename)
	}
	defer file.Close()
	rot := Rotation_EQJ_GAL()
	time := TimeFromCalendar(2000, 1, 1, 0, 0, 0.0)
	maxDiff := 0.0
	ncases := 0
	foundBegin := false
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		if !foundBegin {
			if len(line) >= 5 && line[:5] == "$$SOE" {
				foundBegin = true
			}
			continue
		}
		if len(line) >= 5 && line[:5] == "$$EOE" {
			break
		}
		if len(line) < 18 {
			continue
		}
		tokens := strings.Fields(line[18:])
		if len(tokens) < 4 {
			continue
		}
		ra := atof(tokens[0])
		dec := atof(tokens[1])
		glon := atof(tokens[2])
		glat := atof(tokens[3])

		eqjVec := VectorFromSphere(Spherical{Lat: dec, Lon: ra, Dist: 1.0}, time)
		galVec := RotateVector(rot, eqjVec)
		galSphere, err := SphereFromVector(galVec)
		if err != nil {
			t.Fatal(err)
		}
		dlat := galSphere.Lat - glat
		dlon := math.Cos(RadiansFromDegrees(glat)) * (galSphere.Lon - glon)
		diff := 3600.0 * math.Sqrt(dlon*dlon+dlat*dlat)
		if diff > thresholdArcsec {
			t.Errorf("EQJ->GAL excessive error = %0.3f arcseconds", diff)
		}
		if diff > maxDiff {
			maxDiff = diff
		}
		ncases++
	}
	t.Logf("EQJ->GAL tested %d cases, max error = %0.3f arcseconds", ncases, maxDiff)
	if ncases == 0 {
		t.Error("no test cases found")
	}
}
