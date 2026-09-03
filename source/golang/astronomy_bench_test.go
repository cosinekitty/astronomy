package astronomy

import "testing"

func benchTime() AstroTime {
	return TimeFromCalendar(2019, 6, 24, 15, 45, 37.0)
}

func BenchmarkGeoMoon(b *testing.B) {
	time := benchTime()
	for i := 0; i < b.N; i++ {
		GeoMoon(time)
	}
}

func BenchmarkGeoMoonState(b *testing.B) {
	time := benchTime()
	for i := 0; i < b.N; i++ {
		GeoMoonState(time)
	}
}

func BenchmarkHelioVector(b *testing.B) {
	time := benchTime()
	for i := 0; i < b.N; i++ {
		HelioVector(Mars, time)
	}
}

func BenchmarkHelioVectorPluto(b *testing.B) {
	time := benchTime()
	for i := 0; i < b.N; i++ {
		HelioVector(Pluto, time)
	}
}

func BenchmarkGeoVector(b *testing.B) {
	time := benchTime()
	for i := 0; i < b.N; i++ {
		GeoVector(Mars, time, WithAberration)
	}
}

func BenchmarkSearchMoonPhase(b *testing.B) {
	time := benchTime()
	for i := 0; i < b.N; i++ {
		SearchMoonPhase(180.0, time, 40.0)
	}
}

func BenchmarkSearchRiseSet(b *testing.B) {
	time := benchTime()
	observer := Observer{Latitude: 45.0, Longitude: -90.0}
	for i := 0; i < b.N; i++ {
		SearchRiseSet(Sun, observer, +1, time, 1.0, 0.0)
	}
}

func BenchmarkSeasons(b *testing.B) {
	for i := 0; i < b.N; i++ {
		Seasons(2023)
	}
}

func BenchmarkJupiterMoons(b *testing.B) {
	time := benchTime()
	for i := 0; i < b.N; i++ {
		JupiterMoons(time)
	}
}
