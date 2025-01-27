program Demo;

{$mode objfpc}{$H+}

uses
  Classes,
  SysUtils,
  Astronomy;

const
  Months: array[1..12] of string = (
  		  'January', 'February', 'March', 'April',
  		 'May', 'June', 'July', 'August',
         'September', 'October', 'November', 'December');

  // Observer site: Slovakia/Bratislava
  Latitude: Double = 48.143889;
  Longitude: Double = 17.109722;
  Elevation: Double = 124;

  // ICRS equatorial coordinates of the star HD1 (Henry Draper catalog, star #1)
  RA_ICRS: Double = 1.2961219444444443; // degrees
  DE_ICRS: Double = 67.84007388888888; // degrees

var
  at1: AstroTime;
  st1: Double;
  ptr: PAstroTime;
  mph: AstroAngleResult;
  MoonPhase: string;
  lap: AstroApsis;
  sss: AstroSeasons;
  obs: AstroObserver;
  coo: AstroEquatorial;

function ToDateString(Event: AstroTime): string;
var
  t: AstroUTC;
begin
  t := Astronomy_UtcFromTime(Event);
  ToDateString := '' + IntToStr(t.Year) + '. ' + Months[t.Month] + '. ' + IntToStr(t.Day) + '. ';
end;

BEGIN
  at1 := Astronomy_MakeTime(2025, 1, 20, 22, 20, 0.5);
  WriteLn('Astro time for 2025-Jan-20 at 22:20:30 -> ', at1.UT:20:6, #10#13);

  // Define an observer site
  obs.Latitude := Latitude;
  obs.Longitude := Longitude;
  obs.Height := Elevation;

  at1 := Astronomy_CurrentTime();
  WriteLn('Current astro time: ', at1.UT:20:6, #10#13);

  WriteLn('Body name: ', Astronomy_BodyName(AstroBody.BODY_SUN), #10#13);

  WriteLn('Body code for "Earth": [', Astronomy_BodyCode(PChar('Earth')), ']', #10#13);

  // Calculate Greenwich sidereal time using current time
  ptr := New(PAstroTime);
  ptr^.UT := at1.UT;
  ptr^.TT := at1.TT;
  ptr^.ST := 0 / 0; // gives NaN, must be set, otherwise GAST calculation doesn't work
  st1 := Astronomy_SiderealTime(ptr);
  WriteLn('Greenwich apparent sidereal time: ', st1:20:6, ' hours', #10#13);

  // Get the current moon phase
  mph := Astronomy_MoonPhase(at1);
  if mph.Status <> AstroStatus.ASTRO_SUCCESS then
  	 WriteLn('An error occured during the Moon phase calculation!')
  else begin
  	 if (mph.Angle = 0) or (mph.Angle = 360)    then MoonPhase := 'New Moon';
     if (mph.Angle > 0) and (mph.Angle < 90)    then MoonPhase := 'Waxing Crescent';
     if mph.Angle = 90                          then MoonPhase := 'First quarter';
     if (mph.Angle > 90) and (mph.Angle < 180)  then MoonPhase := 'Waxing Gibbous';
     if mph.Angle = 180                         then MoonPhase := 'Full Moon';
     if (mph.Angle > 180) and (mph.Angle < 270) then MoonPhase := 'Waning Gibbous';
     if mph.Angle = 270                         then MoonPhase := 'Last quarter';
     if (mph.Angle > 270) and (mph.Angle < 360) then MoonPhase := 'Waning Crescent';
     WriteLn('Current Moon phase: ' + MoonPhase + #10#13);
  end;

  // Search for the first lunar apsis event since now
  lap := Astronomy_SearchLunarApsis(at1);
  if lap.Status <> AstroStatus.ASTRO_SUCCESS then
     WriteLn('An error occured during the Moon apsis calculation!')
  else
     WriteLn('Next Lunar ', lap.Kind, ' occurs at ', ToDateString(lap.Time), ' at the distance of ', lap.DistKM:12:3, ' km from Earth.');

  // Search for next lunar apsis event
  lap := Astronomy_NextLunarApsis(lap);
  if lap.Status <> AstroStatus.ASTRO_SUCCESS then
     WriteLn('An error occured during the Moon apsis calculation!')
  else
     WriteLn('Next Lunar ', lap.Kind, ' occurs at ', ToDateString(lap.Time), ' at the distance of ', lap.DistKM:12:3, ' km from Earth.', #10#13);

  // Search for seasons in 2025
  sss := Astronomy_Seasons(2025);
  if sss.Status <> AstroStatus.ASTRO_SUCCESS then
     WriteLn('An error occured during the 2025 seasons calculation!')
  else begin
     WriteLn('2025 march equinox occurs at: ', ToDateString(sss.MarEquinox));
     WriteLn('2025 june solstice occurs at: ', ToDateString(sss.JunSolstice));
     WriteLn('2025 sept equinox occurs at: ', ToDateString(sss.SepEquinox));
     WriteLn('2025 dec solstice occurs at: ', ToDateString(sss.DecSolstice), #10#13);
  end;

  // Find equatorial coordinates of date for the star HD1 as seen from Slovakia/Bratislava
  if Astronomy_DefineStar(AstroBody.BODY_STAR1, RA_ICRS / 15.0, DE_ICRS, 1000.0) <> AstroStatus.ASTRO_SUCCESS then
     WriteLn('An error occured during defining a custom star!')
  else begin
     coo := Astronomy_Equator(AstroBody.BODY_STAR1, ptr, obs, AstroEquatorDate.EQUATOR_OF_DATE, AstroAberration.NO_ABERRATION);
     if coo.Status <> AstroStatus.ASTRO_SUCCESS then
        WriteLn('An error occured during the coordinate conversion!')
     else begin
        WriteLn('Equatorial coordinates of the HD1 star at ', ToDateString(at1));
        WriteLn('Right ascension: ', coo.RA:12:3, ' hours');
        WriteLn('Declination:     ', coo.Dec:12:3, ' degrees');
     end;
  end;
  
  // Release the dynamically allocated object used in the previous calculations
  Dispose(ptr);
END.

