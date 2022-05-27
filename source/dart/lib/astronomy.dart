library astronomy;

import 'dart:core';
import 'dart:math';

import 'package:meta/meta.dart';

/// The number of minutes in a day.
const minutesPerDay = 60.0 * 24.0;

/// The number of seconds in a day.
const secondsPerDay = 60.0 * minutesPerDay;

/// The number of milliseconds in a day.
const millisecondsPerDay = 1000.0 * secondsPerDay;

const daysPerTropicalYear = 365.24217;
const dayPerMillennium = 365250.0;

@visibleForTesting
@immutable
class TerseVector {
  final double x;
  final double y;
  final double z;

  const TerseVector({
    required this.x,
    required this.y,
    required this.z,
  });

  static final TerseVector zero = TerseVector(x: 0.0, y: 0.0, z: 0.0);

  TerseVector operator +(TerseVector other) {
    return TerseVector(x: x + other.x, y: y + other.y, z: z + other.z);
  }

  TerseVector operator -(TerseVector other) {
    return TerseVector(x: x - other.x, y: y - other.y, z: z - other.z);
  }

  TerseVector operator *(double s) {
    return TerseVector(x: s * x, y: s * y, z: s * z);
  }

  TerseVector operator /(double s) {
    return TerseVector(x: x / s, y: y / s, z: z / s);
  }

  double quadrature() {
    return x * x + y * y + z * z;
  }

  double magnitude() {
    return sqrt(quadrature());
  }

//<editor-fold desc="Data Methods">

  @override
  bool operator ==(Object other) =>
      identical(this, other) ||
      (other is TerseVector &&
          runtimeType == other.runtimeType &&
          x == other.x &&
          y == other.y &&
          z == other.z);

  @override
  int get hashCode => x.hashCode ^ y.hashCode ^ z.hashCode;

  @override
  String toString() {
    return 'TerseVector{ x: $x, y: $y, z: $z,}';
  }

  TerseVector copyWith({
    double? x,
    double? y,
    double? z,
  }) {
    return TerseVector(
      x: x ?? this.x,
      y: y ?? this.y,
      z: z ?? this.z,
    );
  }

  Map<String, dynamic> toMap() {
    return {
      'x': x,
      'y': y,
      'z': z,
    };
  }

  factory TerseVector.fromMap(Map<String, dynamic> map) {
    return TerseVector(
      x: map['x'] as double,
      y: map['y'] as double,
      z: map['z'] as double,
    );
  }

//</editor-fold>
}

/// A date and time used for astronomical calculations.
class Time implements Comparable<Time> {
  /// UT1/UTC number of days since noon on January 1, 2000.
  ///
  /// The floating point number of days of Universal Time since noon UTC January 1, 2000.
  /// Astronomy Engine approximates UTC and UT1 as being the same thing, although they are
  /// not exactly equivalent; UTC and UT1 can disagree by up to plus or minus 0.9 seconds.
  /// This approximation is sufficient for the accuracy requirements of Astronomy Engine.
  ///
  /// Universal Time Coordinate (UTC) is the international standard for legal and civil
  /// timekeeping and replaces the older Greenwich Mean Time (GMT) standard.
  /// UTC is kept in sync with unpredictable observed changes in the Earth's rotation
  /// by occasionally adding leap seconds as needed.
  ///
  /// UT1 is an idealized time scale based on observed rotation of the Earth, which
  /// gradually slows down in an unpredictable way over time, due to tidal drag by the Moon and Sun,
  /// large scale weather events like hurricanes, and internal seismic and convection effects.
  /// Conceptually, UT1 drifts from atomic time continuously and erratically, whereas UTC
  /// is adjusted by a scheduled whole number of leap seconds as needed.
  ///
  /// The value in `ut` is appropriate for any calculation involving the Earth's rotation,
  /// such as calculating rise/set times, culumination, and anything involving apparent
  /// sidereal time.
  ///
  /// Before the era of atomic timekeeping, days based on the Earth's rotation
  /// were often known as *mean solar days*.
  final double ut;

  /// Terrestrial Time days since noon on January 1, 2000.
  ///
  /// Terrestrial Time is an atomic time scale defined as a number of days since noon on January 1, 2000.
  /// In this system, days are not based on Earth rotations, but instead by
  /// the number of elapsed [SI seconds](https://physics.nist.gov/cuu/Units/second.html)
  /// divided by 86400. Unlike `ut`, `tt` increases uniformly without adjustments
  /// for changes in the Earth's rotation.
  ///
  /// The value in `tt` is used for calculations of movements not involving the Earth's rotation,
  /// such as the orbits of planets around the Sun, or the Moon around the Earth.
  ///
  /// Historically, Terrestrial Time has also been known by the term *Ephemeris Time* (ET).
  final double tt;

  Time(this.ut, this.tt);

  // For internal use only. Used to optimize Earth tilt calculations.
  //double _psi = double.nan;

  // For internal use only. Used to optimize Earth tilt calculations.
  //double _eps = double.nan

  // For internal use only. Lazy-caches sidereal time (Earth rotation).
  //double st = double.nan

  static Time fromUniversalTime(double ut) =>
      Time(ut, ut /* turn this to terrestrialTime(ut) */);

  /// Creates a `Time` object from a UTC year, month, day, hour, minute and second.
  ///
  /// @param year The UTC year value.
  /// @param month The UTC month value 1..12.
  /// @param day The UTC day of the month 1..31.
  /// @param hour The UTC hour value 0..23.
  /// @param minute The UTC minute value 0..59.
  /// @param second The UTC second in the half-open range [0, 60).
  static Time fromTime(
    int year,
    int month,
    int day,
    int hour,
    int minute,
    double second,
  ) =>
      Time.fromUniversalTime(0);

  static final origin = DateTime(2000, 0, 1, 12, 0, 0).millisecondsSinceEpoch;

  /// Creates an `AstroTime` object from a `Date` object.
  ///
  /// @param d The date and time to be converted to AstroTime format.
  static Time fromDateTime(DateTime dateTime) {
    return Time.fromUniversalTime(
      (dateTime.millisecondsSinceEpoch - origin).toDouble() /
          millisecondsPerDay,
    );
  }

  /// Resolves this `Time` into year, month, day, hour, minute, second.
  DateTime toDateTime() {
    return DateTime.fromMillisecondsSinceEpoch(
      (ut * millisecondsPerDay + origin).round(),
    );
  }

  /// Converts this `Time` to the integer number of millseconds since 1970.
  int toMillisecondsSince1970() =>
      ((ut + 10957.5) * millisecondsPerDay).round();

  /// Converts this `Time` to ISO 8601 format, expressed in UTC with millisecond resolution.
  ///
  /// @return Example: "2019-08-30T17:45:22.763Z".
  @override
  String toString() => toDateTime().toIso8601String();

  /// Calculates the sum or difference of an [Time] with a specified floating point number of days.
  ///
  /// Sometimes we need to adjust a given [Time] value by a certain amount of time.
  /// This function adds the given real number of days in `days` to the date and time in this object.
  ///
  /// More precisely, the result's Universal Time field `ut` is exactly adjusted by `days` and
  /// the Terrestrial Time field `tt` is adjusted for the resulting UTC date and time,
  /// using a best-fit piecewise polynomial model devised by
  /// [Espenak and Meeus](https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html).
  ///
  /// @param days A floating point number of days by which to adjust `time`. May be negative, 0, or positive.
  /// @return A date and time that is conceptually equal to `time + days`.
  Time addDays(double days) => Time.fromUniversalTime(ut + days);

  double _julianCenturies() => tt / 36525.0;

  double _julianMillennia() => tt / dayPerMillennium;

  /// Compares the chronological order of two `Time` values.
  ///
  /// Two instances of `Time` can be compared for chronological order
  /// using the usual operators like `t1 < t2` or `t1 == t2`.
  @override
  int compareTo(Time other) => ut.compareTo(other.ut);

  /// Creates a `Time` object from a Terrestrial Time day value.
  ///
  /// This function can be used in rare cases where a time must be based
  /// on Terrestrial Time (TT) rather than Universal Time (UT).
  /// Most developers will want to use `Time(ut)` with a universal time
  /// instead of this function, because usually time is based on civil time adjusted
  /// by leap seconds to match the Earth's rotation, rather than the uniformly
  /// flowing TT used to calculate solar system dynamics. In rare cases
  /// where the caller already knows TT, this function is provided to create
  /// a `Time` value that can be passed to Astronomy Engine functions.
  ///
  /// @param tt The number of days after the J2000 epoch.
  static Time fromTerrestrialTime(double tt) => Time(/*universalTime(*/ tt, tt);
}
