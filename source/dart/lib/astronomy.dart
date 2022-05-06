library astronomy;

import 'dart:math';

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
    return 'TerseVector{' + ' x: $x,' + ' y: $y,' + ' z: $z,' + '}';
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
      'x': this.x,
      'y': this.y,
      'z': this.z,
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
