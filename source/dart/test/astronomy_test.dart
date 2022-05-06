import 'package:astronomy/astronomy.dart';
import 'package:test/test.dart';

void main() {
  test('TerseVector methods', () {
    final ones = TerseVector(x: 1.0, y: 1.0, z: 1.0);
    expect(ones, ones + TerseVector.zero);
    expect(ones, ones + TerseVector.zero);
    expect(TerseVector(x: 6.0, y: 8.0, z: 4.0), TerseVector(x: 3.0, y: 4.0, z: 2.0) * 2.0);
    expect(TerseVector(x: -1.5, y: 2.0, z: -1.0), TerseVector(x: -3.0, y: 4.0, z: -2.0) / 2.0);
    expect(29.0, TerseVector(x: -3.0, y: 4.0, z: -2.0).quadrature());
    expect(5.744562646538029, TerseVector(x: -2.0, y: -2.0, z: 5.0).magnitude());
  });
}
