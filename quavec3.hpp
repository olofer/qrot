#pragma once

struct quaternion {

  double x;
  double y;
  double z;
  double s;

  quaternion conj() const {
    return {-x, -y, -z, s};
  }

  void conjugate() {
    x *= -1.0;
    y *= -1.0;
    z *= -1.0;
  }

  quaternion neg() const {
    return {-x, -y, -z, -s};
  }

  void negate() {
    (*this) *= -1.0;
  }

  void multiply_right(const quaternion& q) {
    *this *= q;
  }

  void multiply_left(const quaternion& q) {
    *this = q * (*this);
  }

  quaternion operator*(const quaternion& q) const {
    quaternion p = {x, y, z, s};
    p *= q;
    return p;
  }

  void operator*=(const quaternion& q) {
    const double a = y * q.z - z * q.y + s * q.x + q.s * x;
    const double b = - x * q.z + z * q.x + s * q.y + q.s * y;
    const double c = x * q.y - y * q.x + s * q.z + q.s * z;
    const double d = s * q.s - x * q.x - y * q.y - z * q.z;
    x = a;
    y = b;
    z = c;
    s = d;
  }

  double dotself() const {
    return (x * x + y * y + z * z + s * s);
  }

  void operator*=(double a) {
    x *= a;
    y *= a;
    z *= a;
    s *= a;
  }

  void operator/=(double a) {
    x /= a;
    y /= a;
    z /= a;
    s /= a;
  }

  void set_rotation(double nx, double ny, double nz, double theta) {
    const double thh = theta / 2.0;
    const double sinthh = std::sin(thh);
    const double costhh = std::cos(thh);
    x = nx * sinthh;
    y = ny * sinthh;
    z = nz * sinthh;
    s = costhh;
  }

  void renormalize() {
    (*this) /= std::sqrt(dotself());
  }

  void set_unit() {
    x = 0.0;
    y = 0.0;
    z = 0.0;
    s = 1.0;
  }

};

struct vec3 {

  double x;
  double y;
  double z;

  void zero() {
    x = 0;
    y = 0;
    z = 0;
  }

  void operator+=(const vec3& v) {
    x += v.x;
    y += v.y;
    z += v.z;
  }

  void operator-=(const vec3& v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
  }

  void operator*=(const vec3& v) {
    const double a = y * v.z - z * v.y;
    const double b = - x * v.z + z * v.x;
    const double c = x * v.y - y * v.x;
    x = a;
    y = b;
    z = c;
  }

  void operator*=(double a) {
    x *= a;
    y *= a;
    z *= a;
  }

  void operator/=(double a) {
    x /= a;
    y /= a;
    z /= a;
  }

  vec3 cross(const vec3& v) const {
    vec3 p = {x, y, z};
    p *= v;
    return p;
  }

  vec3 operator*(const vec3& v) const {
    return cross(v);
  }

  vec3 operator+(const vec3& v) const {
    vec3 p = {x, y, z};
    p += v;
    return p;
  }

  vec3 operator-(const vec3& v) const {
    vec3 p = {x, y, z};
    p -= v;
    return p;
  }

  vec3 operator*(double a) const {
    return {x * a, y * a, z * a};
  }

  vec3 operator/(double a) const {
    return {x / a, y / a, z / a};
  }

  double dot(const vec3& v) const {
    return (x * v.x + y * v.y + z * v.z);
  }

  double dotself() const {
    return (x * x + y * y + z * z);
  }

};

struct ten3 {

  double xx;
  double xy;
  double xz;
  double yy;
  double yz;
  double zz;

  double l11, l21, l31, l22, l32, l33;

  void clear() {
    xx = 0.0;
    xy = 0.0;
    xz = 0.0;
    yy = 0.0;
    yz = 0.0;
    zz = 0.0;
  }

  void clear_factor() {
    l11 = 0.0;
    l21 = 0.0;
    l31 = 0.0;
    l22 = 0.0;
    l32 = 0.0;
    l33 = 0.0;
  }

  void recalc(const std::vector<vec3>& dr, 
              const std::vector<double>& m)
  {
    clear();
    for (size_t i = 0; i < dr.size(); i++) {
      xx += m[i] * (dr[i].y * dr[i].y + dr[i].z * dr[i].z);
      yy += m[i] * (dr[i].x * dr[i].x + dr[i].z * dr[i].z);
      zz += m[i] * (dr[i].x * dr[i].x + dr[i].y * dr[i].y);
      xy -= m[i] * (dr[i].x * dr[i].y);
      xz -= m[i] * (dr[i].x * dr[i].z);
      yz -= m[i] * (dr[i].y * dr[i].z);
    }
  }

  // Assume the original calc is done in CM frame; and sum(m[i]) = M
  void translate(const vec3& R, double M) {
    xx += M * (R.y * R.y + R.z * R.z);
    yy += M * (R.x * R.x + R.z * R.z);
    zz += M * (R.x * R.x + R.y * R.y);
    xy -= M * (R.x * R.y);
    xz -= M * (R.x * R.z);
    yz -= M * (R.y * R.z);
  }

  vec3 operator*(const vec3& v) const {
    const double a = xx * v.x + xy * v.y + xz * v.z;
    const double b = xy * v.x + yy * v.y + yz * v.z;
    const double c = xz * v.x + yz * v.y + zz * v.z;
    return {a, b, c};
  }

  double trace() const {
    return (xx + yy + zz);
  }

  double logdet() const {
    return 2.0 * (std::log(l11) + std::log(l22) + std::log(l33));
  }

  bool factorize() {
    if (xx <= 0) return false;
    l11 = std::sqrt(xx);
    l21 = xy / l11;
    l31 = xz / l11;
    const double a = yy - l21 * l21;
    if (a <= 0) return false;
    l22 = std::sqrt(a);
    l32 = (yz - l21 * l31) / l22;
    const double b = zz - l31 * l31 - l32 * l32;
    if (b <= 0) return false;
    l33 = std::sqrt(b);
    return (l33 > 0);
  }

  void chol_solve_rhs(vec3& rhs) const {
    rhs.x /= l11;
    rhs.y = (rhs.y - l21 * rhs.x) / l22;
    rhs.z = (rhs.z - l31 * rhs.x - l32 * rhs.y) / l33;
  }

  void chol_transposed_solve_rhs(vec3& rhs) const {
    rhs.z /= l33;
    rhs.y = (rhs.y - l32 * rhs.z) / l22;
    rhs.x = (rhs.x - l21 * rhs.y - l31 * rhs.z) / l11;
  }

  void inplace_solve(vec3& x) const {
    // inplace transform x <- inv(I)*x, assuming the factorization has been done
    chol_solve_rhs(x);
    chol_transposed_solve_rhs(x);
  }

  // solve for x such that I*x = rhs
  vec3 solve(const vec3& rhs) const {
    vec3 x = rhs;
    inplace_solve(x);
    return x;
  }

};
