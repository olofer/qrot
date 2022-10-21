//
// Demo of intermediate axis theorem (aka. Dzhanibekov effect, "tennis racket theorem")
// using quaternions to apply incremental rotations, and Cholesky factorization of the
// intertial matrix in the lab frame.
// 
// COMPILE: clang++ -o qrot -O2 -Wall qrot.cpp
//          g++ -o qrot -O2 -Wall qrot.cpp
// 
// USAGE: ./qrot ...options...
//
// List all available options by calling ./qrot without any option (the below list is incomplete).
//
// --in-file=FN1 
// --out-file=FN2
// --trace-file=FN3 
// --steps=S 
// --dt=T 
// --verbosity=L
// --vcm=VX,VY,VZ
// --omega=OX,OY,OZ
//
// FN1 is a filename, an existing textfile where each row is 4 numbers: mass, x, y, z (comma separated).
// FN2 is a filename, a file to which the final state will be written (same format as FN1).
// FN3 is a filename to which traces (invariants and some state variables) will be written.
// T is a timestep (positive floating point number).
// S is the integer number of integration steps to take.
// VX..VZ are components of the velocity of the center of mass
// OX..OZ are components of the angular velocity initial state
//
// The order of the options on the command line does not matter.
// The trace/log file will contain timestamp, kinetic energy, angular momentum, 
// tensor trace and determinant, angular velocity, and an orientation indicator.
//

// TODO: mechanical variational integrator? basic predictor-corrector integrator? RK45?

// EXAMPLE: 
// ./qrot --dt=.1e-5 --steps=25e6 --omega=0,10,0.01 --verbosity=1 --x3d-step=10000 --x3d-file=qrot-anim.html --x3d-template=example.html 

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
//#include <random>
//#include <chrono>
#include <map>
#include <string>
#include <cmath>

#include "quavec3.hpp"
#include "argsmap.hpp"
#include "csvreader.hpp"
#include "x3dframes.hpp"

struct sim_state {
  vec3 rcm;
  vec3 vcm;
  vec3 omega;
  vec3 omegadot;
  std::vector<vec3> r;
  std::vector<vec3> dr;
  std::vector<double> m;
  std::vector<vec3> v;
  double M;
  ten3 Icm;

  double t;
  vec3 p;
  vec3 L;
  double K;

  vec3 cmtot;      // cumulative cm translation since t = 0
  quaternion qtot; // cumulative rotation since t = 0

  void clear() {
    r.clear();
    dr.clear();
    m.clear();
    v.clear();
  }

  void zero_velocities() {
    vcm = {0.0, 0.0, 0.0};
    omegadot = {0.0, 0.0, 0.0};
  }

  void reset_accumulators() {
    cmtot.zero();
    qtot.set_unit();
  }

  int initialize(const std::vector<vec3>& coord, 
                 const std::vector<double>& mass) 
  {
    clear();

    if (coord.size() != mass.size())
      return r.size();

    for (size_t i = 0; i < coord.size(); i++) {
      if (mass[i] <= 0.0)
        continue;
      r.push_back(coord[i]);
      m.push_back(mass[i]);
    }

    const vec3 zero_vector = {0.0, 0.0, 0.0};
    dr.resize(r.size(), zero_vector);
    v.resize(r.size(), zero_vector);

    zero_velocities();
    calc_cm_from_r();
    zero_velocities();
    reset_accumulators();

    return r.size();
  }

  void calc_cm_from_r() {
    M = 0.0;
    vec3 sumw = {0.0, 0.0, 0.0};
    for (size_t i = 0; i < r.size(); i++) {
      vec3 w = r[i];
      w *= m[i];
      sumw += w;
      M += m[i];
    }
    rcm = sumw;
    rcm /= M;
  }

  void calc_dr() {
    for (size_t i = 0; i < r.size(); i++) {
      dr[i] = r[i];
      dr[i] -= rcm;
    }
  }

  void calc_r() {
    for (size_t i = 0; i < dr.size(); i++) {
      r[i] = dr[i];
      r[i] += rcm;
    }
  }

  void calc_v() {
    for (size_t i = 0; i < dr.size(); i++) {
      v[i] = omega;
      v[i] *= dr[i];
      v[i] += vcm;
    }
  }

  void calc_Ic() {
    this->Icm.recalc(this->dr, this->m);
  }

  vec3 calc_linear_momentum_from_v() const {
    vec3 p = {0.0, 0.0, 0.0};
    for (size_t i = 0; i < v.size(); i++) {
      p += (v[i] * m[i]);
    }
    return p;
  }

  vec3 calc_linear_momentum_from_cm() const {
    return (vcm * M);
  }

  vec3 calc_angular_momentum_from_cm() const {
    return (Icm * omega);
  }

  vec3 calc_angular_momentum_from_dr() const {
    vec3 L = {0.0, 0.0, 0.0};
    for (size_t i = 0; i < dr.size(); i++) {
      L += omega.cross(dr[i]) * m[i];
    }
    return L;
  }

  double calc_kinetic_energy_from_v() const {
    double Ek = 0.0;
    for (size_t i = 0; i < v.size(); i++)
      Ek += m[i] * v[i].dotself() / 2.0;
    return Ek;
  }

  double calc_kinetic_energy_from_cm() const {
    return (M * vcm.dotself() / 2.0) + omega.dot(Icm * omega) / 2.0; 
  }

  bool calc_omegadot() {
    if (!Icm.factorize())
      return false;
    vec3 y = (Icm * omega) * (-1.0);
    vec3 x = omega * y;
    omegadot = Icm.solve(x); // uses Cholesky factorization
    vec3 err = x - Icm * omegadot;
    const double errnorm = std::sqrt(err.dotself());
    if (errnorm > 1.0e-10) {
      std::cout << "|err| = " << std::sqrt(err.dotself()) << std::endl;
    }
    return true;
  }

  void cm_frame_prepare() {
    calc_cm_from_r();
    calc_dr();
    calc_Ic();
    p = calc_linear_momentum_from_cm();
    L = calc_angular_momentum_from_cm();
    K = calc_kinetic_energy_from_cm();
  }

  /* this is NOT an efficient integrator */
  void qevolve(double dt) {
    quaternion q, qconj;
    if (!calc_omegadot()) {
      std::cout << "calc_omegadot() failed" << std::endl;
      omegadot.zero();
    }
    const double norm_omega = std::sqrt(omega.dotself());
    if (norm_omega > 0.0) {
      const vec3 omegahat = omega / norm_omega;
      q.set_rotation(omegahat.x, omegahat.y, omegahat.z, dt * norm_omega);
      qconj = q.conj();
      for (size_t i = 0; i < dr.size(); i++) {
        quaternion xbar = q;
        xbar *= {dr[i].x, dr[i].y, dr[i].z, 0.0};
        xbar *= qconj;
        dr[i] = {xbar.x, xbar.y, xbar.z};
      }
      qtot.multiply_left(q);
    }
    t += dt;
    omega += (omegadot * dt);
    const vec3 drcm = vcm * dt; 
    rcm += drcm;
    calc_r();
    cmtot += drcm;
  }

  // ...
  // TODO: evolve which takes external forces as parameters
  // ...

  // backtransform angular velocity state to the original frame
  // (defined as the state where the accumulator was reset)
  vec3 omega_in_og_frame() const {
    quaternion obar = qtot.conj();
    obar *= {omega.x, omega.y, omega.z, 0.0};
    obar *= qtot;
    return {obar.x, obar.y, obar.z};
  }

  void undo_accumulated_transform() {
    for (size_t i = 0; i < r.size(); i++) {
      r[i] -= cmtot;
      quaternion xbar = qtot.conj();
      xbar *= {r[i].x, r[i].y, r[i].z, 0.0};
      xbar *= qtot;
      r[i] = {xbar.x, xbar.y, xbar.z};
    }
  }

  void print_invariants() const {
    std::cout << "K     = " << K << " J" << std::endl;
    std::cout << "L     = (" << L.x << ", " << L.y << ", " << L.z << ") kg*m*m/s" << std::endl;
    std::cout << "P     = (" << p.x << ", " << p.y << ", " << p.z << ") kg*m/s" << std::endl;
    std::cout << "tr(I) = " << Icm.trace() << " kg*m*m" << std::endl;
  }

  void print_cm_state(bool include_rcm = false, 
                      bool include_vcm = false) const 
  {
    std::cout << "omega = (" << omega.x << ", " << omega.y << ", " << omega.z << ") rad/s" << std::endl;
    if (include_rcm) {
      std::cout << "rcm   = (" << rcm.x << ", " << rcm.y << ", " << rcm.z << ") m" << std::endl;
    }
    if (include_vcm) {
      std::cout << "vcm   = (" << vcm.x << ", " << vcm.y << ", " << vcm.z << ") m/s" << std::endl;
    }
  }

  void print_points_and_masses() const {
    std::cout << "(x, y, z; m)" << std::endl;
    for (size_t i = 0; i < r.size(); i++) {
      std::cout << "(" << r[i].x << ", " << r[i].y << ", " << r[i].z << "; " << m[i] << ")" << std::endl;
    }
  }

  bool write_csv(const std::string& filename, 
                 int digits) const 
  {
    std::ofstream thefile(filename);
    if (!thefile.is_open())
      return false;
    thefile << std::setprecision(digits);
    for (size_t i = 0; i < r.size(); i++) {
      thefile << m[i] << ", " << r[i].x << ", " << r[i].y << ", " << r[i].z << std::endl;
    }
    return true;
  }

  bool read_csv(const std::string& filename) {
    csvreader f(filename);
    if (f.ncol() != 4)
      return false;
    std::vector<vec3> csv_xyz;
    std::vector<double> csv_mass;
    for (size_t r = 0; r < f.nrow(); r++) {
      csv_mass.push_back(f.get_element(r, 0));
      csv_xyz.push_back({f.get_element(r, 1), f.get_element(r, 2), f.get_element(r, 3)});
    }
    return (initialize(csv_xyz, csv_mass) == static_cast<int>(csv_mass.size()));
  }

};

struct sim_params
{
  double dt;
  int steps;
  vec3 omega_nought;
  vec3 vcm_nought;
  std::string infile;
  std::string outfile;
  std::string tracefile;
  int verbosity;
  int tracestep;
  int verbosestep;
  int csvdigits;
  std::string x3dfile;
  std::string x3dtemplate;
  int x3dstep;
  int x3ddigits;
  double x3dcycle;

  bool check_sanity() const {
    return (dt > 0.0 && 
            steps >= 0 && 
            tracestep > 0 && 
            verbosestep > 0 && 
            csvdigits > 0 &&
            x3ddigits > 0 &&
            x3dstep > 0 &&
            x3dcycle >= 0.0);
  }

  bool setup_from_argsmap(const argsmap& args) {

    int int_cast_steps = -1;
    std::vector<double> dummy;
    long parsed_long = -1;

    if (!args.value_as_scalar_double("--dt", dt))
      return false;

    if (!args.value_as_scalar_cast_from_double<int>("--steps", int_cast_steps))
      return false;
    steps = int_cast_steps;

    if (!args.value_as_scalar_cast_from_double<int>("--verbosity", int_cast_steps))
      return false;
    verbosity = int_cast_steps;

    if (!args.value_as_vector_double("--vcm", dummy))
      return false;
    if (dummy.size() != 3)
      return false;
    vcm_nought = {dummy[0], dummy[1], dummy[2]};

    if (!args.value_as_vector_double("--omega", dummy))
      return false;
    if (dummy.size() != 3)
      return false;
    omega_nought = {dummy[0], dummy[1], dummy[2]};

    if (!args.value_as_string("--in-file", infile))
      return false;

    if (!args.value_as_string("--out-file", outfile))
      return false;

    if (!args.value_as_string("--trace-file", tracefile))
      return false;

    if (!args.value_as_scalar_long("--trace-step", parsed_long))
      return false;
    tracestep = static_cast<int>(parsed_long);

    if (!args.value_as_scalar_long("--verbose-step", parsed_long))
      return false;
    verbosestep = static_cast<int>(parsed_long);

    if (!args.value_as_scalar_long("--csv-digits", parsed_long))
      return false;
    csvdigits = static_cast<int>(parsed_long);

    if (!args.value_as_string("--x3d-file", x3dfile))
      return false;

    if (!args.value_as_string("--x3d-template", x3dtemplate))
      return false;

    if (!args.value_as_scalar_long("--x3d-step", parsed_long))
      return false;
    x3dstep = static_cast<int>(parsed_long);

    if (!args.value_as_scalar_long("--x3d-digits", parsed_long))
      return false;
    x3ddigits = static_cast<int>(parsed_long);

    if (!args.value_as_scalar_double("--x3d-cycle", x3dcycle))
      return false;

    return check_sanity();
  }

  void argsmap_defaults(std::map<std::string, std::string>& defaults) const {
    defaults = {{"--steps",        "0"}, 
                {"--dt",           "1.0e-3"},
                {"--in-file",      "#example#"},
                {"--out-file",     "#none#"},
                {"--trace-file",   "#none#"},
                {"--verbosity",    "1"},
                {"--omega",        "0.0,0.0,1.0"},
                {"--vcm",          "0.0,0.0,0.0"},
                {"--trace-step",   "10"},
                {"--verbose-step", "100"},
                {"--csv-digits",   "16"},
                {"--x3d-file",     "#none#"},
                {"--x3d-template", "#none#"},
                {"--x3d-digits",   "8"},
                {"--x3d-step",     "1000"},
                {"--x3d-cycle",    "0.0"}};
  }

};

int main(int argc, const char **argv)
{
  std::map<std::string, std::string> defaults;
  sim_params P;

  P.argsmap_defaults(defaults);

  if (argc == 1) {
    std::cout << "usage: " << argv[0] << " param1=value1 ..." << std::endl;
    std::cout << "recognized parameters:" << std::endl;
    for (auto x : defaults) {
      std::cout << x.first << std::endl;
    }
    return 0;
  }

  argsmap arg(argc - 1, &argv[1]);

  if (arg.size_ignored() != 0) {
    std::cout << "*** arguments with incorrect pattern ***" << std::endl;
    arg.print_ignored();
    return 1;
  }

  arg.insert_if_not_set(defaults);

  if (!arg.all_recognized(defaults)) {
    std::cout << "there is at least one non-recognized parameter name" << std::endl;
    return 1;
  }

  if (!P.setup_from_argsmap(arg)) {
    std::cout << "*** key,value ***" << std::endl;
    arg.print();
    std::cout << "failed to setup parameter struct" << std::endl;
    return 1;
  }

  sim_state S;

  if (P.infile == "#example#") {
    const std::vector<vec3> coord = {{0.0, 1.0, 0.0}, 
                                     {0.0, -1.0, 0.0}, 
                                     {1.0, 0.0, 0.0}, 
                                     {-1.0, 0.0, 0.0}};
    const std::vector<double> mass = {1.0, 1.0, 2.0, 2.0};
    int n = S.initialize(coord, mass);
    if (n != static_cast<int>(coord.size())) {
      std::cout << "builtin example initialization failed" << std::endl;
      return 1;
    }
  } else {
    if (P.infile == "#none#") {
      std::cout << "infile is #none# - nothing to do" << std::endl;
      return 1; 
    }
    if (!S.read_csv(P.infile)) {
      std::cout << "failed to load CSV file (" << P.infile << ")" << std::endl;
      return 1;
    }
  }

  if (S.r.size() == 0) {
    std::cout << "the system of particles has not been loaded" << std::endl;
    return 1;
  }

  std::cout << "num. particles in system = " << S.r.size() << std::endl;
  std::cout << "total mass = " << S.M << " kg" << std::endl;

  const bool use_trace_file = (P.tracefile != "#none#");
  const bool use_out_file = (P.outfile != "#none#");
  const bool use_x3d_file = (P.x3dfile != "#none#");
  const bool use_x3d_template = (P.x3dtemplate != "#none#");

  x3dframes X3D(S.r.size());

  std::ofstream logfile;
  if (use_trace_file) {
    logfile.open(P.tracefile);
    if (!logfile.is_open()) {
      std::cout << "failed to open trace file (" << P.tracefile << ") for writing" << std::endl;
      return 1;
    }
    logfile << "step, time, K, Lx, Ly, Lz, omegax, omegay, omegaz, trI, logdetI, costheta, Kcheck" << std::endl;
    logfile << std::setprecision(P.csvdigits);
  }

  // Set initial condition velocity states
  S.omega = P.omega_nought;
  S.vcm = P.vcm_nought;
  S.t = 0.0;

  if (P.verbosity > 0) {
    std::cout << "*** initial state *** time = " << S.t << " s" << std::endl;
    S.cm_frame_prepare();
    S.print_invariants();
    S.print_cm_state();
    if (P.verbosity > 1) {
      S.print_points_and_masses();
    }
  }

  int64_t trace_counter = 1;
  int64_t verbose_counter = 1;
  int64_t x3d_counter = 1;

  for (int k = 0; k < P.steps; k++) {
    S.cm_frame_prepare();

    vec3 og_omega = S.omega_in_og_frame();
    vec3 ref_omega = P.omega_nought;
    const double cos_theta = og_omega.dot(ref_omega) / std::sqrt(og_omega.dotself() * ref_omega.dotself());

    if (--verbose_counter == 0) {

      if (P.verbosity > 3) {
        std::cout << "k (t) = " << k << " (" << S.t << " s) : ";
        std::cout << "K = " << S.K << " J" << std::endl;
      }

      if (P.verbosity > 2) {
        std::cout << "cos(th) = " << cos_theta << std::endl;
      }

      verbose_counter = P.verbosestep;
    }

    if (use_trace_file && --trace_counter == 0) {

      logfile << k << ", " << S.t << ", " << S.K;
      logfile << ", " << S.L.x << ", " << S.L.y << ", " << S.L.z; 
      logfile << ", " << S.omega.x << ", " << S.omega.y << ", " << S.omega.z; 
      logfile << ", " << S.Icm.trace() << ", " << S.Icm.logdet() << ", " << cos_theta;

      S.calc_v();
      const double Kcheck = S.calc_kinetic_energy_from_v();
      logfile << ", " << Kcheck << std::endl;

      trace_counter = P.tracestep;
    }

    if (use_x3d_file && --x3d_counter == 0) {
      if (!X3D.add(S.t, S.r)) {
        std::cout << "x3d add error; breaking simulation" << std::endl;
        break;
      }
      x3d_counter = P.x3dstep;
    }

    S.qevolve(P.dt);
  }

  if (use_trace_file) {
    logfile.close();
  }

  if (use_out_file) {
    if (!S.write_csv(P.outfile, P.csvdigits)) {
      std::cout << "failed to write CSV file (" << P.outfile << ")" << std::endl;
    }
  }

  if (use_x3d_file) {
    if (P.verbosity > 0) {
      std::cout << "saving X3D object with " << X3D.frames() << " snapshots to file (" << P.x3dfile << ")" << std::endl;
    }
    std::vector<double> radii;
    for (size_t i = 0; i < X3D.points(); i++) {
      const double refmass = 1.0;
      const double refradius = 0.150;
      radii.push_back(refradius * std::pow(S.m[i] / refmass, 1.0 / 3.0));
    }
    bool write_ok = false;
    X3D.set_digits(P.x3ddigits);
    std::string cmd = std::string(argv[0]);
    for (int i = 1; i < argc; i++)
      cmd += " " + std::string(argv[i]);
    X3D.set_caption("generated by: " + cmd);
    const std::string subttl = "n = " 
                               + std::to_string(S.r.size()) 
                               + " point masses (" 
                               + std::to_string(X3D.frames()) 
                               + " snapshots)";
    X3D.set_titles("qrot-x3dom", "qrot-x3dom", subttl);
    const double cycleInterval = (P.x3dcycle == 0.0 ? S.t : P.x3dcycle);
    if (use_x3d_template) {
      write_ok = X3D.write_with_template(P.x3dtemplate, P.x3dfile, radii, cycleInterval);
    } else {
      write_ok = X3D.write_without_template(P.x3dfile, radii, cycleInterval);
    }
    if (!write_ok) {
      std::cout << "failed to write X3D data to file (" << P.x3dfile << ")" << std::endl;
    }
  }

  if (P.verbosity > 0) {
    std::cout << "*** final state *** time = " << S.t << " s" << std::endl;
    S.cm_frame_prepare();
    S.print_invariants();
    S.print_cm_state();
    if (P.verbosity > 1) {
      S.print_points_and_masses();
      S.undo_accumulated_transform();
      S.print_points_and_masses();
    }
  }

  return 0;
}

