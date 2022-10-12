#pragma once

class x3dframes
{
public:
  x3dframes(int npts) : digits_(6) {
    for (int i = 0; i < npts; i++)
      pos_.emplace_back();
  }

  ~x3dframes() {
  }

  bool check() const {
    size_t l = key_.size();
    for (size_t i = 0; i < pos_.size(); i++)
      if (pos_[i].size() != l)
        return false;
    return true;
  }

  bool add(double key, 
           const std::vector<vec3>& pos,
           bool posAndIncreasingKeys = true)
  {
    if (pos_.size() != pos.size())
      return false;
    if (posAndIncreasingKeys) {
      if (key_.size() > 0 && key <= key_[key_.size() - 1])
        return false;
      if (key < 0.0)
        return false;
    }
    key_.push_back(key);
    for (size_t i = 0; i < pos_.size(); i++)
      pos_[i].push_back(pos[i]);
    return check();
  }

  size_t frames() const {
    return key_.size();
  }

  size_t points() const {
    return pos_.size();
  }

  void set_digits(int d) {
    if (d > 0)
      digits_ = d;
  }

  // TODO: similarly have set_title(), set_subtitle(), set_caption()

  void write_sphere(std::ostream& os, 
                    const std::string& name, 
                    double radius) const 
  {
    const double cr = 0.5;
    const double cg = 0.5;
    const double cb = 1.0;
    os << "<transform DEF=\"" << name << "\">" << std::endl;
    os << "<shape>" << std::endl;
    os << "<appearance>" << std::endl;
    os << "<material diffuseColor=\"" << cr << " " << cg << " " << cb << "\"></material>" << std::endl;
    os << "</appearance>" << std::endl;
    os << "<sphere radius=\"" << radius << "\"></sphere>" << std::endl;
    os << "</shape>" << std::endl;
    os << "</transform>" << std::endl;
  }

  void write_timeSensor(std::ostream& os, 
                        const std::string& timename,
                        double cycleInterval, 
                        bool loop = true) const 
  {
    os << "<timeSensor DEF=\"" << timename << "\" cycleInterval=\"" 
       << cycleInterval << "\" loop=\"" 
       << (loop ? "true" : "false") << "\">" << "</timeSensor>" << std::endl;
  }

  void write_PositionInterpolator(std::ostream& os,
                                  const std::string& name, 
                                  int i) const 
  {
    // <PositionInterpolator DEF="move1" key="0 0.5 1" keyValue="0 0 0  0 3 0  0 0 0"></PositionInterpolator>
    const double keymax = (frames() == 0 ? 1.0 : key_[frames() - 1]);
    os << std::setprecision(digits_);
    os << "<PositionInterpolator DEF=\"" << name << "\" key=\"";
    for (size_t r = 0; r < frames(); r++)
      os << " " << key_[r] / keymax;
    os << "\" keyValue=\"";
    for (size_t r = 0; r < frames(); r++)
      os << "  " << pos_[i][r].x << " " << pos_[i][r].y << " " << pos_[i][r].z;
    os << "\"></PositionInterpolator>" << std::endl;
  }

  void write_Route(std::ostream& os, 
                   const std::string& timename,
                   const std::string& interpname, 
                   const std::string& spherename) const 
  {
    os << "<Route fromNode=\"" << timename << "\" fromField=\"fraction_changed\" toNode=\"" 
       << interpname << "\" toField=\"set_fraction\">" << "</Route>" << std::endl;
    os << "<Route fromNode=\"" << interpname << "\" fromField=\"value_changed\" toNode=\"" 
       << spherename << "\" toField=\"translation\">" << "</Route>" << std::endl;
  }

  void write(std::ostream& os, 
             const std::vector<double>& radii,
             double cycleInterval) const 
  {
    const std::string spherenamebase = "point";
    const double default_radius = 0.150;
    for (size_t i = 0; i < points(); i++) {
      const double radi = (radii.size() == points() ? radii[i] : default_radius);
      write_sphere(os, spherenamebase + std::to_string(i), radi);
      os << std::endl;
    }

    const std::string timername = "time";
    write_timeSensor(os, timername, cycleInterval, true);
    os << std::endl;

    const std::string movenamebase = "move";
    for (size_t i = 0; i < points(); i++) {
      write_PositionInterpolator(os, movenamebase + std::to_string(i), i);
      os << std::endl;
      write_Route(os, timername, movenamebase + std::to_string(i), spherenamebase + std::to_string(i));
      os << std::endl;
    }
  }

  bool write_without_template(const std::string& x3dfile,
                              const std::vector<double>& radii, 
                              double cycleInterval) const 
  {
    std::ofstream outfile(x3dfile);
    if (!outfile.is_open())
      return false;
    write(outfile, radii, cycleInterval);
    outfile.close();
    return true;
  }

  // TODO: add option to replace title + subtitle + caption in the template
  bool write_with_template(const std::string& templateFile, 
                           const std::string& htmlFile,
                           const std::vector<double>& radii,
                           double cycleInterval) const 
  {
    std::ifstream infile(templateFile);
    if (!infile.is_open()) {
      //std::cout << "failed to open \"" << templateFile << "\" for reading" << std::endl;
      return false;
    }

    std::ofstream outfile(htmlFile);
    if (!outfile.is_open()) {
      //std::cout << "failed to open \"" << htmlFile << "\" for writing" << std::endl;
      return false;
    }

    std::string therow;
    int state = 0;

    while (std::getline(infile, therow)) {
      if (state == 0 && therow.find("<!-- PASTE BEGIN -->") != std::string::npos) {
        state = 1;
        continue;
      }
      if (state == 1 && therow.find("<!-- PASTE END -->") != std::string::npos) {
        //outfile << "<!-- PASTE BEGIN -->" << std::endl;
        write(outfile, radii, cycleInterval);
        //outfile << "<!-- PASTE END -->" << std::endl;
        state = 0;
        continue;
      }
      if (state == 1)
        continue;
      outfile << therow << std::endl;
    }

    infile.close();
    outfile.close();
    return true;
  }

private:
  int digits_;
  std::vector<double> key_;
  std::vector<std::vector<vec3>> pos_;
};
