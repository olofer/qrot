#pragma once

class csvreader
{
  public:
    csvreader(const std::string& csvname) : name_("") {
      std::ifstream thefile(csvname);
      if (!thefile.is_open()) return;
      int r = 0;
      std::string rowtext;
      if (!std::getline(thefile, rowtext))
        return;
      std::vector<double> rowdata;
      if (!parse_single_line_(rowtext, rowdata))
        return;
      for (size_t c = 0; c < rowdata.size(); c++) {
        cols_.emplace_back();
        cols_[cols_.size() - 1].push_back(rowdata[c]);
      }
      r++;
      while (std::getline(thefile, rowtext)) {
        if (!parse_single_line_(rowtext, rowdata))
          return;
        if (rowdata.size() != cols_.size())
          return;
        for (size_t c = 0; c < rowdata.size(); c++)
          cols_[c].push_back(rowdata[c]);
        r++;
      }
      name_ = csvname;
    }

    ~csvreader() {
    }

    bool isok() const {
      return (name_.size() > 0 && ncol() > 0);
    }

    bool isempty() const {
      return (cols_.size() == 0);
    }

    size_t ncol() const {
      return cols_.size();
    }

    size_t nrow() const {
      if (ncol() == 0)
        return 0;
      return cols_[0].size();
    }

    double get_element(size_t r, size_t c) const {
      return cols_[c][r];
    }

  private:

    bool parse_single_line_(const std::string& txt, 
                            std::vector<double>& val) const 
    {
      val.clear();
      size_t at = 0;
      while (at < txt.size()) {
        size_t csep = txt.find(',', at);
        if (csep == std::string::npos)
          csep = txt.size();
        std::string elemstr(txt.substr(at, csep - at));
        char* endptr = nullptr;
        const double elem = std::strtod(elemstr.c_str(), &endptr);
        if (elemstr.c_str() == endptr)
          return false;
        val.push_back(elem);
        at = csep + 1;
      }
      return true;
    }

    std::string name_;
    std::vector<std::vector<double>> cols_;
};
