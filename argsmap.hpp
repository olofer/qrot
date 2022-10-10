#pragma once

class argsmap {
  public:
    argsmap(int argc, const char** argv) {
      for (int i = 0; i < argc; i++) {
        std::string str(argv[i]);
        auto nfirst = str.find('=');
        auto nlast = str.rfind('=');
        // require pattern: key=value (no spaces) with key and value both non-empty
        if (nfirst != std::string::npos && nlast == nfirst && nfirst != 0 && nfirst != str.size() - 1) {
          std::string key(str.substr(0, nfirst));
          std::string val(str.substr(nfirst + 1, str.size() - nfirst - 1));
          auto pos = kv.find(key);
          if (pos != kv.end()) {
            ignored.push_back(str); // ignored since key already exists
            continue;
          }
          kv[key] = val;
          continue;
        }
        ignored.push_back(str); // ignored due to pattern issues
      }
    }

    ~argsmap() {
    }

    void print() const {
      for (auto& item : kv) {
        std::cout << item.first << "=" << item.second << std::endl;
      }
    }

    void print_ignored() const {
      for (std::string item : ignored) {
        std::cout << item << std::endl;
      }
    }

    int size() const {
      return kv.size();
    }

    int size_ignored() const {
      return ignored.size();
    }

    void insert_if_not_set(const std::map<std::string, std::string>& defs) {
      for (auto& x : defs) {
        auto s = kv.find(x.first);
        if (s != kv.end())
          continue;
        kv[x.first] = x.second;
      }
    }

    bool all_recognized(const std::map<std::string, std::string>& recognized) const {
      for (auto& x : kv) {
        auto s = recognized.find(x.first);
        if (s == recognized.end())
          return false;
      }
      return true;
    }

    bool value_as_string(const std::string& key, std::string& value) const {
      auto s = kv.find(key);
      if (s == kv.end())
        return false;
      value = s->second;
      return true;
    }

    bool value_as_scalar_double(const std::string& key, double& result) const {
      auto s = kv.find(key);
      if (s == kv.end())
        return false;
      const char* str = s->second.c_str();
      char* endptr = nullptr;
      result = std::strtod(str, &endptr);
      return (endptr != str);
    }

    bool value_as_scalar_long(const std::string& key, long& result, int base = 0) const {
      auto s = kv.find(key);
      if (s == kv.end())
        return false;
      const char* str = s->second.c_str();
      char* endptr = nullptr;
      result = std::strtol(str, &endptr, base);
      return (endptr != str);
    }

    template <typename TV>
    bool value_as_scalar_cast_from_double(const std::string& key, TV& result) const {
      auto s = kv.find(key);
      if (s == kv.end())
        return false;
      const char* str = s->second.c_str();
      char* endptr = nullptr;
      const double r = std::strtod(str, &endptr);
      if (endptr == str)
        return false; 
      result = static_cast<TV>(r);
      return true;
    }

    // kv[key] should have the format "n1,n2,n3,..,nk" where n1 and others are interpretable as doubles
    bool value_as_vector_double(const std::string& key, std::vector<double>& result) const {
      result.clear();
      auto s = kv.find(key);
      if (s == kv.end())
        return false;
      size_t at = 0;
      while (at < s->second.size()) {
        size_t csep = s->second.find(',', at);
        if (csep == std::string::npos)
          csep = s->second.size();
        std::string elemstr(s->second.substr(at, csep - at));
        //std::cout << at << "-" << csep << ":" << elemstr << std::endl;
        char* endptr = nullptr;
        const double elem = std::strtod(elemstr.c_str(), &endptr);
        if (elemstr.c_str() == endptr)
          return false;
        result.push_back(elem);
        at = csep + 1;
      }
      return true;
    }

  private:
    std::map<std::string, std::string> kv;
    std::vector<std::string> ignored;
};
